""" Module astropak.catalogs:
    Catalog parsing and presentation. One class per catalog.
    As of March 2022, only ATLAS refcat2 catalog is handled.
"""

__author__ = "Eric Dose, Albuquerque"

# Python core packages:
import os
from datetime import datetime, timezone
from math import floor, ceil

# External packages:
import pandas as pd
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import Angle
from astropy.stats import circmean

# EVD packages:
from .reference import DAYS_PER_YEAR_NOMINAL


__________ATLAS_REFCAT2_NEW________________________________________________________ = 0

from .util import count_files_immediate

ATLAS_REFCAT2_DIRECTORY = 'D:/Astro/Catalogs/ATLAS-refcat2/mag-0-16/'
ATLAS_REFCAT2_EPOCH_UTC = (datetime(2015, 1, 1) +
                           (datetime(2016, 1, 1) - datetime(2015, 1, 1)) / 2.0) \
    .replace(tzinfo=timezone.utc)  # refcat2 proper-motion catalog_epoch is 2015.5


class AtlasRefcat2:
    """Data container for stars found in ATLAS refcat2 catalog, within a RA,Dec rectangle.
    """
    def __init__(self, atlas_top_directory, subdir_listing='subdirs.txt',
                 ra_deg_range=None, dec_deg_range=None, target_epoch=None,
                 overlap_distance=10, sort_by='ra'):
        """Constructor.
        :param atlas_top_directory: All Atlas subdirectories (e.g., '/mag-0-16/'), each of which contains
        one magnitude range ('chunk') of the catalog.
        :param subdir_listing: name of file holding list of ATLAS subdirectories, each line as:
        min_gri_mag : subdir name, e.g., '16 : mag-0-16'. [string]
        ~64800 text (.csv) files, are found directly under this catalog top directory.
        :param ra_deg_range:
        :param dec_deg_range:
        :param target_epoch:
        :param overlap_distance:
        :param sort_by: must be 'ra', 'dec', 'r', or None. [string or None]
        """
        if not (os.path.exists(atlas_top_directory) and os.path.isdir(atlas_top_directory)):
            raise NotADirectoryError('ATLAS refcat2 top directory \'' +
                                     atlas_top_directory + '\' not found.')
        if isinstance(target_epoch, Time):
            self.target_epoch = target_epoch
        elif isinstance(target_epoch, datetime):
            self.target_epoch = Time(target_epoch, scale='utc')
        else:
            raise ValueError('Parameter \'target_epoch\' is required and '
                             'must be a py datetime or astropy Time object.')
        if sort_by is not None:
            if sort_by.lower() not in ['ra', 'dec', 'r']:
                raise ValueError('Parameter \'sort_by\' must be None, \'ra\', \'dec\', or \'r\'.')

        # Load requested stars from catalog and pre-process:
        self.df_subcat = self._locate_subdirs(atlas_top_directory, subdir_listing)
        df = self._get_stars_by_index(self.df_subcat, ra_deg_range, dec_deg_range)
        df = self._trim_to_ra_dec_range(df, ra_deg_range, dec_deg_range)
        df = self._remove_overlapping(df, overlap_distance)
        df = self._add_new_columns(df)
        df = self._update_epoch(df, self.target_epoch)
        self.df_all = self._sort_by(df, sort_by)
        self.df_selected = self.df_all.copy()  # for later star selections by user.

    @staticmethod
    def _locate_subdirs(atlas_top_directory, subdir_listing):
        # Locate subcatalog directories, make dataframe.
        fullpath = os.path.join(atlas_top_directory, subdir_listing)
        df_subcat = pd.read_csv(fullpath, sep=':', engine='python', header=None,
                                skip_blank_lines=True, on_bad_lines='error')
        df_subcat.columns = ['gri_max', 'Subdir']
        df_subcat['gri_max'] = [int(gri) for gri in df_subcat['gri_max']]
        df_subcat['Subdir'] = [subdir.strip() for subdir in df_subcat['Subdir']]
        df_subcat['Path'] = [os.path.join(atlas_top_directory, subdir)
                             for subdir in df_subcat['Subdir']]
        df_subcat['Nfiles'] = [count_files_immediate(subdir) for subdir in df_subcat['Path']]
        df_subcat.index = list(df_subcat['gri_max'].values)

        # Verify that each listed subdirectory appears to be a valid ATLAS refcat2 directory:
        for _, row in df_subcat.iterrows():
            subdir_path = os.path.join(atlas_top_directory, row['Subdir'])
            if not (64798 <= row['Nfiles'] <= 64802):
                raise ValueError('Subdirectory \'' + subdir_path +
                                 '\' is not a valid ATLAS refcat2 directory.')
        return df_subcat

    @staticmethod
    def _get_stars_by_index(df_subcat, ra_deg_range, dec_deg_range):
        """Get star data from catalog within RA and Dec range, return as a dataframe.
        :param df_subcat:
        :param ra_deg_range: (min RA, max RA) as tuple, in degrees. Graceful handling of near-zero RA
        (RA crossover) is guaranteed if the RA range is less than 180 degrees. [2-tuple of floats]
        :param dec_deg_range: (min Dec, max Dec) as tuple, in degrees. [2-tuple of floats]
        correct handling of proper motions. Typically, this should be the approximate date and time
        at which the relevant images were taken. [astropy Time or py datetime object]
        is counted toward overlap detection. Typically, the minimum practical distance between reference
        stars. Default of 10 arcseconds. [float]
        [string, or None]
        :return: dataframe of ATLAS refcat2 data for relevant stars. [pandas DataFrame]
        """
        # Make lists of all integer RA and Dec indices to retrieve from catalog:
        ra_index_list = _make_ra_index_list(ra_deg_range, include_max_bound=False)
        dec_index_list = _make_dec_index_list(dec_deg_range, include_max_bound=False)

        # Make pandas DataFrame with catalog data over superset of requested RA, Dec range:
        df_list = []
        for _, row in df_subcat.iterrows():
            for ra_index in ra_index_list:
                for dec_index in dec_index_list:
                    df_degsq = AtlasRefcat2._read_one_subdir_one_file(row['gri_max'], row['Path'],
                                                                      ra_index, dec_index)
                    df_list.append(df_degsq)
        df = pd.DataFrame(pd.concat(df_list, ignore_index=True))  # new index of unique integers
        return df

    @staticmethod
    def _read_one_subdir_one_file(gri_max, subdir_path, ra_index, dec_index, max_stars=None):
        """
        :param gri_max:
        :param subdir_path:
        :param ra_index:
        :param dec_index:
        :param max_stars:
        :return:
        """
        filename = '{:03d}'.format(ra_index) + '{:+03d}'.format(dec_index) + '.rc2'
        fullpath = os.path.join(subdir_path, filename)
        df = pd.read_csv(fullpath, sep=',', engine='python', header=None,
                         skip_blank_lines=True, error_bad_lines=False, nrows=max_stars,
                         usecols=[0, 1, 4, 5, 6, 7,
                                  8, 9, 10, 11, 12, 13,
                                  14, 16, 18, 19, 20,
                                  21, 22, 25, 26, 29, 30, 33, 34], prefix='col')  # col nums: zero-origin.
        df.columns = ['RA_deg', 'Dec_deg', 'PM_ra', 'dPM_ra', 'PM_dec', 'dPM_dec',
                      'G_gaia', 'dG_gaia', 'BP_gaia', 'dBP_gaia', 'RP_gaia', 'dRP_gaia',
                      'T_eff', 'dupvar', 'RP1', 'R1', 'R10',
                      'g', 'dg', 'r', 'dr', 'i', 'di', 'z', 'dz']
        df['RA_deg'] *= 0.00000001
        df['Dec_deg'] *= 0.00000001
        df['PM_ra'] *= 0.00001    # proper motion in arcsec/year
        df['dPM_ra'] *= 0.00001   # uncert in PM, arcsec/year
        df['PM_dec'] *= 0.00001   # proper motion in arcsec/year
        df['dPM_dec'] *= 0.00001  # uncert in PM, arcsec/year
        df['G_gaia'] *= 0.001  # in magnitudes; dG_gaia remains in millimagnitudes
        df['BP_gaia'] *= 0.001  # in magnitudes; dBP_gaia remains in millimagnitudes
        df['RP_gaia'] *= 0.001  # in magnitudes; dRP_gaia remains in millimagnitudes
        df['RP1'] = [None if rp1 == 999 else rp1 / 10.0 for rp1 in df['RP1']]  # radius in arcseconds
        df['R1'] = [None if r1 == 999 else r1 / 10.0 for r1 in df['R1']]       # "
        df['R10'] = [None if r10 == 999 else r10 / 10.0 for r10 in df['R10']]  # "
        df['g'] *= 0.001  # in magnitudes; dg remains in millimagnitudes
        df['r'] *= 0.001  # in magnitudes; dr remains in millimagnitudes
        df['i'] *= 0.001  # in magnitudes; di remains in millimagnitudes
        df['z'] *= 0.001  # in magnitudes; dz remains in millimagnitudes

        # Add CatalogID, unique to each star:
        gri_max_prefix = '{:02d}'.format(gri_max) + ':'
        id_prefix = '{:03d}'.format(ra_index) + '{:+03d}'.format(dec_index) + '_'
        id_list = [gri_max_prefix + id_prefix + '{:0>6d}'.format(i + 1)
                   for i in range(len(df))]  # unique in entire catalog.
        df.insert(0, 'CatalogID', id_list)
        return df

    @staticmethod
    def _trim_to_ra_dec_range(df, ra_deg_range, dec_deg_range):
        """Remove rows of dataframe that lie outside RA range or Dec range. Respects RA zero-crossings.
        :param df: input dataframe of stars (1 per row). [pandas Dataframe]
        :param ra_deg_range: (min RA, max RA) as tuple, in degrees. Graceful handling of near-zero RA
        (RA crossover) is guaranteed if the RA range is less than 180 degrees. [2-tuple of floats]
        :param dec_deg_range: (min Dec, max Dec) as tuple, in degrees. [2-tuple of floats]
        :return: dataframe of stars (1 per row) all being inside RA and Dec ranges. [pandas Dataframe]
        """
        # Detect RA out of range, managing gracefully any RA zero-crossing:
        ra_range = Angle(ra_deg_range * u.deg)
        ra_center = circmean(ra_range)
        wrap_angle = Angle(180.0 * u.deg + ra_center)
        ra_wrapped = ra_range.wrap_at(wrap_angle)
        ra_deg_wrapped = ra_wrapped.degree
        ra_cat_wrapped = Angle(df.loc[:, 'RA_deg'] * u.deg).wrap_at(wrap_angle).degree
        ra_too_low = ra_cat_wrapped < ra_deg_wrapped[0]
        ra_too_high = ra_cat_wrapped > ra_deg_wrapped[1]

        # Detect Dec out of range:
        dec_cat = df.loc[:, 'Dec_deg']
        dec_too_low = dec_cat < dec_deg_range[0]
        dec_too_high = dec_cat > dec_deg_range[1]

        # Trim:
        is_outside_range = ra_too_low | ra_too_high | dec_too_low | dec_too_high
        df_trimmed = df.loc[~is_outside_range, :]
        return df_trimmed

    @staticmethod
    def _remove_overlapping(df, overlap_distance):
        """Remove rows (stars) which have other nearby stars (which will interfere with photometry).
        Current implementation uses only criterion based on 'RP1' = 'radius point-one', or
        the minimum radius inside which other light sources comprise 0.1% of the target's own flux.
        Rows (stars) with missing/null rp1 entries in catalog are assumed to have no nearby flux at all.
            :param df: input dataframe. [pandas Dataframe]
            :param overlap_distance: minimum allowable distance to the nearest other flux,
            in arcseconds. If None (or not specified), then ANY mention of nearby flux
            in this catalog will remove that star (very strict but common case). [float]
            :return: dataframe with stars removed that are indicated to overlap with other light sources.
            [pandas Dataframe]
        """
        rp1_too_close = pd.Series([False if pd.isnull(rp1)
                                   else ((overlap_distance is None) or (rp1 < overlap_distance))
                                   for rp1 in df['RP1']])
        return df.loc[list(~rp1_too_close), :].copy()

    @staticmethod
    def _add_new_columns(df):
        """Add new columns with derived data. Currently, new columns are 'BminusV' and 'ri_color'.
        :param df: input dataframe. [pandas Dataframe]
        """
        df.loc[:, 'BminusV'] = [(0.830 * g - 0.803 * r) for (g, r) in zip(df['g'], df['r'])]
        df.loc[:, 'APASS_R'] = [(0.950 * r + 0.05 * i) for (r, i) in zip(df['r'], df['i'])]
        df.loc[:, 'ri_color'] = [(r - i) for (r, i) in zip(df['r'], df['i'])]
        return df.copy()

    @staticmethod
    def _update_epoch(df, target_epoch):
        """Update RA and Dec of each star for proper motion between catalog epoch and target epoch.
        :param df: input dataframe. [pandas Dataframe]
        :param target_epoch: time for which this catalog retrieval is effective, to ensure
        correct handling of proper motions. Typically, this should be the approximate date and time
        at which the relevant images were taken. [astropy Time or py datetime object]
        :return:
        """
        delta_years = (target_epoch - Time(ATLAS_REFCAT2_EPOCH_UTC)).jd / DAYS_PER_YEAR_NOMINAL
        ra_date = [(ra_epoch + delta_years * pm_ra / 3600.0) % 360
                   for (ra_epoch, pm_ra) in zip(df['RA_deg'], df['PM_ra'])]
        dec_date = [dec_epoch + delta_years * pm_dec / 3600.0
                    for (dec_epoch, pm_dec) in zip(df['Dec_deg'], df['PM_dec'])]
        df.loc[:, 'RA_deg'] = ra_date
        df.loc[:, 'Dec_deg'] = dec_date
        return df

    @staticmethod
    def _sort_by(df, sort_by):
        """Sort stars (df rows) by one of: 'ra', 'dec', or 'r'; or None to skip sorting."""
        if sort_by is None:
            return df
        if sort_by.lower() == 'ra':
            ra_angles = Angle(df.loc[:, 'RA_deg'] * u.deg)
            ra_center = Angle(circmean(ra_angles))
            wrap_angle = Angle(180.0 * u.deg + ra_center)
            df['RA_deg_wrapped'] = ra_angles.wrap_at(wrap_angle).degree
            df = df.copy().sort_values(by='RA_deg_wrapped')
            del df['RA_deg_wrapped']
            return df
        if sort_by.lower() == 'dec':
            return df.copy().sort_values(by='Dec_deg')
        if sort_by.lower() == 'r':
            return df.copy().sort_values(by='r')
        raise ValueError('Parameter \'sort by\' must be \'ra\', \'dec\', \'r\', or None.')

    def _select_on(self, column_name, minimum=None, maximum=None):
        """Select rows (stars) on minimum and/or maximum value in column_name.
           Modifies self.df_selected in place; no return value.
           Either minimum or maximum, or both, should be given.
           Most users will prefer convenience functions (e.g., .select_on_g_mag()) instead.
           :param column_name: name of column on which to make selection. [string]
           :param minimum: minimum value to cause row (star) to be retained. [float]
           :param maximum: maximum value to cause row (star) to be retained. [float]
           :return: [None] Modifies df_selected in place.
           """
        at_least_minimum = len(self.df_selected) * [True]  # default value
        at_most_maximum = at_least_minimum.copy()          # default value
        if minimum is not None:
            at_least_minimum = [(x >= minimum) for x in self.df_selected[column_name]]
        if maximum is not None:
            at_most_maximum = [(x <= maximum) for x in self.df_selected[column_name]]
        rows_to_keep = [(mn and mx) for (mn, mx) in zip(at_least_minimum, at_most_maximum)]
        self.df_selected = self.df_selected.loc[rows_to_keep, :]

    def select_on_g_mag(self, min_g_mag=None, max_g_mag=None):
        """Select rows (stars) on minimum and/or maximum Sloan g magnitude."""
        self._select_on(column_name='g', minimum=min_g_mag, maximum=max_g_mag)

    def select_on_r_mag(self, min_r_mag=None, max_r_mag=None):
        """Select rows (stars) on minimum and/or maximum Sloan r magnitude."""
        self._select_on(column_name='r', minimum=min_r_mag, maximum=max_r_mag)

    def select_on_i_mag(self, min_i_mag=None, max_i_mag=None):
        """Select rows (stars) on minimum and/or maximum Sloan i magnitude."""
        self._select_on(column_name='i', minimum=min_i_mag, maximum=max_i_mag)

    def select_on_g_uncert(self, min_g_uncert=None, max_g_uncert=None):
        """Select rows (stars) on (minimum and/or) maximum Sloan g uncertainty in millimagnitude."""
        self._select_on(column_name='dg', minimum=min_g_uncert, maximum=max_g_uncert)

    def select_on_r_uncert(self, min_r_uncert=None, max_r_uncert=None):
        """Select rows (stars) on (minimum and/or) maximum Sloan r uncertainty in millimagnitude."""
        self._select_on(column_name='dr', minimum=min_r_uncert, maximum=max_r_uncert)

    def select_on_i_uncert(self, min_i_uncert=None, max_i_uncert=None):
        """Select rows (stars) on (minimum and/or) maximum Sloan u uncertainty in millimagnitude."""
        self._select_on(column_name='di', minimum=min_i_uncert, maximum=max_i_uncert)

    def select_on_bv_color(self, min_bv_color=None, max_bv_color=None):
        """Select rows (stars) on minimum and/or maximum B-V color, in magnitudes."""
        self._select_on(column_name='BminusV', minimum=min_bv_color, maximum=max_bv_color)

    def select_on_ri_color(self, min_ri_color=None, max_ri_color=None):
        """Select rows (stars) on minimum and/or maximum Sloan r-i color, in magnitudes."""
        self._select_on(column_name='ri_color', minimum=min_ri_color, maximum=max_ri_color)


__________GENERAL_CATALOG_FUNCTIONS_____________________________________________ = 0


def _make_ra_index_list(ra_deg_range, include_max_bound=False):
    """Accept tuple (RA_min, RA_max), return all integer RA values covering the range.
       Handles RA zero-crossing gracefully. Assumes RA range less than 180 deg.
    :param ra_deg_range: (RA_min, RA_max) in degrees. [tuple of 2 floats]
    :param include_max_bound: if True, also include the first RA integer greater than input RA_min. [bool]
    :return: all integer RA values, in degrees, that cover RA range. [list of integers]
    """
    # Generate list of integer RA values needed to cover the entire RA range:
    ra_range = Angle(ra_deg_range * u.deg)
    ra_center = circmean(ra_range)
    wrap_angle = Angle(180.0 * u.deg + ra_center)
    ra_deg_wrapped = ra_range.wrap_at(wrap_angle).degree
    ra_index_first = floor(ra_deg_wrapped[0])
    ra_index_last = ceil(ra_deg_wrapped[1]) if include_max_bound else floor(ra_deg_wrapped[1])
    ra_index_list = [ra % 360 for ra in range(ra_index_first, ra_index_last + 1)]
    return ra_index_list


def _make_dec_index_list(dec_deg_range, include_max_bound=False):
    """Accept tuple (Dec_min, Dec_max), return all integer Declination values covering the range.
    :param dec_deg_range: (Dec_min, Dec_max) in degrees. [tuple of 2 floats]
    :param include_max_bound: if True, also include the first Dec integer greater than input Dec_min. [bool]
    :return: all integer Dec values, in degrees, that cover Dec range. [list of integers]
    """
    dec_deg_first = floor(dec_deg_range[0])
    dec_deg_last = ceil(dec_deg_range[1]) if include_max_bound else floor(dec_deg_range[1])
    dec_index_list = [dec for dec in range(dec_deg_first, dec_deg_last + 1)]
    return dec_index_list
