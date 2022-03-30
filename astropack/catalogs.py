"""Catalog parsing and presentation.

    As of March 2022, only ATLAS refcat2 catalog is handled."""

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


__all__ = ['AtlasRefcat2']


__________ATLAS_REFCAT2_NEW________________________________________________________ = 0

from .util import count_files_immediate

# ATLAS_REFCAT2_DIRECTORY = 'D:/Astro/Catalogs/ATLAS-refcat2/mag-0-16/'
ATLAS_REFCAT2_EPOCH_UTC = (datetime(2015, 1, 1) +
                           (datetime(2016, 1, 1) - datetime(2015, 1, 1)) / 2.0) \
    .replace(tzinfo=timezone.utc)  # refcat2 proper-motion catalog_epoch is 2015.5


class AtlasRefcat2:
    """Represents data for stars found in ATLAS refcat2 catalog,
    within a user-specified rectangle in RA and Dec.

    Parameters
    ----------

    atlas_top_directory : str
        Full path of the top directory where the local copy of ATLAS refcat2 catalog
        resides.

    index_filename : str
        Name of text file housing a minimal index to the catalog's component
        subdirectories. The user must usually construct this file by hand.
        Default 'subdirs.txt'.

    ra_deg_range : tuple of 2 float
        Minimum and maximum Right Ascension of range from which to select catalog stars.
        Tuple is (ra_min, ra_max), in degrees.

        RA zero-crossing is handled gracefully, e.g.,
        ``ra_deg_range`` of (359, 2) will include stars with RA in range [359, 360]
        as well as those with RA in range [0, 2].

    dec_deg_range : tuple of 2 float
        Minimum and maximum Declination of range from which to select catalog stars.
        Tuple is (ra_min, dec_max), in degrees.

    target_epoch : |py.datetime| or |Time|
        Date for which catalog RA and Dec should be updated for proper motion.

        Typically, this will be the date when images were taken (images to which
        catalog data will be applied).

    overlap_distance : float, optional
        Minimum distance from a target light source (star, asteroid) that a nearby
        lightsource can be considered not to interfere photometrically. In arcseconds.
        Default 10.

    sort_by : str {'ra', 'dec', 'r'} or None, optional
        DataFrame column on which to sort, ascending. Default None.

    Attributes
    ----------
    target_epoch : |Time|
        Date to which catalog RA and Dec should be updated for proper motion,
        from input parameter `target_epoch`.

    df_all : |DataFrame|
        Pandas Dataframe of stars from ATLAS refcat2 catalog within RA and Dec ranges.
        Intended to be a persistent record of the originally extracted data.

    df_selected : |DataFrame|
        Independent copy of df_all.
        Intended to be further selected by user to a desired subset of ``df_all``
        after class instance construction, via user calls to ``select_on_*()`` methods.
    """

    def __init__(self, atlas_top_directory, index_filename='subdirs.txt',
                 ra_deg_range=None, dec_deg_range=None, target_epoch=None,
                 overlap_distance=10, sort_by='ra'):
        if not (os.path.exists(atlas_top_directory) and
                os.path.isdir(atlas_top_directory)):
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
                raise ValueError('Parameter \'sort_by\' must be None, '
                                 '\'ra\', \'dec\', or \'r\'.')

        # Load requested stars from catalog and pre-process:
        self.df_subcat = self._locate_subdirs(atlas_top_directory, index_filename)
        df = self._get_stars_by_index(self.df_subcat, ra_deg_range, dec_deg_range)
        df = self._trim_to_ra_dec_range(df, ra_deg_range, dec_deg_range)
        df = self._remove_overlapping(df, overlap_distance)
        df = self._add_new_columns(df)
        df = self._update_epoch(df, self.target_epoch)
        self.df_all = self._sort_by(df, sort_by)
        self.df_selected = self.df_all.copy()  # for later star selections by user.

    @staticmethod
    def _locate_subdirs(atlas_top_directory, subdir_listing):
        """Locate subcatalog directories, make dataframe."""
        fullpath = os.path.join(atlas_top_directory, subdir_listing)
        df_subcat = pd.read_csv(fullpath, sep=':', engine='python', header=None,
                                skip_blank_lines=True, on_bad_lines='error')
        df_subcat.columns = ['gri_max', 'Subdir']
        df_subcat['gri_max'] = [int(gri) for gri in df_subcat['gri_max']]
        df_subcat['Subdir'] = [subdir.strip() for subdir in df_subcat['Subdir']]
        df_subcat['Path'] = [os.path.join(atlas_top_directory, subdir)
                             for subdir in df_subcat['Subdir']]
        df_subcat['Nfiles'] = [count_files_immediate(subdir)
                               for subdir in df_subcat['Path']]
        df_subcat.index = list(df_subcat['gri_max'].values)

        # Verify listed subdirectories appears to be valid ATLAS refcat2 directories:
        for _, row in df_subcat.iterrows():
            subdir_path = os.path.join(atlas_top_directory, row['Subdir'])
            if not (64798 <= row['Nfiles'] <= 64802):
                raise ValueError('Subdirectory \'' + subdir_path +
                                 '\' is not a valid ATLAS refcat2 directory.')
        return df_subcat

    @staticmethod
    def _get_stars_by_index(df_subcat, ra_deg_range, dec_deg_range):
        """Get catalog star data within RA and Dec range, return as a dataframe."""
        # Make lists of all integer RA and Dec indices to retrieve from catalog:
        ra_index_list = _make_ra_index_list(ra_deg_range, include_max_bound=False)
        dec_index_list = _make_dec_index_list(dec_deg_range, include_max_bound=False)

        # Make DataFrame with catalog data over superset of requested RA, Dec range:
        df_list = []
        for _, row in df_subcat.iterrows():
            for ra_index in ra_index_list:
                for dec_index in dec_index_list:
                    df_degsq = AtlasRefcat2._read_one_subdir_one_file(row['gri_max'],
                                                                      row['Path'],
                                                                      ra_index,
                                                                      dec_index)
                    df_list.append(df_degsq)
        df = pd.DataFrame(pd.concat(df_list, ignore_index=True))  # new index
        return df

    @staticmethod
    def _read_one_subdir_one_file(gri_max, subdir_path,
                                  ra_index, dec_index, max_stars=None):
        """Read one file (square degree) from one catalog subdirectory."""
        filename = '{:03d}'.format(ra_index) + '{:+03d}'.format(dec_index) + '.rc2'
        fullpath = os.path.join(subdir_path, filename)
        df = pd.read_csv(fullpath, sep=',', engine='python', header=None,
                         skip_blank_lines=True, error_bad_lines=False, nrows=max_stars,
                         # Column numbers have origin zero:
                         usecols=[0, 1, 4, 5, 6, 7,
                                  8, 9, 10, 11, 12, 13,
                                  14, 16, 18, 19, 20,
                                  21, 22, 25, 26, 29, 30, 33, 34], prefix='col')
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
        # Radii are in arcseconds:
        df['RP1'] = [None if rp1 == 999 else rp1 / 10.0 for rp1 in df['RP1']]
        df['R1'] = [None if r1 == 999 else r1 / 10.0 for r1 in df['R1']]
        df['R10'] = [None if r10 == 999 else r10 / 10.0 for r10 in df['R10']]
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
        """Remove rows of dataframe that lie outside RA range or Dec range.
        Respects RA zero-crossings."""
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
        """Remove rows (stars) which have other nearby stars (i.e., which will
        interfere with photometry).
        Current implementation uses only criterion based on 'RP1' = 'radius point-one',
        or the minimum radius inside which other light sources comprise 0.1% of
        the target's own flux. Rows (stars) with missing/null rp1 entries in catalog
        are assumed to have no nearby flux."""
        rp1_too_close = pd.Series([False if pd.isnull(rp1)
                                   else ((overlap_distance is None) or
                                         (rp1 < overlap_distance))
                                   for rp1 in df['RP1']])
        return df.loc[list(~rp1_too_close), :].copy()

    @staticmethod
    def _add_new_columns(df):
        """Add new columns 'BminusV', 'APASS_R', and 'ri_color' with derived data."""
        df.loc[:, 'BminusV'] = [(0.830 * g - 0.803 * r)
                                for (g, r) in zip(df['g'], df['r'])]
        df.loc[:, 'APASS_R'] = [(0.950 * r + 0.05 * i)
                                for (r, i) in zip(df['r'], df['i'])]
        df.loc[:, 'ri_color'] = [(r - i)
                                 for (r, i) in zip(df['r'], df['i'])]
        return df.copy()

    @staticmethod
    def _update_epoch(df, target_epoch):
        """Update RA and Dec of each star for proper motion between catalog epoch
        and target epoch."""
        delta_years = (target_epoch - Time(ATLAS_REFCAT2_EPOCH_UTC)).jd /\
            DAYS_PER_YEAR_NOMINAL
        ra_date = [(ra_epoch + delta_years * pm_ra / 3600.0) % 360
                   for (ra_epoch, pm_ra) in zip(df['RA_deg'], df['PM_ra'])]
        dec_date = [dec_epoch + delta_years * pm_dec / 3600.0
                    for (dec_epoch, pm_dec) in zip(df['Dec_deg'], df['PM_dec'])]
        df.loc[:, 'RA_deg'] = ra_date
        df.loc[:, 'Dec_deg'] = dec_date
        return df

    @staticmethod
    def _sort_by(df, sort_by):
        """Sort stars (Dataframe rows) by one of: 'ra', 'dec', or 'r';
        or None to skip sorting."""
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
        raise ValueError('Parameter \'sort by\' must be '
                         '\'ra\', \'dec\', \'r\', or None.')

    def _select_on(self, column_name, minimum=None, maximum=None):
        """Select rows (stars) on minimum and/or maximum value in column_name."""
        at_least_minimum = len(self.df_selected) * [True]  # default value
        at_most_maximum = at_least_minimum.copy()          # default value
        if minimum is not None:
            at_least_minimum = [(x >= minimum) for x in self.df_selected[column_name]]
        if maximum is not None:
            at_most_maximum = [(x <= maximum) for x in self.df_selected[column_name]]
        rows_to_keep = [(mn and mx)
                        for (mn, mx) in zip(at_least_minimum, at_most_maximum)]
        self.df_selected = self.df_selected.loc[rows_to_keep, :]

    def select_on_g_mag(self, min_g_mag=None, max_g_mag=None):
        """Select rows (stars) on minimum and/or maximum Sloan g' magnitude.

        Parameters
        ----------
        min_g_mag : float, or None
            Minimum Sloan g' magnitude a star (Dataframe row) may have to be retained.
            If None, criterion is not applied.
        max_g_mag : float, or None
            Maximum Sloan g' magnitude a star (Dataframe row) may have to be retained.
            If None, criterion is not applied.

        Returns
        -------
        None.
            Modifies ``df_selected`` in place.
        """
        self._select_on(column_name='g', minimum=min_g_mag, maximum=max_g_mag)

    def select_on_r_mag(self, min_r_mag=None, max_r_mag=None):
        """Select rows (stars) on minimum and/or maximum Sloan r' magnitude.

        Parameters
        ----------
        min_r_mag : float, or None
            Minimum Sloan r' magnitude a star (Dataframe row) may have to be retained.
            If None, criterion is not applied.
        max_r_mag : float, or None
            Maximum Sloan r' magnitude a star (Dataframe row) may have to be retained.
            If None, criterion is not applied.

        Returns
        -------
        None.
            Modifies ``df_selected`` in place.
        """
        self._select_on(column_name='r', minimum=min_r_mag, maximum=max_r_mag)

    def select_on_i_mag(self, min_i_mag=None, max_i_mag=None):
        """Select rows (stars) on minimum and/or maximum Sloan i' magnitude.

        Parameters
        ----------
        min_i_mag : float, or None
            Minimum Sloan i' magnitude a star (Dataframe row) may have to be retained.
            If None, criterion is not applied.
        max_i_mag : float, or None
            Maximum Sloan i' magnitude a star (Dataframe row) may have to be retained.
            If None, criterion is not applied.

        Returns
        -------
        None.
            Modifies ``df_selected`` in place.
        """
        self._select_on(column_name='i', minimum=min_i_mag, maximum=max_i_mag)

    def select_on_g_uncert(self, min_g_uncert=None, max_g_uncert=None):
        """Select rows (stars) on minimum and/or maximum catalog uncertainty in
        Sloan g' magnitude.

        Parameters
        ----------
        min_g_uncert : float, or None
            Minimum Sloan g' magnitude catalog uncertainty a star (Dataframe row)
            may have to be retained. In millimagnitudes.
            If None, criterion is not applied.
        max_g_uncert : float, or None
            Maximum Sloan g' magnitude catalog uncertainty a star (Dataframe row)
            may have to be retained. In millimagnitudes.
            If None, criterion is not applied.

        Returns
        -------
        None.
            Modifies ``df_selected`` in place.
        """
        self._select_on(column_name='dg', minimum=min_g_uncert, maximum=max_g_uncert)

    def select_on_r_uncert(self, min_r_uncert=None, max_r_uncert=None):
        """Select rows (stars) on minimum and/or maximum catalog uncertainty in
        Sloan r' magnitude.

        Parameters
        ----------
        min_r_uncert : float, or None
            Minimum Sloan r' magnitude catalog uncertainty a star (Dataframe row)
            may have to be retained. In millimagnitudes.
            If None, criterion is not applied.
        max_r_uncert : float, or None
            Maximum Sloan r' magnitude catalog uncertainty a star (Dataframe row)
            may have to be retained. In millimagnitudes.
            If None, criterion is not applied.

        Returns
        -------
        None.
            Modifies ``df_selected`` in place.
        """
        self._select_on(column_name='dr', minimum=min_r_uncert, maximum=max_r_uncert)

    def select_on_i_uncert(self, min_i_uncert=None, max_i_uncert=None):
        """Select rows (stars) on minimum and/or maximum catalog uncertainty in
        Sloan i' magnitude.

        Parameters
        ----------
        min_i_uncert : float, or None
            Minimum Sloan i' magnitude catalog uncertainty a star (Dataframe row)
            may have to be retained. In millimagnitudes.
            If None, criterion is not applied.
        max_i_uncert : float, or None
            Maximum Sloan i' magnitude catalog uncertainty a star (Dataframe row)
            may have to be retained. In millimagnitudes.
            If None, criterion is not applied.

        Returns
        -------
        None.
            Modifies ``df_selected`` in place.
        """
        self._select_on(column_name='di', minimum=min_i_uncert, maximum=max_i_uncert)

    def select_on_bv_color(self, min_bv_color=None, max_bv_color=None):
        """Select rows (stars) on minimum and/or maximum Johnson B-V color index.

        Parameters
        ----------
        min_bv_color : float, or None
            Minimum catalog Johnson B-V color index a star (Dataframe row) may
            have to be retained. In magnitudes. If None, criterion is not applied.
        max_bv_color : float, or None
            Maximum catalog Johnson B-V color index a star (Dataframe row) may
            have to be retained. In magnitudes. If None, criterion is not applied.

        Returns
        -------
        None.
            Modifies ``df_selected`` in place.
        """
        self._select_on(column_name='BminusV',
                        minimum=min_bv_color, maximum=max_bv_color)

    def select_on_ri_color(self, min_ri_color=None, max_ri_color=None):
        """Select rows (stars) on minimum and/or maximum Sloan r'-i' color index.

        Parameters
        ----------
        min_ri_color : float, or None
            Minimum catalog Sloan r'-i' color index a star (Dataframe row) may
            have to be retained. In magnitudes. If None, criterion is not applied.
        max_ri_color : float, or None
            Maximum catalog Sloan r'-i' color index a star (Dataframe row) may
            have to be retained. In magnitudes. If None, criterion is not applied.

        Returns
        -------
        None.
            Modifies ``df_selected`` in place.
        """
        self._select_on(column_name='ri_color',
                        minimum=min_ri_color, maximum=max_ri_color)


__________GENERAL_CATALOG_FUNCTIONS_____________________________________________ = 0


def _make_ra_index_list(ra_deg_range, include_max_bound=False):
    """Accept tuple (RA_min, RA_max), return all integer RA values covering the range.
       Handles RA zero-crossing gracefully. Assumes RA range less than 180 deg."""
    # Generate list of integer RA values needed to cover the entire RA range:
    ra_range = Angle(ra_deg_range * u.deg)
    ra_center = circmean(ra_range)
    wrap_angle = Angle(180.0 * u.deg + ra_center)
    ra_deg_wrapped = ra_range.wrap_at(wrap_angle).degree
    ra_index_first = floor(ra_deg_wrapped[0])
    ra_index_last = ceil(ra_deg_wrapped[1]) if include_max_bound \
        else floor(ra_deg_wrapped[1])
    ra_index_list = [ra % 360 for ra in range(ra_index_first, ra_index_last + 1)]
    return ra_index_list


def _make_dec_index_list(dec_deg_range, include_max_bound=False):
    """Accept tuple (Dec_min, Dec_max), return all integer Declination values
    covering the range."""
    dec_deg_first = floor(dec_deg_range[0])
    dec_deg_last = ceil(dec_deg_range[1]) if include_max_bound \
        else floor(dec_deg_range[1])
    dec_index_list = [dec for dec in range(dec_deg_first, dec_deg_last + 1)]
    return dec_index_list
