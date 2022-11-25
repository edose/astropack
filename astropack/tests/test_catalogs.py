"""test_catalogs.py"""

__author__ = "Eric Dose, Albuquerque"

# Python core packages:
import os
from random import shuffle

# External packages:
import pytest
import pandas as pd
from astropy.time import Time

# Test target:
from astropack import catalogs

THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname
                                              (os.path.abspath(__file__)))
TEST_TOP_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, "test")

ATLAS_CATALOG_TOP_DIRECTORY = 'D:/Astro/Catalogs/ATLAS-refcat2'
INDEX_FILENAME = 'index.txt'
TARGET_EPOCH = Time('2022-02-28 00:00:00', scale='utc')


__________TEST_MODULE_FUNCTIONS________________________________________________ = 0


def test__make_ra_index_list():
    fn = catalogs._make_ra_index_list
    assert fn((11.5, 13.9)) == [11, 12, 13]
    assert fn((11.5, 13.9), include_max_bound=True) == [11, 12, 13, 14]
    assert fn((12, 14)) == [12, 13, 14]
    assert fn((358.75, 1.44)) == [358, 359, 0, 1]
    assert fn((358.75, 1.44), include_max_bound=True) == [358, 359, 0, 1, 2]


def test__make_dec_index_list():
    fn = catalogs._make_dec_index_list
    assert fn((11.5, 13.9)) == [11, 12, 13]
    assert fn((11.5, 13.9), include_max_bound=True) == [11, 12, 13, 14]
    assert fn((-14, -12)) == [-14, -13, -12]
    assert fn((-0.1, 0.7)) == [-1, 0]
    assert fn((-0.1, 0.7), include_max_bound=True) == [-1, 0, 1]


__________TEST_CLASS_ATLASREFCAT2_____________________________________________ = 0
_____constructor_support_methods_____ = 0


def test_class_atlasrefcat2__validate_input_values():
    fn = catalogs.AtlasRefcat2._validate_input_values
    top_dir = ATLAS_CATALOG_TOP_DIRECTORY
    index = INDEX_FILENAME
    epoch = Time('2022-11-22 03:04:05')

    # Normal case passes:
    fn(atlas_top_directory=top_dir, index_filename=index,
       ra_deg_range=(34, 35), dec_deg_range=(64, 66), target_epoch=epoch,
       overlap_distance=10, sort_by='ra')  # passes.

    # Test Exception cases in order, so that first failing condition only operates.

    # NotADirectoryError:
    with pytest.raises(NotADirectoryError):
        fn((top_dir + 'XXX'), index, (34, 35), (64, 66), epoch, 10, 'ra')

    with pytest.raises(catalogs.MissingIndexFileError):
        fn(top_dir, index + 'XXX', (34, 35), (64, 66), epoch, 10, 'ra')

    with pytest.raises(catalogs.InvalidTargetEpochError):
        fn(top_dir, index, (34, 35), (64, 66), 3.14159, 10, 'ra')

    with pytest.raises(catalogs.Invalid_RA_RangeError):
        fn(top_dir, index, (34, 435), (64, 66), epoch, 10, 'ra')

    with pytest.raises(catalogs.Invalid_Dec_RangeError):
        fn(top_dir, index, (34, 35), (64, 96), epoch, 10, 'ra')

    with pytest.raises(catalogs.InvalidOverlapDistanceError):
        fn(top_dir, index, (34, 35), (64, 66), epoch, -5, 'ra')
    with pytest.raises(catalogs.InvalidOverlapDistanceError):
        fn(top_dir, index, (34, 35), (64, 66), epoch, 65, 'ra')

    with pytest.raises(catalogs.InvalidSortTypeError):
        fn(top_dir, index, (34, 35), (64, 66), epoch, 10, 'raXXX')


def test_class_atlasrefcat2__locate_subdirs():
    df_subcat = catalogs.AtlasRefcat2._locate_subdirs(ATLAS_CATALOG_TOP_DIRECTORY,
                                                      'index.txt')
    assert isinstance(df_subcat, pd.DataFrame)
    assert list(df_subcat['gri_max']) == [16, 17, 18, 19, 20]
    assert list(df_subcat['Path'])[0] == os.path.join(ATLAS_CATALOG_TOP_DIRECTORY,
                                                      'mag-0-16')
    assert list(df_subcat['Path'])[4] == os.path.join(ATLAS_CATALOG_TOP_DIRECTORY,
                                                      'mag-19-20')
    assert list(df_subcat['Nfiles']) == [64800, 64799, 64800, 64800, 64800]


def test_class_atlasrefcat2__read_one_subdir_one_file():
    fn = catalogs.AtlasRefcat2._read_one_subdir_one_file
    gri_max = 16
    subdir_path = os.path.join(ATLAS_CATALOG_TOP_DIRECTORY, 'mag-0-16')
    ra_index, dec_index = 240, 13

    # Read all stars:
    df = fn(gri_max, subdir_path, ra_index, dec_index)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 716
    assert df.loc[2, 'RA_deg'] == pytest.approx(240.00214, abs=0.00002)

    # Read limited number of stars:
    df = fn(gri_max, subdir_path, ra_index, dec_index, max_stars=22)
    assert len(df) == 22
    assert df.loc[2, 'RA_deg'] == pytest.approx(240.00214, abs=0.00002)


def test_class_atlasrefcat2__get_stars_by_index():
    df_subcat = catalogs.AtlasRefcat2._locate_subdirs(ATLAS_CATALOG_TOP_DIRECTORY,
                                                      'index.txt')
    # Case: normal, away from RA=0:
    df = catalogs.AtlasRefcat2._get_stars_by_index(df_subcat,
                                                   ra_deg_range=(23.75, 24.125),
                                                   dec_deg_range=(-23.125, -22.625))
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 5691
    assert all((ra >= 23) and (ra <= 25) for ra in df['RA_deg'])
    assert all((dec >= -24) and (dec <= -22) for dec in df['Dec_deg'])
    assert not all((ra >= 23.75) and (ra <= 24.125) for ra in df['RA_deg'])
    assert not all((dec >= -23.125) and (dec <= -22.625) for dec in df['Dec_deg'])

    # Case: normal, crossing RA=0:
    df = catalogs.AtlasRefcat2._get_stars_by_index(df_subcat,
                                                   ra_deg_range=(359.5, 0.125),
                                                   dec_deg_range=(-0.375, 0.25))
    assert len(df) == 8135
    assert all((ra >= 359) or (ra <= 1) for ra in df['RA_deg'])
    assert all((dec >= -1) and (dec <= 1) for dec in df['Dec_deg'])
    assert not all((ra >= 359.5) or (ra <= 0.125) for ra in df['RA_deg'])
    assert not all((dec >= -0.375) and (dec <= 0.25) for dec in df['Dec_deg'])

    # Case: entirely within one square degree:
    df = catalogs.AtlasRefcat2._get_stars_by_index(df_subcat,
                                                   ra_deg_range=(120.125, 120.75),
                                                   dec_deg_range=(0.25, 0.5))
    assert len(df) == 11026
    assert all((ra >= 120) and (ra <= 121) for ra in df['RA_deg'])
    assert all((dec >= 0) and (dec <= 1) for dec in df['Dec_deg'])
    assert not all((ra >= 120.125) and (ra <= 120.75) for ra in df['RA_deg'])
    assert not all((dec >= 0.25) and (dec <= 0.5) for dec in df['Dec_deg'])


def test_class_atlasrefcat2__trim_to_ra_dec_range():
    df_subcat = catalogs.AtlasRefcat2._locate_subdirs(ATLAS_CATALOG_TOP_DIRECTORY,
                                                      INDEX_FILENAME)
    # Case: normal, away from RA=0:
    df_all = catalogs.AtlasRefcat2._get_stars_by_index(df_subcat,
                                                       ra_deg_range=(23.75, 24.125),
                                                       dec_deg_range=(-23.125, -22.625))
    assert not all(ra >= 23.75 for ra in df_all['RA_deg'])
    assert not all(ra <= 24.125 for ra in df_all['RA_deg'])
    assert not all(dec >= -23.125 for dec in df_all['Dec_deg'])
    assert not all(dec <= -22.625 for dec in df_all['Dec_deg'])

    df = catalogs.AtlasRefcat2._trim_to_ra_dec_range(df_all,
                                                     ra_deg_range=(23.75, 24.125),
                                                     dec_deg_range=(-23.125, -22.625))
    assert len(df) < len(df_all)
    assert len(df) == 268
    assert all((ra >= 23.75) and (ra <= 24.125) for ra in df['RA_deg'])
    assert all((dec >= -23.125) and (dec <= -22.625) for dec in df['Dec_deg'])

    # Case: normal, crossing RA=0:
    df_all = catalogs.AtlasRefcat2._get_stars_by_index(df_subcat,
                                                       ra_deg_range=(359.5, 0.125),
                                                       dec_deg_range=(-0.375, 0.25))
    df = catalogs.AtlasRefcat2._trim_to_ra_dec_range(df_all,
                                                     ra_deg_range=(359.5, 0.125),
                                                     dec_deg_range=(-0.375, 0.25))
    assert len(df) == 815
    assert all((ra >= 359.5) or (ra <= 0.125) for ra in df['RA_deg'])
    assert all((dec >= -0.375) and (dec <= 0.25) for dec in df['Dec_deg'])

    # Case: entirely within one square degree:
    df_all = catalogs.AtlasRefcat2._get_stars_by_index(df_subcat,
                                                       ra_deg_range=(120.125, 120.75),
                                                       dec_deg_range=(0.25, 0.5))
    df = catalogs.AtlasRefcat2._trim_to_ra_dec_range(df_all,
                                                     ra_deg_range=(120.125, 120.75),
                                                     dec_deg_range=(0.25, 0.5))
    assert len(df) == 1655
    assert all((ra >= 120.125) and (ra <= 120.75) for ra in df['RA_deg'])
    assert all((dec >= 0.25) and (dec <= 0.5) for dec in df['Dec_deg'])


def test_class_atlasrefcat2__remove_overlapping():
    df_subcat = catalogs.AtlasRefcat2._locate_subdirs(ATLAS_CATALOG_TOP_DIRECTORY,
                                                      'index.txt')
    df_all = catalogs.AtlasRefcat2._get_stars_by_index(df_subcat,
                                                       ra_deg_range=(23.75, 24.125),
                                                       dec_deg_range=(-23.125, -22.625))
    n_overlaps_all = sum(1 for rp1 in df_all['RP1'] if (rp1 is not None and rp1 < 10))
    assert len(df_all) == 5691
    assert n_overlaps_all == 356
    df = catalogs.AtlasRefcat2._remove_overlapping(df_all, overlap_distance=10)
    assert len(df) == len(df_all) - n_overlaps_all
    n_overlaps = sum(1 for rp1 in df['RP1'] if (rp1 is not None and rp1 < 10))
    assert n_overlaps == 0


def test_class_atlasrefcat2__add_new_columns():
    df = make_atlasrefcat2_test_object().df_all
    assert all(colname in df.columns
               for colname in ['BminusV', 'APASS_R', 'ri_color'])
    assert all(bv == pytest.approx(0.830 * g - 0.803 * r)
               for (bv, g, r) in zip(df['BminusV'], df['g'], df['r']))
    assert all(ar == pytest.approx(0.950 * r + 0.05 * i)
               for (ar, r, i) in zip(df['APASS_R'], df['r'], df['i']))
    assert all(ri == pytest.approx(r - i)
               for (ri, r, i) in zip(df['ri_color'], df['r'], df['i']))


def test_class_atlasrefcat2__update_epoch():
    df_subcat = catalogs.AtlasRefcat2._locate_subdirs(ATLAS_CATALOG_TOP_DIRECTORY,
                                                      'index.txt')
    df_raw = catalogs.AtlasRefcat2._get_stars_by_index(df_subcat,
                                                       ra_deg_range=(120.125, 120.75),
                                                       dec_deg_range=(0.25, 0.5))
    df_raw_copy = df_raw.copy()
    df_updated = catalogs.AtlasRefcat2._update_epoch(df_raw_copy, TARGET_EPOCH)
    delta_years = (TARGET_EPOCH.jd - Time(catalogs.ATLAS_REFCAT2_EPOCH_UTC).jd) / \
        catalogs.DAYS_PER_YEAR_NOMINAL
    assert all(ra_updated == pytest.approx(ra_raw + (ra_motion / 3600.0) *
                                           delta_years, abs=0.000001)
               for (ra_raw, ra_updated, ra_motion)
               in zip(df_raw['RA_deg'], df_updated['RA_deg'], df_raw['PM_ra']))
    assert all(dec_updated == pytest.approx(dec_raw + (dec_motion / 3600.0) *
                                            delta_years, abs=0.000001)
               for (dec_raw, dec_updated, dec_motion)
               in zip(df_raw['Dec_deg'], df_updated['Dec_deg'], df_raw['PM_dec']))


def test_class_atlasrefcat2__sort_by():
    df_subcat = catalogs.AtlasRefcat2._locate_subdirs(ATLAS_CATALOG_TOP_DIRECTORY,
                                                      'index.txt')
    df = catalogs.AtlasRefcat2._get_stars_by_index(df_subcat,
                                                   ra_deg_range=(120.125, 120.75),
                                                   dec_deg_range=(0.25, 0.5))

    # Test for 'ra':
    df_ra = df.copy()
    ra_list = list(df_ra['RA_deg'].copy())
    shuffle(ra_list)
    df_ra['RA_deg'] = ra_list
    assert not all(df_ra['RA_deg'].diff().iloc[1:] >= 0)  # i.e., not sorted.
    df_ra_sorted = catalogs.AtlasRefcat2._sort_by(df_ra, 'RA')
    assert all(df_ra_sorted['RA_deg'].diff().iloc[1:] >= 0)  # i.e., sorted in RA.
    assert all(df_ra_sorted.loc[i, 'CatalogID'] == df_ra.loc[i, 'CatalogID']
               for i in range(len(df_ra)))  # other columns ordered, too.
    del ra_list, df_ra, df_ra_sorted

    # Test for 'dec':
    df_dec = df.copy()
    dec_list = list(df_dec['Dec_deg'].copy())
    shuffle(dec_list)
    df_dec['Dec_deg'] = dec_list
    assert not all(df_dec['Dec_deg'].diff().iloc[1:] >= 0)
    df_dec_sorted = catalogs.AtlasRefcat2._sort_by(df_dec, 'dec')
    assert all(df_dec_sorted['Dec_deg'].diff().iloc[1:] >= 0)
    assert all(df_dec_sorted.loc[i, 'CatalogID'] == df_dec.loc[i, 'CatalogID']
               for i in range(len(df_dec)))  # other columns ordered, too.
    del dec_list, df_dec, df_dec_sorted

    # Test for 'r' (Sloan r magnitude):
    df_r = df.copy()
    r_list = list(df_r['r'].copy())
    shuffle(r_list)
    df_r['r'] = r_list
    assert not all(df_r['r'].diff().iloc[1:] >= 0)  # i.e., not sorted.
    df_r_sorted = catalogs.AtlasRefcat2._sort_by(df_r, 'r')
    assert all(df_r_sorted['r'].diff().iloc[1:] >= 0)  # is sorted.
    assert all(df_r_sorted.loc[i, 'CatalogID'] == df_r.loc[i, 'CatalogID']
               for i in range(len(df_r)))  # other columns ordered, too.
    del r_list, df_r, df_r_sorted


_____constructor_____ = 0


def test_class_atlasrefcat2_constructor():
    # Case: normal, multi-sq-degree range away from RA-zero:
    cat = catalogs.AtlasRefcat2(ATLAS_CATALOG_TOP_DIRECTORY, 'index.txt',
                                ra_deg_range=(120.125, 120.75),
                                dec_deg_range=(0.25, 0.5),
                                target_epoch=TARGET_EPOCH)
    assert len(cat.df_selected) == len(cat.df_all) == 1222
    assert all((ra >= 120.125) and (ra <= 120.75) for ra in cat.df_selected['RA_deg'])
    assert all((dec >= 0.25) and (dec <= 0.5) for dec in cat.df_selected['Dec_deg'])
    assert all(cat.df_selected['RA_deg'].diff().iloc[1:] >= 0)  # i.e., sorted in RA.

    # Case: range crossing RA=zero:
    cat = catalogs.AtlasRefcat2(ATLAS_CATALOG_TOP_DIRECTORY, 'index.txt',
                                ra_deg_range=(358.5, 0.25),
                                dec_deg_range=(-0.5, +0.125),
                                target_epoch=TARGET_EPOCH, sort_by='dec')
    assert len(cat.df_selected) == len(cat.df_all) == 2154
    assert all((ra >= 358.5) or (ra <= 0.25) for ra in cat.df_selected['RA_deg'])
    assert all((dec >= -0.5) and (dec <= 0.125) for dec in cat.df_selected['Dec_deg'])
    assert all(cat.df_selected['Dec_deg'].diff().iloc[1:] >= 0)  # i.e., sorted in Dec.

    # Case: range within one single sq degree:
    cat = catalogs.AtlasRefcat2(ATLAS_CATALOG_TOP_DIRECTORY, 'index.txt',
                                ra_deg_range=(240.2, 240.3), dec_deg_range=(13.1, 13.2),
                                target_epoch=TARGET_EPOCH, sort_by='r')
    assert len(cat.df_selected) == len(cat.df_all) == 32
    assert all((ra >= 240.2) and (ra <= 240.3) for ra in cat.df_selected['RA_deg'])
    assert all((dec >= 13.1) and (dec <= 13.2) for dec in cat.df_selected['Dec_deg'])
    assert all(cat.df_selected['r'].diff().iloc[1:] >= 0)  # i.e., sorted in r mag.


_____select_on_methods_____ = 0


def test_class_atlasrefcat2__select_on():
    # Test both minimum and maximum limits:
    cat = make_atlasrefcat2_test_object()
    assert len(cat.df_all) == len(cat.df_selected) == 1865
    assert not all([13 <= r <= 14 for r in cat.df_all['r']])
    cat._select_on('r', minimum=13, maximum=14)
    assert len(cat.df_all) == 1865
    assert len(cat.df_selected) == 53
    assert all([13 <= r <= 14 for r in cat.df_selected['r']])

    # Test minimum-only selection:
    cat = make_atlasrefcat2_test_object()
    cat._select_on('r', minimum=13)
    assert len(cat.df_selected) == 1795
    assert all([r >= 13 for r in cat.df_selected['r']])

    # Test maximum-only selection:
    cat = make_atlasrefcat2_test_object()
    cat._select_on('r', maximum=13)
    assert len(cat.df_selected) == 70
    assert all([r <= 13 for r in cat.df_selected['r']])


def test_class_atlasrefcat2_select_on_g_mag():
    cat = make_atlasrefcat2_test_object()
    assert len(cat.df_all) == len(cat.df_selected) == 1865
    cat.select_on_g_mag(14, 14.5)
    assert len(cat.df_selected) == 31
    assert not all([14 <= g <= 14.5 for g in cat.df_all['g']])
    assert all([14 <= g <= 14.5 for g in cat.df_selected['g']])


def test_class_atlasrefcat2_select_on_r_mag():
    cat = make_atlasrefcat2_test_object()
    assert len(cat.df_all) == len(cat.df_selected) == 1865
    cat.select_on_r_mag(14.25, 14.75)
    assert len(cat.df_selected) == 47
    assert not all([14.25 <= r <= 14.75 for r in cat.df_all['r']])
    assert all([14.25 <= r <= 14.75 for r in cat.df_selected['r']])


def test_class_atlasrefcat2_select_on_i_mag():
    cat = make_atlasrefcat2_test_object()
    assert len(cat.df_all) == len(cat.df_selected) == 1865
    cat.select_on_i_mag(14, 15.25)
    assert len(cat.df_selected) == 133
    assert not all([14 <= i <= 15.25 for i in cat.df_all['i']])
    assert all([14 <= i <= 15.25 for i in cat.df_selected['i']])


def test_class_atlasrefcat2_select_on_g_uncert():
    cat = make_atlasrefcat2_test_object()
    assert len(cat.df_all) == len(cat.df_selected) == 1865
    cat.select_on_g_uncert(max_g_uncert=15)
    assert len(cat.df_selected) == 1034
    assert not all([dg <= 15 for dg in cat.df_all['dg']])
    assert all([dg <= 15 for dg in cat.df_selected['dg']])


def test_class_atlasrefcat2_select_on_r_uncert():
    cat = make_atlasrefcat2_test_object()
    assert len(cat.df_all) == len(cat.df_selected) == 1865
    cat.select_on_r_uncert(max_r_uncert=14)
    assert len(cat.df_selected) == 1275
    assert not all([dr <= 14 for dr in cat.df_all['dr']])
    assert all([dr <= 14 for dr in cat.df_selected['dr']])


def test_class_atlasrefcat2_select_on_i_uncert():
    cat = make_atlasrefcat2_test_object()
    assert len(cat.df_all) == len(cat.df_selected) == 1865
    cat.select_on_i_uncert(max_i_uncert=12)
    assert len(cat.df_selected) == 1318
    assert not all([di <= 12 for di in cat.df_all['di']])
    assert all([di <= 12 for di in cat.df_selected['di']])


def test_class_atlasrefcat2_select_on_bv_color():
    cat = make_atlasrefcat2_test_object()
    assert len(cat.df_all) == len(cat.df_selected) == 1865
    cat.select_on_bv_color(min_bv_color=0.5, max_bv_color=0.6)
    assert len(cat.df_selected) == 19
    assert not all([0.5 <= bv <= 0.6 for bv in cat.df_all['BminusV']])
    assert all([0.5 <= bv <= 0.6 for bv in cat.df_selected['BminusV']])


def test_class_atlasrefcat2_select_on_ri_color():
    cat = make_atlasrefcat2_test_object()
    assert len(cat.df_all) == len(cat.df_selected) == 1865
    cat.select_on_ri_color(min_ri_color=0.25, max_ri_color=0.325)
    assert len(cat.df_selected) == 176
    assert not all([0.25 <= ri <= 0.325 for ri in cat.df_all['ri_color']])
    assert all([0.25 <= ri <= 0.325 for ri in cat.df_selected['ri_color']])


__________HELPER_FUNCTIONS_____________________________________________ = 0


def make_atlasrefcat2_test_object():
    """Get an ATLAS test object for use in testing. RA spans zero."""
    cat = catalogs.AtlasRefcat2(ATLAS_CATALOG_TOP_DIRECTORY, 'index.txt',
                                ra_deg_range=(358.5, 0.75),
                                dec_deg_range=(13.875, 14.25),
                                target_epoch=TARGET_EPOCH, sort_by='r')
    return cat
