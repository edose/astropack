"""test_geometry.py"""

__author__ = "Eric Dose, Albuquerque"

from math import pi, sqrt

import pytest
# import numpy as np

# Test target:
import astropack.geometry as geom


class Test_Class_XY:
    @staticmethod
    @pytest.fixture(scope='function')
    def setup_xy():
        """ Helper function for testing XY. """
        xy1a = geom.XY(4, 5)
        xy1b = geom.XY(4, 5)
        xy2 = geom.XY(7, 9)
        dxy = geom.DXY(3, 4)
        return xy1a, xy1b, xy2, dxy

    def test_xy_constructor(self, setup_xy):
        xy1a, xy1b, xy2, dxy = setup_xy
        assert xy1a == geom.XY(4, 5)
        # assert geom.XY(xy1a) == xy1a

    def test_from_tuple(self):
        tupa = (44, 56)
        xya = geom.XY.from_tuple(tupa)
        assert isinstance(xya, geom.XY)
        assert xya.x == tupa[0] and xya.y == tupa[1]

    def test_xy_eq_ne_is(self, setup_xy):
        xy1a, xy1b, xy2, dxy = setup_xy
        assert xy1a == xy1b
        assert not xy1a != xy1b
        assert xy1a is not xy1b

    def test_xy_add(self, setup_xy):
        xy1a, xy1b, xy2, dxy = setup_xy
        assert xy1a + dxy == xy2
        assert isinstance(xy1a + dxy, geom.XY)
        with pytest.raises(TypeError):
            _ = xy1a + xy2
            _ = xy1a + (5, 6)
            _ = xy1a + 5

    def test_xy_subtract(self, setup_xy):
        xy1a, xy1b, xy2, dxy = setup_xy
        assert xy2 - xy1a == dxy
        assert xy2 - dxy == xy1a
        assert isinstance(xy2 - xy1a, geom.DXY)
        assert isinstance(xy2 - dxy, geom.XY)
        with pytest.raises(TypeError):
            _ = xy2 - (5, 6)
            _ = xy2 - 5

    def test_vector_to(self, setup_xy):
        xy1a, xy1b, xy2, dxy = setup_xy
        assert xy1a.vector_to(xy2) == dxy
        assert xy1a.vector_to(xy1b) == geom.DXY(0, 0)
        assert isinstance(xy1a.vector_to(xy2), geom.DXY)
        assert isinstance(xy1a.vector_to(xy1b), geom.DXY)
        with pytest.raises(TypeError):
            _ = xy1a.vector_to(dxy)
            _ = xy1a.vector_to((5, 6))

    def test_xy_operators_not_implemented(self, setup_xy):
        xy1a, xy1b, xy2, dxy = setup_xy
        with pytest.raises(Exception):
            _ = xy1a.__and__(xy2)
            _ = xy1a is True
            _ = xy1a.__or__(xy2)
            _ = xy1a / xy2
            _ = xy1a / dxy


class Test_Class_DXY:
    @staticmethod
    @pytest.fixture(scope='function')
    def setup_dxy():
        xy1 = geom.XY(4, 5)
        xy2 = geom.XY(7, 9)
        dxy1 = geom.DXY(3, 4)
        dxy2a = geom.DXY(55, 66)
        dxy2b = geom.DXY(55, 66)
        return xy1, xy2, dxy1, dxy2a, dxy2b

    def test_dxy_constructor(self, setup_dxy):
        xy1, xy2, dxy1, dxy2a, dxy2b = setup_dxy
        assert (dxy1.dx, dxy1.dy) == (3, 4)

    def test_from_tuple(self):
        tupa = (44, 56)
        dxya = geom.DXY.from_tuple(tupa)
        assert isinstance(dxya, geom.DXY)
        assert dxya.dx == tupa[0] and dxya.dy == tupa[1]

    def test_eq_ne_is(self, setup_dxy):
        xy1, xy2, dxy1, dxy2a, dxy2b = setup_dxy
        assert dxy2a == dxy2b
        assert dxy2a != dxy1
        assert dxy2a is not dxy2b

    def test_dxy_add(self, setup_dxy):
        xy1, xy2, dxy1, dxy2a, dxy2b = setup_dxy
        assert dxy1 + dxy2a == geom.DXY(58, 70)
        assert isinstance(dxy1 + dxy2a, geom.DXY)
        assert dxy1 + xy1 == xy2
        with pytest.raises(TypeError):
            _ = dxy2a + (5, 6)
            _ = dxy2a + 5

    def test_dxy_subtract(self, setup_dxy):
        xy1, xy2, dxy1, dxy2a, dxy2b = setup_dxy
        assert dxy2a - dxy1 == geom.DXY(52, 62)
        assert isinstance(dxy2a - dxy1, geom.DXY)
        with pytest.raises(TypeError):
            _ = dxy1 - xy1
            _ = dxy2a - (5, 6)
            _ = dxy2a - 5

    def test_dxy_multiply_fns(self, setup_dxy):
        xy1, xy2, dxy1, dxy2a, dxy2b = setup_dxy
        assert dxy1 * 3.5 == 3.5 * dxy1 == geom.DXY(10.5, 14)
        assert isinstance(dxy1 * 3.5, geom.DXY)
        assert isinstance(3.5 * dxy1, geom.DXY)
        with pytest.raises(TypeError):
            _ = dxy1 * xy1
            _ = dxy1 * (5, 6)

    def test_dxy_division(self, setup_dxy):
        xy1, xy2, dxy1, dxy2a, dxy2b = setup_dxy
        assert dxy1 / 2 == geom.DXY(1.5, 2)
        assert isinstance(dxy1 / 2, geom.DXY)
        with pytest.raises(ZeroDivisionError):
            _ = dxy1 / 0.0
        with pytest.raises(TypeError):
            _ = dxy1 / dxy2a
            _ = dxy1 / xy1
            _ = dxy1 / (5, 6)

    def test_dxy_angle_with(self, setup_dxy):
        xy1, _, dxy1, dxy2a, dxy2b = setup_dxy
        assert dxy1.angle_with(dxy1) == 0.0
        assert geom.DXY(5, 2).angle_with(geom.DXY(-2, 5)) == \
               pytest.approx(pi / 2, abs=0.000001)
        assert geom.DXY(5, 2).angle_with(geom.DXY(2, -5)) == \
               pytest.approx((3 / 2) * pi, abs=0.000001)
        with pytest.raises(TypeError):
            _ = dxy1.angle_with((5, 5))
            _ = dxy1.angle_with(xy1)
            _ = dxy1.angle_with(5)

    def test_dxy_dot(self, setup_dxy):
        xy1, xy2, dxy1, dxy2a, dxy2b = setup_dxy
        assert dxy1.dot(dxy2a) == 429
        assert dxy1.dot(dxy1) == 25
        with pytest.raises(TypeError):
            _ = dxy1.dot((5, 5))
            _ = dxy1.dot(xy1)
            _ = dxy1.dot(5)

    def test_dxy_properties(self, setup_dxy):
        _, _, dxy1, dxy2a, dxy2b = setup_dxy
        assert dxy1.length2 == 25
        assert dxy2a.length2 == 7381
        assert dxy1.length == 5
        assert dxy2a.length == pytest.approx(sqrt(7381), abs=0.000001)
        assert geom.DXY(5, 0).direction == 0
        assert geom.DXY(0, 5).direction == pytest.approx(pi / 2, abs=0.000001)
        assert geom.DXY(5, 5).direction == pytest.approx(pi / 4, abs=0.000001)


class Test_Class_Rectangle_in_2D:
    @staticmethod
    @pytest.fixture(scope='function')
    def setup_rectangle_in_2d():
        rect_pt_a = geom.XY(5, 0)
        rect_pt_b = geom.XY(0, 2)
        rect_pt_c = geom.XY(1, 4.5)
        return geom.Rectangle_in_2D(rect_pt_a, rect_pt_b, rect_pt_c)

    def test_constructor(self, setup_rectangle_in_2d):
        rect = setup_rectangle_in_2d
        assert (rect.a, rect.b, rect.c) == \
               (geom.XY(5, 0), geom.XY(0, 2), geom.XY(1, 4.5))
        assert rect.ab == geom.DXY(-5, 2)
        assert rect.bc == geom.DXY(1, 2.5)
        bad_pt_c = geom.XY(1, 4.49)
        with pytest.raises(ValueError):
            _ = geom.Rectangle_in_2D(rect.a, rect.b, bad_pt_c)

    def test_properties(self, setup_rectangle_in_2d):
        rect = setup_rectangle_in_2d
        assert rect.area == pytest.approx(sqrt(29) * sqrt(7.25), abs=0.000001)

    def test_contains_point(self, setup_rectangle_in_2d):
        rect = setup_rectangle_in_2d
        assert rect.contains_point(geom.XY(5.49, 1.25)) is True
        assert rect.contains_point(geom.XY(5.51, 1.25)) is False
        assert rect.contains_point(geom.XY(3.5, 3.49)) is True
        assert rect.contains_point(geom.XY(3.5, 3.51)) is False
        assert rect.contains_point(geom.XY(1, 4.49)) is True
        assert rect.contains_point(geom.XY(1, 4.51)) is False
        assert rect.contains_point(geom.XY(3, 2.25)) is True
        assert rect.contains_point(geom.XY(100, 100)) is False
        pts_on_edge = [rect.a, rect.b, rect.c,
                       geom.XY(2.5, 1),
                       geom.XY(0.5, 3.25),
                       geom.XY(3.5, 3.5),
                       geom.XY(5.5, 1.25)]
        assert all([rect.contains_point(pt, include_edges=True)
                    for pt in pts_on_edge])
        assert not any([rect.contains_point(pt, include_edges=False)
                        for pt in pts_on_edge])

    def test_contains_points(self, setup_rectangle_in_2d):
        rect = setup_rectangle_in_2d
        # Just inside and outside of rectangle:
        xy_array = (geom.XY(5.49, 1.25),
                    geom.XY(5.51, 1.25),
                    geom.XY(3.5, 3.49),
                    geom.XY(3.5, 3.51),
                    geom.XY(1, 4.49),
                    geom.XY(1, 4.51),
                    geom.XY(3, 2.25),
                    geom.XY(100, 100))
        contains = rect.contains_points(xy_array, include_edges=True)
        assert contains == 4 * [True, False]
        # Edges and vertices:
        xy_array = (geom.XY(0, 2),
                    geom.XY(1, 4.5),
                    geom.XY(5.5, 1.25))
        contains = rect.contains_points(xy_array, include_edges=True)
        assert contains == 3 * [True]
        contains = rect.contains_points(xy_array, include_edges=False)
        assert contains == 3 * [False]

    def test_contains_points_unitgrid(self):
        rect_pt_a = geom.XY(10, -2)  # all (x,y)
        rect_pt_b = geom.XY(-2, 1)
        rect_pt_c = geom.XY(0.5, 11)  # fourth pt is thus (12.5, 10)
        rect = geom.Rectangle_in_2D(rect_pt_a, rect_pt_b, rect_pt_c)
        with pytest.raises(TypeError):
            _ = rect.contains_points_unitgrid(-3, 15.5, -3, 12)
        with pytest.raises(ValueError):
            _ = rect.contains_points_unitgrid(11, 8, -3, 12)
            _ = rect.contains_points_unitgrid(-3, 15, 12, 3)
        result_grid = rect.contains_points_unitgrid(-3, 15, -3, 12, include_edges=True)
        assert result_grid.shape == (19, 16)  # (n x grid pts, n y grid pts), not numpy
        assert result_grid.dtype == bool
        assert result_grid.sum() == 132
        assert result_grid[0, :].sum() == 0
        assert result_grid[:, 0].sum() == 0

        assert result_grid[14, :].sum() == 7
        assert result_grid[10, :].sum() == 11
        assert result_grid[:, 10].sum() == 13

        assert result_grid[4, :].sum() == 10
        assert result_grid[:, 3].sum() == 9
        result_grid2 = rect.contains_points_unitgrid(-3, 15, -3, 12,
                                                     include_edges=False)
        assert result_grid2.shape == (19, 16)
        assert result_grid2.dtype == bool
        assert result_grid2.sum() == 124


class Test_Class_Circle_in_2D:
    @staticmethod
    @pytest.fixture(scope='function')
    def setup_circle_in_2d():
        return geom.Circle_in_2D(xy_origin=geom.XY(4, 6), radius=5)

    def test_constructor(self, setup_circle_in_2d):
        circle = setup_circle_in_2d
        assert circle.origin == geom.XY(4, 6)
        assert circle.radius == 5.0
        assert circle.x == 4.0
        assert circle.y == 6.0
        with pytest.raises(TypeError):
            _ = geom.Circle_in_2D(55, 5)
            _ = geom.Circle_in_2D((6, 7), 5)

    def test_properties(self, setup_circle_in_2d):
        circle = setup_circle_in_2d
        assert circle.area == pytest.approx(pi * circle.radius ** 2, abs=1e-10)

    def test_contains_point(self, setup_circle_in_2d):
        circle = setup_circle_in_2d
        assert circle.contains_point(geom.XY(0, 6), include_edges=True) is True
        assert circle.contains_point(geom.XY(-1, 6), include_edges=True) is True
        assert circle.contains_point(geom.XY(8.9, 6), include_edges=True) is True
        assert circle.contains_point(geom.XY(9.1, 6), include_edges=True) is False
        assert circle.contains_point(geom.XY(100, 100), include_edges=True) is False
        assert circle.contains_point(geom.XY(0, 6), include_edges=False) is True
        assert circle.contains_point(geom.XY(-1, 6), include_edges=False) is False
        assert circle.contains_point(geom.XY(8.9, 6), include_edges=False) is True
        assert circle.contains_point(geom.XY(9.1, 6), include_edges=False) is False
        assert circle.contains_point(geom.XY(100, 100), include_edges=False) is False

    def test_contains_points(self, setup_circle_in_2d):
        circle = setup_circle_in_2d
        points = (geom.XY(0, 6),
                  geom.XY(-1, 6),
                  geom.XY(8.9, 6),
                  geom.XY(9.1, 6),
                  geom.XY(100, 100))
        assert circle.contains_points(points, include_edges=True) == \
               [True, True, True, False, False]
        assert circle.contains_points(points, include_edges=False) == \
               [True, False, True, False, False]

    def test_contains_points_unitgrid(self):
        circle_origin = geom.XY(6, 5)
        circle_radius = 7
        circle = geom.Circle_in_2D(circle_origin, circle_radius)
        with pytest.raises(TypeError):
            _ = circle.contains_points_unitgrid(-3, 15.5, -3, 12)
        with pytest.raises(ValueError):
            _ = circle.contains_points_unitgrid(11, 8, -3, 12)
            _ = circle.contains_points_unitgrid(-3, 15, 12, 3)
        result_grid = circle.contains_points_unitgrid(-3, 15, -3, 12,
                                                      include_edges=True)
        assert result_grid.shape == (19, 16)  # (x,y), not numpy's convention
        assert result_grid.dtype == bool
        assert result_grid.sum() == 149
        assert result_grid[0, :].sum() == 0
        assert result_grid[15, :].sum() == 7
        assert result_grid[:, :2].sum() == 1
        assert result_grid[:, 17:].sum() == 0
        assert result_grid[4, :].sum() == 9
        assert result_grid[:, 9].sum() == 13
        result_grid2 = circle.contains_points_unitgrid(-3, 15, -3, 12,
                                                       include_edges=False)
        assert result_grid2.shape == (19, 16)
        assert result_grid2.dtype == bool
        assert result_grid2.sum() == 145


__________FUNCTIONS_____________________________________________________________ = 0


def test_distance_to_line():
    pta = geom.XY(1, 6)
    ptb = geom.XY(5, 9)
    # Case: line points are coincident:
    assert geom.distance_to_line((5, 6), pta, pta) == 4.0
    assert geom.distance_to_line((-1, -1), pta, pta) == pytest.approx(sqrt(53.0))
    # Normal case:
    assert geom.distance_to_line((5, 6), pta, ptb) == pytest.approx(2.4)
    assert geom.distance_to_line((21, 18), pta, ptb) == pytest.approx(2.4)
    assert geom.distance_to_line(((pta.x + ptb.x) / 2.0,
                                  (pta.y + ptb.y) / 2.0), pta, ptb) == 0.0
    # Case: dist_12 is passed in:
    assert geom.distance_to_line((10, 0), (0, 3), (0, 7), dist_12=4) == 10.0
    # Exception cases:
    with pytest.raises(TypeError):
        _ = geom.distance_to_line((2, 3), 5, (2, 3))
        _ = geom.distance_to_line([2, 3], (4, 5), (2, 3))


def test_make_golden_spiral():
    n = 225
    points = geom.make_golden_spiral(n)
    assert isinstance(points, list)
    assert len(points) == n
    assert all([isinstance(p, tuple) for p in points])
    assert points[0].az == pytest.approx(291.25, abs=0.01)
    assert points[0].alt == pytest.approx(84.60, abs=0.01)
    assert points[-1].az == pytest.approx(89.51, abs=0.01)
    assert points[-1].alt == pytest.approx(-84.60, abs=0.01)
