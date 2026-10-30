import pytest

import geomcore
from geomcore import curves, surfaces

SUBCLASS_PARENTS = {
    "ConstructionError": "GeomcoreError",
    "FrameConstructionError": "ConstructionError",
    "TrimmedCurveConstructionError": "ConstructionError",
    "KnotInsertionError": "GeomcoreError",
    "KnotRemovalError": "GeomcoreError",
    "DegreeElevationError": "GeomcoreError",
    "BSplineSplitError": "GeomcoreError",
    "BSplineConversionError": "GeomcoreError",
    "Curve3DToUvError": "GeomcoreError",
    "UvCurveTo3DError": "GeomcoreError",
    "SurfaceIsoError": "GeomcoreError",
    "DiscretizationCriteriaError": "GeomcoreError",
}


def test_base_is_value_error():
    assert issubclass(geomcore.GeomcoreError, ValueError)


@pytest.mark.parametrize("name, parent", SUBCLASS_PARENTS.items())
def test_hierarchy(name, parent):
    cls = getattr(geomcore, name)
    assert issubclass(cls, getattr(geomcore, parent))
    assert issubclass(cls, geomcore.GeomcoreError)
    assert issubclass(cls, ValueError)
    assert cls.__doc__


@pytest.mark.parametrize(
    "catch", [geomcore.ConstructionError, geomcore.GeomcoreError, ValueError]
)
def test_negative_radius_circle(catch):
    with pytest.raises(catch):
        curves.Circle(geomcore.Point(0, 0, 0), geomcore.Vector(0, 0, 1), -1.0)


def test_negative_radius_sphere_is_construction_error():
    with pytest.raises(geomcore.ConstructionError):
        surfaces.Sphere(geomcore.Point(0, 0, 0), -1.0)


def test_discretization_criteria_error():
    circle = curves.Circle(geomcore.Point(0, 0, 0), geomcore.Vector(0, 0, 1), 1.0)
    with pytest.raises(geomcore.DiscretizationCriteriaError):
        circle.discretize(-1.0)


def test_frame_construction_error():
    with pytest.raises(geomcore.FrameConstructionError):
        geomcore.Frame(
            geomcore.Point(0, 0, 0),
            geomcore.Vector(0, 0, 0),
            geomcore.Vector(1, 0, 0),
        )


def test_trimmed_curve_construction_error():
    line = curves.Line.from_two_points(
        geomcore.Point(0, 0, 0), geomcore.Point(1, 0, 0)
    )
    with pytest.raises(geomcore.TrimmedCurveConstructionError):
        curves.TrimmedCurve(line, 1.0, 0.0)
    with pytest.raises(geomcore.ConstructionError):
        curves.TrimmedCurve(line, 1.0, 0.0)


def _straight_bspline():
    return curves.BSplineCurve.interpolate(
        [
            geomcore.Point(0, 0, 0),
            geomcore.Point(1, 1, 0),
            geomcore.Point(2, 0, 0),
        ],
        2,
    )


def test_knot_insertion_error():
    spline = _straight_bspline()
    with pytest.raises(geomcore.KnotInsertionError):
        spline.insert_knot(1e9, 1)


def test_bspline_split_error():
    spline = _straight_bspline()
    with pytest.raises(geomcore.BSplineSplitError):
        spline.split_at(1e9)


def _cubic_spline():
    return curves.BSplineCurve.interpolate(
        [
            geomcore.Point(0, 0, 0),
            geomcore.Point(1, 1, 0),
            geomcore.Point(2, 0, 0),
            geomcore.Point(3, 1, 0),
        ],
        3,
    )


def _bilinear_surface():
    poles = [
        [geomcore.Point(0, 0, 0), geomcore.Point(0, 1, 0)],
        [geomcore.Point(1, 0, 0), geomcore.Point(1, 1, 0)],
    ]
    return surfaces.BSplineSurface(
        1, 1, poles, [0.0, 1.0], [2, 2], [0.0, 1.0], [2, 2], False, False
    )


def test_knot_removal_error():
    spline = _cubic_spline()
    with pytest.raises(geomcore.KnotRemovalError):
        spline.remove_knot(0.123456, 1, 1e-6)
    with pytest.raises(geomcore.KnotRemovalError):
        spline.remove_knots([0.5], -1.0)


def test_degree_elevation_error():
    spline = _cubic_spline()
    with pytest.raises(geomcore.DegreeElevationError):
        spline.elevate_degree(100)
    with pytest.raises(geomcore.DegreeElevationError):
        _bilinear_surface().elevate_u_degree(100)


def test_bspline_conversion_error():
    hyperbola = curves.Hyperbola(
        geomcore.Point(0, 0, 0),
        geomcore.Vector(0, 0, 1),
        geomcore.Vector(1, 0, 0),
        1.0,
        1.0,
    )
    far = curves.TrimmedCurve(hyperbola, -800.0, 800.0)
    with pytest.raises(geomcore.BSplineConversionError):
        far.to_bspline()
    plane = surfaces.Plane(geomcore.Point(0, 0, 0), geomcore.Vector(0, 0, 1))
    with pytest.raises(geomcore.BSplineConversionError):
        plane.to_bspline_patch((1.0, 0.0), (0.0, 1.0))


def test_curve_3d_to_uv_error():
    sphere = surfaces.Sphere(geomcore.Point(0, 0, 0), 1.0)
    far_line = curves.Line.from_two_points(
        geomcore.Point(5, 0, 0), geomcore.Point(5, 1, 0)
    )
    with pytest.raises(geomcore.Curve3DToUvError):
        sphere.curve_3d_to_uv(far_line)


def test_uv_curve_to_3d_error():
    cylinder = surfaces.Cylinder(geomcore.Point(0, 0, 0), geomcore.Vector(0, 0, 1), 2.0)
    helix = curves.Line2D(geomcore.Point2D(0, 0), geomcore.Vector2D(0.6, 0.8))
    with pytest.raises(geomcore.UvCurveTo3DError):
        cylinder.uv_curve_to_3d(helix)


def test_surface_iso_error():
    surface = _bilinear_surface()
    with pytest.raises(geomcore.SurfaceIsoError):
        surface.u_iso(float("nan"))
    with pytest.raises(geomcore.SurfaceIsoError):
        surface.v_iso(5.0)
