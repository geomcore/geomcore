"""Surface tessellation: numpy arrays within a deviation (#252)."""

import math

import numpy
import pytest

import geomcore
from geomcore import Cylinder, Plane, Point, Sphere, TessellationError, Torus, Vector


def _deviation_at_triangle_middles(surface, tessellation):
    """The largest distance from the middle of a triangle to the surface."""
    vertices = tessellation.vertices
    middles = vertices[tessellation.triangles].mean(axis=1)
    return max(surface.project_point(Point(*middle)).distance for middle in middles)


def test_a_sphere_gives_numpy_arrays_of_matching_shapes():
    sphere = Sphere(Point(0, 0, 0), 2.0)
    tessellation = sphere.tessellate(1e-2)
    vertices, normals, uvs, triangles = (
        tessellation.vertices,
        tessellation.normals,
        tessellation.uvs,
        tessellation.triangles,
    )
    count = len(vertices)
    assert vertices.shape == (count, 3) and vertices.dtype == numpy.float64
    assert normals.shape == (count, 3) and normals.dtype == numpy.float64
    assert uvs.shape == (count, 2) and uvs.dtype == numpy.float64
    assert triangles.shape == (len(tessellation), 3) and triangles.dtype == numpy.uint32
    assert triangles.max() < count
    assert tessellation.is_proven
    assert 0 < tessellation.deviation <= 1e-2
    assert repr(tessellation) == f"Tessellation({count} vertices, {len(tessellation)} triangles)"


def test_the_vertices_lie_on_the_surface_at_their_parameters():
    torus = Torus(Point(0, 0, 0), Vector(0, 0, 1), 3.0, 1.0)
    tessellation = torus.tessellate(1e-2)
    numpy.testing.assert_allclose(
        torus.points_at_array(tessellation.uvs), tessellation.vertices, atol=1e-12
    )
    numpy.testing.assert_allclose(numpy.linalg.norm(tessellation.normals, axis=1), 1.0)
    assert _deviation_at_triangle_middles(torus, tessellation) <= tessellation.deviation


def test_the_normals_point_outward_from_a_sphere():
    tessellation = Sphere(Point(0, 0, 0), 1.0).tessellate(1e-1)
    numpy.testing.assert_allclose(tessellation.normals, tessellation.vertices, atol=1e-12)


def test_an_unbounded_direction_needs_a_range():
    cylinder = Cylinder(Point(0, 0, 0), Vector(0, 0, 1), 1.0)
    with pytest.raises(TessellationError):
        cylinder.tessellate(1e-3)
    tessellation = cylinder.tessellate(1e-3, v_range=(0.0, 2.0))
    # The fewest equal steps around with 1 - cos(π / n) within 1e-3.
    assert len(tessellation) == 71 * 2
    assert tessellation.deviation == pytest.approx(1 - math.cos(math.pi / 71))


def test_a_part_of_a_plane_is_two_triangles():
    plane = Plane(Point(0, 0, 0), Vector(0, 0, 1))
    tessellation = plane.tessellate(1e-3, u_range=(0.0, 2.0), v_range=(-1.0, 1.0))
    assert len(tessellation) == 2
    assert tessellation.deviation == 0.0


def test_a_trimmed_surface_tessellates_its_own_ranges():
    sphere = Sphere(Point(0, 0, 0), 1.0)
    trimmed = sphere.trimmed((0.0, 1.0), (0.0, 0.5))
    uvs = trimmed.tessellate(1e-2).uvs
    assert uvs[:, 0].min() == 0.0 and uvs[:, 0].max() == 1.0
    assert uvs[:, 1].min() == 0.0 and uvs[:, 1].max() == 0.5


@pytest.mark.parametrize(
    "arguments",
    [
        {"max_deviation": 0.0},
        {"max_deviation": -1.0},
        {"max_deviation": math.nan},
        {"max_deviation": 1e-2, "u_range": (0.0, 10.0)},
        {"max_deviation": 1e-2, "v_range": (0.0, 3.0)},
        {"max_deviation": 1e-2, "u_range": (1.0, 1.0)},
        {"max_deviation": 1e-300},
    ],
)
def test_invalid_arguments_raise_a_tessellation_error(arguments):
    with pytest.raises(TessellationError):
        Sphere(Point(0, 0, 0), 1.0).tessellate(**arguments)
    assert issubclass(TessellationError, geomcore.GeomcoreError)
