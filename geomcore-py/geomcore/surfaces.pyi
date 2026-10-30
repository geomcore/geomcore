from typing import (
    Any,
    ClassVar,
    Generic,
    List,
    Optional,
    Tuple,
    TypeVar,
    Union,
    final,
    overload,
)

from geomcore import (
    Axis,
    BoundingBox,
    CurveSurfaceIntersection,
    Frame,
    Point,
    Side,
    SurfaceCurvatures,
    SurfaceProjection,
    SurfaceSurfaceIntersection,
    Tolerance,
    Transform,
    Vector,
)
from geomcore.curves import (
    BezierCurve,
    BSplineCurve,
    Circle,
    Line,
    _Curve,
    _Curve2D,
    _WholeCurve,
)

_Surface = Union[
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    ExtrusionSurface,
    RevolutionSurface,
    BSplineSurface,
    BezierSurface,
]

_UIsoCurveType = TypeVar("_UIsoCurveType")
_VIsoCurveType = TypeVar("_VIsoCurveType")
_Self = TypeVar("_Self")

# The method sets of geomcore-py/src/method_sets.rs, one stub-only class
# each, as in curves.pyi. The type arguments of `_SurfaceMethods` are the
# classes of the `u` and `v` isoparametric curves.

class _SurfaceMethods(Generic[_UIsoCurveType, _VIsoCurveType]):
    def transformed(self: _Self, transform: Transform) -> _Self:
        """The image of the surface under `transform`, parametrized by
        the image of its parametrization: a parameter that is a
        length scales with the transform. Raises
        `ConstructionError` when the transform collapses the surface
        (a zero scale) or overflows.
        """
    def reversed_u(self: _Self) -> _Self:
        """The same surface with the `u` direction reversed:
        `reversed_u().eval_point(u_reversed_parameter(u), v)` is
        `eval_point(u, v)`.
        """
    def reversed_v(self: _Self) -> _Self:
        """The same surface with the `v` direction reversed, as
        `reversed_u`.
        """
    def u_reversed_parameter(self, u: float) -> float:
        """Maps `u` onto the `reversed_u()` surface: the `u` parameter
        of the same point there.
        """
    def v_reversed_parameter(self, v: float) -> float:
        """Maps `v` onto the `reversed_v()` surface: the `v` parameter
        of the same point there.
        """
    def eval_point(self, u: float, v: float) -> Point:
        """Evaluate the point at surface parameters (u, v)."""
    def eval_points(self, uvs: List[Tuple[float, float]]) -> List[Point]:
        """Evaluate many (u, v) pairs at once."""
    def eval_derivative(self, u: float, v: float, du: int, dv: int) -> Vector:
        """Evaluate the partial derivative of order (du, dv), with
        1 <= du+dv <= 2.
        """
    def eval_d1(self, u: float, v: float) -> Tuple[Point, Tuple[Vector, Vector]]:
        """The point and the first partial derivatives `(du, dv)` at
        `(u, v)`.
        """
    def eval_d2(
        self, u: float, v: float
    ) -> Tuple[Point, Tuple[Vector, Vector], Tuple[Vector, Vector, Vector]]:
        """The point, the first partial derivatives `(du, dv)` and the
        second partial derivatives `(du du, du dv, dv dv)` at
        `(u, v)`.
        """
    def u_bounds(self) -> Tuple[float, float]:
        """The (first, last) `u` parameter bounds."""
    def v_bounds(self) -> Tuple[float, float]:
        """The (first, last) `v` parameter bounds."""
    def u_period(self) -> Optional[float]:
        """The period of the u parameter, or `None` when the surface is
        not periodic in u.
        """
    def v_period(self) -> Optional[float]:
        """The period of the v parameter, or `None` when the surface is
        not periodic in v.
        """
    def is_u_closed(self, tol: Optional[Tolerance] = ...) -> bool:
        """Whether travelling in `u` returns to its start: the boundary
        curves at the two `u` bounds coincide (`Tolerance.default()`
        when omitted). Periodic in `u` implies closed in `u`, not
        conversely.
        """
    def is_v_closed(self, tol: Optional[Tolerance] = ...) -> bool:
        """Whether travelling in `v` returns to its start, as
        `is_u_closed`.
        """
    def continuity_u(self, u: float, tol: Optional[Tolerance] = ...) -> Union[int, float]:
        """Continuity order `n` of `C^n` in `u` at parameter `u` as an
        `int`, or `math.inf` where smooth. A B-spline reports
        `u_degree - multiplicity` at a `u` knot (within
        `tol.parametric`); the ends of a non-periodic direction
        report smooth, and a periodic direction's seam counts as a
        knot. Use `BSplineSurface.eval_d1_sided` and its siblings to
        evaluate each side of a knot reporting a low order.
        """
    def continuity_v(self, v: float, tol: Optional[Tolerance] = ...) -> Union[int, float]:
        """Continuity order in `v` at parameter `v`, as `continuity_u`."""
    def min_continuity_u(self) -> Union[int, float]:
        """The lowest `continuity_u` anywhere inside the `u` range, ends
        excluded (`math.inf` when smooth everywhere).
        """
    def min_continuity_v(self) -> Union[int, float]:
        """The lowest `continuity_v` anywhere inside the `v` range, ends
        excluded (`math.inf` when smooth everywhere).
        """
    def curvatures(
        self, u: float, v: float, tol: Optional[Tolerance] = ...
    ) -> Optional[SurfaceCurvatures]:
        """Principal curvatures with their directions, and the mean and
        Gaussian curvatures, at `(u, v)` (`Tolerance.default()` when
        omitted), as a `SurfaceCurvatures`. `None` where the unit
        normal degenerates, as at a sphere's pole or a cone's apex,
        or a derivative is not finite. Curvatures are measured along
        the unit normal: negative bends away from it.
        """
    def project_point(
        self, point: Point, tol: Optional[Tolerance] = ...
    ) -> SurfaceProjection:
        """Projects `point` onto the surface: the closest point, as a
        `SurfaceProjection`.
        """
    def project_points(
        self, points: List[Point], tol: Optional[Tolerance] = ...
    ) -> List[SurfaceProjection]:
        """Projects each point in `points` onto the surface, returning a
        `SurfaceProjection` per point.
        """
    def contains(self, point: Point, tol: Optional[Tolerance] = ...) -> bool:
        """Returns whether `point` lies on the surface within `tol`
        (`Tolerance.default()` when omitted).
        """
    def parameters_of(
        self, point: Point, tol: Optional[Tolerance] = ...
    ) -> Tuple[float, float]:
        """The `(u, v)` parameters of `point`, a point on the surface,
        within `tol` (`Tolerance.default()` when omitted). Off the
        surface, a B-spline or Bézier surface gives the parameters of
        the closest point, and an elementary surface those of its
        closed-form inversion.
        """
    def bounding_box(self) -> Optional[BoundingBox]:
        """The axis-aligned box around the shape (`None` when unbounded)."""
    def u_iso(self, u: float) -> _UIsoCurveType:
        """The `u` isoparametric curve at `u`, parametrized by `v` like
        the surface: `u_iso(u).eval_point(t)` is `eval_point(u, t)`.
        Raises `SurfaceIsoError` for a non-finite parameter or one
        outside a non-periodic direction (a periodic one wraps).
        """
    def v_iso(self, v: float) -> _VIsoCurveType:
        """The `v` isoparametric curve at `v`, parametrized by `u` like
        the surface, as `u_iso`.
        """
    def pcurve_to_curve(
        self, pcurve: _Curve2D, tol: Optional[Tolerance] = ...
    ) -> _WholeCurve:
        """The exact 3D curve that `pcurve`, a curve in this surface's (u, v)
        parameter space, traces on the surface, in the same parameter: its
        point at `t` is the surface's point at the pcurve's (u, v) at `t`.

        Closed forms: any curve on a plane, as the same kind of curve; any
        curve on the extrusion of a line, where the affine image keeps its
        kind in the same parameter; and a line along u or v on any surface,
        as the iso curve run from the line's origin. Raises
        `PcurveToCurveError` for any other pair, when the line runs
        along an iso curve that is a point (a sphere's pole, a cone's apex),
        or when the surface has no iso curve there. `tol` is
        `Tolerance.default()` when omitted.
        """
    @overload
    def intersect(
        self, other: _Curve, tol: Optional[Tolerance] = ...
    ) -> CurveSurfaceIntersection:
        """Intersects with any curve or surface, in either order.

        Returns a `CurveSurfaceIntersection` or
        `SurfaceSurfaceIntersection`, by pair family, holding this
        object as `first` and `other` as `second`. A pair with a
        closed form is solved exactly (`method` is `EXACT`), any
        other pair numerically (`NUMERIC`), whose `completeness`
        tells whether it proved that nothing was missed. Raises
        `TypeError` for a 2D curve, and for an object that is not a
        geomcore curve or surface. Other Python threads run
        meanwhile.
        """
    @overload
    def intersect(
        self, other: _Surface, tol: Optional[Tolerance] = ...
    ) -> SurfaceSurfaceIntersection: ...

class _SurfaceNormalMethods:
    def normal(
        self, u: float, v: float, tol: Optional[Tolerance] = ...
    ) -> Optional[Vector]:
        """Unit normal at surface parameters `(u, v)`
        (`Tolerance.default()` when omitted). Singular points such as
        poles give their limiting normal; returns `None` at an apex,
        whose meridians each have their own, and where no normal can
        be recovered.
        """

class _SidedSurfaceEvaluationMethods:
    def eval_d1_sided(
        self, u: float, v: float, u_side: Side, v_side: Side
    ) -> Tuple[Point, Tuple[Vector, Vector]]:
        """The point and the first partial derivatives `(du, dv)` at
        `(u, v)`, each direction from its own side: `Side.BELOW`
        takes the span ending at the parameter, `Side.ABOVE` the span
        starting there. Both sides agree strictly inside a span, so
        everywhere on a Bézier surface. Plain evaluation takes the
        above side, except at the upper bound of a non-periodic
        direction. A periodic direction's seam is a knot like any
        other: at either bound, `BELOW` takes the last span of the
        period and `ABOVE` the first.
        """
    def eval_d2_sided(
        self, u: float, v: float, u_side: Side, v_side: Side
    ) -> Tuple[Point, Tuple[Vector, Vector], Tuple[Vector, Vector, Vector]]:
        """The point, the first partial derivatives `(du, dv)` and the
        second partial derivatives `(du du, du dv, dv dv)` at
        `(u, v)`, each direction from its own side, as
        `eval_d1_sided`.
        """
    def eval_derivative_sided(
        self, u: float, v: float, du: int, dv: int, u_side: Side, v_side: Side
    ) -> Vector:
        """The partial derivative of order `(du, dv)`, with
        `1 <= du + dv <= 2`, at `(u, v)`, each direction from its own
        side, as `eval_d1_sided`.
        """

class _SurfaceExtremaMethods:
    def extrema(
        self, point: Point, tol: Optional[Tolerance] = ...
    ) -> List[SurfaceProjection]:
        """The local minima of the distance from `point` to the
        surface, as `SurfaceProjection` objects ordered by ascending
        distance.
        """


@final
class Plane(_SurfaceMethods[Line, Line]):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, point: Point, normal: Vector) -> Plane: ...
    @staticmethod
    def from_frame(frame: Frame) -> Plane: ...
    @staticmethod
    def from_three_points(p1: Point, p2: Point, p3: Point) -> Plane: ...
    @staticmethod
    def from_coefficients(a: float, b: float, c: float, d: float) -> Plane: ...
    def normal(self) -> Vector: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...
    def frame(self) -> Frame: ...

@final
class Cylinder(_SurfaceMethods[Line, Circle], _SurfaceNormalMethods):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, center: Point, axis_direction: Vector, radius: float) -> Cylinder: ...
    @staticmethod
    def from_frame(frame: Frame, radius: float) -> Cylinder: ...
    @staticmethod
    def from_axis(axis: Axis, radius: float) -> Cylinder: ...
    @staticmethod
    def from_circle(circle: Circle) -> Cylinder: ...
    def radius(self) -> float: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...
    def frame(self) -> Frame: ...
    def axis(self) -> Axis: ...

@final
class Cone(_SurfaceMethods[Line, Circle], _SurfaceNormalMethods):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, center: Point, axis_direction: Vector, semi_angle: float, ref_radius: float
    ) -> Cone: ...
    @staticmethod
    def from_frame(frame: Frame, semi_angle: float, ref_radius: float) -> Cone: ...
    @staticmethod
    def from_two_points_and_radii(
        p1: Point, p2: Point, r1: float, r2: float
    ) -> Cone: ...
    def apex(self) -> Point: ...
    def semi_angle(self) -> float: ...
    def ref_radius(self) -> float: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...
    def frame(self) -> Frame: ...

@final
class Sphere(_SurfaceMethods[Circle, Circle], _SurfaceNormalMethods):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, center: Point, radius: float) -> Sphere: ...
    @staticmethod
    def from_frame(frame: Frame, radius: float) -> Sphere: ...
    def center(self) -> Point: ...
    def radius(self) -> float: ...
    def area(self) -> float: ...
    def volume(self) -> float: ...
    def to_bspline(self) -> BSplineSurface: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...
    def frame(self) -> Frame: ...

@final
class Torus(_SurfaceMethods[Circle, Circle], _SurfaceNormalMethods):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, center: Point, normal: Vector, major_radius: float, minor_radius: float
    ) -> Torus: ...
    @staticmethod
    def from_frame(frame: Frame, major_radius: float, minor_radius: float) -> Torus: ...
    def major_radius(self) -> float: ...
    def minor_radius(self) -> float: ...
    def to_bspline(self) -> BSplineSurface: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...
    def frame(self) -> Frame: ...
    def center(self) -> Point: ...

@final
class ExtrusionSurface(_SurfaceMethods[Line, _Curve], _SurfaceNormalMethods):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, profile: _Curve, direction: Vector) -> ExtrusionSurface: ...
    def profile(self) -> _Curve: ...
    def direction(self) -> Vector: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...

@final
class RevolutionSurface(_SurfaceMethods[_Curve, Circle], _SurfaceNormalMethods):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, profile: _Curve, axis: Axis) -> RevolutionSurface: ...
    def profile(self) -> _Curve: ...
    def axis(self) -> Axis: ...
    def to_bspline(self) -> BSplineSurface: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...

@final
class BSplineSurface(
    _SurfaceMethods[BSplineCurve, BSplineCurve],
    _SurfaceNormalMethods,
    _SidedSurfaceEvaluationMethods,
    _SurfaceExtremaMethods,
):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        u_degree: int,
        v_degree: int,
        poles: List[List[Point]],
        u_knots: List[float],
        u_multiplicities: List[int],
        v_knots: List[float],
        v_multiplicities: List[int],
        u_periodic: bool,
        v_periodic: bool,
    ) -> BSplineSurface: ...
    @staticmethod
    def new_rational(
        u_degree: int,
        v_degree: int,
        poles: List[List[Point]],
        weights: List[List[float]],
        u_knots: List[float],
        u_multiplicities: List[int],
        v_knots: List[float],
        v_multiplicities: List[int],
        u_periodic: bool,
        v_periodic: bool,
    ) -> BSplineSurface: ...
    def u_degree(self) -> int: ...
    def v_degree(self) -> int: ...
    def is_u_periodic(self) -> bool: ...
    def is_v_periodic(self) -> bool: ...
    def is_rational(self) -> bool: ...
    def u_pole_count(self) -> int: ...
    def v_pole_count(self) -> int: ...
    def poles(self) -> List[List[Point]]: ...
    def weights(self) -> Optional[List[List[float]]]: ...
    def u_knots(self) -> List[float]: ...
    def u_multiplicities(self) -> List[int]: ...
    def v_knots(self) -> List[float]: ...
    def v_multiplicities(self) -> List[int]: ...
    @staticmethod
    def approximate(
        points: List[List[Point]],
        u_degree: int,
        v_degree: int,
        u_pole_count: int,
        v_pole_count: int,
        parametrization: str = ...,
    ) -> BSplineSurface: ...
    @staticmethod
    def approximate_within(
        points: List[List[Point]],
        u_degree: int,
        v_degree: int,
        max_deviation: float,
        parametrization: str = ...,
    ) -> BSplineSurface: ...
    def insert_u_knot(self, u: float, times: int) -> BSplineSurface: ...
    def insert_v_knot(self, v: float, times: int) -> BSplineSurface: ...
    def refine_u_knots(self, knots: List[float]) -> BSplineSurface: ...
    def refine_v_knots(self, knots: List[float]) -> BSplineSurface: ...
    def elevate_u_degree(self, amount: int) -> BSplineSurface: ...
    def elevate_v_degree(self, amount: int) -> BSplineSurface: ...
    def remove_u_knot(
        self, u: float, times: int, max_deviation: float
    ) -> Tuple[BSplineSurface, int]: ...
    def remove_v_knot(
        self, v: float, times: int, max_deviation: float
    ) -> Tuple[BSplineSurface, int]: ...
    def remove_u_knots(
        self, knots: List[float], max_deviation: float
    ) -> Tuple[BSplineSurface, int]: ...
    def remove_v_knots(
        self, knots: List[float], max_deviation: float
    ) -> Tuple[BSplineSurface, int]: ...
    def to_bezier_patches(self) -> List[List[BezierSurface]]: ...
    def extract(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...
    def extract_u(self, first: float, last: float) -> BSplineSurface: ...
    def extract_v(self, first: float, last: float) -> BSplineSurface: ...
    def split_u_at(self, u: float) -> Tuple[BSplineSurface, BSplineSurface]: ...
    def split_v_at(self, v: float) -> Tuple[BSplineSurface, BSplineSurface]: ...

@final
class BezierSurface(
    _SurfaceMethods[BezierCurve, BezierCurve],
    _SurfaceNormalMethods,
    _SidedSurfaceEvaluationMethods,
    _SurfaceExtremaMethods,
):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, u_degree: int, v_degree: int, poles: List[List[Point]]
    ) -> BezierSurface: ...
    @staticmethod
    def new_rational(
        u_degree: int,
        v_degree: int,
        poles: List[List[Point]],
        weights: List[List[float]],
    ) -> BezierSurface: ...
    @staticmethod
    def from_bspline(surface: BSplineSurface) -> BezierSurface: ...
    def to_bspline(self) -> BSplineSurface: ...
    def u_degree(self) -> int: ...
    def v_degree(self) -> int: ...
    def is_rational(self) -> bool: ...
    def u_pole_count(self) -> int: ...
    def v_pole_count(self) -> int: ...
    def poles(self) -> List[List[Point]]: ...
    def weights(self) -> Optional[List[List[float]]]: ...
    def split_u_at(self, u: float) -> Tuple[BezierSurface, BezierSurface]: ...
    def split_v_at(self, v: float) -> Tuple[BezierSurface, BezierSurface]: ...
    def elevate_u_degree(self, amount: int) -> BezierSurface: ...
    def elevate_v_degree(self, amount: int) -> BezierSurface: ...

__all__ = [
    "BSplineSurface",
    "BezierSurface",
    "Cone",
    "Cylinder",
    "ExtrusionSurface",
    "Plane",
    "RevolutionSurface",
    "Sphere",
    "Torus",
]
