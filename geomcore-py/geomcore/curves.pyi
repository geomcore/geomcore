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
    Axis2D,
    Axis,
    BoundingBox,
    BoundingBox2D,
    CurveCurveIntersection,
    CurveCurveIntersection2D,
    CurveProjection,
    CurveProjection2D,
    CurveSurfaceIntersection,
    Discretization2D,
    Discretization,
    Frame2D,
    Frame,
    Point2D,
    Point,
    Side,
    Tolerance,
    Transform,
    Vector2D,
    Vector,
)
from geomcore.surfaces import _Surface

# The whole curves are the classes of the `Curve` enum's variants.
_WholeCurve = Union[
    Line, Circle, Ellipse, Parabola, Hyperbola, BSplineCurve, BezierCurve
]
_Curve = Union[_WholeCurve, TrimmedCurve]
_Curve2D = Union[
    Line2D,
    Circle2D,
    Ellipse2D,
    Parabola2D,
    Hyperbola2D,
    BSplineCurve2D,
    BezierCurve2D,
    TrimmedCurve2D,
]

_PointType = TypeVar("_PointType")
_VectorType = TypeVar("_VectorType")
_BoundingBoxType = TypeVar("_BoundingBoxType")
_ProjectionType = TypeVar("_ProjectionType")
_DiscretizationType = TypeVar("_DiscretizationType")
_Curve2DType = TypeVar("_Curve2DType")
_Self = TypeVar("_Self")

# The method sets of geomcore-py/src/method_sets.rs, one stub-only class
# each. A curve class inherits the sets its `shape_methods!` line names, so
# each shared method is typed and documented here once; the type arguments
# carry what differs by class or dimension.

class _CurveMethodsInEitherDimension(
    Generic[_PointType, _VectorType, _BoundingBoxType, _ProjectionType]
):
    def eval_point(self, u: float) -> _PointType:
        """Evaluate the point at parameter `u`."""
    def eval_points(self, us: List[float]) -> List[_PointType]:
        """Evaluate many parameters at once."""
    def eval_derivative(self, u: float, order: int) -> _VectorType:
        """Evaluate the derivative of the given order (>= 1) at `u`."""
    def eval_d1(self, u: float) -> Tuple[_PointType, _VectorType]:
        """The point and the first derivative at parameter `u`."""
    def eval_d2(self, u: float) -> Tuple[_PointType, _VectorType, _VectorType]:
        """The point and the first two derivatives at parameter `u`."""
    def tangent(self, u: float, tol: Optional[Tolerance] = ...) -> Optional[_VectorType]:
        """Unit tangent at parameter `u` (`Tolerance.default()` when
        omitted), or `None` where the speed vanishes, as at a cusp.
        The speed counts as zero when the distance it would cover
        across the parameter range is within `tol.confusion`, so
        rescaling the parameter does not change the answer.
        """
    def curvature(self, u: float, tol: Optional[Tolerance] = ...) -> Optional[float]:
        """Curvature at parameter `u`: the turn rate of the unit tangent
        per unit arc length, zero on a straight run. `None` where the
        speed vanishes, by the test of `tangent`.
        """
    def curvature_center(
        self, u: float, tol: Optional[Tolerance] = ...
    ) -> Optional[_PointType]:
        """Centre of the osculating circle at parameter `u`, or `None`
        where the speed vanishes or the curve runs straight (the
        centre is then at infinity).
        """
    def bounds(self) -> Tuple[float, float]:
        """The (first, last) parameter bounds."""
    def period(self) -> Optional[float]:
        """The period of the parameter, or `None` when the curve is not
        periodic.
        """
    def is_periodic(self) -> bool:
        """Whether the curve is periodic."""
    def is_closed(self, tol: Optional[Tolerance] = ...) -> bool:
        """Whether travelling along the curve returns to its start
        (`Tolerance.default()` when omitted). A periodic curve is
        closed; a clamped B-spline whose ends coincide is closed
        without being periodic.
        """
    def continuity(self, u: float, tol: Optional[Tolerance] = ...) -> Union[int, float]:
        """Continuity order `n` of `C^n` at parameter `u` as an `int`,
        or `math.inf` where no knot constrains the curve. A B-spline
        reports `degree - multiplicity` at a knot (within
        `tol.parametric`); the ends of a non-periodic curve report
        smooth, and a periodic curve's seam counts as a knot. Use
        the `*_sided` evaluations of a B-spline to evaluate each
        side of a knot reporting a low order.
        """
    def min_continuity(self) -> Union[int, float]:
        """The lowest `continuity` anywhere inside the bounds, ends
        excluded (`math.inf` when smooth everywhere).
        """
    def extrema(
        self, point: _PointType, tol: Optional[Tolerance] = ...
    ) -> List[_ProjectionType]:
        """The stationary points of the distance from `point` to the
        curve (local minima plus the ends of a bounded curve), as
        projection objects ordered by ascending distance.
        """
    def project_point(
        self, point: _PointType, tol: Optional[Tolerance] = ...
    ) -> _ProjectionType:
        """Projects `point` onto the curve: the closest point, as a
        projection object.
        """
    def project_points(
        self, points: List[_PointType], tol: Optional[Tolerance] = ...
    ) -> List[_ProjectionType]:
        """Projects each point in `points` onto the curve, returning a
        projection object per point.
        """
    def contains(self, point: _PointType, tol: Optional[Tolerance] = ...) -> bool:
        """Returns whether `point` lies on the curve within `tol`
        (`Tolerance.default()` when omitted).
        """
    def parameter_of(self, point: _PointType, tol: Optional[Tolerance] = ...) -> float:
        """The parameter of `point`, a point on the curve, within `tol`
        (`Tolerance.default()` when omitted). Off the curve, a
        B-spline or Bézier curve gives the parameter of the closest
        point, and a line or conic that of its closed-form
        inversion.
        """
    def reversed(self: _Self) -> _Self:
        """The same curve traversed the other way:
        `reversed().eval_point(reversed_parameter(u))` is
        `eval_point(u)`.
        """
    def reversed_parameter(self, u: float) -> float:
        """Maps `u` onto the `reversed()` curve: the parameter of the
        same point there. It is `first + last - u` over bounded
        parameters (`2π - u` on a circle or ellipse), `-u` on a
        line, parabola or hyperbola, and a trimmed curve maps as
        its basis.
        """
    def bounding_box(self) -> Optional[_BoundingBoxType]:
        """The axis-aligned box around the shape (`None` when unbounded)."""

class _CurveMethods(
    _CurveMethodsInEitherDimension[Point, Vector, BoundingBox, CurveProjection]
):
    def transformed(self: _Self, transform: Transform) -> _Self:
        """The image of the curve under `transform`, parametrized by
        the image of its parametrization: a line's or parabola's
        parameter, a length, scales with the transform. Raises
        `ConstructionError` when the transform collapses the
        curve (a zero scale) or overflows.
        """
    @overload
    def intersect(
        self, other: _Curve, tol: Optional[Tolerance] = ...
    ) -> CurveCurveIntersection:
        """Intersects with any curve or surface, in either order.

        Returns a `CurveCurveIntersection`, `CurveCurveIntersection2D`,
        `CurveSurfaceIntersection` or `SurfaceSurfaceIntersection`, by
        pair family, holding this object as `first` and `other` as
        `second`. A pair with a closed form is solved exactly
        (`method` is `EXACT`), any other pair by a numeric search
        (`NUMERIC`), whose `completeness` tells whether it proved that
        nothing was missed. Raises `TypeError` for a 2D curve against
        a 3D curve or a surface, and for an object that is not a
        geomcore curve or surface. Other Python threads run meanwhile.
        """
    @overload
    def intersect(
        self, other: _Surface, tol: Optional[Tolerance] = ...
    ) -> CurveSurfaceIntersection: ...

class _Curve2DMethods(
    _CurveMethodsInEitherDimension[Point2D, Vector2D, BoundingBox2D, CurveProjection2D]
):
    def intersect(
        self, other: _Curve2D, tol: Optional[Tolerance] = ...
    ) -> CurveCurveIntersection2D:
        """Intersects with any curve or surface, in either order.

        Returns a `CurveCurveIntersection`, `CurveCurveIntersection2D`,
        `CurveSurfaceIntersection` or `SurfaceSurfaceIntersection`, by
        pair family, holding this object as `first` and `other` as
        `second`. A pair with a closed form is solved exactly
        (`method` is `EXACT`), any other pair by a numeric search
        (`NUMERIC`), whose `completeness` tells whether it proved that
        nothing was missed. Raises `TypeError` for a 2D curve against
        a 3D curve or a surface, and for an object that is not a
        geomcore curve or surface. Other Python threads run meanwhile.
        """

class _DiscretizeMethods(Generic[_DiscretizationType]):
    def discretize(
        self, max_chord_deviation: float, max_angle_deviation: Optional[float] = ...
    ) -> _DiscretizationType:
        """Discretize the curve into points joined by straight
        segments, each segment's chord deviation at most
        `max_chord_deviation` and, if given, its angle deviation at
        most `max_angle_deviation` (radians).
        """
    def point_count(
        self, max_chord_deviation: float, max_angle_deviation: Optional[float] = ...
    ) -> int:
        """How many points `discretize` returns for the same arguments,
        without computing them.
        """

class _SidedEvaluationMethods(Generic[_PointType, _VectorType]):
    def eval_d1_sided(self, u: float, side: Side) -> Tuple[_PointType, _VectorType]:
        """The point and the first derivative at parameter `u`, from
        one `side`: `Side.BELOW` takes the span ending at `u`,
        `Side.ABOVE` the span starting there. Both sides agree
        strictly inside a span, so everywhere on a Bézier curve.
        Plain evaluation takes the above side, except at the upper
        bound of a non-periodic curve. A periodic curve's seam is a
        knot like any other: at either bound, `BELOW` takes the last
        span of the period and `ABOVE` the first.
        """
    def eval_d2_sided(
        self, u: float, side: Side
    ) -> Tuple[_PointType, _VectorType, _VectorType]:
        """The point and the first two derivatives at parameter `u`,
        from one `side`, as `eval_d1_sided`.
        """
    def eval_derivative_sided(self, u: float, order: int, side: Side) -> _VectorType:
        """The derivative of the given order (>= 1) at parameter `u`,
        from one `side`, as `eval_d1_sided`.
        """

class _ParametrizeOnMethods(Generic[_Curve2DType]):
    def parametrize_on(
        self, surface: Any, tol: Optional[Tolerance] = ...
    ) -> _Curve2DType:
        """Compute this curve's 2D representation in a surface's
        (u, v) space.

        Raises `ValueError` if no closed-form representation exists
        for the pair, or if the curve does not lie on the surface
        within `tol` (`Tolerance.default()` when omitted).
        """


@final
class Line(_CurveMethods, _ParametrizeOnMethods[Any]):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, origin: Point, direction: Vector) -> Line: ...
    @staticmethod
    def from_axis(axis: Axis) -> Line: ...
    @staticmethod
    def from_two_points(p1: Point, p2: Point) -> Line: ...
    def origin(self) -> Point: ...
    def direction(self) -> Vector: ...
    def axis(self) -> Axis: ...

@final
class Circle(
    _CurveMethods,
    _DiscretizeMethods[Discretization],
    _ParametrizeOnMethods[Any],
):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, center: Point, normal: Vector, radius: float) -> Circle: ...
    @staticmethod
    def from_axis(axis: Axis, radius: float) -> Circle: ...
    @staticmethod
    def from_frame(frame: Frame, radius: float) -> Circle: ...
    @staticmethod
    def from_three_points(p1: Point, p2: Point, p3: Point) -> Circle: ...
    def center(self) -> Point: ...
    def radius(self) -> float: ...
    def circumference(self) -> float: ...
    def disk_area(self) -> float: ...
    def normal(self) -> Vector: ...
    def to_bspline(self) -> BSplineCurve: ...
    def frame(self) -> Frame: ...

@final
class Ellipse(
    _CurveMethods,
    _DiscretizeMethods[Discretization],
    _ParametrizeOnMethods[Ellipse2D],
):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        center: Point,
        normal: Vector,
        x_direction: Vector,
        major_radius: float,
        minor_radius: float,
    ) -> Ellipse: ...
    @staticmethod
    def from_frame(frame: Frame, major_radius: float, minor_radius: float) -> Ellipse: ...
    @staticmethod
    def from_center_and_points(center: Point, s1: Point, s2: Point) -> Ellipse: ...
    def center(self) -> Point: ...
    def major_radius(self) -> float: ...
    def minor_radius(self) -> float: ...
    def area(self) -> float: ...
    def to_bspline(self) -> BSplineCurve: ...
    def frame(self) -> Frame: ...

@final
class Parabola(_CurveMethods, _ParametrizeOnMethods[Parabola2D]):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, apex: Point, normal: Vector, x_direction: Vector, focal: float
    ) -> Parabola: ...
    @staticmethod
    def from_frame(frame: Frame, focal: float) -> Parabola: ...
    def focal(self) -> float: ...
    def frame(self) -> Frame: ...
    def apex(self) -> Point: ...

@final
class Hyperbola(_CurveMethods, _ParametrizeOnMethods[Hyperbola2D]):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        center: Point,
        normal: Vector,
        x_direction: Vector,
        major_radius: float,
        minor_radius: float,
    ) -> Hyperbola: ...
    @staticmethod
    def from_frame(frame: Frame, major_radius: float, minor_radius: float) -> Hyperbola: ...
    @staticmethod
    def from_center_and_points(center: Point, s1: Point, s2: Point) -> Hyperbola: ...
    def center(self) -> Point: ...
    def major_radius(self) -> float: ...
    def minor_radius(self) -> float: ...
    def frame(self) -> Frame: ...

@final
class BSplineCurve(
    _CurveMethods,
    _DiscretizeMethods[Discretization],
    _SidedEvaluationMethods[Point, Vector],
    _ParametrizeOnMethods[BSplineCurve2D],
):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        degree: int,
        poles: List[Point],
        knots: List[float],
        multiplicities: List[int],
        periodic: bool,
    ) -> BSplineCurve: ...
    @staticmethod
    def new_rational(
        degree: int,
        poles: List[Point],
        weights: List[float],
        knots: List[float],
        multiplicities: List[int],
        periodic: bool,
    ) -> BSplineCurve: ...
    @staticmethod
    def interpolate(
        points: List[Point], degree: int, parametrization: str = ...
    ) -> BSplineCurve: ...
    @staticmethod
    def approximate(
        points: List[Point], degree: int, num_poles: int, parametrization: str = ...
    ) -> BSplineCurve: ...
    @staticmethod
    def approximate_within(
        points: List[Point],
        degree: int,
        max_deviation: float,
        parametrization: str = ...,
    ) -> BSplineCurve: ...
    def degree(self) -> int: ...
    def is_rational(self) -> bool: ...
    def poles(self) -> List[Point]: ...
    def weights(self) -> Optional[List[float]]: ...
    def knots(self) -> List[float]: ...
    def multiplicities(self) -> List[int]: ...
    def insert_knot(self, u: float, times: int) -> BSplineCurve: ...
    def refine_knots(self, knots: List[float]) -> BSplineCurve: ...
    def elevate_degree(self, amount: int) -> BSplineCurve: ...
    def remove_knot(
        self, u: float, times: int, max_deviation: float
    ) -> Tuple[BSplineCurve, int]: ...
    def remove_knots(
        self, knots: List[float], max_deviation: float
    ) -> Tuple[BSplineCurve, int]: ...
    def to_bezier_pieces(self) -> List[BezierCurve]: ...
    def extract(self, first: float, last: float) -> BSplineCurve: ...
    def split_at(self, u: float) -> Tuple[BSplineCurve, BSplineCurve]: ...
    def length(self, max_error: float) -> float: ...

@final
class BSplineCurve2D(
    _Curve2DMethods,
    _DiscretizeMethods[Discretization2D],
    _SidedEvaluationMethods[Point2D, Vector2D],
):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        degree: int,
        poles: List[Point2D],
        knots: List[float],
        multiplicities: List[int],
        periodic: bool,
    ) -> BSplineCurve2D: ...
    @staticmethod
    def new_rational(
        degree: int,
        poles: List[Point2D],
        weights: List[float],
        knots: List[float],
        multiplicities: List[int],
        periodic: bool,
    ) -> BSplineCurve2D: ...
    @staticmethod
    def interpolate(
        points: List[Point2D], degree: int, parametrization: str = ...
    ) -> BSplineCurve2D: ...
    @staticmethod
    def approximate(
        points: List[Point2D], degree: int, num_poles: int, parametrization: str = ...
    ) -> BSplineCurve2D: ...
    @staticmethod
    def approximate_within(
        points: List[Point2D],
        degree: int,
        max_deviation: float,
        parametrization: str = ...,
    ) -> BSplineCurve2D: ...
    def degree(self) -> int: ...
    def is_rational(self) -> bool: ...
    def poles(self) -> List[Point2D]: ...
    def weights(self) -> Optional[List[float]]: ...
    def knots(self) -> List[float]: ...
    def multiplicities(self) -> List[int]: ...
    def insert_knot(self, u: float, times: int) -> BSplineCurve2D: ...
    def refine_knots(self, knots: List[float]) -> BSplineCurve2D: ...
    def elevate_degree(self, amount: int) -> BSplineCurve2D: ...
    def remove_knot(
        self, u: float, times: int, max_deviation: float
    ) -> Tuple[BSplineCurve2D, int]: ...
    def remove_knots(
        self, knots: List[float], max_deviation: float
    ) -> Tuple[BSplineCurve2D, int]: ...
    def to_bezier_pieces(self) -> List[BezierCurve2D]: ...
    def extract(self, first: float, last: float) -> BSplineCurve2D: ...
    def split_at(self, u: float) -> Tuple[BSplineCurve2D, BSplineCurve2D]: ...
    def length(self, max_error: float) -> float: ...

@final
class BezierCurve(
    _CurveMethods,
    _DiscretizeMethods[Discretization],
    _SidedEvaluationMethods[Point, Vector],
):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, degree: int, poles: List[Point]) -> BezierCurve: ...
    @staticmethod
    def new_rational(
        degree: int, poles: List[Point], weights: List[float]
    ) -> BezierCurve: ...
    @staticmethod
    def from_bspline(curve: BSplineCurve) -> BezierCurve: ...
    def to_bspline(self) -> BSplineCurve: ...
    def degree(self) -> int: ...
    def is_rational(self) -> bool: ...
    def poles(self) -> List[Point]: ...
    def weights(self) -> Optional[List[float]]: ...
    def split_at(self, u: float) -> Tuple[BezierCurve, BezierCurve]: ...
    def elevate_degree(self, amount: int) -> BezierCurve: ...
    def length(self, max_error: float) -> float: ...

@final
class BezierCurve2D(
    _Curve2DMethods,
    _DiscretizeMethods[Discretization2D],
    _SidedEvaluationMethods[Point2D, Vector2D],
):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, degree: int, poles: List[Point2D]) -> BezierCurve2D: ...
    @staticmethod
    def new_rational(
        degree: int, poles: List[Point2D], weights: List[float]
    ) -> BezierCurve2D: ...
    @staticmethod
    def from_bspline(curve: BSplineCurve2D) -> BezierCurve2D: ...
    def to_bspline(self) -> BSplineCurve2D: ...
    def degree(self) -> int: ...
    def is_rational(self) -> bool: ...
    def poles(self) -> List[Point2D]: ...
    def weights(self) -> Optional[List[float]]: ...
    def split_at(self, u: float) -> Tuple[BezierCurve2D, BezierCurve2D]: ...
    def elevate_degree(self, amount: int) -> BezierCurve2D: ...
    def length(self, max_error: float) -> float: ...

@final
class Line2D(_Curve2DMethods):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, origin: Point2D, direction: Vector2D) -> Line2D: ...
    def origin(self) -> Point2D: ...
    def direction(self) -> Vector2D: ...
    def __repr__(self) -> str: ...
    @staticmethod
    def from_axis(axis: Axis2D) -> Line2D: ...
    def axis(self) -> Axis2D: ...
    @staticmethod
    def from_two_points(p1: Point2D, p2: Point2D) -> Line2D: ...

@final
class Circle2D(_Curve2DMethods, _DiscretizeMethods[Discretization2D]):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, center: Point2D, radius: float) -> Circle2D: ...
    def center(self) -> Point2D: ...
    def radius(self) -> float: ...
    def to_bspline(self) -> BSplineCurve2D: ...
    def __repr__(self) -> str: ...
    def frame(self) -> Frame2D: ...

@final
class Ellipse2D(_Curve2DMethods, _DiscretizeMethods[Discretization2D]):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        center: Point2D,
        x_direction: Vector2D,
        major_radius: float,
        minor_radius: float,
    ) -> Ellipse2D: ...
    @staticmethod
    def from_frame(
        frame: Frame2D, major_radius: float, minor_radius: float
    ) -> Ellipse2D: ...
    @staticmethod
    def from_center_and_points(
        center: Point2D, s1: Point2D, s2: Point2D
    ) -> Ellipse2D: ...
    def frame(self) -> Frame2D: ...
    def center(self) -> Point2D: ...
    def major_radius(self) -> float: ...
    def minor_radius(self) -> float: ...
    def area(self) -> float: ...
    def to_bspline(self) -> BSplineCurve2D: ...
    def __repr__(self) -> str: ...

@final
class Parabola2D(_Curve2DMethods):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, apex: Point2D, x_direction: Vector2D, focal: float) -> Parabola2D: ...
    @staticmethod
    def from_frame(frame: Frame2D, focal: float) -> Parabola2D: ...
    def frame(self) -> Frame2D: ...
    def apex(self) -> Point2D: ...
    def focal(self) -> float: ...
    def __repr__(self) -> str: ...

@final
class Hyperbola2D(_Curve2DMethods):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        center: Point2D,
        x_direction: Vector2D,
        major_radius: float,
        minor_radius: float,
    ) -> Hyperbola2D: ...
    @staticmethod
    def from_frame(
        frame: Frame2D, major_radius: float, minor_radius: float
    ) -> Hyperbola2D: ...
    @staticmethod
    def from_center_and_points(
        center: Point2D, s1: Point2D, s2: Point2D
    ) -> Hyperbola2D: ...
    def frame(self) -> Frame2D: ...
    def center(self) -> Point2D: ...
    def major_radius(self) -> float: ...
    def minor_radius(self) -> float: ...
    def __repr__(self) -> str: ...

@final
class TrimmedCurve(_CurveMethods, _DiscretizeMethods[Discretization]):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, basis: Any, first: float, last: float) -> TrimmedCurve: ...
    @staticmethod
    def segment(start: Point, end: Point) -> TrimmedCurve: ...
    def basis(self) -> Any: ...
    def first(self) -> float: ...
    def last(self) -> float: ...
    def to_bspline(self) -> BSplineCurve: ...
    def start_point(self) -> Point: ...
    def end_point(self) -> Point: ...
    def __repr__(self) -> str: ...

@final
class TrimmedCurve2D(_Curve2DMethods, _DiscretizeMethods[Discretization2D]):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, basis: Any, first: float, last: float) -> TrimmedCurve2D: ...
    @staticmethod
    def segment(start: Point2D, end: Point2D) -> TrimmedCurve2D: ...
    def basis(self) -> Any: ...
    def to_bspline(self) -> BSplineCurve2D: ...
    def first(self) -> float: ...
    def last(self) -> float: ...
    def start_point(self) -> Point2D: ...
    def end_point(self) -> Point2D: ...
    def __repr__(self) -> str: ...

__all__ = [
    "BSplineCurve2D",
    "BSplineCurve",
    "BezierCurve2D",
    "BezierCurve",
    "Circle2D",
    "Circle",
    "Ellipse2D",
    "Ellipse",
    "Hyperbola2D",
    "Hyperbola",
    "Line2D",
    "Line",
    "Parabola2D",
    "Parabola",
    "TrimmedCurve2D",
    "TrimmedCurve",
]
