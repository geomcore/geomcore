import numpy
import numpy.typing

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
    _Parametrization,
    Axis2D,
    Axis,
    BoundingBox,
    BoundingBox2D,
    CurveApproximation,
    CurveApproximation2D,
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
_WholeCurve2D = Union[
    Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D, BSplineCurve2D, BezierCurve2D
]
_Curve2D = Union[_WholeCurve2D, TrimmedCurve2D]

_FloatArray = numpy.typing.NDArray[numpy.float64]

_PointType = TypeVar("_PointType")
_VectorType = TypeVar("_VectorType")
_BoundingBoxType = TypeVar("_BoundingBoxType")
_ProjectionType = TypeVar("_ProjectionType")
_DiscretizationType = TypeVar("_DiscretizationType")
_Self = TypeVar("_Self")

# The method sets of geomcore-py/src/method_sets.rs, one stub-only class
# each. A curve class inherits the sets its `shape_methods!` line names, so
# each shared method is typed and documented here once; the type arguments
# carry what differs by class or dimension.

class _CurveMethodsInEitherDimension(
    Generic[_PointType, _VectorType, _BoundingBoxType, _ProjectionType]
):
    def point_at(self, u: float) -> _PointType:
        """Evaluate the point at parameter `u`.

        Every `u` evaluates. Past the bounds, a line, conic or
        periodic curve continues along its own equation, and a
        trimmed curve evaluates its basis curve. A non-periodic
        B-spline or Bézier curve clamps `u` to the nearer bound, so
        the point and every derivative are those at the bound. An
        operation that builds from a parameter, such as a B-spline's
        `split_at`, raises on one outside the bounds instead.
        """
    def points_at(self, us: List[float]) -> List[_PointType]:
        """Evaluate many parameters at once."""
    def points_at_array(self, us: numpy.typing.ArrayLike) -> _FloatArray:
        """The points at each parameter of the array-like `us`, as a
        numpy array with one row of coordinates per parameter.
        """
    def derivative(self, u: float, order: int) -> _VectorType:
        """Evaluate the derivative of the given order (>= 1) at `u`."""
    def first_derivatives(self, u: float) -> Tuple[_PointType, _VectorType]:
        """The point and the first derivative at parameter `u`."""
    def second_derivatives(self, u: float) -> Tuple[_PointType, _VectorType, _VectorType]:
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

        A distance within `tol.confusion` (`1e-7` by default) snaps
        to `0.0`, so projection agrees with `contains`; pass a
        smaller `Tolerance` to read a smaller distance.
        """
    def project_points(
        self, points: List[_PointType], tol: Optional[Tolerance] = ...
    ) -> List[_ProjectionType]:
        """Projects each point in `points` onto the curve, returning a
        projection object per point.
        """
    def project_points_array(
        self, points: numpy.typing.ArrayLike, tol: Optional[Tolerance] = ...
    ) -> Tuple[_FloatArray, _FloatArray, _FloatArray]:
        """Projects each row of the array-like `points` onto the curve,
        returning numpy arrays `(parameters, points, distances)`:
        the parameter, the closest point and its distance per row.
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
        `reversed().point_at(reversed_parameter(u))` is
        `point_at(u)`.
        """
    def reversed_parameter(self, u: float) -> float:
        """Maps `u` onto the `reversed()` curve: the parameter of the
        same point there. It is `first + last - u` over bounded
        parameters (`2π - u` on a circle or ellipse), `-u` on a
        line, parabola or hyperbola, and a trimmed curve maps as
        its basis.
        """
    def bounding_box(self) -> Optional[_BoundingBoxType]:
        """The axis-aligned box around the shape (`None` when unbounded).
        It touches an elementary shape; for a B-spline it is the box of the
        control points, which may be wider than the shape, as
        `tight_bounding_box` is not."""
    def tight_bounding_box(
        self, tol: Optional[Tolerance] = ...
    ) -> Optional[_BoundingBoxType]:
        """The axis-aligned box around the shape as tightly as
        `tol.confusion` (`None` when unbounded): it holds the shape, and each
        of its sides is within `tol.confusion` of a point of it. The same as
        `bounding_box` on an elementary shape; a B-spline's comes from halving
        its Bézier pieces, a revolution surface's from its profile's.
        `tol` is `Tolerance.default()` when omitted."""

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
    def first_derivatives_sided(self, u: float, side: Side) -> Tuple[_PointType, _VectorType]:
        """The point and the first derivative at parameter `u`, from
        one `side`: `Side.BELOW` takes the span ending at `u`,
        `Side.ABOVE` the span starting there. Both sides agree
        strictly inside a span, so everywhere on a Bézier curve.
        Plain evaluation takes the above side, except at the upper
        bound of a non-periodic curve. A periodic curve's seam is a
        knot like any other: at either bound, `BELOW` takes the last
        span of the period and `ABOVE` the first.
        """
    def second_derivatives_sided(
        self, u: float, side: Side
    ) -> Tuple[_PointType, _VectorType, _VectorType]:
        """The point and the first two derivatives at parameter `u`,
        from one `side`, as `first_derivatives_sided`.
        """
    def derivative_sided(self, u: float, order: int, side: Side) -> _VectorType:
        """The derivative of the given order (>= 1) at parameter `u`,
        from one `side`, as `first_derivatives_sided`.
        """

@final
class Line(_CurveMethods):
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
    def length(self) -> float: ...
    def enclosed_area(self) -> float: ...
    def normal(self) -> Vector: ...
    def to_bspline(self) -> BSplineCurve: ...
    def frame(self) -> Frame: ...

@final
class Ellipse(
    _CurveMethods,
    _DiscretizeMethods[Discretization],
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
    def enclosed_area(self) -> float: ...
    def length(self, max_error: float) -> float: ...
    def eccentricity(self) -> float: ...
    def foci(self) -> Tuple[Point, Point]: ...
    def directrices(self) -> Optional[Tuple[Line, Line]]: ...
    def to_bspline(self) -> BSplineCurve: ...
    def frame(self) -> Frame: ...

@final
class Parabola(_CurveMethods):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, apex: Point, normal: Vector, x_direction: Vector, focal_distance: float
    ) -> Parabola: ...
    @staticmethod
    def from_frame(frame: Frame, focal_distance: float) -> Parabola: ...
    @staticmethod
    def from_focus_and_directrix(focus: Point, directrix: Line) -> Parabola: ...
    def focal_distance(self) -> float: ...
    def eccentricity(self) -> float: ...
    def focus(self) -> Point: ...
    def directrix(self) -> Line: ...
    def frame(self) -> Frame: ...
    def apex(self) -> Point: ...

@final
class Hyperbola(_CurveMethods):
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
    def eccentricity(self) -> float: ...
    def foci(self) -> Tuple[Point, Point]: ...
    def directrices(self) -> Tuple[Line, Line]: ...
    def other_branch(self) -> Hyperbola: ...
    def frame(self) -> Frame: ...

@final
class BSplineCurve(
    _CurveMethods,
    _DiscretizeMethods[Discretization],
    _SidedEvaluationMethods[Point, Vector],
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
        points: Union[List[Point], numpy.typing.ArrayLike],
        degree: int,
        parametrization: _Parametrization = ...,
    ) -> BSplineCurve: ...
    @staticmethod
    def approximate(
        points: Union[List[Point], numpy.typing.ArrayLike],
        degree: int,
        num_poles: int,
        parametrization: _Parametrization = ...,
    ) -> BSplineCurve: ...
    @staticmethod
    def approximate_within(
        points: Union[List[Point], numpy.typing.ArrayLike],
        degree: int,
        max_deviation: float,
        parametrization: _Parametrization = ...,
    ) -> CurveApproximation: ...
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
        points: Union[List[Point2D], numpy.typing.ArrayLike],
        degree: int,
        parametrization: _Parametrization = ...,
    ) -> BSplineCurve2D: ...
    @staticmethod
    def approximate(
        points: Union[List[Point2D], numpy.typing.ArrayLike],
        degree: int,
        num_poles: int,
        parametrization: _Parametrization = ...,
    ) -> BSplineCurve2D: ...
    @staticmethod
    def approximate_within(
        points: Union[List[Point2D], numpy.typing.ArrayLike],
        degree: int,
        max_deviation: float,
        parametrization: _Parametrization = ...,
    ) -> CurveApproximation2D: ...
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
    def length(self) -> float: ...
    def is_counterclockwise(self) -> bool: ...
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
    def enclosed_area(self) -> float: ...
    def length(self, max_error: float) -> float: ...
    def eccentricity(self) -> float: ...
    def foci(self) -> Tuple[Point2D, Point2D]: ...
    def directrices(self) -> Optional[Tuple[Line2D, Line2D]]: ...
    def to_bspline(self) -> BSplineCurve2D: ...
    def __repr__(self) -> str: ...

@final
class Parabola2D(_Curve2DMethods):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, apex: Point2D, x_direction: Vector2D, focal_distance: float) -> Parabola2D: ...
    @staticmethod
    def from_frame(frame: Frame2D, focal_distance: float) -> Parabola2D: ...
    @staticmethod
    def from_focus_and_directrix(focus: Point2D, directrix: Line2D) -> Parabola2D: ...
    def frame(self) -> Frame2D: ...
    def apex(self) -> Point2D: ...
    def focal_distance(self) -> float: ...
    def eccentricity(self) -> float: ...
    def focus(self) -> Point2D: ...
    def directrix(self) -> Line2D: ...
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
    def eccentricity(self) -> float: ...
    def foci(self) -> Tuple[Point2D, Point2D]: ...
    def directrices(self) -> Tuple[Line2D, Line2D]: ...
    def other_branch(self) -> Hyperbola2D: ...
    def __repr__(self) -> str: ...

@final
class TrimmedCurve(_CurveMethods, _DiscretizeMethods[Discretization]):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, basis: _WholeCurve, first: float, last: float) -> TrimmedCurve: ...
    @staticmethod
    def segment(start: Point, end: Point) -> TrimmedCurve: ...
    @staticmethod
    def arc_through(start: Point, middle: Point, end: Point) -> TrimmedCurve: ...
    @staticmethod
    def arc(
        center: Point,
        start: Point,
        end: Point,
        normal: Vector,
        tol: Optional[Tolerance] = ...,
    ) -> TrimmedCurve: ...
    def sweep_angle(self) -> float: ...
    def length(self, max_error: float) -> float: ...
    def basis(self) -> _WholeCurve: ...
    def first(self) -> float: ...
    def last(self) -> float: ...
    def to_bspline(self) -> BSplineCurve: ...
    def start_point(self) -> Point: ...
    def end_point(self) -> Point: ...
    def __repr__(self) -> str: ...

@final
class TrimmedCurve2D(_Curve2DMethods, _DiscretizeMethods[Discretization2D]):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, basis: _WholeCurve2D, first: float, last: float) -> TrimmedCurve2D: ...
    @staticmethod
    def segment(start: Point2D, end: Point2D) -> TrimmedCurve2D: ...
    @staticmethod
    def arc_through(start: Point2D, middle: Point2D, end: Point2D) -> TrimmedCurve2D: ...
    @staticmethod
    def arc(
        center: Point2D,
        start: Point2D,
        end: Point2D,
        counterclockwise: bool = ...,
        tol: Optional[Tolerance] = ...,
    ) -> TrimmedCurve2D: ...
    def sweep_angle(self) -> float: ...
    def length(self, max_error: float) -> float: ...
    def basis(self) -> _WholeCurve2D: ...
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
