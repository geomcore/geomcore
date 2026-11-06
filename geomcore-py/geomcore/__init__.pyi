import enum
from typing import (
    Any,
    ClassVar,
    Generic,
    Iterator,
    List,
    Literal,
    Optional,
    Sequence,
    Tuple,
    TypeVar,
    Union,
    final,
    overload,
)

import numpy
import numpy.typing
from typing_extensions import deprecated

_FloatArray = numpy.typing.NDArray[numpy.float64]
_Self = TypeVar("_Self")

__version__: str

class GeomcoreError(ValueError):
    """Base class of every error raised by geomcore. Subclasses `ValueError`."""

class ConstructionError(GeomcoreError):
    """A geometric object could not be built from the given arguments (for example a negative radius or a zero-length direction)."""

class FrameConstructionError(ConstructionError):
    """A frame could not be built, for example from parallel or zero directions."""

class TrimmedCurveConstructionError(ConstructionError):
    """A trimmed curve could not be built, for example from an invalid parameter range."""

class KnotInsertionError(GeomcoreError):
    """A knot could not be inserted into a B-spline."""

class KnotRemovalError(GeomcoreError):
    """A knot could not be removed from a B-spline."""

class DegreeElevationError(GeomcoreError):
    """The degree of a B-spline could not be elevated."""

class BSplineSplitError(GeomcoreError):
    """A B-spline could not be split at the given parameter."""

class BSplineConversionError(GeomcoreError):
    """A curve or surface patch could not be converted to a B-spline."""

class Curve3DToUvError(GeomcoreError):
    """A 3D curve has no exact UV curve, a curve in a surface's parameter space, in the same parameter."""

class Curve3DToUvApproximationError(GeomcoreError):
    """A 3D curve has no UV curve, a curve in a surface's parameter space, within the asked deviation."""

class UvCurveTo3DError(GeomcoreError):
    """A UV curve, a curve in a surface's parameter space, has no exact 3D curve in the same parameter."""

class UvCurveTo3DApproximationError(GeomcoreError):
    """A UV curve, a curve in a surface's parameter space, has no 3D curve within the asked deviation."""

class SurfaceIsoError(GeomcoreError):
    """An iso-parametric curve of a surface could not be extracted."""

class DiscretizationCriteriaError(GeomcoreError):
    """Discretization criteria were invalid, for example a non-positive deviation."""

@final
class Point:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, x: float, y: float, z: float) -> Point: ...
    @staticmethod
    def origin() -> Point: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @property
    def z(self) -> float: ...
    def distance(self, other: Point) -> float: ...
    def coincides(self, other: Point, tol: Optional[_ToleranceLike] = ...) -> bool: ...
    def components(self) -> Tuple[float, float, float]: ...
    def midpoint(self, other: Point) -> Point: ...
    def transformed(self, transform: Transform) -> Point: ...
    def __repr__(self) -> str: ...
    def is_finite(self) -> bool: ...
    def __add__(self, offset: Vector, /) -> Point: ...
    @overload
    def __sub__(self, other: Point, /) -> Vector: ...
    @overload
    def __sub__(self, other: Vector, /) -> Point: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> float: ...
    def __iter__(self) -> Iterator[float]: ...
    def __array__(
        self, dtype: Optional[numpy.typing.DTypeLike] = None, copy: Optional[bool] = None
    ) -> numpy.typing.NDArray[Any]: ...

@final
class Point2D:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, x: float, y: float) -> Point2D: ...
    @staticmethod
    def origin() -> Point2D: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    def distance(self, other: Point2D) -> float: ...
    def coincides(self, other: Point2D, tol: Optional[_ToleranceLike] = ...) -> bool: ...
    def components(self) -> Tuple[float, float]: ...
    def midpoint(self, other: Point2D) -> Point2D: ...
    def is_finite(self) -> bool: ...
    def __add__(self, offset: Vector2D, /) -> Point2D: ...
    @overload
    def __sub__(self, other: Point2D, /) -> Vector2D: ...
    @overload
    def __sub__(self, other: Vector2D, /) -> Point2D: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> float: ...
    def __iter__(self) -> Iterator[float]: ...
    def __array__(
        self, dtype: Optional[numpy.typing.DTypeLike] = None, copy: Optional[bool] = None
    ) -> numpy.typing.NDArray[Any]: ...

@final
class Vector:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, x: float, y: float, z: float) -> Vector: ...
    @staticmethod
    def zero() -> Vector: ...
    @staticmethod
    def unit_x() -> Vector: ...
    @staticmethod
    def unit_y() -> Vector: ...
    @staticmethod
    def unit_z() -> Vector: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @property
    def z(self) -> float: ...
    def components(self) -> Tuple[float, float, float]: ...
    def dot(self, other: Vector) -> float: ...
    def cross(self, other: Vector) -> Vector: ...
    def magnitude(self) -> float: ...
    def square_magnitude(self) -> float: ...
    def is_finite(self) -> bool: ...
    def normalized(self) -> Optional[Vector]: ...
    def angle_with_reference(self, other: Vector, reference: Vector) -> float: ...
    def angle(self, other: Vector) -> float: ...
    def is_parallel(self, other: Vector, angular: Optional[float] = ...) -> bool: ...
    def is_perpendicular(self, other: Vector, angular: Optional[float] = ...) -> bool: ...
    def component_along(self, direction: Vector) -> Vector: ...
    def component_perpendicular_to(self, direction: Vector) -> Vector: ...
    def __add__(self, other: Vector, /) -> Vector: ...
    def __sub__(self, other: Vector, /) -> Vector: ...
    def __neg__(self) -> Vector: ...
    def __mul__(self, scalar: float, /) -> Vector: ...
    def __rmul__(self, scalar: float, /) -> Vector: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> float: ...
    def __iter__(self) -> Iterator[float]: ...
    def __array__(
        self, dtype: Optional[numpy.typing.DTypeLike] = None, copy: Optional[bool] = None
    ) -> numpy.typing.NDArray[Any]: ...
    def __truediv__(self, scalar: float, /) -> Vector: ...
    def __abs__(self) -> float: ...
    def __repr__(self) -> str: ...

@final
class Vector2D:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, x: float, y: float) -> Vector2D: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @staticmethod
    def zero() -> Vector2D: ...
    @staticmethod
    def unit_x() -> Vector2D: ...
    @staticmethod
    def unit_y() -> Vector2D: ...
    def components(self) -> Tuple[float, float]: ...
    def dot(self, other: Vector2D) -> float: ...
    def cross(self, other: Vector2D) -> float: ...
    def perpendicular(self) -> Vector2D: ...
    def angle(self, other: Vector2D) -> float: ...
    def signed_angle(self, other: Vector2D) -> float: ...
    def rotated(self, angle: float) -> Vector2D: ...
    def is_parallel(self, other: Vector2D, angular: Optional[float] = ...) -> bool: ...
    def is_perpendicular(self, other: Vector2D, angular: Optional[float] = ...) -> bool: ...
    def component_along(self, direction: Vector2D) -> Vector2D: ...
    def component_perpendicular_to(self, direction: Vector2D) -> Vector2D: ...
    def magnitude(self) -> float: ...
    def square_magnitude(self) -> float: ...
    def is_finite(self) -> bool: ...
    def normalized(self) -> Optional[Vector2D]: ...
    def __add__(self, other: Vector2D, /) -> Vector2D: ...
    def __sub__(self, other: Vector2D, /) -> Vector2D: ...
    def __neg__(self) -> Vector2D: ...
    def __mul__(self, scalar: float, /) -> Vector2D: ...
    def __rmul__(self, scalar: float, /) -> Vector2D: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> float: ...
    def __iter__(self) -> Iterator[float]: ...
    def __array__(
        self, dtype: Optional[numpy.typing.DTypeLike] = None, copy: Optional[bool] = None
    ) -> numpy.typing.NDArray[Any]: ...
    def __truediv__(self, scalar: float, /) -> Vector2D: ...
    def __abs__(self) -> float: ...
    def __repr__(self) -> str: ...

@final
class Axis2D:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, origin: Point2D, direction: Vector2D) -> Axis2D: ...
    def origin(self) -> Point2D: ...
    def direction(self) -> Vector2D: ...
    def point_at(self, t: float) -> Point2D: ...
    def parameter_of(self, point: Point2D) -> float: ...
    def closest_point(self, point: Point2D) -> Point2D: ...
    def distance(self, point: Point2D) -> float: ...
    def contains(self, point: Point2D, tol: Optional[_ToleranceLike] = ...) -> bool: ...
    def is_parallel(self, other: Axis2D, angular: Optional[float] = ...) -> bool: ...
    def reversed(self) -> Axis2D: ...
    def __repr__(self) -> str: ...

@final
class Frame2D:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, origin: Point2D, x_direction: Vector2D, y_direction: Vector2D
    ) -> Frame2D: ...
    @staticmethod
    def from_x(origin: Point2D, x_direction: Vector2D) -> Frame2D: ...
    @staticmethod
    def world() -> Frame2D: ...
    def origin(self) -> Point2D: ...
    def x_direction(self) -> Vector2D: ...
    def y_direction(self) -> Vector2D: ...
    def is_direct(self) -> bool: ...
    def x_axis(self) -> Axis2D: ...
    def y_axis(self) -> Axis2D: ...
    def local_coordinates(self, point: Point2D) -> Tuple[float, float]: ...
    def point_at(self, u: float, v: float) -> Point2D: ...
    def vector_at(self, u: float, v: float) -> Vector2D: ...
    def local_components(self, vector: Vector2D) -> Tuple[float, float]: ...
    def __repr__(self) -> str: ...

@final
class Axis:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, origin: Point, direction: Vector) -> Axis: ...
    def origin(self) -> Point: ...
    def direction(self) -> Vector: ...
    def point_at(self, t: float) -> Point: ...
    def parameter_of(self, point: Point) -> float: ...
    def closest_point(self, point: Point) -> Point: ...
    def distance(self, point: Point) -> float: ...
    def contains(self, point: Point, tol: Optional[_ToleranceLike] = ...) -> bool: ...
    def is_parallel(self, other: Axis, angular: Optional[float] = ...) -> bool: ...
    def is_coaxial(self, other: Axis, tol: Optional[_ToleranceLike] = ...) -> bool: ...
    def reversed(self) -> Axis: ...

@final
class Frame:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, origin: Point, z_direction: Vector, x_hint: Vector) -> Frame: ...
    @staticmethod
    def from_z(origin: Point, z_direction: Vector) -> Frame: ...
    @staticmethod
    def world() -> Frame: ...
    @staticmethod
    def from_axes(
        origin: Point, x_direction: Vector, y_direction: Vector, z_direction: Vector
    ) -> Frame: ...
    def is_direct(self) -> bool: ...
    def flipped_x(self) -> Frame: ...
    def flipped_y(self) -> Frame: ...
    def flipped_z(self) -> Frame: ...
    def origin(self) -> Point: ...
    def x_direction(self) -> Vector: ...
    def y_direction(self) -> Vector: ...
    def z_direction(self) -> Vector: ...
    def axis(self) -> Axis: ...
    def local_coordinates(self, point: Point) -> Tuple[float, float, float]: ...
    def point_at(self, u: float, v: float, w: float) -> Point: ...
    def vector_at(self, u: float, v: float, w: float) -> Vector: ...
    def local_components(self, vector: Vector) -> Tuple[float, float, float]: ...
    def transformed(self, transform: Transform) -> Frame: ...

@final
class Tolerance:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, confusion: float = ..., angular: float = ..., parametric: float = ...
    ) -> Tolerance: ...
    @staticmethod
    def default() -> Tolerance: ...
    @property
    def confusion(self) -> float: ...
    @property
    def angular(self) -> float: ...
    @property
    def parametric(self) -> float: ...
    MODEL_SIZE: ClassVar[float]
    def angular_for_extent(self, extent: float) -> float: ...
    def __repr__(self) -> str: ...

# A `tol` argument: a `Tolerance`, or a number taken as `Tolerance(confusion=number)`.
_ToleranceLike = Union[Tolerance, float]

@final
class Transform:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @staticmethod
    def identity() -> Transform: ...
    @staticmethod
    def translation(offset: Vector) -> Transform: ...
    @staticmethod
    def rotation(axis: Axis, angle: float) -> Transform: ...
    @staticmethod
    def scaling(center: Point, factor: float) -> Transform: ...
    @staticmethod
    def mirror_point(center: Point) -> Transform: ...
    @staticmethod
    def mirror_axis(axis: Axis) -> Transform: ...
    @staticmethod
    def mirror_plane(frame: Frame) -> Transform: ...
    @staticmethod
    def between_frames(source: Frame, target: Frame) -> Transform: ...
    def scale_factor(self) -> float: ...
    def preserves_orientation(self) -> bool: ...
    def is_identity(self, tol: Optional[_ToleranceLike] = ...) -> bool: ...
    def then(self, other: Transform) -> Transform: ...
    def apply_point(self, p: Point) -> Point: ...
    def apply_vector(self, v: Vector) -> Vector: ...
    @staticmethod
    def from_matrix(matrix: Sequence[Sequence[float]]) -> Transform: ...
    def to_matrix(self) -> List[List[float]]: ...
    def inverse(self) -> Transform: ...
    def __matmul__(self, other: Transform, /) -> Transform: ...
    def __repr__(self) -> str: ...

@final
class BoundingBox:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, a: Point, b: Point) -> BoundingBox: ...
    @staticmethod
    def from_points(points: List[Point]) -> BoundingBox: ...
    def min_corner(self) -> Point: ...
    def max_corner(self) -> Point: ...
    def center(self) -> Point: ...
    def diagonal(self) -> Vector: ...
    def union(self, other: BoundingBox) -> BoundingBox: ...
    def enlarged(self, distance: float) -> BoundingBox: ...
    def encloses(self, point: Point) -> bool: ...
    def __contains__(self, point: Point, /) -> bool: ...
    def including(self, point: Point) -> BoundingBox: ...
    def intersection(self, other: BoundingBox) -> Optional[BoundingBox]: ...
    def distance(self, point: Point) -> float: ...
    def overlaps(self, other: BoundingBox) -> bool: ...
    def __repr__(self) -> str: ...

@final
class BoundingBox2D:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, a: Point2D, b: Point2D) -> BoundingBox2D: ...
    @staticmethod
    def from_points(points: List[Point2D]) -> BoundingBox2D: ...
    def min_corner(self) -> Point2D: ...
    def max_corner(self) -> Point2D: ...
    def center(self) -> Point2D: ...
    def diagonal(self) -> Vector2D: ...
    def union(self, other: BoundingBox2D) -> BoundingBox2D: ...
    def enlarged(self, distance: float) -> BoundingBox2D: ...
    def encloses(self, point: Point2D) -> bool: ...
    def __contains__(self, point: Point2D, /) -> bool: ...
    def including(self, point: Point2D) -> BoundingBox2D: ...
    def intersection(self, other: BoundingBox2D) -> Optional[BoundingBox2D]: ...
    def distance(self, point: Point2D) -> float: ...
    def overlaps(self, other: BoundingBox2D) -> bool: ...
    def __repr__(self) -> str: ...

@final
class Discretization2D:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[Point2D]: ...
    @overload
    def __getitem__(self, index: int, /) -> Point2D: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[Point2D]: ...
    @property
    def parameters(self) -> List[float]: ...
    @property
    def points(self) -> List[Point2D]: ...
    @property
    def chord_deviations(self) -> List[float]: ...
    @property
    def angle_deviations(self) -> List[float]: ...
    def __repr__(self) -> str: ...

@final
class Discretization:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[Point]: ...
    @overload
    def __getitem__(self, index: int, /) -> Point: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[Point]: ...
    @property
    def parameters(self) -> List[float]: ...
    @property
    def points(self) -> List[Point]: ...
    @property
    def chord_deviations(self) -> List[float]: ...
    @property
    def angle_deviations(self) -> List[float]: ...
    def __repr__(self) -> str: ...


class Side(enum.IntEnum):
    BELOW = 0
    ABOVE = 1

class PointParametrization(str, enum.Enum):
    CHORDAL = "chordal"
    CENTRIPETAL = "centripetal"

_Parametrization = Union[PointParametrization, Literal["chordal", "centripetal"]]

@final
class SurfaceCurvatures:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def min_curvature(self) -> float: ...
    @property
    def max_curvature(self) -> float: ...
    @property
    def min_direction(self) -> Vector: ...
    @property
    def max_direction(self) -> Vector: ...
    @property
    def mean(self) -> float: ...
    @property
    def gaussian(self) -> float: ...

@final
class CurveApproximation:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve(self) -> BSplineCurve: ...
    @property
    def deviation(self) -> float: ...
    def is_proven(self) -> bool: ...

@final
class CurveApproximation2D:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve(self) -> BSplineCurve2D: ...
    @property
    def deviation(self) -> float: ...
    def is_proven(self) -> bool: ...

@final
class SurfaceApproximation:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def surface(self) -> BSplineSurface: ...
    @property
    def deviation(self) -> float: ...
    def is_proven(self) -> bool: ...

class IntersectionMethod(enum.IntEnum):
    EXACT = 0
    NUMERIC = 1

class CompletenessStatus(enum.IntEnum):
    COMPLETE = 0
    INCOMPLETE = 1
    SAMPLED = 2

class IncompleteReason(enum.IntEnum):
    BUDGET_EXHAUSTED = 0
    UNSEPARATED = 1
    UNBOUNDED = 2

@final
class CurveCurveRegion:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_range(self) -> Tuple[float, float]: ...
    @property
    def second_range(self) -> Tuple[float, float]: ...
    @property
    def min_corner(self) -> Point: ...
    @property
    def max_corner(self) -> Point: ...

@final
class CurveSurfaceRegion:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve_range(self) -> Tuple[float, float]: ...
    @property
    def u_range(self) -> Tuple[float, float]: ...
    @property
    def v_range(self) -> Tuple[float, float]: ...
    @property
    def min_corner(self) -> Point: ...
    @property
    def max_corner(self) -> Point: ...

@final
class SurfaceSurfaceRegion:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_u_range(self) -> Tuple[float, float]: ...
    @property
    def first_v_range(self) -> Tuple[float, float]: ...
    @property
    def second_u_range(self) -> Tuple[float, float]: ...
    @property
    def second_v_range(self) -> Tuple[float, float]: ...
    @property
    def min_corner(self) -> Point: ...
    @property
    def max_corner(self) -> Point: ...

@final
class Completeness:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def status(self) -> CompletenessStatus: ...
    @property
    def reason(self) -> Optional[IncompleteReason]: ...
    @property
    def regions(
        self,
    ) -> Union[List[CurveCurveRegion], List[CurveSurfaceRegion], List[SurfaceSurfaceRegion]]: ...
    def is_complete(self) -> bool: ...

@final
class CurveCurveHit:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_parameter(self) -> float: ...
    @property
    def second_parameter(self) -> float: ...
    @property
    def point(self) -> Point: ...
    @property
    def multiplicity(self) -> int: ...

@final
class CurveCurveHit2D:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_parameter(self) -> float: ...
    @property
    def second_parameter(self) -> float: ...
    @property
    def point(self) -> Point2D: ...
    @property
    def multiplicity(self) -> int: ...

@final
class CurveSurfaceHit:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve_parameter(self) -> float: ...
    @property
    def u(self) -> float: ...
    @property
    def v(self) -> float: ...
    @property
    def point(self) -> Point: ...
    @property
    def multiplicity(self) -> int: ...

@final
class CurveCurveOverlap:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_range(self) -> Tuple[float, float]: ...
    @property
    def second_range(self) -> Tuple[float, float]: ...
    @property
    def same_direction(self) -> bool: ...

@final
class CurveSurfaceOverlap:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve_range(self) -> Tuple[float, float]: ...

@final
class IntersectionCurve:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve(self) -> _WholeCurve: ...
    @property
    def first_uv_curve(self) -> Optional[_WholeCurve2D]: ...
    @property
    def second_uv_curve(self) -> Optional[_WholeCurve2D]: ...
    @property
    def tangent(self) -> bool: ...
    @property
    def deviation(self) -> float: ...

@final
class CurveCurveIntersection:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[CurveCurveHit]: ...
    @overload
    def __getitem__(self, index: int, /) -> CurveCurveHit: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[CurveCurveHit]: ...
    @property
    def points(self) -> List[Point]: ...
    @property
    def hits(self) -> List[CurveCurveHit]: ...
    @property
    def overlaps(self) -> List[CurveCurveOverlap]: ...
    @property
    def method(self) -> IntersectionMethod: ...
    @property
    def completeness(self) -> Completeness: ...
    @property
    def first(self) -> Curve: ...
    @property
    def second(self) -> Curve: ...
    def is_complete(self) -> bool: ...
    def is_empty(self) -> bool: ...
    def __bool__(self) -> bool: ...

@final
class CurveCurveIntersection2D:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[CurveCurveHit2D]: ...
    @overload
    def __getitem__(self, index: int, /) -> CurveCurveHit2D: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[CurveCurveHit2D]: ...
    @property
    def points(self) -> List[Point2D]: ...
    @property
    def hits(self) -> List[CurveCurveHit2D]: ...
    @property
    def overlaps(self) -> List[CurveCurveOverlap]: ...
    @property
    def method(self) -> IntersectionMethod: ...
    @property
    def completeness(self) -> Completeness: ...
    @property
    def first(self) -> Curve2D: ...
    @property
    def second(self) -> Curve2D: ...
    def is_complete(self) -> bool: ...
    def is_empty(self) -> bool: ...
    def __bool__(self) -> bool: ...

@final
class CurveSurfaceIntersection:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[CurveSurfaceHit]: ...
    @overload
    def __getitem__(self, index: int, /) -> CurveSurfaceHit: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[CurveSurfaceHit]: ...
    @property
    def points(self) -> List[Point]: ...
    @property
    def hits(self) -> List[CurveSurfaceHit]: ...
    @property
    def overlaps(self) -> List[CurveSurfaceOverlap]: ...
    @property
    def method(self) -> IntersectionMethod: ...
    @property
    def completeness(self) -> Completeness: ...
    @property
    def first(self) -> Union[Curve, Surface]: ...
    @property
    def second(self) -> Union[Curve, Surface]: ...
    @property
    def curve(self) -> Curve: ...
    @property
    def surface(self) -> Surface: ...
    def is_complete(self) -> bool: ...
    def is_empty(self) -> bool: ...
    def __bool__(self) -> bool: ...

@final
class SurfaceSurfaceIntersection:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[IntersectionCurve]: ...
    @overload
    def __getitem__(self, index: int, /) -> IntersectionCurve: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[IntersectionCurve]: ...
    @property
    def curves(self) -> List[IntersectionCurve]: ...
    @property
    def contact_points(self) -> List[Point]: ...
    @property
    def method(self) -> IntersectionMethod: ...
    @property
    def completeness(self) -> Completeness: ...
    @property
    def first(self) -> Surface: ...
    @property
    def second(self) -> Surface: ...
    def is_coincident(self) -> bool: ...
    def is_complete(self) -> bool: ...
    def is_empty(self) -> bool: ...
    def __bool__(self) -> bool: ...

@final
class CurveProjection:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def parameter(self) -> float: ...
    @property
    def point(self) -> Point: ...
    @property
    def distance(self) -> float: ...
    @property
    def target(self) -> Curve: ...
    @property
    def query_point(self) -> Point: ...

@final
class CurveProjection2D:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def parameter(self) -> float: ...
    @property
    def point(self) -> Point2D: ...
    @property
    def distance(self) -> float: ...
    @property
    def target(self) -> Curve2D: ...
    @property
    def query_point(self) -> Point2D: ...

@final
class SurfaceProjection:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def u(self) -> float: ...
    @property
    def v(self) -> float: ...
    @property
    def point(self) -> Point: ...
    @property
    def distance(self) -> float: ...
    @property
    def target(self) -> Surface: ...
    @property
    def query_point(self) -> Point: ...

@final
class CurveCurveExtremum:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_parameter(self) -> float: ...
    @property
    def second_parameter(self) -> float: ...
    @property
    def first_point(self) -> Point: ...
    @property
    def second_point(self) -> Point: ...
    @property
    def distance(self) -> float: ...

@final
class CurveCurveExtrema:
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def closest(self) -> Optional[CurveCurveExtremum]: ...
    def __iter__(self) -> Iterator[CurveCurveExtremum]: ...
    @overload
    def __getitem__(self, index: int, /) -> CurveCurveExtremum: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[CurveCurveExtremum]: ...
    @property
    def minima(self) -> List[CurveCurveExtremum]: ...
    @property
    def completeness(self) -> Completeness: ...
    @property
    def first(self) -> Curve: ...
    @property
    def second(self) -> Curve: ...
    def is_complete(self) -> bool: ...

@deprecated("use surface.curve_3d_to_uv_within(curve, max_deviation, tol)")
def curve_3d_to_uv_within(
    surface: Surface,
    curve: Curve,
    max_deviation: float,
    tol: Optional[_ToleranceLike] = ...,
) -> CurveApproximation2D: ...
@deprecated("use surface.uv_curve_to_3d_within(uv_curve, max_deviation, tol)")
def uv_curve_to_3d_within(
    surface: Surface,
    uv_curve: Curve2D,
    max_deviation: float,
    tol: Optional[_ToleranceLike] = ...,
) -> CurveApproximation: ...
@deprecated("use first.extrema(second, tol)")
def curve_curve_extrema(
    first: Curve, second: Curve, tol: Optional[_ToleranceLike] = ...
) -> CurveCurveExtrema: ...

# ---------------------------------------------------------------- curves

# Every class of each kind of shape, as the union to annotate with; at
# runtime each name is the tuple of the classes, for `isinstance`.
# The whole curves are the classes of the `Curve` enum's variants.
_WholeCurve = Union[
    Line, Circle, Ellipse, Parabola, Hyperbola, BSplineCurve, BezierCurve
]
Curve = (
    Line
    | Circle
    | Ellipse
    | Parabola
    | Hyperbola
    | BSplineCurve
    | BezierCurve
    | TrimmedCurve
)
_WholeCurve2D = Union[
    Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D, BSplineCurve2D, BezierCurve2D
]
Curve2D = (
    Line2D
    | Circle2D
    | Ellipse2D
    | Parabola2D
    | Hyperbola2D
    | BSplineCurve2D
    | BezierCurve2D
    | TrimmedCurve2D
)

_PointType = TypeVar("_PointType")
_VectorType = TypeVar("_VectorType")
_BoundingBoxType = TypeVar("_BoundingBoxType")
_ProjectionType = TypeVar("_ProjectionType")
_DiscretizationType = TypeVar("_DiscretizationType")

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
    def tangent(self, u: float, tol: Optional[_ToleranceLike] = ...) -> Optional[_VectorType]:
        """Unit tangent at parameter `u` (`Tolerance.default()` when
        omitted), or `None` where the speed vanishes, as at a cusp.
        The speed counts as zero when the distance it would cover
        across the parameter range is within `tol.confusion`, so
        rescaling the parameter does not change the answer.
        """
    def curvature(self, u: float, tol: Optional[_ToleranceLike] = ...) -> Optional[float]:
        """Curvature at parameter `u`: the turn rate of the unit tangent
        per unit arc length, zero on a straight run. `None` where the
        speed vanishes, by the test of `tangent`.
        """
    def curvature_center(
        self, u: float, tol: Optional[_ToleranceLike] = ...
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
    def is_closed(self, tol: Optional[_ToleranceLike] = ...) -> bool:
        """Whether travelling along the curve returns to its start
        (`Tolerance.default()` when omitted). A periodic curve is
        closed; a clamped B-spline whose ends coincide is closed
        without being periodic.
        """
    def continuity(self, u: float, tol: Optional[_ToleranceLike] = ...) -> Union[int, float]:
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
    def project_point(
        self, point: _PointType, tol: Optional[_ToleranceLike] = ...
    ) -> _ProjectionType:
        """Projects `point` onto the curve: the closest point, as a
        projection object.

        A distance within `tol.confusion` (`1e-7` by default) snaps
        to `0.0`, so projection agrees with `contains`; pass a
        smaller `Tolerance` to read a smaller distance.
        """
    def project_points(
        self, points: List[_PointType], tol: Optional[_ToleranceLike] = ...
    ) -> List[_ProjectionType]:
        """Projects each point in `points` onto the curve, returning a
        projection object per point.
        """
    def project_points_array(
        self, points: numpy.typing.ArrayLike, tol: Optional[_ToleranceLike] = ...
    ) -> Tuple[_FloatArray, _FloatArray, _FloatArray]:
        """Projects each row of the array-like `points` onto the curve,
        returning numpy arrays `(parameters, points, distances)`:
        the parameter, the closest point and its distance per row.
        """
    def contains(self, point: _PointType, tol: Optional[_ToleranceLike] = ...) -> bool:
        """Returns whether `point` lies on the curve within `tol`
        (`Tolerance.default()` when omitted).
        """
    def parameter_of(self, point: _PointType, tol: Optional[_ToleranceLike] = ...) -> float:
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
        self, tol: Optional[_ToleranceLike] = ...
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
        self, other: Curve, tol: Optional[_ToleranceLike] = ...
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
        self, other: Surface, tol: Optional[_ToleranceLike] = ...
    ) -> CurveSurfaceIntersection: ...
    @overload
    def extrema(
        self, other: Point, tol: Optional[_ToleranceLike] = ...
    ) -> List[CurveProjection]:
        """The extrema of the distance to `other`, a point or a curve.

        From a point: its stationary points on the curve (local
        minima plus the ends of a bounded curve), as
        `CurveProjection` objects ordered by ascending distance.

        From a curve or `TrimmedCurve`: a `CurveCurveExtrema` of the
        local minima of the distance between the two curves, closest
        first, with whether the search proved there is no other.
        Parallel lines, concentric circles and other continuous
        families of equal minima report the distance once, with a
        representative pair. The result holds this curve as `first`
        and `other` as `second`. Other Python threads run during the
        search. Raises `TypeError` for anything else, a 2D curve
        included.
        """
    @overload
    def extrema(
        self, other: Curve, tol: Optional[_ToleranceLike] = ...
    ) -> CurveCurveExtrema: ...

class _Curve2DMethods(
    _CurveMethodsInEitherDimension[Point2D, Vector2D, BoundingBox2D, CurveProjection2D]
):
    def intersect(
        self, other: Curve2D, tol: Optional[_ToleranceLike] = ...
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
    def extrema(
        self, point: Point2D, tol: Optional[_ToleranceLike] = ...
    ) -> List[CurveProjection2D]:
        """The stationary points of the distance from `point` to the
        curve (local minima plus the ends of a bounded curve), as
        `CurveProjection2D` objects ordered by ascending distance.
        """

class _BoundedCurveMethods(Generic[_PointType, _DiscretizationType]):
    def start_point(self) -> _PointType:
        """The point at the first bound: on a circle or an ellipse, the
        seam, where the end point is too.
        """
    def end_point(self) -> _PointType:
        """The point at the last bound."""
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
    _BoundedCurveMethods[Point, Discretization],
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
    def axis(self) -> Axis: ...
    def to_bspline(self) -> BSplineCurve: ...
    def frame(self) -> Frame: ...

@final
class Ellipse(
    _CurveMethods,
    _BoundedCurveMethods[Point, Discretization],
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
    def length(self, max_length_error: float) -> float: ...
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
    _BoundedCurveMethods[Point, Discretization],
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
        pole_count: int,
        parametrization: _Parametrization = ...,
    ) -> BSplineCurve: ...
    @staticmethod
    def approximate_within(
        points: Union[List[Point], numpy.typing.ArrayLike],
        degree: int,
        max_deviation: float,
        parametrization: _Parametrization = ...,
    ) -> CurveApproximation: ...
    @staticmethod
    def from_flat_knots(
        degree: int,
        poles: Union[List[Point], numpy.typing.ArrayLike],
        flat_knots: List[float],
        weights: Optional[List[float]] = ...,
        periodic: bool = ...,
    ) -> BSplineCurve:
        """Build a curve from its flat knot vector, the form textbooks use:
        each distinct knot repeated by its multiplicity, so a non-periodic
        curve of `n` poles takes `n + degree + 1` knots. `weights`, one per
        pole, make it rational. The flat knots of a periodic curve are as
        `flat_knots()` gives them, with no padding beyond the period.
        """
    @staticmethod
    def clamped_uniform(
        degree: int, poles: Union[List[Point], numpy.typing.ArrayLike]
    ) -> BSplineCurve:
        """Build the clamped curve of `degree` with uniform knots over
        `[0, 1]`: both end knots of multiplicity `degree + 1`, and
        `len(poles) - degree - 1` simple knots between them, evenly spaced.
        The curve starts at the first pole and ends at the last.
        """
    def degree(self) -> int: ...
    def is_rational(self) -> bool: ...
    def poles(self) -> List[Point]: ...
    def weights(self) -> Optional[List[float]]: ...
    def knots(self) -> List[float]: ...
    def multiplicities(self) -> List[int]: ...
    def flat_knots(self) -> List[float]:
        """The flat knot vector: each distinct knot repeated by its
        multiplicity, `len(poles) + degree + 1` knots for a non-periodic
        curve.
        """
    def with_pole(self, index: int, pole: Point) -> BSplineCurve:
        """The curve with pole `index` (negative counting from the end) moved
        to `pole`, its degree, knots and weights kept.
        """
    def with_weight(self, index: int, weight: float) -> BSplineCurve:
        """The curve with the weight of pole `index` (negative counting from
        the end) set to `weight`, a positive finite number. A non-rational
        curve becomes rational, every other pole weighing 1.
        """
    def geometric_continuity(
        self, u: float, tol: Optional[_ToleranceLike] = ...
    ) -> Union[int, float]:
        """The geometric continuity order at parameter `u` (`math.inf` where
        no knot constrains the shape): the shape's own smoothness,
        whatever the parametrization. `continuity` reads the knot's
        multiplicity, so a knot of multiplicity `degree` reports 0 even
        where the poles around it keep the curve smooth.

        At a knot where `continuity` is below 2, the two sides are
        compared: 1 when their unit tangents agree within `tol.angular`, 2
        when their curvature vectors also agree, the osculating circles
        parting by no more than `tol.confusion` across the control
        polygon's size (the largest distance from the first pole). Where
        the speed of a side vanishes, there is no tangent to compare and
        `continuity` is returned; anywhere else too, and the geometric
        order never falls below it. `tol` is `Tolerance.default()` when
        omitted.
        """
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
    def to_bspline(self) -> BSplineCurve: ...
    def extract(self, start: float, end: float) -> BSplineCurve: ...
    def split_at(self, u: float) -> Tuple[BSplineCurve, BSplineCurve]: ...
    def length(self, max_length_error: float) -> float: ...

@final
class BSplineCurve2D(
    _Curve2DMethods,
    _BoundedCurveMethods[Point2D, Discretization2D],
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
        pole_count: int,
        parametrization: _Parametrization = ...,
    ) -> BSplineCurve2D: ...
    @staticmethod
    def approximate_within(
        points: Union[List[Point2D], numpy.typing.ArrayLike],
        degree: int,
        max_deviation: float,
        parametrization: _Parametrization = ...,
    ) -> CurveApproximation2D: ...
    @staticmethod
    def from_flat_knots(
        degree: int,
        poles: Union[List[Point2D], numpy.typing.ArrayLike],
        flat_knots: List[float],
        weights: Optional[List[float]] = ...,
        periodic: bool = ...,
    ) -> BSplineCurve2D:
        """Build a curve from its flat knot vector, the form textbooks use:
        each distinct knot repeated by its multiplicity, so a non-periodic
        curve of `n` poles takes `n + degree + 1` knots. `weights`, one per
        pole, make it rational. The flat knots of a periodic curve are as
        `flat_knots()` gives them, with no padding beyond the period.
        """
    @staticmethod
    def clamped_uniform(
        degree: int, poles: Union[List[Point2D], numpy.typing.ArrayLike]
    ) -> BSplineCurve2D:
        """Build the clamped curve of `degree` with uniform knots over
        `[0, 1]`: both end knots of multiplicity `degree + 1`, and
        `len(poles) - degree - 1` simple knots between them, evenly spaced.
        The curve starts at the first pole and ends at the last.
        """
    def degree(self) -> int: ...
    def is_rational(self) -> bool: ...
    def poles(self) -> List[Point2D]: ...
    def weights(self) -> Optional[List[float]]: ...
    def knots(self) -> List[float]: ...
    def multiplicities(self) -> List[int]: ...
    def flat_knots(self) -> List[float]:
        """The flat knot vector: each distinct knot repeated by its
        multiplicity, `len(poles) + degree + 1` knots for a non-periodic
        curve.
        """
    def with_pole(self, index: int, pole: Point2D) -> BSplineCurve2D:
        """The curve with pole `index` (negative counting from the end) moved
        to `pole`, its degree, knots and weights kept.
        """
    def with_weight(self, index: int, weight: float) -> BSplineCurve2D:
        """The curve with the weight of pole `index` (negative counting from
        the end) set to `weight`, a positive finite number. A non-rational
        curve becomes rational, every other pole weighing 1.
        """
    def geometric_continuity(
        self, u: float, tol: Optional[_ToleranceLike] = ...
    ) -> Union[int, float]:
        """The geometric continuity order at parameter `u` (`math.inf` where
        no knot constrains the shape): the shape's own smoothness,
        whatever the parametrization. `continuity` reads the knot's
        multiplicity, so a knot of multiplicity `degree` reports 0 even
        where the poles around it keep the curve smooth.

        At a knot where `continuity` is below 2, the two sides are
        compared: 1 when their unit tangents agree within `tol.angular`, 2
        when their curvature vectors also agree, the osculating circles
        parting by no more than `tol.confusion` across the control
        polygon's size (the largest distance from the first pole). Where
        the speed of a side vanishes, there is no tangent to compare and
        `continuity` is returned; anywhere else too, and the geometric
        order never falls below it. `tol` is `Tolerance.default()` when
        omitted.
        """
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
    def to_bspline(self) -> BSplineCurve2D: ...
    def extract(self, start: float, end: float) -> BSplineCurve2D: ...
    def split_at(self, u: float) -> Tuple[BSplineCurve2D, BSplineCurve2D]: ...
    def length(self, max_length_error: float) -> float: ...

@final
class BezierCurve(
    _CurveMethods,
    _BoundedCurveMethods[Point, Discretization],
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
    def length(self, max_length_error: float) -> float: ...

@final
class BezierCurve2D(
    _Curve2DMethods,
    _BoundedCurveMethods[Point2D, Discretization2D],
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
    def length(self, max_length_error: float) -> float: ...

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
class Circle2D(_Curve2DMethods, _BoundedCurveMethods[Point2D, Discretization2D]):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, center: Point2D, radius: float) -> Circle2D: ...
    @staticmethod
    def from_frame(frame: Frame2D, radius: float) -> Circle2D: ...
    def center(self) -> Point2D: ...
    def radius(self) -> float: ...
    def length(self) -> float: ...
    def is_counterclockwise(self) -> bool: ...
    def to_bspline(self) -> BSplineCurve2D: ...
    def __repr__(self) -> str: ...
    def frame(self) -> Frame2D: ...

@final
class Ellipse2D(_Curve2DMethods, _BoundedCurveMethods[Point2D, Discretization2D]):
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
    def length(self, max_length_error: float) -> float: ...
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
class TrimmedCurve(_CurveMethods, _BoundedCurveMethods[Point, Discretization]):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, basis: _WholeCurve, start_parameter: float, end_parameter: float
    ) -> TrimmedCurve: ...
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
        tol: Optional[_ToleranceLike] = ...,
    ) -> TrimmedCurve: ...
    def sweep_angle(self) -> float: ...
    def length(self, max_length_error: float) -> float: ...
    def basis(self) -> _WholeCurve: ...
    def start_parameter(self) -> float: ...
    def end_parameter(self) -> float: ...
    def to_bspline(self) -> BSplineCurve: ...
    def __repr__(self) -> str: ...

@final
class TrimmedCurve2D(_Curve2DMethods, _BoundedCurveMethods[Point2D, Discretization2D]):
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, basis: _WholeCurve2D, start_parameter: float, end_parameter: float
    ) -> TrimmedCurve2D: ...
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
        tol: Optional[_ToleranceLike] = ...,
    ) -> TrimmedCurve2D: ...
    def sweep_angle(self) -> float: ...
    def length(self, max_length_error: float) -> float: ...
    def basis(self) -> _WholeCurve2D: ...
    def to_bspline(self) -> BSplineCurve2D: ...
    def start_parameter(self) -> float: ...
    def end_parameter(self) -> float: ...
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

# -------------------------------------------------------------- surfaces

Surface = (
    Plane
    | Cylinder
    | Cone
    | Sphere
    | Torus
    | ExtrusionSurface
    | RevolutionSurface
    | BSplineSurface
    | BezierSurface
)

_UIsoCurveType = TypeVar("_UIsoCurveType")
_VIsoCurveType = TypeVar("_VIsoCurveType")

# The method sets of geomcore-py/src/method_sets.rs, one stub-only class
# each, as for the curves. The type arguments of `_SurfaceMethods` are the
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
        `reversed_u().point_at(u_reversed_parameter(u), v)` is
        `point_at(u, v)`.
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
    def point_at(self, u: float, v: float) -> Point:
        """Evaluate the point at surface parameters (u, v).

        Every `(u, v)` evaluates. Past the bounds, an elementary
        surface continues along its own equation, and a swept surface
        along its sweep, with its profile evaluated as that curve is
        past its own bounds. A non-periodic direction of
        a B-spline or Bézier surface clamps its parameter to the
        nearer bound, so the point and every partial are those at the
        bound. An operation that builds from a parameter, such as
        `u_iso`, raises on one outside the bounds instead.
        """
    def points_at(self, uvs: List[Tuple[float, float]]) -> List[Point]:
        """Evaluate many (u, v) pairs at once."""
    def points_at_array(self, uvs: numpy.typing.ArrayLike) -> _FloatArray:
        """The points at each `(u, v)` row of the `(N, 2)` array-like
        `uvs`, as an `(N, 3)` numpy array.
        """
    def derivative(self, u: float, v: float, du: int, dv: int) -> Vector:
        """Evaluate the partial derivative of order (du, dv), with
        du+dv >= 1: of every order on the plane, cylinder, cone,
        sphere and torus, and du+dv <= 2 on the other surfaces.
        """
    def first_derivatives(self, u: float, v: float) -> Tuple[Point, Tuple[Vector, Vector]]:
        """The point and the first partial derivatives `(du, dv)` at
        `(u, v)`.
        """
    def second_derivatives(
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
    def is_u_periodic(self) -> bool:
        """Whether the surface is periodic in u."""
    def is_v_periodic(self) -> bool:
        """Whether the surface is periodic in v."""
    def is_u_closed(self, tol: Optional[_ToleranceLike] = ...) -> bool:
        """Whether travelling in `u` returns to its start: the boundary
        curves at the two `u` bounds coincide (`Tolerance.default()`
        when omitted). Periodic in `u` implies closed in `u`, not
        conversely.
        """
    def is_v_closed(self, tol: Optional[_ToleranceLike] = ...) -> bool:
        """Whether travelling in `v` returns to its start, as
        `is_u_closed`.
        """
    def continuity_u(self, u: float, tol: Optional[_ToleranceLike] = ...) -> Union[int, float]:
        """Continuity order `n` of `C^n` in `u` at parameter `u` as an
        `int`, or `math.inf` where smooth. A B-spline reports
        `u_degree - multiplicity` at a `u` knot (within
        `tol.parametric`); the ends of a non-periodic direction
        report smooth, and a periodic direction's seam counts as a
        knot. Use `BSplineSurface.first_derivatives_sided` and its siblings to
        evaluate each side of a knot reporting a low order.
        """
    def continuity_v(self, v: float, tol: Optional[_ToleranceLike] = ...) -> Union[int, float]:
        """Continuity order in `v` at parameter `v`, as `continuity_u`."""
    def min_continuity_u(self) -> Union[int, float]:
        """The lowest `continuity_u` anywhere inside the `u` range, ends
        excluded (`math.inf` when smooth everywhere).
        """
    def min_continuity_v(self) -> Union[int, float]:
        """The lowest `continuity_v` anywhere inside the `v` range, ends
        excluded (`math.inf` when smooth everywhere).
        """
    def normal(
        self, u: float, v: float, tol: Optional[_ToleranceLike] = ...
    ) -> Optional[Vector]:
        """Unit normal at surface parameters `(u, v)`
        (`Tolerance.default()` when omitted). Singular points such as
        poles give their limiting normal; returns `None` at an apex,
        whose meridians each have their own, and where no normal can
        be recovered.
        """
    def curvatures(
        self, u: float, v: float, tol: Optional[_ToleranceLike] = ...
    ) -> Optional[SurfaceCurvatures]:
        """Principal curvatures with their directions, and the mean and
        Gaussian curvatures, at `(u, v)` (`Tolerance.default()` when
        omitted), as a `SurfaceCurvatures`. `None` where the unit
        normal degenerates, as at a sphere's pole or a cone's apex,
        or a derivative is not finite. Curvatures are measured along
        the unit normal: negative bends away from it.
        """
    def project_point(
        self, point: Point, tol: Optional[_ToleranceLike] = ...
    ) -> SurfaceProjection:
        """Projects `point` onto the surface: the closest point, as a
        `SurfaceProjection`.

        A distance within `tol.confusion` (`1e-7` by default) snaps
        to `0.0`, so projection agrees with `contains`; pass a
        smaller `Tolerance` to read a smaller distance.
        """
    def project_points(
        self, points: List[Point], tol: Optional[_ToleranceLike] = ...
    ) -> List[SurfaceProjection]:
        """Projects each point in `points` onto the surface, returning a
        `SurfaceProjection` per point.
        """
    def project_points_array(
        self, points: numpy.typing.ArrayLike, tol: Optional[_ToleranceLike] = ...
    ) -> Tuple[_FloatArray, _FloatArray, _FloatArray]:
        """Projects each row of the `(N, 3)` array-like `points` onto
        the surface, returning numpy arrays `(parameters, points,
        distances)`: the `(u, v)` row, the closest point and its
        distance per row.
        """
    def contains(self, point: Point, tol: Optional[_ToleranceLike] = ...) -> bool:
        """Returns whether `point` lies on the surface within `tol`
        (`Tolerance.default()` when omitted).
        """
    def parameters_of(
        self, point: Point, tol: Optional[_ToleranceLike] = ...
    ) -> Tuple[float, float]:
        """The `(u, v)` parameters of `point`, a point on the surface,
        within `tol` (`Tolerance.default()` when omitted). Off the
        surface, a B-spline or Bézier surface gives the parameters of
        the closest point, and an elementary surface those of its
        closed-form inversion.
        """
    def bounding_box(self) -> Optional[BoundingBox]:
        """The axis-aligned box around the shape (`None` when unbounded).
        It touches an elementary shape; for a B-spline it is the box of the
        control points, which may be wider than the shape, as
        `tight_bounding_box` is not."""
    def tight_bounding_box(self, tol: Optional[_ToleranceLike] = ...) -> Optional[BoundingBox]:
        """The axis-aligned box around the shape as tightly as
        `tol.confusion` (`None` when unbounded): it holds the shape, and each
        of its sides is within `tol.confusion` of a point of it. The same as
        `bounding_box` on an elementary shape; a B-spline's comes from halving
        its Bézier pieces, a revolution surface's from its profile's.
        `tol` is `Tolerance.default()` when omitted."""
    def u_iso(self, u: float) -> _UIsoCurveType:
        """The `u` isoparametric curve at `u`, parametrized by `v` like
        the surface: `u_iso(u).point_at(t)` is `point_at(u, t)`.
        Raises `SurfaceIsoError` for a non-finite parameter or one
        outside a non-periodic direction (a periodic one wraps).
        """
    def v_iso(self, v: float) -> _VIsoCurveType:
        """The `v` isoparametric curve at `v`, parametrized by `u` like
        the surface, as `u_iso`.
        """
    @overload
    def uv_curve_to_3d(
        self, uv_curve: _WholeCurve2D, tol: Optional[_ToleranceLike] = ...
    ) -> _WholeCurve:
        """The exact 3D curve that `uv_curve`, a curve in this surface's (u, v)
        parameter space, traces on the surface, in the same parameter: its
        point at `t` is the surface's point at the UV curve's (u, v) at `t`.

        Closed forms: any curve on a plane, as the same kind of curve; any
        curve on the extrusion of a line, where the affine image keeps its
        kind in the same parameter; a line along u or v on any surface, as
        the iso curve run from the line's origin, as is a degree-1 B-spline
        or Bezier curve with two poles; and any other line, parabola,
        B-spline or Bezier curve on a B-spline or Bezier surface, or on the
        extrusion of a B-spline or Bezier profile, composed exactly into a
        B-spline. Raises `UvCurveTo3DError` for any other pair, when the
        line runs along an iso curve that is a point (a sphere's pole, a
        cone's apex), when the surface has no iso curve there, when a
        composed curve leaves the surface's domain, or when its degree would
        be above 25. `tol` is `Tolerance.default()` when omitted.

        A `TrimmedCurve2D` traces the `TrimmedCurve` over the same range
        of its basis's 3D curve. Raises `UvCurveTo3DError` when the piece
        leaves the surface's domain, and `TrimmedCurveConstructionError`
        when it wraps more than once around a periodic direction.
        """
    @overload
    def uv_curve_to_3d(
        self, uv_curve: TrimmedCurve2D, tol: Optional[_ToleranceLike] = ...
    ) -> TrimmedCurve: ...
    @overload
    def curve_3d_to_uv(
        self, curve: _WholeCurve, tol: Optional[_ToleranceLike] = ...
    ) -> _WholeCurve2D:
        """The exact UV curve of `curve`, a 3D curve lying on this surface: the
        curve in the surface's (u, v) parameter space, in the same parameter,
        whose point at `t` is the (u, v) of the curve's point at `t`.

        Closed forms: any line, circle, ellipse, parabola, hyperbola, B-spline
        or Bezier curve on a plane, as the 2D curve of the same kind; a line
        along the axis of a cylinder or a generator of a cone, as a line along
        v; a circle about the axis of a cylinder or cone, a meridian or
        parallel of a sphere, or a poloidal or toroidal circle of a torus, as a
        line along u or v; and on an extrusion or revolution, a generator line,
        a parallel circle, or the profile translated or turned in its own
        parameter, as a line along u or v. A `TrimmedCurve` maps to the
        `TrimmedCurve2D` of its basis's image over the same range. Raises
        `Curve3DToUvError` for any other pair, or when the curve does not
        lie on the surface within `tol` (`Tolerance.default()` when
        omitted).
        """
    @overload
    def curve_3d_to_uv(
        self, curve: TrimmedCurve, tol: Optional[_ToleranceLike] = ...
    ) -> TrimmedCurve2D: ...
    def curve_3d_to_uv_within(
        self,
        curve: Curve,
        max_deviation: float,
        tol: Optional[_ToleranceLike] = ...,
    ) -> CurveApproximation2D:
        """The UV curve of `curve` on this surface, within `max_deviation`
        at the same parameter, over the curve's bounds: one B-spline
        whose image `point_at(u, v)` at every `t` lies within
        `max_deviation` of `curve.point_at(t)`, continuous across the
        surface's seams.

        The UV curve is one non-rational, non-periodic `BSplineCurve2D`
        of degree 5 in the 3D curve's parameter, C2 except at a knot
        where the 3D curve is less smooth. Its deviation is the largest
        sampled distance between the surface's point at its point at `t`
        and `curve.point_at(t)`: each span is checked at nine points
        against half `max_deviation`, so it is never proven.

        The counterpart of `curve_3d_to_uv`, which gives the exact UV
        curve and refuses the pairs that have none, such as an ellipse
        on a cylinder or any curve on a spline surface, and the mirror
        of `uv_curve_to_3d_within`. `curve` is any 3D curve or
        `TrimmedCurve`; a line must be trimmed. Raises
        `Curve3DToUvApproximationError` when `max_deviation` is not a
        positive finite number, when the curve is unbounded or a point,
        when it is farther from the surface than half `max_deviation`,
        when it lingers at a sphere's pole or a cone's apex, or when the
        deviation is not reached. A curve that starts or ends at such a
        point takes the limit of the UV curve from inside: a meridian
        keeps its own `u` up to the pole. `tol` is `Tolerance.default()`
        when omitted. Other Python threads run during the
        approximation.
        """
    def uv_curve_to_3d_within(
        self,
        uv_curve: Curve2D,
        max_deviation: float,
        tol: Optional[_ToleranceLike] = ...,
    ) -> CurveApproximation:
        """The 3D curve that `uv_curve` traces on this surface, within
        `max_deviation` at the same parameter, over the UV curve's
        bounds: one B-spline whose point at every `t` lies within
        `max_deviation` of the surface's point at the UV curve's (u, v)
        at `t`.

        The 3D curve is one non-rational, non-periodic `BSplineCurve` in
        the UV curve's parameter. It is of degree 9 and C4 on an
        elementary surface along a UV curve with no knot, unless degree
        5 takes no more poles, and of degree 5 and C2 wherever the UV
        curve is otherwise. Its deviation is a proven bound for a
        straight UV curve on a plane, cylinder, cone, sphere or torus;
        for every other pair it is sampled, each span checked at nine
        points against half `max_deviation`.

        The counterpart of `uv_curve_to_3d`, which gives the exact 3D
        curve and refuses the pairs that have none, such as a slanted
        line on a cylinder. `uv_curve` is any 2D curve or
        `TrimmedCurve2D`; a line must be trimmed. Raises
        `UvCurveTo3DApproximationError` when `max_deviation` is not a
        positive finite number, when the UV curve is unbounded, when the
        3D curve collapses to a point, when the UV curve leaves the
        surface's domain, or when the deviation is not reached. `tol` is
        `Tolerance.default()` when omitted. Other Python threads run
        during the approximation.
        """
    @overload
    def intersect(
        self, other: Curve, tol: Optional[_ToleranceLike] = ...
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
        self, other: Surface, tol: Optional[_ToleranceLike] = ...
    ) -> SurfaceSurfaceIntersection: ...

class _SidedSurfaceEvaluationMethods:
    def first_derivatives_sided(
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
    def second_derivatives_sided(
        self, u: float, v: float, u_side: Side, v_side: Side
    ) -> Tuple[Point, Tuple[Vector, Vector], Tuple[Vector, Vector, Vector]]:
        """The point, the first partial derivatives `(du, dv)` and the
        second partial derivatives `(du du, du dv, dv dv)` at
        `(u, v)`, each direction from its own side, as
        `first_derivatives_sided`.
        """
    def derivative_sided(
        self, u: float, v: float, du: int, dv: int, u_side: Side, v_side: Side
    ) -> Vector:
        """The partial derivative of order `(du, dv)`, with
        `1 <= du + dv <= 2`, at `(u, v)`, each direction from its own
        side, as `first_derivatives_sided`.
        """

class _SurfaceExtremaMethods:
    def extrema(
        self, point: Point, tol: Optional[_ToleranceLike] = ...
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
    def axis(self) -> Axis: ...
    def signed_distance(self, point: Point) -> float: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...
    def frame(self) -> Frame: ...

@final
class Cylinder(_SurfaceMethods[Line, Circle]):
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
class Cone(_SurfaceMethods[Line, Circle]):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        center: Point,
        axis_direction: Vector,
        half_angle: float,
        reference_radius: float,
    ) -> Cone: ...
    @staticmethod
    def from_frame(frame: Frame, half_angle: float, reference_radius: float) -> Cone: ...
    @staticmethod
    def from_two_points_and_radii(
        p1: Point, p2: Point, r1: float, r2: float
    ) -> Cone: ...
    @staticmethod
    def from_apex(apex: Point, axis_direction: Vector, half_angle: float) -> Cone: ...
    def axis(self) -> Axis: ...
    def apex(self) -> Point: ...
    def half_angle(self) -> float: ...
    def reference_radius(self) -> float: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...
    def frame(self) -> Frame: ...

@final
class Sphere(_SurfaceMethods[Circle, Circle]):
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
class Torus(_SurfaceMethods[Circle, Circle]):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, center: Point, axis_direction: Vector, major_radius: float, minor_radius: float
    ) -> Torus: ...
    @staticmethod
    def from_frame(frame: Frame, major_radius: float, minor_radius: float) -> Torus: ...
    def major_radius(self) -> float: ...
    def minor_radius(self) -> float: ...
    def to_bspline(self) -> BSplineSurface: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...
    def axis(self) -> Axis: ...
    def frame(self) -> Frame: ...
    def center(self) -> Point: ...

@final
class ExtrusionSurface(_SurfaceMethods[Line, _WholeCurve]):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, profile: _WholeCurve, direction: Vector) -> ExtrusionSurface: ...
    def profile(self) -> _WholeCurve: ...
    def direction(self) -> Vector: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...

@final
class RevolutionSurface(_SurfaceMethods[_WholeCurve, Circle]):
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, profile: _WholeCurve, axis: Axis) -> RevolutionSurface: ...
    def profile(self) -> _WholeCurve: ...
    def axis(self) -> Axis: ...
    def to_bspline(self) -> BSplineSurface: ...
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...

@final
class BSplineSurface(
    _SurfaceMethods[BSplineCurve, BSplineCurve],
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
        points: Union[List[List[Point]], numpy.typing.ArrayLike],
        u_degree: int,
        v_degree: int,
        u_pole_count: int,
        v_pole_count: int,
        parametrization: _Parametrization = ...,
    ) -> BSplineSurface: ...
    @staticmethod
    def approximate_within(
        points: Union[List[List[Point]], numpy.typing.ArrayLike],
        u_degree: int,
        v_degree: int,
        max_deviation: float,
        parametrization: _Parametrization = ...,
    ) -> SurfaceApproximation: ...
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
    def to_bspline(self) -> BSplineSurface: ...
    def extract(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface: ...
    def extract_u(self, start: float, end: float) -> BSplineSurface: ...
    def extract_v(self, start: float, end: float) -> BSplineSurface: ...
    def split_u_at(self, u: float) -> Tuple[BSplineSurface, BSplineSurface]: ...
    def split_v_at(self, v: float) -> Tuple[BSplineSurface, BSplineSurface]: ...

@final
class BezierSurface(
    _SurfaceMethods[BezierCurve, BezierCurve],
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

__all__ = [
    "GeomcoreError",
    "ConstructionError",
    "FrameConstructionError",
    "TrimmedCurveConstructionError",
    "KnotInsertionError",
    "KnotRemovalError",
    "DegreeElevationError",
    "BSplineSplitError",
    "BSplineConversionError",
    "Curve3DToUvError",
    "Curve3DToUvApproximationError",
    "UvCurveTo3DError",
    "UvCurveTo3DApproximationError",
    "SurfaceIsoError",
    "DiscretizationCriteriaError",
    "Point",
    "Point2D",
    "Vector",
    "Vector2D",
    "Axis2D",
    "Axis",
    "Frame2D",
    "Frame",
    "Tolerance",
    "Transform",
    "BoundingBox",
    "BoundingBox2D",
    "Discretization2D",
    "Discretization",
    "Side",
    "PointParametrization",
    "IntersectionMethod",
    "CompletenessStatus",
    "IncompleteReason",
    "Completeness",
    "CurveCurveRegion",
    "CurveSurfaceRegion",
    "SurfaceSurfaceRegion",
    "CurveCurveHit",
    "CurveCurveHit2D",
    "CurveSurfaceHit",
    "CurveCurveOverlap",
    "CurveSurfaceOverlap",
    "IntersectionCurve",
    "CurveCurveIntersection",
    "CurveCurveIntersection2D",
    "CurveSurfaceIntersection",
    "SurfaceSurfaceIntersection",
    "CurveProjection",
    "CurveProjection2D",
    "SurfaceProjection",
    "CurveCurveExtremum",
    "CurveCurveExtrema",
    "SurfaceCurvatures",
    "CurveApproximation",
    "CurveApproximation2D",
    "SurfaceApproximation",
    "Line",
    "Circle",
    "Ellipse",
    "Parabola",
    "Hyperbola",
    "BSplineCurve",
    "BezierCurve",
    "Line2D",
    "Circle2D",
    "Ellipse2D",
    "Parabola2D",
    "Hyperbola2D",
    "BSplineCurve2D",
    "BezierCurve2D",
    "TrimmedCurve",
    "TrimmedCurve2D",
    "Plane",
    "Cylinder",
    "Cone",
    "Sphere",
    "Torus",
    "ExtrusionSurface",
    "RevolutionSurface",
    "BSplineSurface",
    "BezierSurface",
    "curve_3d_to_uv_within",
    "uv_curve_to_3d_within",
    "curve_curve_extrema",
    "Curve",
    "Curve2D",
    "Surface",
]
