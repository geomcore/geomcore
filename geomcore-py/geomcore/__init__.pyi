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

import sys

import numpy
import numpy.typing

if sys.version_info >= (3, 13):
    from warnings import deprecated
else:
    from typing_extensions import deprecated

_FloatArray = numpy.typing.NDArray[numpy.float64]
# Parameters take any sequence (a list, a tuple) or a numpy array.
_Floats = Union[Sequence[float], numpy.typing.NDArray[Any]]
_Ints = Union[Sequence[int], numpy.typing.NDArray[Any]]
# A point argument takes a point, or a sequence of its coordinates: a tuple,
# a list or a row of an array (#322). A list of points takes a list of
# those or an `(N, 3)` (`(N, 2)`) array, empty included; a grid, rows of
# them or an `(N, M, 3)` array.
_PointLike = Union["Point", Sequence[float], numpy.typing.NDArray[Any]]
_Point2DLike = Union["Point2D", Sequence[float], numpy.typing.NDArray[Any]]
_Points = Union[Sequence[_PointLike], numpy.typing.NDArray[Any]]
_Points2D = Union[Sequence[_Point2DLike], numpy.typing.NDArray[Any]]
_PointGrid = Union[Sequence[_Points], numpy.typing.NDArray[Any]]
# `(u, v)` pairs: a list of pairs or an `(N, 2)` array.
_Pairs = Union[Sequence[Sequence[float]], numpy.typing.NDArray[Any]]
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

class TrimmedSurfaceConstructionError(ConstructionError):
    """A trimmed surface could not be built, for example from an invalid parameter range."""

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
    """A point in 3D space with `x`, `y` and `z` coordinates.

    An argument that takes a point takes a `Point` or a sequence of three
    numbers, such as a tuple or a row of an array (#322).
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, x: float, y: float, z: float) -> Point: ...
    @staticmethod
    def origin() -> Point:
        """The origin (0, 0, 0)."""
    @property
    def x(self) -> float:
        """The x coordinate."""
    @property
    def y(self) -> float:
        """The y coordinate."""
    @property
    def z(self) -> float:
        """The z coordinate."""
    def distance(self, other: _PointLike) -> float:
        """Euclidean distance to another point."""
    def coincides(self, other: _PointLike, tol: Optional[_ToleranceLike] = None) -> bool:
        """Returns whether `other` is within `tol.confusion` of this point
        (`Tolerance.default()` when omitted). `==` compares coordinates
        exactly instead.
        """
    def components(self) -> Tuple[float, float, float]:
        """The coordinates as a tuple `(x, y, z)`."""
    def midpoint(self, other: _PointLike) -> Point:
        """The point halfway between this point and `other`."""
    def transformed(self, transform: Transform) -> Point:
        """The image of this point under `transform`, the same as
        `transform.apply_point(point)`.
        """
    def __repr__(self) -> str: ...
    def is_finite(self) -> bool:
        """Whether every coordinate is finite."""
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
    ) -> numpy.typing.NDArray[Any]:
        """The coordinates as a numpy array of shape `(3,)`, so
        `numpy.array(p)` and `numpy.array(list_of_points)` give floats.
        """
    __array_ufunc__: ClassVar[None]

@final
class Point2D:
    """A point in 2D space (used by curves living in a surface's (u, v) domain).

    An argument that takes a 2D point takes a `Point2D` or a sequence of
    two numbers (#322).
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, x: float, y: float) -> Point2D: ...
    @staticmethod
    def origin() -> Point2D:
        """The origin (0, 0)."""
    @property
    def x(self) -> float:
        """The x coordinate."""
    @property
    def y(self) -> float:
        """The y coordinate."""
    def distance(self, other: _Point2DLike) -> float:
        """Euclidean distance to another point."""
    def coincides(self, other: _Point2DLike, tol: Optional[_ToleranceLike] = None) -> bool:
        """Returns whether `other` is within `tol.confusion` of this point
        (`Tolerance.default()` when omitted). `==` compares coordinates
        exactly instead.
        """
    def components(self) -> Tuple[float, float]:
        """The coordinates as a tuple `(x, y)`."""
    def midpoint(self, other: _Point2DLike) -> Point2D:
        """The point halfway between this point and `other`."""
    def is_finite(self) -> bool:
        """Whether every coordinate is finite."""
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
    ) -> numpy.typing.NDArray[Any]:
        """The coordinates as a numpy array of shape `(2,)`, so
        `numpy.array(p)` and `numpy.array(list_of_points)` give floats.
        """
    __array_ufunc__: ClassVar[None]

@final
class Vector:
    """A vector in 3D space.

    Components are read with the properties `x`, `y` and `z` or with
    `components()`; the unit vectors along the world axes are the static
    methods `unit_x()`, `unit_y()` and `unit_z()`. Vectors add, subtract,
    negate and scale by a number with the usual operators.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, x: float, y: float, z: float) -> Vector: ...
    @staticmethod
    def zero() -> Vector:
        """The zero vector."""
    @staticmethod
    def unit_x() -> Vector:
        """The unit vector along the world x axis."""
    @staticmethod
    def unit_y() -> Vector:
        """The unit vector along the world y axis."""
    @staticmethod
    def unit_z() -> Vector:
        """The unit vector along the world z axis."""
    @property
    def x(self) -> float:
        """The x coordinate."""
    @property
    def y(self) -> float:
        """The y coordinate."""
    @property
    def z(self) -> float:
        """The z coordinate."""
    def components(self) -> Tuple[float, float, float]:
        """The `(x, y, z)` components as a tuple."""
    def dot(self, other: Vector) -> float:
        """Dot product with another vector."""
    def cross(self, other: Vector) -> Vector:
        """Cross product with another vector."""
    def magnitude(self) -> float:
        """Euclidean length of the vector."""
    def square_magnitude(self) -> float:
        """The squared length of the vector."""
    def is_finite(self) -> bool:
        """Whether every component is finite."""
    def normalized(self) -> Optional[Vector]:
        """The unit vector with the same direction, or `None` for a vector too
        short to normalize: check the result before passing it on, as a
        constructor taking a vector refuses `None`.
        """
    def angle_with_reference(self, other: Vector, reference: Vector) -> float:
        """The signed angle from this vector to `other` around `reference`
        (right-hand rule), in (-pi, pi]. A zero `reference` leaves the sign
        unknown: the result is then 0 or pi; `angle(other)` is the unsigned
        angle.
        """
    def angle(self, other: Vector) -> float:
        """The unsigned angle between the two vectors, in [0, pi]; 0 when
        either is zero.
        """
    def is_parallel(self, other: Vector, angular: Optional[float] = None) -> bool:
        """Whether the vectors are parallel, in either orientation: the sine
        of the angle between them is at most `angular`
        (`Tolerance.default().angular` when omitted). A zero vector is
        parallel to nothing.
        """
    def is_perpendicular(self, other: Vector, angular: Optional[float] = None) -> bool:
        """Whether the vectors are perpendicular: the cosine of the angle
        between them is at most `angular` (`Tolerance.default().angular`
        when omitted). A zero vector is perpendicular to nothing.
        """
    def component_along(self, direction: Vector) -> Vector:
        """The part of this vector along `direction`, of any length; the zero
        vector when `direction` is zero.
        """
    def component_perpendicular_to(self, direction: Vector) -> Vector:
        """The part of this vector across `direction`: the vector less its
        component along `direction`.
        """
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
    ) -> numpy.typing.NDArray[Any]:
        """The coordinates as a numpy array of shape `(3,)`, so
        `numpy.array(p)` and `numpy.array(list_of_points)` give floats.
        """
    __array_ufunc__: ClassVar[None]
    def __truediv__(self, scalar: float, /) -> Vector: ...
    def __abs__(self) -> float: ...
    def __repr__(self) -> str: ...

@final
class Vector2D:
    """A vector in 2D space."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, x: float, y: float) -> Vector2D: ...
    @property
    def x(self) -> float:
        """The x coordinate."""
    @property
    def y(self) -> float:
        """The y coordinate."""
    @staticmethod
    def zero() -> Vector2D:
        """The zero vector."""
    @staticmethod
    def unit_x() -> Vector2D:
        """The unit vector along the x axis."""
    @staticmethod
    def unit_y() -> Vector2D:
        """The unit vector along the y axis."""
    def components(self) -> Tuple[float, float]:
        """The `(x, y)` components as a tuple."""
    def dot(self, other: Vector2D) -> float:
        """Dot product with another vector."""
    def cross(self, other: Vector2D) -> float:
        """The z component of the 3D cross product: positive when `other` lies
        counter-clockwise from this vector.
        """
    def perpendicular(self) -> Vector2D:
        """This vector turned a quarter turn counter-clockwise."""
    def angle(self, other: Vector2D) -> float:
        """The unsigned angle between the two vectors, in [0, pi]; 0 when
        either is zero.
        """
    def signed_angle(self, other: Vector2D) -> float:
        """The signed angle from this vector to `other`, in (-pi, pi]:
        positive when `other` is counter-clockwise from this one.
        """
    def rotated(self, angle: float) -> Vector2D:
        """This vector turned counter-clockwise by `angle` radians."""
    def is_parallel(self, other: Vector2D, angular: Optional[float] = None) -> bool:
        """Whether the vectors are parallel, in either orientation: the sine
        of the angle between them is at most `angular`
        (`Tolerance.default().angular` when omitted). A zero vector is
        parallel to nothing.
        """
    def is_perpendicular(self, other: Vector2D, angular: Optional[float] = None) -> bool:
        """Whether the vectors are perpendicular: the cosine of the angle
        between them is at most `angular` (`Tolerance.default().angular`
        when omitted). A zero vector is perpendicular to nothing.
        """
    def component_along(self, direction: Vector2D) -> Vector2D:
        """The part of this vector along `direction`, of any length; the zero
        vector when `direction` is zero.
        """
    def component_perpendicular_to(self, direction: Vector2D) -> Vector2D:
        """The part of this vector across `direction`: the vector less its
        component along `direction`.
        """
    def magnitude(self) -> float:
        """Euclidean length of the vector."""
    def square_magnitude(self) -> float:
        """The squared length of the vector."""
    def is_finite(self) -> bool:
        """Whether every component is finite."""
    def normalized(self) -> Optional[Vector2D]:
        """The unit vector with the same direction, or `None` for a vector too
        short to normalize: check the result before passing it on, as a
        constructor taking a vector refuses `None`.
        """
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
    ) -> numpy.typing.NDArray[Any]:
        """The coordinates as a numpy array of shape `(2,)`, so
        `numpy.array(p)` and `numpy.array(list_of_points)` give floats.
        """
    __array_ufunc__: ClassVar[None]
    def __truediv__(self, scalar: float, /) -> Vector2D: ...
    def __abs__(self) -> float: ...
    def __repr__(self) -> str: ...

@final
class Axis2D:
    """An axis in the plane: a point plus a unit direction."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, origin: _Point2DLike, direction: Vector2D) -> Axis2D: ...
    def origin(self) -> Point2D:
        """The axis origin."""
    def direction(self) -> Vector2D:
        """The unit direction of the axis."""
    def point_at(self, t: float) -> Point2D:
        """The point at `t` along the axis, `origin + t * direction`."""
    def parameter_of(self, point: _Point2DLike) -> float:
        """The parameter of the foot of `point` on the axis, so that
        `point_at` gives the `closest_point`.
        """
    def closest_point(self, point: _Point2DLike) -> Point2D:
        """The point of the axis closest to `point`."""
    def distance(self, point: _Point2DLike) -> float:
        """The distance from `point` to the axis."""
    def contains(self, point: _Point2DLike, tol: Optional[_ToleranceLike] = None) -> bool:
        """Whether `point` lies on the axis, within `tol.confusion`
        (`Tolerance.default()` when omitted).
        """
    def is_parallel(self, other: Axis2D, angular: Optional[float] = None) -> bool:
        """Whether the axes have parallel directions, in either orientation,
        within `angular` (`Tolerance.default().angular` when omitted).
        """
    def reversed(self) -> Axis2D:
        """The axis with its direction negated, through the same origin."""
    def __repr__(self) -> str: ...

@final
class Frame2D:
    """An orthonormal placement frame in the plane (origin + x/y directions).
    The frame is direct when y is x turned a quarter turn counter-clockwise.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, origin: _Point2DLike, x_direction: Vector2D, y_direction: Vector2D
    ) -> Frame2D: ...
    @staticmethod
    def from_x(origin: _Point2DLike, x_direction: Vector2D) -> Frame2D:
        """Build a direct frame from an origin and an x direction."""
    @staticmethod
    def world() -> Frame2D:
        """The world frame at the origin."""
    def origin(self) -> Point2D:
        """The frame origin."""
    def x_direction(self) -> Vector2D:
        """The unit x direction."""
    def y_direction(self) -> Vector2D:
        """The unit y direction."""
    def is_direct(self) -> bool:
        """Whether y is x turned a quarter turn counter-clockwise."""
    def x_axis(self) -> Axis2D:
        """The axis through the origin along the x direction."""
    def y_axis(self) -> Axis2D:
        """The axis through the origin along the y direction."""
    def local_coordinates(self, point: _Point2DLike) -> Tuple[float, float]:
        """The coordinates `(u, v)` of `point` along the x and y directions,
        measured from the origin.
        """
    def point_at(self, u: float, v: float) -> Point2D:
        """The point at `u` along x and `v` along y from the origin."""
    def vector_at(self, u: float, v: float) -> Vector2D:
        """The vector `u * x_direction + v * y_direction`: `point_at` for a
        direction, which ignores the origin.
        """
    def local_components(self, vector: Vector2D) -> Tuple[float, float]:
        """The components `(u, v)` of `vector` along the x and y directions:
        `local_coordinates` for a direction, which ignores the origin.
        """
    def __repr__(self) -> str: ...

@final
class Axis:
    """An axis: a point plus a unit direction."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, origin: _PointLike, direction: Vector) -> Axis: ...
    def origin(self) -> Point:
        """The axis origin."""
    def direction(self) -> Vector:
        """The unit direction of the axis."""
    def point_at(self, t: float) -> Point:
        """The point at `t` along the axis, `origin + t * direction`."""
    def parameter_of(self, point: _PointLike) -> float:
        """The parameter of the foot of `point` on the axis, so that
        `point_at` gives the `closest_point`.
        """
    def closest_point(self, point: _PointLike) -> Point:
        """The point of the axis closest to `point`."""
    def distance(self, point: _PointLike) -> float:
        """The distance from `point` to the axis."""
    def contains(self, point: _PointLike, tol: Optional[_ToleranceLike] = None) -> bool:
        """Whether `point` lies on the axis, within `tol.confusion`
        (`Tolerance.default()` when omitted).
        """
    def is_parallel(self, other: Axis, angular: Optional[float] = None) -> bool:
        """Whether the axes have parallel directions, in either orientation,
        within `angular` (`Tolerance.default().angular` when omitted).
        """
    def is_coaxial(self, other: Axis, tol: Optional[_ToleranceLike] = None) -> bool:
        """Whether the axes are the same line, in either orientation: parallel
        within `tol.angular`, the other origin on this axis within
        `tol.confusion` (`Tolerance.default()` when omitted).
        """
    def reversed(self) -> Axis:
        """The axis with its direction negated, through the same origin."""

@final
class Frame:
    """An orthonormal placement frame (origin + x/y/z directions).

    The constructors build direct (right-handed, `x × y = z`) frames;
    `from_axes` accepts either handedness, and reversing a surface's
    direction or mirroring it turns its frame indirect.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, origin: _PointLike, z_direction: Vector, x_hint: Vector) -> Frame: ...
    @staticmethod
    def from_z(origin: _PointLike, z_direction: Vector) -> Frame:
        """Build a frame from an origin and main (z) direction, deriving an
        arbitrary perpendicular x direction.
        """
    @staticmethod
    def world() -> Frame:
        """The world frame at the origin."""
    @staticmethod
    def from_axes(
        origin: _PointLike, x_direction: Vector, y_direction: Vector, z_direction: Vector
    ) -> Frame:
        """Build a frame from an origin and its three directions (each
        normalized internally; they must be mutually orthogonal), of
        either handedness. Raises `FrameConstructionError` for a
        non-finite input, a zero direction or non-orthogonal directions.
        """
    def is_direct(self) -> bool:
        """Whether the frame is right-handed (`x × y = z`) rather than
        left-handed (`x × y = -z`).
        """
    def flipped_x(self) -> Frame:
        """The frame with its x direction negated (handedness flips)."""
    def flipped_y(self) -> Frame:
        """The frame with its y direction negated (handedness flips)."""
    def flipped_z(self) -> Frame:
        """The frame with its z direction negated (handedness flips)."""
    def origin(self) -> Point:
        """The frame origin."""
    def x_direction(self) -> Vector:
        """The unit x direction."""
    def y_direction(self) -> Vector:
        """The unit y direction."""
    def z_direction(self) -> Vector:
        """The unit z direction."""
    def axis(self) -> Axis:
        """The axis through the origin along the z direction."""
    def local_coordinates(self, point: _PointLike) -> Tuple[float, float, float]:
        """The coordinates `(u, v, w)` of `point` along the x, y and z
        directions, measured from the origin.
        """
    def point_at(self, u: float, v: float, w: float) -> Point:
        """The point at `u` along x, `v` along y and `w` along z from the origin."""
    def vector_at(self, u: float, v: float, w: float) -> Vector:
        """The vector `u * x_direction + v * y_direction + w * z_direction`:
        `point_at` for a direction, which ignores the origin.
        """
    def local_components(self, vector: Vector) -> Tuple[float, float, float]:
        """The components `(u, v, w)` of `vector` along the x, y and z
        directions: `local_coordinates` for a direction, which ignores the
        origin.
        """
    def transformed(self, transform: Transform) -> Frame:
        """The image of the frame under `transform`: a point at local
        coordinates `(u, v, w)` maps to the point at the same coordinates,
        scaled by `transform.scale_factor()`, in the image. A mirror
        changes the handedness. Raises `ConstructionError` when the
        transform collapses space or the image overflows.
        """

@final
class Tolerance:
    """Distance and angle tolerances shared by all geometric queries:
    containment, projection, intersection and fitting.

    `confusion` is a distance in model units (default 1e-7), `angular` an
    angle in radians (default 1e-12, the default confusion over
    `MODEL_SIZE`) and `parametric` a distance between curve or surface
    parameters (default 1e-9). Scale `confusion` with the model: models
    are assumed to fit in a cube of side `confusion / angular`.

    `Tolerance()` with no arguments equals `Tolerance.default()`. Every
    `tol=` argument also takes a number, read as
    `Tolerance(confusion=number)`, so `tol=1e-3` loosens the distance
    tolerance alone.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, confusion: float = 1e-07, angular: float = 1e-12, parametric: float = 1e-09
    ) -> Tolerance: ...
    @staticmethod
    def default() -> Tolerance:
        """The default tolerances (confusion 1e-7, angular 1e-12, parametric 1e-9)."""
    @property
    def confusion(self) -> float:
        """Distance below which two points are considered coincident, in
        model units (default 1e-7).
        """
    @property
    def angular(self) -> float:
        """Angle below which two directions are parallel or orthogonal, in
        radians (default 1e-12).
        """
    @property
    def parametric(self) -> float:
        """Distance between two parameter values below which they are the
        same parameter (default 1e-9).
        """
    MODEL_SIZE: ClassVar[float]
    def angular_for_extent(self, extent: float) -> float:
        """The angle (radians) below which two directions are parallel for
        geometry reaching `extent` from the pivot of the turn: turning by it
        moves that geometry by at most `confusion`. Never below `angular`.
        """
    def __repr__(self) -> str: ...

# A `tol` argument: a `Tolerance`, or a number taken as `Tolerance(confusion=number)`.
_ToleranceLike = Union[Tolerance, float]

@final
class Transform:
    """A rigid (or mirrored/scaled) affine transformation."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @staticmethod
    def identity() -> Transform:
        """The identity transformation."""
    @staticmethod
    def translation(offset: Vector) -> Transform:
        """Translation by an offset vector."""
    @staticmethod
    def rotation(axis: Axis, angle: float) -> Transform:
        """Rotation by `angle` radians about `axis` (right-hand rule)."""
    @staticmethod
    def scaling(center: _PointLike, factor: float) -> Transform:
        """Uniform scaling about `center`."""
    @staticmethod
    def mirror_point(center: _PointLike) -> Transform:
        """Point reflection through `center`."""
    @staticmethod
    def mirror_axis(axis: Axis) -> Transform:
        """Reflection about an axis (rotation by pi around the line)."""
    @staticmethod
    def mirror_plane(frame: Frame) -> Transform:
        """Reflection across the plane through the frame origin with normal
        equal to the frame's z direction.
        """
    @staticmethod
    def between_frames(source: Frame, target: Frame) -> Transform:
        """The rigid motion taking the `source` placement to the `target` one:
        the point at local coordinates `(u, v, w)` in `source` maps to the
        point at the same coordinates in `target`. It includes a mirror when
        the frames have opposite handedness.
        """
    def scale_factor(self) -> float:
        """The factor by which the transformation scales lengths: 1 for a
        rotation, a translation or a mirror, 0 for a collapse onto a point.
        """
    def preserves_orientation(self) -> bool:
        """Whether the transformation keeps the handedness of space: false for
        mirrors (whose images of surfaces have their normals flipped) and
        for a collapse onto a point.
        """
    def is_identity(self, tol: Optional[_ToleranceLike] = None) -> bool:
        """Whether the transformation leaves space unchanged within `tol`
        (`Tolerance.default()` when omitted): its linear part within
        `tol.angular` of the identity's, entry by entry, and its translation
        shorter than `tol.confusion`.
        """
    def then(self, other: Transform) -> Transform:
        """Compose transformations: `self` is applied first, then `other`."""
    def apply_point(self, p: _PointLike) -> Point:
        """Apply the transformation to a point."""
    def apply_vector(self, v: Vector) -> Vector:
        """Apply the transformation to a vector (ignores the translation part)."""
    @staticmethod
    def from_matrix(matrix: Union[Sequence[Sequence[float]], numpy.typing.NDArray[Any]]) -> Transform:
        """The transformation of a 4×4 homogeneous `matrix`, four rows of
        four numbers (nested lists, or a numpy array): a point `p` maps to
        the upper 3×3 block times `p` plus the last column. Raises
        `ConstructionError` unless the last row is `[0, 0, 0, 1]` and the
        3×3 block is a uniform scale of a rotation or a mirror, its columns
        of one length and orthogonal within `1e-12`. `to_matrix()` gives
        the matrix back unchanged.
        """
    def to_matrix(self) -> List[List[float]]:
        """The 4×4 homogeneous matrix, as four rows of four floats: the linear
        part in the upper 3×3 block, the translation in the last column,
        and `[0, 0, 0, 1]` below. `numpy.array(t.to_matrix())` makes it an
        array.
        """
    def inverse(self) -> Transform:
        """The transformation that undoes this one: `t.then(t.inverse())`
        maps every point back to itself, within rounding. Raises
        `ConstructionError` for a transformation that collapses space onto
        a point (a scaling by 0).
        """
    def __matmul__(self, other: Transform, /) -> Transform: ...
    def __repr__(self) -> str: ...

@final
class BoundingBox:
    """The axis-aligned box around a curve or surface, holding the shape but
    possibly loose for B-splines (their control points' box); `None`
    wherever `bounding_box` reports no box for an unbounded shape.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, a: _PointLike, b: _PointLike) -> BoundingBox: ...
    @staticmethod
    def from_points(points: _Points) -> BoundingBox:
        """The smallest box holding every point of `points` (`ValueError`
        when there are none).
        """
    def min_corner(self) -> Point:
        """The corner with the smallest coordinates."""
    def max_corner(self) -> Point:
        """The corner with the largest coordinates."""
    def center(self) -> Point:
        """The middle of the box."""
    def diagonal(self) -> Vector:
        """The vector from the `min` corner to the `max` corner."""
    def union(self, other: BoundingBox) -> BoundingBox:
        """The smallest box holding both boxes."""
    def enlarged(self, distance: float) -> BoundingBox:
        """The box grown by `distance` on every side (negative shrinks it)."""
    def encloses(self, point: _PointLike) -> bool:
        """Whether the box holds `point` on every axis."""
    def __contains__(self, point: _PointLike, /) -> bool: ...
    def including(self, point: _PointLike) -> BoundingBox:
        """The smallest box holding this box and `point`."""
    def intersection(self, other: BoundingBox) -> Optional[BoundingBox]:
        """The box of the points the two boxes share, or `None` when they are
        apart. Boxes that only touch share a flat box.
        """
    def distance(self, point: _PointLike) -> float:
        """The distance from `point` to the box: 0 for a point in it, the
        distance to the nearest point of the box otherwise.
        """
    def overlaps(self, other: BoundingBox) -> bool:
        """Whether the boxes touch or overlap on every axis."""
    def __repr__(self) -> str: ...

@final
class BoundingBox2D:
    """The axis-aligned box around a 2D curve, holding the curve but possibly
    loose for B-splines (their control points' box); `None` wherever
    `bounding_box` reports no box for an unbounded curve.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, a: _Point2DLike, b: _Point2DLike) -> BoundingBox2D: ...
    @staticmethod
    def from_points(points: _Points2D) -> BoundingBox2D:
        """The smallest box holding every point of `points` (`ValueError`
        when there are none).
        """
    def min_corner(self) -> Point2D:
        """The corner with the smallest coordinates."""
    def max_corner(self) -> Point2D:
        """The corner with the largest coordinates."""
    def center(self) -> Point2D:
        """The middle of the box."""
    def diagonal(self) -> Vector2D:
        """The vector from the `min` corner to the `max` corner."""
    def union(self, other: BoundingBox2D) -> BoundingBox2D:
        """The smallest box holding both boxes."""
    def enlarged(self, distance: float) -> BoundingBox2D:
        """The box grown by `distance` on every side (negative shrinks it)."""
    def encloses(self, point: _Point2DLike) -> bool:
        """Whether the box holds `point` on both axes."""
    def __contains__(self, point: _Point2DLike, /) -> bool: ...
    def including(self, point: _Point2DLike) -> BoundingBox2D:
        """The smallest box holding this box and `point`."""
    def intersection(self, other: BoundingBox2D) -> Optional[BoundingBox2D]:
        """The box of the points the two boxes share, or `None` when they are
        apart. Boxes that only touch share a flat box.
        """
    def distance(self, point: _Point2DLike) -> float:
        """The distance from `point` to the box: 0 for a point in it, the
        distance to the nearest point of the box otherwise.
        """
    def overlaps(self, other: BoundingBox2D) -> bool:
        """Whether the boxes touch or overlap on both axes."""
    def __repr__(self) -> str: ...

@final
class Discretization2D:
    """A discretized 2D curve, such as a curve in a surface's parameter
    plane: the 2D counterpart of `Discretization`, its chord deviations
    being distances in the plane.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[Point2D]: ...
    def __reversed__(self) -> Iterator[Point2D]:
        """Iterates over the points from the last to the first."""
    def __array__(
        self, dtype: Optional[numpy.typing.DTypeLike] = None, copy: Optional[bool] = None
    ) -> numpy.typing.NDArray[Any]:
        """The points as a numpy array of shape `(N, 2)`, so
        `numpy.array(discretization)` gives floats, not an array of objects.
        """
    @overload
    def __getitem__(self, index: int, /) -> Point2D: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[Point2D]: ...
    @property
    def parameters(self) -> List[float]:
        """The parameter of each point, increasing."""
    @property
    def points(self) -> List[Point2D]:
        """The points on the curve."""
    @property
    def chord_deviations(self) -> List[float]:
        """The chord deviation of each segment."""
    @property
    def angle_deviations(self) -> List[float]:
        """The angle deviation of each segment, in radians."""
    def __repr__(self) -> str: ...

@final
class Discretization:
    """A discretized 3D curve: points on the curve in order of parameter,
    joined by straight segments. `parameters[i]` is the parameter of
    `points[i]`; segment `i` joins `points[i]` and `points[i + 1]` and
    carries a proven upper bound on its chord deviation (the largest
    distance from the curve to the segment) and angle deviation (the
    largest angle between the curve's tangent and the segment, π where no
    bound exists).
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[Point]: ...
    def __reversed__(self) -> Iterator[Point]:
        """Iterates over the points from the last to the first."""
    def __array__(
        self, dtype: Optional[numpy.typing.DTypeLike] = None, copy: Optional[bool] = None
    ) -> numpy.typing.NDArray[Any]:
        """The points as a numpy array of shape `(N, 3)`, so
        `numpy.array(discretization)` gives floats, not an array of objects.
        """
    @overload
    def __getitem__(self, index: int, /) -> Point: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[Point]: ...
    @property
    def parameters(self) -> List[float]:
        """The parameter of each point, increasing."""
    @property
    def points(self) -> List[Point]:
        """The points on the curve."""
    @property
    def chord_deviations(self) -> List[float]:
        """The chord deviation of each segment."""
    @property
    def angle_deviations(self) -> List[float]:
        """The angle deviation of each segment, in radians."""
    def __repr__(self) -> str: ...


class Side(enum.IntEnum):
    """Which side of a parameter to evaluate from, for one-sided limits at
    B-spline knots and periodic seams.

    `BELOW` approaches from the span ending at the parameter, `ABOVE`
    from the span starting there (the side plain evaluation takes,
    except at the upper bound). Away from knots and seams both sides
    agree.
    """
    BELOW = 0
    ABOVE = 1

class PointParametrization(str, enum.Enum):
    """How the parameters of fitted points are spaced: `CHORDAL` by chord
    length, `CENTRIPETAL` by its square root (tamer at sharp turns, the
    default). A member equals its lowercase name, which is accepted in
    its place.
    """
    CHORDAL = "chordal"
    CENTRIPETAL = "centripetal"

_Parametrization = Union[PointParametrization, Literal["chordal", "centripetal"]]

@final
class SurfaceCurvatures:
    """Principal, mean and Gaussian curvatures at a surface point, with the
    principal directions, as `curvatures(u, v)` returns them.

    Curvatures are measured along the unit normal: negative bends away
    from it, so a sphere with an outward normal has `-1 / radius`. The two
    principal directions are unit tangents, always perpendicular, even at
    an umbilical point (a sphere, a plane) where every direction is
    principal. Their signs are not specified.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def min_curvature(self) -> float:
        """The smaller principal curvature."""
    @property
    def max_curvature(self) -> float:
        """The larger principal curvature."""
    @property
    def min_direction(self) -> Vector:
        """The unit tangent direction of `min_curvature`."""
    @property
    def max_direction(self) -> Vector:
        """The unit tangent direction of `max_curvature`, perpendicular to
        `min_direction`.
        """
    @property
    def mean(self) -> float:
        """The mean curvature, `(min_curvature + max_curvature) / 2`."""
    @property
    def gaussian(self) -> float:
        """The Gaussian curvature, `min_curvature * max_curvature`."""

@final
class CurveApproximation:
    """A 3D B-spline curve within a deviation of what it approximates, from
    `uv_curve_to_3d_within` and `BSplineCurve.approximate_within`. The
    operation says what the curve is and how its deviation is measured.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve(self) -> BSplineCurve:
        """The curve."""
    @property
    def deviation(self) -> float:
        """The largest distance between the curve and what it approximates,
        as the operation measures it: a bound when `is_proven()`, the
        largest sampled distance otherwise. It is within the asked
        deviation, up to rounding.
        """
    def is_proven(self) -> bool:
        """Whether `deviation` is a proven bound, up to rounding, rather than
        sampled.
        """

@final
class CurveApproximation2D:
    """A 2D B-spline curve within a deviation of what it approximates, from
    `curve_3d_to_uv_within` and `BSplineCurve2D.approximate_within`: as
    `CurveApproximation` is for a 3D curve.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve(self) -> BSplineCurve2D:
        """The curve."""
    @property
    def deviation(self) -> float:
        """The largest distance between the curve and what it approximates,
        as the operation measures it: a bound when `is_proven()`, the
        largest sampled distance otherwise. It is within the asked
        deviation, up to rounding.
        """
    def is_proven(self) -> bool:
        """Whether `deviation` is a proven bound, up to rounding, rather than
        sampled.
        """

@final
class SurfaceApproximation:
    """A B-spline surface within a deviation of what it approximates, from
    `BSplineSurface.approximate_within`: as `CurveApproximation` is for a
    curve.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def surface(self) -> BSplineSurface:
        """The surface."""
    @property
    def deviation(self) -> float:
        """The largest distance between the surface and what it approximates,
        as the operation measures it: a bound when `is_proven()`, the
        largest sampled distance otherwise. It is within the asked
        deviation, up to rounding.
        """
    def is_proven(self) -> bool:
        """Whether `deviation` is a proven bound, up to rounding, rather than
        sampled.
        """

class IntersectionMethod(enum.IntEnum):
    """How an intersection was computed.

    `EXACT` is a closed form, always complete. `NUMERIC` is a search,
    whose `completeness` tells whether it proved that nothing was
    missed.
    """
    EXACT = 0
    NUMERIC = 1

class CompletenessStatus(enum.IntEnum):
    """Whether a numeric search proved that it found everything.

    - `COMPLETE`: nothing was missed (always the case for an exact
      result).
    - `INCOMPLETE`: the search may have missed something inside its
      `regions`, and found everything outside them.
    - `SAMPLED`: the result comes from samples; nothing is proved.
    """
    COMPLETE = 0
    INCOMPLETE = 1
    SAMPLED = 2

class IncompleteReason(enum.IntEnum):
    """Why an incomplete search stopped.

    - `BUDGET_EXHAUSTED`: it ran out of subdivisions.
    - `UNSEPARATED`: pieces of the operands stayed within tolerance of
      each other down to the smallest piece, as along an overlap or at a
      tangent contact of a curve with a surface.
    - `UNBOUNDED`: a branch of two unbounded surfaces runs beyond the
      box it was traced in, and continues past the regions reported.
    """
    BUDGET_EXHAUSTED = 0
    UNSEPARATED = 1
    UNBOUNDED = 2

@final
class CurveCurveRegion:
    """Where an incomplete curve–curve search may have missed something: a
    parameter interval on each curve, and a box in space holding both
    curves over their intervals.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_range(self) -> Tuple[float, float]:
        """The first curve's parameter interval."""
    @property
    def second_range(self) -> Tuple[float, float]:
        """The second curve's parameter interval."""
    @property
    def min_corner(self) -> Point:
        """The corner of the box with the smallest coordinates."""
    @property
    def max_corner(self) -> Point:
        """The corner of the box with the largest coordinates."""

@final
class CurveSurfaceRegion:
    """Where an incomplete curve–surface search may have missed something: a
    parameter interval on the curve, a `(u, v)` rectangle on the surface,
    and a box in space holding both over them.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve_range(self) -> Tuple[float, float]:
        """The curve parameter interval."""
    @property
    def u_range(self) -> Tuple[float, float]:
        """The surface's `u` parameter interval."""
    @property
    def v_range(self) -> Tuple[float, float]:
        """The surface's `v` parameter interval."""
    @property
    def min_corner(self) -> Point:
        """The corner of the box with the smallest coordinates."""
    @property
    def max_corner(self) -> Point:
        """The corner of the box with the largest coordinates."""

@final
class SurfaceSurfaceRegion:
    """Where an incomplete surface–surface search may have missed something:
    a `(u, v)` rectangle on each surface, and a box in space holding both
    over them.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_u_range(self) -> Tuple[float, float]:
        """The first surface's `u` parameter interval."""
    @property
    def first_v_range(self) -> Tuple[float, float]:
        """The first surface's `v` parameter interval."""
    @property
    def second_u_range(self) -> Tuple[float, float]:
        """The second surface's `u` parameter interval."""
    @property
    def second_v_range(self) -> Tuple[float, float]:
        """The second surface's `v` parameter interval."""
    @property
    def min_corner(self) -> Point:
        """The corner of the box with the smallest coordinates."""
    @property
    def max_corner(self) -> Point:
        """The corner of the box with the largest coordinates."""

@final
class Completeness:
    """Whether a result is proved complete and, when it is not, why and where
    something may have been missed.

    `regions` is empty unless `status` is `INCOMPLETE`; `reason` is `None`
    unless it is.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def status(self) -> CompletenessStatus:
        """`COMPLETE`, `INCOMPLETE` or `SAMPLED`."""
    @property
    def reason(self) -> Optional[IncompleteReason]:
        """Why the search stopped short, when `status` is `INCOMPLETE`;
        `None` otherwise.
        """
    @property
    def regions(
        self,
    ) -> Union[List[CurveCurveRegion], List[CurveSurfaceRegion], List[SurfaceSurfaceRegion]]:
        """Where something may have been missed, when `status` is
        `INCOMPLETE`: `CurveCurveRegion`, `CurveSurfaceRegion` or
        `SurfaceSurfaceRegion` objects, by pair family. Empty otherwise.
        """
    def is_complete(self) -> bool:
        """Whether the result is proved to hold everything."""

@final
class CurveCurveHit:
    """Where two 3D curves meet, with the parameter of the point on each curve."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_parameter(self) -> float:
        """The parameter on the first curve."""
    @property
    def second_parameter(self) -> float:
        """The parameter on the second curve."""
    @property
    def point(self) -> Point:
        """The meeting point."""
    @property
    def multiplicity(self) -> int:
        """1 for a crossing, 2 or more where the curves are tangent."""

@final
class CurveCurveHit2D:
    """Where two 2D curves meet, with the parameter of the point on each curve."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_parameter(self) -> float:
        """The parameter on the first curve."""
    @property
    def second_parameter(self) -> float:
        """The parameter on the second curve."""
    @property
    def point(self) -> Point2D:
        """The meeting point."""
    @property
    def multiplicity(self) -> int:
        """1 for a crossing, 2 or more where the curves are tangent."""

@final
class CurveSurfaceHit:
    """Where a curve meets a surface, with the parameter of the point on the
    curve and its `(u, v)` parameters on the surface.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve_parameter(self) -> float:
        """The parameter on the curve."""
    @property
    def u(self) -> float:
        """The `u` parameter on the surface."""
    @property
    def v(self) -> float:
        """The `v` parameter on the surface."""
    @property
    def point(self) -> Point:
        """The meeting point."""
    @property
    def multiplicity(self) -> int:
        """1 for a crossing, 2 or more where the curve is tangent to the
        surface.
        """

@final
class CurveCurveOverlap:
    """A range where two curves coincide: a parameter interval on each.

    A whole-curve coincidence is one overlap over both curves' bounds,
    which may be infinite.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_range(self) -> Tuple[float, float]:
        """The first curve's parameter interval."""
    @property
    def second_range(self) -> Tuple[float, float]:
        """The second curve's parameter interval."""
    @property
    def same_direction(self) -> bool:
        """Whether both curves run the same way along the overlap."""

@final
class CurveSurfaceOverlap:
    """A stretch of a curve that lies on a surface: a parameter interval on
    the curve. "The circle lies on the sphere" is one overlap over
    `(0, 2π)`; a line lying in a plane is one over `(-inf, inf)`.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve_range(self) -> Tuple[float, float]:
        """The curve's parameter interval."""

@final
class IntersectionCurve:
    """A curve where two surfaces meet."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def curve(self) -> Curve:
        """The curve: a line, circle, ellipse, parabola or hyperbola from a
        closed form, a `BSplineCurve` traced by marching otherwise.
        """
    @property
    def first_uv_curve(self) -> Optional[Curve2D]:
        """The curve in the first surface's `(u, v)` parameters, in the same
        parameter as `curve`: the surface's point at its `(u, v)` at `t` is
        `curve`'s point at `t`. A marched curve always has one, fitted from
        the parameters the march walked, as does a curve fitted from a
        closed form (a loop where two cylinders cross). Another exact curve
        has the exact one where a closed form exists, and `None` otherwise
        (an ellipse on a cylinder): `curve_3d_to_uv_within` maps it within
        a deviation.
        """
    @property
    def second_uv_curve(self) -> Optional[Curve2D]:
        """The curve in the second surface's parameters, as `first_uv_curve`
        is on the first.
        """
    @property
    def tangent(self) -> bool:
        """Whether the surfaces touch along the curve instead of crossing it."""
    @property
    def deviation(self) -> float:
        """How far `curve` strays from the surfaces through its curves in
        parameters: the largest distance between `curve`'s point at `t` and
        each surface's point at its UV curve's point at `t`, the tolerance an
        edge along the curve needs. A marched curve samples it at the middle
        of each knot span, where its fits stray furthest from the traced
        points; an exact curve is `0.0`.
        """
    @property
    def parameter_range(self) -> Tuple[float, float]:
        """The `(first, last)` range of `curve`'s parameter, and of the curves
        in parameters, along which the surfaces meet: `curve`'s own bounds
        (infinite for a line), or the piece of it inside the ranges of a
        `TrimmedSurface`. A periodic `curve` cut by the ranges keeps a piece
        across its seam whole, its range running past the period.
        """

@final
class CurveCurveIntersection:
    """The intersection of two 3D curves: the points where they meet and the
    ranges where they coincide.

    `first` and `second` are the operands, in the order of the call.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[CurveCurveHit]: ...
    def __reversed__(self) -> Iterator[CurveCurveHit]:
        """Iterates over `hits` from the last to the first."""
    @overload
    def __getitem__(self, index: int, /) -> CurveCurveHit: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[CurveCurveHit]: ...
    @property
    def points(self) -> List[Point]:
        """The point of each hit, in the order of `hits`."""
    @property
    def hits(self) -> List[CurveCurveHit]:
        """The points where the curves meet, by first-curve parameter."""
    @property
    def overlaps(self) -> List[CurveCurveOverlap]:
        """The ranges where the curves coincide."""
    @property
    def method(self) -> IntersectionMethod:
        """`EXACT` for a closed form, `NUMERIC` for a search or marching."""
    @property
    def completeness(self) -> Completeness:
        """Whether the result is proved to hold everything: exact
        results are complete, numeric ones report the search status.
        """
    @property
    def first(self) -> Curve | TrimmedCurve:
        """The first operand."""
    @property
    def second(self) -> Curve | TrimmedCurve:
        """The second operand."""
    def is_complete(self) -> bool:
        """Whether the result is proved to hold everything."""
    def is_empty(self) -> bool:
        """Whether the result holds nothing. An empty result that is not
        complete means "found nothing, may have missed something".
        """
    def __bool__(self) -> bool: ...

@final
class CurveCurveIntersection2D:
    """The intersection of two 2D curves: exact and complete for every pair of
    lines and conics, numeric when one is a B-spline.

    `first` and `second` are the operands, in the order of the call.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[CurveCurveHit2D]: ...
    def __reversed__(self) -> Iterator[CurveCurveHit2D]:
        """Iterates over `hits` from the last to the first."""
    @overload
    def __getitem__(self, index: int, /) -> CurveCurveHit2D: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[CurveCurveHit2D]: ...
    @property
    def points(self) -> List[Point2D]:
        """The point of each hit, in the order of `hits`."""
    @property
    def hits(self) -> List[CurveCurveHit2D]:
        """The points where the curves meet, by first-curve parameter."""
    @property
    def overlaps(self) -> List[CurveCurveOverlap]:
        """The ranges where the curves coincide."""
    @property
    def method(self) -> IntersectionMethod:
        """`EXACT` for a closed form, `NUMERIC` for a search or marching."""
    @property
    def completeness(self) -> Completeness:
        """Whether the result is proved to hold everything: exact
        results are complete, numeric ones report the search status.
        """
    @property
    def first(self) -> Curve2D | TrimmedCurve2D:
        """The first operand."""
    @property
    def second(self) -> Curve2D | TrimmedCurve2D:
        """The second operand."""
    def is_complete(self) -> bool:
        """Whether the result is proved to hold everything."""
    def is_empty(self) -> bool:
        """Whether the result holds nothing. An empty result that is not
        complete means "found nothing, may have missed something".
        """
    def __bool__(self) -> bool: ...

@final
class CurveSurfaceIntersection:
    """The intersection of a curve and a surface: the points where they meet
    and the curve's parameter ranges that lie on the surface.

    `first` and `second` are the operands, in the order of the call;
    whichever order that was, the hit parameters keep their roles, and
    `curve` and `surface` name the operands by role.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[CurveSurfaceHit]: ...
    def __reversed__(self) -> Iterator[CurveSurfaceHit]:
        """Iterates over `hits` from the last to the first."""
    @overload
    def __getitem__(self, index: int, /) -> CurveSurfaceHit: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[CurveSurfaceHit]: ...
    @property
    def points(self) -> List[Point]:
        """The point of each hit, in the order of `hits`."""
    @property
    def hits(self) -> List[CurveSurfaceHit]:
        """The points where the curve meets the surface, by curve parameter."""
    @property
    def overlaps(self) -> List[CurveSurfaceOverlap]:
        """The stretches of the curve that lie on the surface."""
    @property
    def method(self) -> IntersectionMethod:
        """`EXACT` for a closed form, `NUMERIC` for a search or marching."""
    @property
    def completeness(self) -> Completeness:
        """Whether the result is proved to hold everything: exact
        results are complete, numeric ones report the search status.
        """
    @property
    def first(self) -> Union[Curve, TrimmedCurve, Surface]:
        """The first operand."""
    @property
    def second(self) -> Union[Curve, TrimmedCurve, Surface]:
        """The second operand."""
    @property
    def curve(self) -> Curve | TrimmedCurve:
        """The curve operand, whichever of `first` and `second` it is."""
    @property
    def surface(self) -> Surface:
        """The surface operand, whichever of `first` and `second` it is."""
    def is_complete(self) -> bool:
        """Whether the result is proved to hold everything."""
    def is_empty(self) -> bool:
        """Whether the result holds nothing. An empty result that is not
        complete means "found nothing, may have missed something".
        """
    def __bool__(self) -> bool: ...

@final
class SurfaceSurfaceIntersection:
    """The intersection of two surfaces: the curves where they meet, the
    isolated points where they touch, and whether they coincide.

    `first` and `second` are the operands, in the order of the call.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[IntersectionCurve]: ...
    def __reversed__(self) -> Iterator[IntersectionCurve]:
        """Iterates over `curves` from the last to the first."""
    @overload
    def __getitem__(self, index: int, /) -> IntersectionCurve: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[IntersectionCurve]: ...
    @property
    def curves(self) -> List[IntersectionCurve]:
        """The curves where the surfaces meet."""
    @property
    def contact_points(self) -> List[Point]:
        """The points where the surfaces touch, their normals parallel, with
        no curve of contact through them: isolated points, and the points
        where curves cross.
        """
    @property
    def method(self) -> IntersectionMethod:
        """`EXACT` for a closed form, `NUMERIC` for a search or marching."""
    @property
    def completeness(self) -> Completeness:
        """Whether the result is proved to hold everything: exact
        results are complete, numeric ones report the search status.
        """
    @property
    def first(self) -> Surface:
        """The first operand."""
    @property
    def second(self) -> Surface:
        """The second operand."""
    def is_coincident(self) -> bool:
        """Whether the surfaces are the same surface."""
    def is_complete(self) -> bool:
        """Whether the result is proved to hold everything."""
    def is_empty(self) -> bool:
        """Whether the result holds nothing. An empty result that is not
        complete means "found nothing, may have missed something".
        """
    def __bool__(self) -> bool: ...

@final
class CurveProjection:
    """Where a 3D curve is closest to a point.

    `target` is the curve that was projected onto, as the exact object the
    method was called on, and `query_point` is the point that was projected.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def parameter(self) -> float:
        """The parameter of the closest point."""
    @property
    def point(self) -> Point:
        """The closest point: the curve evaluated at `parameter`."""
    @property
    def distance(self) -> float:
        """The distance from `query_point` to the curve."""
    @property
    def target(self) -> Curve | TrimmedCurve:
        """The curve that was projected onto."""
    @property
    def query_point(self) -> Point:
        """The point that was projected."""

@final
class CurveProjection2D:
    """Where a 2D curve is closest to a point.

    `target` is the curve that was projected onto, as the exact object the
    method was called on, and `query_point` is the point that was projected.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def parameter(self) -> float:
        """The parameter of the closest point."""
    @property
    def point(self) -> Point2D:
        """The closest point: the curve evaluated at `parameter`."""
    @property
    def distance(self) -> float:
        """The distance from `query_point` to the curve."""
    @property
    def target(self) -> Curve2D | TrimmedCurve2D:
        """The curve that was projected onto."""
    @property
    def query_point(self) -> Point2D:
        """The point that was projected."""

@final
class SurfaceProjection:
    """Where a surface is closest to a point.

    `target` is the surface that was projected onto, as the exact object the
    method was called on, and `query_point` is the point that was projected.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def u(self) -> float:
        """The `u` parameter of the closest point."""
    @property
    def v(self) -> float:
        """The `v` parameter of the closest point."""
    @property
    def point(self) -> Point:
        """The closest point: the surface evaluated at `(u, v)`."""
    @property
    def distance(self) -> float:
        """The distance from `query_point` to the surface."""
    @property
    def signed_distance(self) -> Optional[float]:
        """The distance with a side: positive when `query_point` lies on the
        side the surface normal at `point` points to, negative on the
        other, `0.0` on the surface. The sign follows the surface's own
        orientation, so `reversed_u()` or `reversed_v()` flips it. `None`
        where no side is defined: the normal at `point` is undefined, as at
        a cone's apex, or `query_point` lies beside the surface, in its
        tangent plane past a boundary.
        """
    @property
    def target(self) -> Surface:
        """The surface that was projected onto."""
    @property
    def query_point(self) -> Point:
        """The point that was projected."""

@final
class CurveCurveExtremum:
    """A local minimum of the distance between two curves."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @property
    def first_parameter(self) -> float:
        """The parameter on the first curve."""
    @property
    def second_parameter(self) -> float:
        """The parameter on the second curve."""
    @property
    def first_point(self) -> Point:
        """The point of the first curve at `first_parameter`."""
    @property
    def second_point(self) -> Point:
        """The point of the second curve at `second_parameter`."""
    @property
    def distance(self) -> float:
        """The distance between `first_point` and `second_point`."""

@final
class CurveCurveExtrema:
    """The local minima of the distance between two curves, closest first,
    with whether the search proved there is no other.

    `first` and `second` are the curves, in the order of the call.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __len__(self) -> int: ...
    def closest(self) -> Optional[CurveCurveExtremum]:
        """The closest minimum, the first of `minima`, or `None` when there
        is none.
        """
    def __iter__(self) -> Iterator[CurveCurveExtremum]: ...
    def __reversed__(self) -> Iterator[CurveCurveExtremum]:
        """Iterates over the minima, farthest first."""
    @overload
    def __getitem__(self, index: int, /) -> CurveCurveExtremum: ...
    @overload
    def __getitem__(self, index: slice, /) -> List[CurveCurveExtremum]: ...
    @property
    def minima(self) -> List[CurveCurveExtremum]:
        """The minima, closest first."""
    @property
    def completeness(self) -> Completeness:
        """Whether the search proved there is no other minimum."""
    @property
    def first(self) -> Curve | TrimmedCurve:
        """The first curve."""
    @property
    def second(self) -> Curve | TrimmedCurve:
        """The second curve."""
    def is_complete(self) -> bool:
        """Whether the search proved there is no other minimum."""

@deprecated("use surface.curve_3d_to_uv_within(curve, max_deviation, tol)")
def curve_3d_to_uv_within(
    surface: _SurfaceMethods[Any, Any],
    curve: _CurveMethods,
    max_deviation: float,
    tol: Optional[_ToleranceLike] = None,
) -> CurveApproximation2D:
    """Deprecated: use `surface.curve_3d_to_uv_within(curve, max_deviation,
    tol)`.
    """
@deprecated("use surface.uv_curve_to_3d_within(uv_curve, max_deviation, tol)")
def uv_curve_to_3d_within(
    surface: _SurfaceMethods[Any, Any],
    uv_curve: _Curve2DMethods,
    max_deviation: float,
    tol: Optional[_ToleranceLike] = None,
) -> CurveApproximation:
    """Deprecated: use `surface.uv_curve_to_3d_within(uv_curve,
    max_deviation, tol)`.
    """
@deprecated("use first.extrema(second, tol)")
def curve_curve_extrema(
    first: _CurveMethods,
    second: _CurveMethods,
    tol: Optional[_ToleranceLike] = None,
) -> CurveCurveExtrema:
    """Deprecated: use `first.extrema(second, tol)`."""

# ---------------------------------------------------------------- curves

# Every class of each kind of shape, as the union to annotate with. A
# curve is one of the classes of the `Curve` enum's variants, as in Rust;
# a `TrimmedCurve` is a piece of one, so a method that takes either says
# `Curve | TrimmedCurve`.
#
# A parameter that takes any curve or surface is annotated with the base
# class every one of them inherits (`_CurveMethods`, `_Curve2DMethods`,
# `_SurfaceMethods`) rather than the union, so a list mixing a circle and
# a line, which mypy types as their common base, can be passed along.
Curve = (
    Line | Circle | Ellipse | Parabola | Hyperbola | BSplineCurve | BezierCurve | QuarticBranch
)
Curve2D = (
    Line2D | Circle2D | Ellipse2D | Parabola2D | Hyperbola2D | BSplineCurve2D | BezierCurve2D
)

_PointType = TypeVar("_PointType")
# A point argument of either dimension: the point, or its coordinates.
_PointArgument = Union[_PointType, Sequence[float], numpy.typing.NDArray[Any]]
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
    def points_at(self, us: _Floats) -> List[_PointType]:
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
    def tangent(self, u: float, tol: Optional[_ToleranceLike] = None) -> Optional[_VectorType]:
        """Unit tangent at parameter `u` (`Tolerance.default()` when
        omitted), or `None` where the speed vanishes, as at a cusp.
        The speed counts as zero when the distance it would cover
        across the parameter range is within `tol.confusion`, so
        rescaling the parameter does not change the answer.
        """
    def curvature(self, u: float, tol: Optional[_ToleranceLike] = None) -> Optional[float]:
        """Curvature at parameter `u`: the turn rate of the unit tangent
        per unit arc length, zero on a straight run. `None` where the
        speed vanishes, by the test of `tangent`.
        """
    def curvature_center(
        self, u: float, tol: Optional[_ToleranceLike] = None
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
    def is_closed(self, tol: Optional[_ToleranceLike] = None) -> bool:
        """Whether travelling along the curve returns to its start
        (`Tolerance.default()` when omitted). A periodic curve is
        closed; a clamped B-spline whose ends coincide is closed
        without being periodic.
        """
    def continuity(self, u: float, tol: Optional[_ToleranceLike] = None) -> Union[int, float]:
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
        self, point: _PointArgument[_PointType], tol: Optional[_ToleranceLike] = None
    ) -> _ProjectionType:
        """Projects `point` onto the curve: the closest point, as a
        projection object.

        A distance within `tol.confusion` (`1e-7` by default) snaps
        to `0.0`, so projection agrees with `contains`; pass a
        smaller `Tolerance` to read a smaller distance.
        """
    def project_points(
        self, points: Union[Sequence[_PointArgument[_PointType]], numpy.typing.NDArray[Any]], tol: Optional[_ToleranceLike] = None
    ) -> List[_ProjectionType]:
        """Projects each point in `points` onto the curve, returning a
        projection object per point.
        """
    def project_points_array(
        self, points: numpy.typing.ArrayLike, tol: Optional[_ToleranceLike] = None
    ) -> Tuple[_FloatArray, _FloatArray, _FloatArray]:
        """Projects each row of the array-like `points` onto the curve,
        returning numpy arrays `(parameters, points, distances)`:
        the parameter, the closest point and its distance per row.
        """
    def contains(self, point: _PointArgument[_PointType], tol: Optional[_ToleranceLike] = None) -> bool:
        """Returns whether `point` lies on the curve within `tol`
        (`Tolerance.default()` when omitted).
        """
    def parameter_of(self, point: _PointArgument[_PointType], tol: Optional[_ToleranceLike] = None) -> float:
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
        self, tol: Optional[_ToleranceLike] = None
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
        self, other: _CurveMethods, tol: Optional[_ToleranceLike] = None
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
        self, other: _SurfaceMethods[Any, Any], tol: Optional[_ToleranceLike] = None
    ) -> CurveSurfaceIntersection: ...
    @overload
    def extrema(
        self, other: _PointLike, tol: Optional[_ToleranceLike] = None
    ) -> List[CurveProjection]:
        """The extrema of the distance to `other`, a point or a curve.

        From a point: the local minima of its distance to the
        curve, an end of a bounded curve among them where the
        distance grows from it into the curve, as `CurveProjection`
        objects ordered by ascending distance.

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
        self, other: _CurveMethods, tol: Optional[_ToleranceLike] = None
    ) -> CurveCurveExtrema: ...

class _Curve2DMethods(
    _CurveMethodsInEitherDimension[Point2D, Vector2D, BoundingBox2D, CurveProjection2D]
):
    def intersect(
        self, other: _Curve2DMethods, tol: Optional[_ToleranceLike] = None
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
        self, point: _Point2DLike, tol: Optional[_ToleranceLike] = None
    ) -> List[CurveProjection2D]:
        """The local minima of the distance from `point` to the
        curve, an end of a bounded curve among them where the
        distance grows from it into the curve, as
        `CurveProjection2D` objects ordered by ascending distance.
        `point` is a `Point2D` or a sequence of its two
        coordinates. Raises `TypeError` for anything else: the
        extrema between two 2D curves are not supported.
        """

class _BoundedCurveMethods(Generic[_PointType, _DiscretizationType]):
    def start_point(self) -> _PointType:
        """The point at the first bound: on a circle or an ellipse, the
        seam, where the end point is too.
        """
    def end_point(self) -> _PointType:
        """The point at the last bound."""
    def discretize(
        self, max_chord_deviation: float, max_angle_deviation: Optional[float] = None
    ) -> _DiscretizationType:
        """Discretize the curve into points joined by straight
        segments, each segment's chord deviation at most
        `max_chord_deviation` and, if given, its angle deviation at
        most `max_angle_deviation` (radians).
        """
    def point_count(
        self, max_chord_deviation: float, max_angle_deviation: Optional[float] = None
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
    """An infinite straight line in 3D, parametrized by arc length from its origin."""
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, origin: _PointLike, direction: Vector) -> Line: ...
    @staticmethod
    def from_axis(axis: Axis) -> Line:
        """Build a line from an axis."""
    @staticmethod
    def from_two_points(p1: _PointLike, p2: _PointLike) -> Line:
        """Build a line through two distinct points."""
    def origin(self) -> Point:
        """The line origin (point at parameter 0)."""
    def direction(self) -> Vector:
        """The unit direction of the line."""
    def axis(self) -> Axis:
        """The axis of the line: its origin and unit direction."""

@final
class Circle(
    _CurveMethods,
    _BoundedCurveMethods[Point, Discretization],
):
    """A circle in 3D, parametrized by the angle from its frame's x direction."""
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, center: _PointLike, normal: Vector, radius: float) -> Circle: ...
    @staticmethod
    def from_axis(axis: Axis, radius: float) -> Circle:
        """Build a circle from an axis (center + normal) and a radius."""
    @staticmethod
    def from_frame(frame: Frame, radius: float) -> Circle:
        """Build a circle from a placement frame and a radius."""
    @staticmethod
    def from_three_points(p1: _PointLike, p2: _PointLike, p3: _PointLike) -> Circle:
        """Build the circle through three non-collinear points; the curve starts
        at the first point.
        """
    def center(self) -> Point:
        """The circle center."""
    def radius(self) -> float:
        """The circle radius."""
    def length(self) -> float:
        """The length of the circle, its circumference (`2*PI*radius`)."""
    def enclosed_area(self) -> float:
        """The area the circle encloses: its disk's (`PI*radius^2`)."""
    def normal(self) -> Vector:
        """The unit normal of the circle's plane, on the side the circle turns
        counterclockwise about.
        """
    def axis(self) -> Axis:
        """The axis of the circle: through the center along the normal."""
    def to_bspline(self) -> BSplineCurve:
        """The circle as an exact periodic rational quadratic B-spline over
        `[0, 2π]`.
        Its rational arcs do not keep the angle as parameter: the same
        parameter gives a different point on each curve between the knots,
        so map a point across with `bspline.parameter_of(point)`.
        """
    def frame(self) -> Frame:
        """The placement frame of the circle."""

@final
class Ellipse(
    _CurveMethods,
    _BoundedCurveMethods[Point, Discretization],
):
    """An ellipse in 3D, evaluated as
    `center + major_radius*cos(t)*x + minor_radius*sin(t)*y`. The major
    axis `x` is `x_direction` projected onto the plane perpendicular to
    `normal`, and `y = normal × x`, so `t` turns counterclockwise seen from
    the side `normal` points to. `major_radius` must be at least
    `minor_radius` (`ConstructionError` otherwise): to lay the long axis
    along another direction, turn `x_direction` rather than swap the radii.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        center: _PointLike,
        normal: Vector,
        x_direction: Vector,
        major_radius: float,
        minor_radius: float,
    ) -> Ellipse: ...
    @staticmethod
    def from_frame(frame: Frame, major_radius: float, minor_radius: float) -> Ellipse:
        """Build an ellipse from a placement frame and the two radii."""
    @staticmethod
    def from_center_and_points(center: _PointLike, s1: _PointLike, s2: _PointLike) -> Ellipse:
        """Build an ellipse from its center, the major-axis end point `s1`
        and a point `s2` whose distance from the major axis is the minor
        radius.
        """
    def center(self) -> Point:
        """The ellipse center."""
    def major_radius(self) -> float:
        """The semi-major radius."""
    def minor_radius(self) -> float:
        """The semi-minor radius."""
    def enclosed_area(self) -> float:
        """The enclosed area (`PI*major*minor`)."""
    def length(self, max_length_error: float = 1e-07) -> float:
        """The perimeter, integrating the speed to within `max_length_error`: it has
        no closed form.
        `max_length_error` defaults to the default confusion, 1e-7.
        """
    def eccentricity(self) -> float:
        """The eccentricity: the ratio of a point's distance to a focus over
        its distance to that focus's directrix.
        """
    def foci(self) -> Tuple[Point, Point]:
        """The two foci on the major axis, the one toward the frame's x
        direction first.
        """
    def directrices(self) -> Optional[Tuple[Line, Line]]:
        """The two directrices across the major axis, in the order of `foci`,
        or None for a circle, whose directrices are at infinity.
        """
    def to_bspline(self) -> BSplineCurve:
        """The ellipse as an exact periodic rational quadratic B-spline over
        `[0, 2π]`.
        Its rational arcs do not keep the angle as parameter: the same
        parameter gives a different point on each curve between the knots,
        so map a point across with `bspline.parameter_of(point)`.
        """
    def frame(self) -> Frame:
        """The placement frame of the ellipse."""

@final
class Parabola(_CurveMethods):
    """A parabola in 3D, parametrized so the apex is at parameter 0."""
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, apex: _PointLike, normal: Vector, x_direction: Vector, focal_distance: float
    ) -> Parabola: ...
    @staticmethod
    def from_frame(frame: Frame, focal_distance: float) -> Parabola:
        """Build a parabola from a placement frame and focal distance."""
    @staticmethod
    def from_focus_and_directrix(focus: _PointLike, directrix: Line) -> Parabola:
        """Build the parabola of every point as far from `focus` as from
        `directrix`; the parameter runs along the directrix's direction.
        """
    def focal_distance(self) -> float:
        """The focal distance."""
    def eccentricity(self) -> float:
        """The eccentricity: the ratio of a point's distance to a focus over
        its distance to that focus's directrix.
        """
    def focus(self) -> Point:
        """The focus, `focal_distance` from the apex along the axis."""
    def directrix(self) -> Line:
        """The directrix, across the axis `focal_distance` behind the apex:
        every point of the curve is as far from it as from the focus.
        """
    def frame(self) -> Frame:
        """The placement frame of the parabola."""
    def apex(self) -> Point:
        """The apex of the parabola: its point of smallest radius, at the frame
        origin.
        """

@final
class Hyperbola(_CurveMethods):
    """One branch of a hyperbola in 3D."""
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        center: _PointLike,
        normal: Vector,
        x_direction: Vector,
        major_radius: float,
        minor_radius: float,
    ) -> Hyperbola: ...
    @staticmethod
    def from_frame(frame: Frame, major_radius: float, minor_radius: float) -> Hyperbola:
        """Build a hyperbola from a placement frame and the two radii."""
    @staticmethod
    def from_center_and_points(center: _PointLike, s1: _PointLike, s2: _PointLike) -> Hyperbola:
        """Build a hyperbola from its center, the vertex point and a second
        point on the curve.
        """
    def center(self) -> Point:
        """The hyperbola center."""
    def major_radius(self) -> float:
        """The semi-major (transverse) radius."""
    def minor_radius(self) -> float:
        """The semi-minor (conjugate) radius."""
    def eccentricity(self) -> float:
        """The eccentricity: the ratio of a point's distance to a focus over
        its distance to that focus's directrix.
        """
    def foci(self) -> Tuple[Point, Point]:
        """The two foci on the major axis, the one toward the frame's x
        direction first.
        """
    def directrices(self) -> Tuple[Line, Line]:
        """The two directrices across the major axis, in the order of `foci`."""
    def other_branch(self) -> Hyperbola:
        """The other branch: the same hyperbola mirrored through its center,
        with the same parametrization direction around the center.
        """
    def frame(self) -> Frame:
        """The placement frame of the hyperbola."""

@final
class BSplineCurve(
    _CurveMethods,
    _BoundedCurveMethods[Point, Discretization],
    _SidedEvaluationMethods[Point, Vector],
):
    """A B-spline (optionally rational, optionally periodic) curve in 3D."""
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        degree: int,
        poles: _Points,
        knots: _Floats,
        multiplicities: _Ints,
        periodic: bool,
    ) -> BSplineCurve: ...
    @staticmethod
    def new_rational(
        degree: int,
        poles: _Points,
        weights: _Floats,
        knots: _Floats,
        multiplicities: _Ints,
        periodic: bool,
    ) -> BSplineCurve:
        """Build a rational B-spline (NURBS) curve with a weight per pole."""
    @staticmethod
    def interpolate(
        points: _Points,
        degree: int,
        parametrization: _Parametrization = "centripetal",
    ) -> BSplineCurve:
        """Interpolate a clamped B-spline through `points`.

        `points` is a list of points or an `(N, 3)` array-like.

        `parametrization` is a `PointParametrization` or its string value,
        `"centripetal"` (the default, tamer at sharp turns) or `"chordal"`.
        """
    @staticmethod
    def approximate(
        points: _Points,
        degree: int,
        pole_count: int,
        parametrization: _Parametrization = "centripetal",
    ) -> BSplineCurve:
        """Least-squares B-spline fit through `points` with `pole_count` poles.

        `points` is a list of points or an `(N, 3)` array-like.

        `parametrization` is a `PointParametrization` or its string value,
        `"centripetal"` (the default) or `"chordal"`.
        """
    @staticmethod
    def approximate_within(
        points: _Points,
        degree: int,
        max_deviation: float,
        parametrization: _Parametrization = "centripetal",
    ) -> CurveApproximation:
        """Least-squares B-spline fit through `points` with as few poles as
        keep every point within `max_deviation` of the curve at its
        parameter.

        Returns the curve with its `deviation`, the largest distance
        between a point and the curve at its parameter, proven.

        `points` is a list of points or an `(N, 3)` array-like.

        `parametrization` is a `PointParametrization` or its string value,
        `"centripetal"` (the default) or `"chordal"`.
        """
    @staticmethod
    def from_flat_knots(
        degree: int,
        poles: _Points,
        flat_knots: _Floats,
        weights: Optional[_Floats] = None,
        periodic: bool = False,
    ) -> BSplineCurve:
        """Build a curve from its flat knot vector, the form textbooks use:
        each distinct knot repeated by its multiplicity, so a non-periodic
        curve of `n` poles takes `n + degree + 1` knots. `weights`, one per
        pole, make it rational. The flat knots of a periodic curve are as
        `flat_knots()` gives them, with no padding beyond the period.
        """
    @staticmethod
    def clamped_uniform(
        degree: int, poles: _Points
    ) -> BSplineCurve:
        """Build the clamped curve of `degree` with uniform knots over
        `[0, 1]`: both end knots of multiplicity `degree + 1`, and
        `len(poles) - degree - 1` simple knots between them, evenly spaced.
        The curve starts at the first pole and ends at the last.
        """
    def degree(self) -> int:
        """The polynomial degree."""
    def is_rational(self) -> bool:
        """Whether the curve carries rational weights."""
    def poles(self) -> List[Point]:
        """The control poles."""
    def weights(self) -> Optional[List[float]]:
        """The per-pole weights, or `None` when the curve is not rational."""
    def knots(self) -> List[float]:
        """The distinct knots."""
    def multiplicities(self) -> List[int]:
        """The multiplicity of each knot."""
    def flat_knots(self) -> List[float]:
        """The flat knot vector: each distinct knot repeated by its
        multiplicity, `len(poles) + degree + 1` knots for a non-periodic
        curve.
        """
    def span_index(self, u: float) -> int:
        """The index `s` of the knot span holding `u`, counted on
        `flat_knots()`: `flat[s] <= u < flat[s + 1]` (`<=` at the upper
        bound). As in evaluation, `u` past the bounds is clamped, or
        wrapped into the period of a periodic curve.
        """
    def basis_functions(self, u: float) -> Tuple[int, List[float]]:
        """The `degree + 1` basis functions non-zero at `u`, as
        `(first, values)`: the point at `u` of a non-rational curve is
        the sum of `values[j]` times pole `first + j` (wrapping past the
        last pole of a periodic curve). They are the polynomial basis,
        before any weights, and sum to 1.
        """
    def with_pole(self, index: int, pole: _PointLike) -> BSplineCurve:
        """The curve with pole `index` (negative counting from the end) moved
        to `pole`, its degree, knots and weights kept.
        """
    def with_weight(self, index: int, weight: float) -> BSplineCurve:
        """The curve with the weight of pole `index` (negative counting from
        the end) set to `weight`, a positive finite number. A non-rational
        curve becomes rational, every other pole weighing 1.
        """
    def geometric_continuity(
        self, u: float, tol: Optional[_ToleranceLike] = None
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
    def insert_knot(self, u: float, times: int) -> BSplineCurve:
        """The curve with the knot `u` inserted `times` times (same shape,
        more poles).
        """
    def refine_knots(self, knots: _Floats) -> BSplineCurve:
        """The curve with each of `knots` inserted once (same shape, more
        poles).
        """
    def elevate_degree(self, amount: int) -> BSplineCurve:
        """The curve of degree raised by `amount` (same shape)."""
    def remove_knot(
        self, u: float, times: int, max_deviation: float
    ) -> Tuple[BSplineCurve, int]:
        """The curve with the knot `u` removed up to `times` times, and how
        many copies were removed: a removal is kept only while every
        point stays within `max_deviation` of this one.
        """
    def remove_knots(
        self, knots: _Floats, max_deviation: float
    ) -> Tuple[BSplineCurve, int]:
        """The curve with each of `knots` removed once, in order, and how
        many copies were removed in all.
        """
    def to_bezier_pieces(self) -> List[BezierCurve]:
        """The polynomial (or rational) piece of each knot span as a
        `BezierCurve`, in order, each keeping its span as its bounds.
        """
    def to_bspline(self) -> BSplineCurve:
        """The curve itself, as the other curves give their exact B-spline
        with `to_bspline()`.
        """
    def extract(self, start: float, end: float) -> BSplineCurve:
        """The piece over `[start, end]` as a clamped, non-periodic B-spline
        (the range may cross the seam of a periodic curve, up to one
        period).
        """
    def split_at(self, u: float) -> Tuple[BSplineCurve, BSplineCurve]:
        """The two halves split at `u`, each a clamped, non-periodic B-spline."""
    def length(self, max_length_error: float = 1e-07) -> float:
        """The arc length, integrating the speed to within `max_length_error`.
        `max_length_error` defaults to the default confusion, 1e-7.
        """

@final
class BSplineCurve2D(
    _Curve2DMethods,
    _BoundedCurveMethods[Point2D, Discretization2D],
    _SidedEvaluationMethods[Point2D, Vector2D],
):
    """A B-spline (optionally rational, optionally periodic) curve in 2D: the
    planar counterpart of `BSplineCurve`.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        degree: int,
        poles: _Points2D,
        knots: _Floats,
        multiplicities: _Ints,
        periodic: bool,
    ) -> BSplineCurve2D: ...
    @staticmethod
    def new_rational(
        degree: int,
        poles: _Points2D,
        weights: _Floats,
        knots: _Floats,
        multiplicities: _Ints,
        periodic: bool,
    ) -> BSplineCurve2D:
        """Build a rational 2D B-spline (NURBS) curve with a weight per pole."""
    @staticmethod
    def interpolate(
        points: _Points2D,
        degree: int,
        parametrization: _Parametrization = "centripetal",
    ) -> BSplineCurve2D:
        """Interpolate a clamped 2D B-spline through `points`.

        `points` is a list of points or an `(N, 2)` array-like.

        `parametrization` is a `PointParametrization` or its string value,
        `"centripetal"` (the default, tamer at sharp turns) or `"chordal"`.
        """
    @staticmethod
    def approximate(
        points: _Points2D,
        degree: int,
        pole_count: int,
        parametrization: _Parametrization = "centripetal",
    ) -> BSplineCurve2D:
        """Least-squares 2D B-spline fit through `points` with `pole_count`
        poles.

        `points` is a list of points or an `(N, 2)` array-like.

        `parametrization` is a `PointParametrization` or its string value,
        `"centripetal"` (the default) or `"chordal"`.
        """
    @staticmethod
    def approximate_within(
        points: _Points2D,
        degree: int,
        max_deviation: float,
        parametrization: _Parametrization = "centripetal",
    ) -> CurveApproximation2D:
        """Least-squares 2D B-spline fit through `points` with as few poles as
        keep every point within `max_deviation` of the curve at its
        parameter.

        Returns the curve with its `deviation`, the largest distance
        between a point and the curve at its parameter, proven.

        `points` is a list of points or an `(N, 2)` array-like.

        `parametrization` is a `PointParametrization` or its string value,
        `"centripetal"` (the default) or `"chordal"`.
        """
    @staticmethod
    def from_flat_knots(
        degree: int,
        poles: _Points2D,
        flat_knots: _Floats,
        weights: Optional[_Floats] = None,
        periodic: bool = False,
    ) -> BSplineCurve2D:
        """Build a curve from its flat knot vector, the form textbooks use:
        each distinct knot repeated by its multiplicity, so a non-periodic
        curve of `n` poles takes `n + degree + 1` knots. `weights`, one per
        pole, make it rational. The flat knots of a periodic curve are as
        `flat_knots()` gives them, with no padding beyond the period.
        """
    @staticmethod
    def clamped_uniform(
        degree: int, poles: _Points2D
    ) -> BSplineCurve2D:
        """Build the clamped curve of `degree` with uniform knots over
        `[0, 1]`: both end knots of multiplicity `degree + 1`, and
        `len(poles) - degree - 1` simple knots between them, evenly spaced.
        The curve starts at the first pole and ends at the last.
        """
    def degree(self) -> int:
        """The polynomial degree."""
    def is_rational(self) -> bool:
        """Whether the curve carries rational weights."""
    def poles(self) -> List[Point2D]:
        """The control poles."""
    def weights(self) -> Optional[List[float]]:
        """The per-pole weights, or `None` when the curve is not rational."""
    def knots(self) -> List[float]:
        """The distinct knots."""
    def multiplicities(self) -> List[int]:
        """The multiplicity of each knot."""
    def flat_knots(self) -> List[float]:
        """The flat knot vector: each distinct knot repeated by its
        multiplicity, `len(poles) + degree + 1` knots for a non-periodic
        curve.
        """
    def span_index(self, u: float) -> int:
        """The index `s` of the knot span holding `u`, counted on
        `flat_knots()`: `flat[s] <= u < flat[s + 1]` (`<=` at the upper
        bound). As in evaluation, `u` past the bounds is clamped, or
        wrapped into the period of a periodic curve.
        """
    def basis_functions(self, u: float) -> Tuple[int, List[float]]:
        """The `degree + 1` basis functions non-zero at `u`, as
        `(first, values)`: the point at `u` of a non-rational curve is
        the sum of `values[j]` times pole `first + j` (wrapping past the
        last pole of a periodic curve). They are the polynomial basis,
        before any weights, and sum to 1.
        """
    def with_pole(self, index: int, pole: _Point2DLike) -> BSplineCurve2D:
        """The curve with pole `index` (negative counting from the end) moved
        to `pole`, its degree, knots and weights kept.
        """
    def with_weight(self, index: int, weight: float) -> BSplineCurve2D:
        """The curve with the weight of pole `index` (negative counting from
        the end) set to `weight`, a positive finite number. A non-rational
        curve becomes rational, every other pole weighing 1.
        """
    def geometric_continuity(
        self, u: float, tol: Optional[_ToleranceLike] = None
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
    def insert_knot(self, u: float, times: int) -> BSplineCurve2D:
        """The curve with the knot `u` inserted `times` times (same shape,
        more poles).
        """
    def refine_knots(self, knots: _Floats) -> BSplineCurve2D:
        """The curve with each of `knots` inserted once (same shape, more
        poles).
        """
    def elevate_degree(self, amount: int) -> BSplineCurve2D:
        """The curve of degree raised by `amount` (same shape)."""
    def remove_knot(
        self, u: float, times: int, max_deviation: float
    ) -> Tuple[BSplineCurve2D, int]:
        """The curve with the knot `u` removed up to `times` times, and how
        many copies were removed: a removal is kept only while every
        point stays within `max_deviation` of this one.
        """
    def remove_knots(
        self, knots: _Floats, max_deviation: float
    ) -> Tuple[BSplineCurve2D, int]:
        """The curve with each of `knots` removed once, in order, and how
        many copies were removed in all.
        """
    def to_bezier_pieces(self) -> List[BezierCurve2D]:
        """The polynomial (or rational) piece of each knot span as a
        `BezierCurve2D`, in order, each keeping its span as its bounds.
        """
    def to_bspline(self) -> BSplineCurve2D:
        """The curve itself, as the other curves give their exact B-spline
        with `to_bspline()`.
        """
    def extract(self, start: float, end: float) -> BSplineCurve2D:
        """The piece over `[start, end]` as a clamped, non-periodic B-spline
        (the range may cross the seam of a periodic curve, up to one
        period).
        """
    def split_at(self, u: float) -> Tuple[BSplineCurve2D, BSplineCurve2D]:
        """The two halves split at `u`, each a clamped, non-periodic B-spline."""
    def length(self, max_length_error: float = 1e-07) -> float:
        """The arc length, integrating the speed to within `max_length_error`.
        `max_length_error` defaults to the default confusion, 1e-7.
        """

@final
class BezierCurve(
    _CurveMethods,
    _BoundedCurveMethods[Point, Discretization],
    _SidedEvaluationMethods[Point, Vector],
):
    """A Bézier curve (optionally rational) in 3D: `degree + 1` control
    poles over one parameter span, `[0, 1]` when built from poles. It has
    no knot vector: `to_bspline` gives the `BSplineCurve`.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, degree: int, poles: _Points) -> BezierCurve: ...
    @staticmethod
    def new_rational(
        degree: int, poles: _Points, weights: _Floats
    ) -> BezierCurve:
        """Build a rational Bézier curve with a positive weight per pole."""
    @staticmethod
    def from_bspline(curve: BSplineCurve) -> BezierCurve:
        """The Bézier curve a `BSplineCurve` is: one clamped, non-periodic
        span. Raises `ConstructionError` for a B-spline with interior
        knots, a periodic one, or one not clamped at its ends.
        """
    def to_bspline(self) -> BSplineCurve:
        """The curve as a `BSplineCurve` over the same bounds, with the same
        point at each parameter.
        """
    def degree(self) -> int:
        """The polynomial degree."""
    def is_rational(self) -> bool:
        """Whether the curve carries rational weights."""
    def poles(self) -> List[Point]:
        """The `degree + 1` control poles."""
    def weights(self) -> Optional[List[float]]:
        """The per-pole weights, or `None` when the curve is not rational."""
    def split_at(self, u: float) -> Tuple[BezierCurve, BezierCurve]:
        """The two pieces either side of `u`, over `[first, u]` and `[u, last]`,
        each a `BezierCurve` of the same degree.
        """
    def elevate_degree(self, amount: int) -> BezierCurve:
        """The curve of degree raised by `amount` (same shape)."""
    def length(self, max_length_error: float = 1e-07) -> float:
        """The arc length, integrating the speed to within `max_length_error`.
        `max_length_error` defaults to the default confusion, 1e-7.
        """

@final
class BezierCurve2D(
    _Curve2DMethods,
    _BoundedCurveMethods[Point2D, Discretization2D],
    _SidedEvaluationMethods[Point2D, Vector2D],
):
    """A Bézier curve (optionally rational) in 2D: `degree + 1` control
    poles over one parameter span, `[0, 1]` when built from poles. It has
    no knot vector: `to_bspline` gives the `BSplineCurve2D`.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, degree: int, poles: _Points2D) -> BezierCurve2D: ...
    @staticmethod
    def new_rational(
        degree: int, poles: _Points2D, weights: _Floats
    ) -> BezierCurve2D:
        """Build a rational Bézier curve with a positive weight per pole."""
    @staticmethod
    def from_bspline(curve: BSplineCurve2D) -> BezierCurve2D:
        """The Bézier curve a `BSplineCurve2D` is: one clamped, non-periodic
        span. Raises `ConstructionError` for a B-spline with interior
        knots, a periodic one, or one not clamped at its ends.
        """
    def to_bspline(self) -> BSplineCurve2D:
        """The curve as a `BSplineCurve2D` over the same bounds, with the same
        point at each parameter.
        """
    def degree(self) -> int:
        """The polynomial degree."""
    def is_rational(self) -> bool:
        """Whether the curve carries rational weights."""
    def poles(self) -> List[Point2D]:
        """The `degree + 1` control poles."""
    def weights(self) -> Optional[List[float]]:
        """The per-pole weights, or `None` when the curve is not rational."""
    def split_at(self, u: float) -> Tuple[BezierCurve2D, BezierCurve2D]:
        """The two pieces either side of `u`, over `[first, u]` and `[u, last]`,
        each a `BezierCurve2D` of the same degree.
        """
    def elevate_degree(self, amount: int) -> BezierCurve2D:
        """The curve of degree raised by `amount` (same shape)."""
    def length(self, max_length_error: float = 1e-07) -> float:
        """The arc length, integrating the speed to within `max_length_error`.
        `max_length_error` defaults to the default confusion, 1e-7.
        """

@final
class QuarticBranch(_CurveMethods):
    """A piece of the curve where two cones, or a cone and a cylinder, meet
    that runs off to infinity at both ends: every point exact, on a
    generator of the first cone or of the cylinder, the parameter running
    over the whole real line.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    @overload
    def __new__(
        cls, first: Cone, second: Cone, index: int, tol: Optional[_ToleranceLike] = None
    ) -> QuarticBranch: ...
    @overload
    def __new__(
        cls, first: Cone, second: Cylinder, index: int, tol: Optional[_ToleranceLike] = None
    ) -> QuarticBranch: ...
    def __setstate__(self, state: Tuple[bool, Transform]) -> None:
        """Restore the orientation and placement a pickle saved."""
    def surfaces(self) -> Tuple[Cone, Union[Cone, Cylinder]]:
        """The two surfaces, as they were before any transform: two cones, the
        one whose generators carry the points first, or a cone and a
        cylinder.
        """
    def index(self) -> int:
        """The branch's place among those of its quartic."""
    def tolerance(self) -> Tolerance:
        """The tolerance its quartic was read with: with `surfaces()` and
        `index()`, it gives the branch back.
        """
    def is_reversed(self) -> bool:
        """Whether the branch runs the other way from the one the constructor
        builds.
        """
    def placement(self) -> Transform:
        """The transform taking the branch the constructor builds to this one."""

@final
class Line2D(_Curve2DMethods):
    """A straight line in a surface's 2D (u, v) parameter space."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, origin: _Point2DLike, direction: Vector2D) -> Line2D: ...
    def origin(self) -> Point2D:
        """The line origin (point at parameter 0)."""
    def direction(self) -> Vector2D:
        """The unit direction of the line."""
    def __repr__(self) -> str: ...
    @staticmethod
    def from_axis(axis: Axis2D) -> Line2D:
        """Build a 2D line from an axis."""
    def axis(self) -> Axis2D:
        """The axis of the line: its origin and unit direction."""
    @staticmethod
    def from_two_points(p1: _Point2DLike, p2: _Point2DLike) -> Line2D:
        """Build the line through two distinct points, with the first as origin
        and the second at positive parameter.
        """

@final
class Circle2D(_Curve2DMethods, _BoundedCurveMethods[Point2D, Discretization2D]):
    """A circle in a surface's 2D (u, v) parameter space."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, center: _Point2DLike, radius: float) -> Circle2D: ...
    @staticmethod
    def from_frame(frame: Frame2D, radius: float) -> Circle2D:
        """Build a 2D circle from its placement frame and radius: centered on
        the frame origin, starting along its x direction, and turning
        clockwise when the frame is indirect.
        """
    def center(self) -> Point2D:
        """The circle center."""
    def radius(self) -> float:
        """The circle radius."""
    def length(self) -> float:
        """The length of the circle, its circumference (`2*PI*radius`)."""
    def is_counterclockwise(self) -> bool:
        """Whether the circle turns counterclockwise as its parameter grows:
        whether its frame is direct.
        """
    def to_bspline(self) -> BSplineCurve2D:
        """The circle as an exact rational B-spline of degree 2, periodic over
        `[0, 2*pi]`; its parameter agrees with the angle at the knots and
        the middle of each quarter turn. In between the same parameter
        gives a different point on each curve, so map a point across with
        `bspline.parameter_of(point)`.
        """
    def __repr__(self) -> str: ...
    def frame(self) -> Frame2D:
        """The placement frame of the circle."""

@final
class Ellipse2D(_Curve2DMethods, _BoundedCurveMethods[Point2D, Discretization2D]):
    """An ellipse in 2D, evaluated as
    `center + major_radius*cos(t)*x + minor_radius*sin(t)*y`. The major
    axis `x` is the normalized `x_direction` and `y` its counterclockwise
    perpendicular. `major_radius` must be at least `minor_radius`
    (`ConstructionError` otherwise): to lay the long axis along another
    direction, turn `x_direction` rather than swap the radii.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        center: _Point2DLike,
        x_direction: Vector2D,
        major_radius: float,
        minor_radius: float,
    ) -> Ellipse2D: ...
    @staticmethod
    def from_frame(
        frame: Frame2D, major_radius: float, minor_radius: float
    ) -> Ellipse2D:
        """Build a 2D ellipse from a placement frame and the two radii."""
    @staticmethod
    def from_center_and_points(
        center: _Point2DLike, s1: _Point2DLike, s2: _Point2DLike
    ) -> Ellipse2D:
        """Build a 2D ellipse from its center, the major-axis end point `s1`
        and a point `s2` whose distance from the major axis is the minor
        radius.
        """
    def frame(self) -> Frame2D:
        """The placement frame of the ellipse."""
    def center(self) -> Point2D:
        """The ellipse center."""
    def major_radius(self) -> float:
        """The semi-major radius."""
    def minor_radius(self) -> float:
        """The semi-minor radius."""
    def enclosed_area(self) -> float:
        """The enclosed area (`PI*major*minor`)."""
    def length(self, max_length_error: float = 1e-07) -> float:
        """The perimeter, integrating the speed to within `max_length_error`: it has
        no closed form.
        `max_length_error` defaults to the default confusion, 1e-7.
        """
    def eccentricity(self) -> float:
        """The eccentricity: the ratio of a point's distance to a focus over
        its distance to that focus's directrix.
        """
    def foci(self) -> Tuple[Point2D, Point2D]:
        """The two foci on the major axis, the one toward the frame's x
        direction first.
        """
    def directrices(self) -> Optional[Tuple[Line2D, Line2D]]:
        """The two directrices across the major axis, in the order of `foci`,
        or None for a circle, whose directrices are at infinity.
        """
    def to_bspline(self) -> BSplineCurve2D:
        """The ellipse as an exact rational B-spline of degree 2, periodic over
        `[0, 2*pi]`; its parameter agrees with the angle at the knots and
        the middle of each quarter turn. In between the same parameter
        gives a different point on each curve, so map a point across with
        `bspline.parameter_of(point)`.
        """
    def __repr__(self) -> str: ...

@final
class Parabola2D(_Curve2DMethods):
    """A parabola in 2D, parametrized so the apex is at parameter 0."""
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, apex: _Point2DLike, x_direction: Vector2D, focal_distance: float) -> Parabola2D: ...
    @staticmethod
    def from_frame(frame: Frame2D, focal_distance: float) -> Parabola2D:
        """Build a 2D parabola from a placement frame and focal distance."""
    @staticmethod
    def from_focus_and_directrix(focus: _Point2DLike, directrix: Line2D) -> Parabola2D:
        """Build the parabola of every point as far from `focus` as from
        `directrix`; the parameter runs along the directrix's direction.
        """
    def frame(self) -> Frame2D:
        """The placement frame of the parabola."""
    def apex(self) -> Point2D:
        """The apex of the parabola, at the frame origin."""
    def focal_distance(self) -> float:
        """The focal distance."""
    def eccentricity(self) -> float:
        """The eccentricity: the ratio of a point's distance to a focus over
        its distance to that focus's directrix.
        """
    def focus(self) -> Point2D:
        """The focus, `focal_distance` from the apex along the axis."""
    def directrix(self) -> Line2D:
        """The directrix, across the axis `focal_distance` behind the apex:
        every point of the curve is as far from it as from the focus.
        """
    def __repr__(self) -> str: ...

@final
class Hyperbola2D(_Curve2DMethods):
    """A hyperbola branch in 2D, parametrized as
    `center + major*cosh(u)*x_dir + minor*sinh(u)*y_dir`.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        center: _Point2DLike,
        x_direction: Vector2D,
        major_radius: float,
        minor_radius: float,
    ) -> Hyperbola2D: ...
    @staticmethod
    def from_frame(
        frame: Frame2D, major_radius: float, minor_radius: float
    ) -> Hyperbola2D:
        """Build a 2D hyperbola from a placement frame and semi-major and
        semi-minor radii.
        """
    @staticmethod
    def from_center_and_points(
        center: _Point2DLike, s1: _Point2DLike, s2: _Point2DLike
    ) -> Hyperbola2D:
        """Build a 2D hyperbola from its center, a major-axis end point `s1`
        and a second point `s2` fixing the minor radius and the side the
        frame's y direction points to.
        """
    def frame(self) -> Frame2D:
        """The placement frame of the hyperbola."""
    def center(self) -> Point2D:
        """The center of the hyperbola, at the frame origin."""
    def major_radius(self) -> float:
        """The semi-major radius."""
    def minor_radius(self) -> float:
        """The semi-minor radius."""
    def eccentricity(self) -> float:
        """The eccentricity: the ratio of a point's distance to a focus over
        its distance to that focus's directrix.
        """
    def foci(self) -> Tuple[Point2D, Point2D]:
        """The two foci on the major axis, the one toward the frame's x
        direction first.
        """
    def directrices(self) -> Tuple[Line2D, Line2D]:
        """The two directrices across the major axis, in the order of `foci`."""
    def other_branch(self) -> Hyperbola2D:
        """The other branch: the same hyperbola mirrored through its center,
        with the same parametrization direction around the center.
        """
    def __repr__(self) -> str: ...

@final
class TrimmedCurve(_CurveMethods, _BoundedCurveMethods[Point, Discretization]):
    """A 3D curve restricted to the parameter range `[start_parameter, end_parameter]` of a
    basis curve: a line segment, an arc of a circle or conic, or a piece of
    a B-spline. Geometry only, with no vertices or orientation.

    The basis keeps its parametrisation and `bounds()` is
    `(start_parameter, end_parameter)`; evaluation does not check a parameter against it. On
    a periodic basis (circle, ellipse, periodic
    B-spline) the range may cross the seam, up to one period. Pass it to
    `intersect`, `curve_curve_extrema` or `curve_3d_to_uv_within` to search
    only the range.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, basis: Curve, start_parameter: float, end_parameter: float
    ) -> TrimmedCurve: ...
    @staticmethod
    def segment(start: _PointLike, end: _PointLike) -> TrimmedCurve:
        """The segment from `start` to `end`, parametrised by arc length from
        `start` over `[0, |end - start|]`.
        """
    @staticmethod
    def arc_through(start: _PointLike, middle: _PointLike, end: _PointLike) -> TrimmedCurve:
        """The arc from `start` through `middle` to `end`, on the circle
        through the three points, whose normal the points turn
        counterclockwise about. Raises `ConstructionError` when two points
        coincide or the three are collinear.
        """
    @staticmethod
    def arc(
        center: _PointLike,
        start: _PointLike,
        end: _PointLike,
        normal: Vector,
        tol: Optional[_ToleranceLike] = None,
    ) -> TrimmedCurve:
        """The arc about `center` from `start` to `end`, counterclockwise about
        `normal`; an `end` confused with `start` makes the full turn. Raises
        `ConstructionError` when `start` or `end` is off the circle (out of
        the plane through `center` normal to `normal`, or `end` at another
        radius) by more than `tol.confusion` (`Tolerance.default()` when
        omitted), when `start` is at `center` or for a zero normal.
        """
    def sweep_angle(self) -> float:
        """The angle an arc sweeps, `last - first`, counterclockwise about its
        circle's `normal()`. Raises `TypeError` when the basis is not a
        `Circle`.
        """
    def length(self, max_length_error: float = 1e-07) -> float:
        """The length of the piece: exact on a segment or a circular arc
        (`max_length_error` is then unused), within `max_length_error` otherwise.
        `max_length_error` defaults to the default confusion, 1e-7.
        """
    def basis(self) -> Curve:
        """The basis curve."""
    def start_parameter(self) -> float:
        """The start of the parameter range, the parameter of `start_point`."""
    def end_parameter(self) -> float:
        """The end of the parameter range, the parameter of `end_point`."""
    def to_bspline(self) -> BSplineCurve:
        """The trimmed piece as an exact B-spline over the same range: a
        clamped B-spline for an elementary basis, the extracted piece for
        a B-spline or Bézier one (raises `ValueError` when the poles would
        overflow `f64`). A line, parabola, B-spline or Bézier basis keeps
        its parameter; the rational arcs of a circle, ellipse or hyperbola
        do not, so map a point across with `bspline.parameter_of(point)`.
        """
    def __repr__(self) -> str: ...

@final
class TrimmedCurve2D(_Curve2DMethods, _BoundedCurveMethods[Point2D, Discretization2D]):
    """A 2D curve restricted to the parameter range `[start_parameter, end_parameter]` of a
    basis curve: a line segment, a circular arc or a B-spline piece,
    typically a trimmed pcurve in a surface's (u, v) parameter space. Mirrors `TrimmedCurve`.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, basis: Curve2D, start_parameter: float, end_parameter: float
    ) -> TrimmedCurve2D: ...
    @staticmethod
    def segment(start: _Point2DLike, end: _Point2DLike) -> TrimmedCurve2D:
        """The segment from `start` to `end`, parametrised by arc length from
        `start` over `[0, |end - start|]`.
        """
    @staticmethod
    def arc_through(start: _Point2DLike, middle: _Point2DLike, end: _Point2DLike) -> TrimmedCurve2D:
        """The arc from `start` through `middle` to `end`, on the circle
        through the three points: counterclockwise when they turn left,
        clockwise otherwise. Raises `ConstructionError` when two points
        coincide or the three are collinear.
        """
    @staticmethod
    def arc(
        center: _Point2DLike,
        start: _Point2DLike,
        end: _Point2DLike,
        counterclockwise: bool = True,
        tol: Optional[_ToleranceLike] = None,
    ) -> TrimmedCurve2D:
        """The arc about `center` from `start` to `end`, counterclockwise or
        clockwise; an `end` confused with `start` makes the full turn.
        Raises `ConstructionError` when `end` is at another distance from
        `center` than `start` by more than `tol.confusion`
        (`Tolerance.default()` when omitted), or `start` is at `center`.
        """
    def sweep_angle(self) -> float:
        """The angle an arc sweeps, `last - first`, in its circle's turning
        sense (`basis().is_counterclockwise()`). Raises `TypeError` when the
        basis is not a `Circle2D`.
        """
    def length(self, max_length_error: float = 1e-07) -> float:
        """The length of the piece: exact on a segment or a circular arc
        (`max_length_error` is then unused), within `max_length_error` otherwise.
        `max_length_error` defaults to the default confusion, 1e-7.
        """
    def basis(self) -> Curve2D:
        """The basis curve."""
    def to_bspline(self) -> BSplineCurve2D:
        """The trimmed piece as an exact B-spline over the same range: the
        extracted piece for a B-spline basis, a clamped B-spline for a line
        or circle (raises `ValueError` when the poles would overflow
        `f64`). A line or B-spline basis keeps its parameter; a circle's
        rational arcs do not, so map a point across with
        `bspline.parameter_of(point)`.
        """
    def start_parameter(self) -> float:
        """The start of the parameter range, the parameter of `start_point`."""
    def end_parameter(self) -> float:
        """The end of the parameter range, the parameter of `end_point`."""
    def __repr__(self) -> str: ...


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

_UIsoCurveType = TypeVar("_UIsoCurveType", covariant=True)
_VIsoCurveType = TypeVar("_VIsoCurveType", covariant=True)

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
        `point_at(u, v)`, and the normal there is flipped.
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
    def points_at(self, uvs: _Pairs) -> List[Point]:
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
    def is_u_closed(self, tol: Optional[_ToleranceLike] = None) -> bool:
        """Whether travelling in `u` returns to its start: the boundary
        curves at the two `u` bounds coincide (`Tolerance.default()`
        when omitted). Periodic in `u` implies closed in `u`, not
        conversely.
        """
    def is_v_closed(self, tol: Optional[_ToleranceLike] = None) -> bool:
        """Whether travelling in `v` returns to its start, as
        `is_u_closed`.
        """
    def continuity_u(self, u: float, tol: Optional[_ToleranceLike] = None) -> Union[int, float]:
        """Continuity order `n` of `C^n` in `u` at parameter `u` as an
        `int`, or `math.inf` where smooth. A B-spline reports
        `u_degree - multiplicity` at a `u` knot (within
        `tol.parametric`); the ends of a non-periodic direction
        report smooth, and a periodic direction's seam counts as a
        knot. Use `BSplineSurface.first_derivatives_sided` and its siblings to
        evaluate each side of a knot reporting a low order.
        """
    def continuity_v(self, v: float, tol: Optional[_ToleranceLike] = None) -> Union[int, float]:
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
        self, u: float, v: float, tol: Optional[_ToleranceLike] = None
    ) -> Optional[Vector]:
        """Unit normal at surface parameters `(u, v)`
        (`Tolerance.default()` when omitted): the direction of
        `du × dv`, the cross product of the `u` and `v` derivatives,
        so reversing either direction flips it. Singular points such
        as poles give their limiting normal; returns `None` at an
        apex, whose meridians each have their own, and where no
        normal can be recovered.
        """
    def curvatures(
        self, u: float, v: float, tol: Optional[_ToleranceLike] = None
    ) -> Optional[SurfaceCurvatures]:
        """Principal curvatures with their directions, and the mean and
        Gaussian curvatures, at `(u, v)` (`Tolerance.default()` when
        omitted), as a `SurfaceCurvatures`. `None` where the unit
        normal degenerates, as at a sphere's pole or a cone's apex,
        or a derivative is not finite. Curvatures are measured along
        the unit normal: negative bends away from it.
        """
    def project_point(
        self, point: _PointLike, tol: Optional[_ToleranceLike] = None
    ) -> SurfaceProjection:
        """Projects `point` onto the surface: the closest point, as a
        `SurfaceProjection`.

        A distance within `tol.confusion` (`1e-7` by default) snaps
        to `0.0`, so projection agrees with `contains`; pass a
        smaller `Tolerance` to read a smaller distance.
        """
    def project_points(
        self, points: _Points, tol: Optional[_ToleranceLike] = None
    ) -> List[SurfaceProjection]:
        """Projects each point in `points` onto the surface, returning a
        `SurfaceProjection` per point.
        """
    def project_points_array(
        self, points: numpy.typing.ArrayLike, tol: Optional[_ToleranceLike] = None
    ) -> Tuple[_FloatArray, _FloatArray, _FloatArray]:
        """Projects each row of the `(N, 3)` array-like `points` onto
        the surface, returning numpy arrays `(parameters, points,
        distances)`: the `(u, v)` row, the closest point and its
        distance per row.
        """
    def contains(self, point: _PointLike, tol: Optional[_ToleranceLike] = None) -> bool:
        """Returns whether `point` lies on the surface within `tol`
        (`Tolerance.default()` when omitted).
        """
    def parameters_of(
        self, point: _PointLike, tol: Optional[_ToleranceLike] = None
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
    def tight_bounding_box(self, tol: Optional[_ToleranceLike] = None) -> Optional[BoundingBox]:
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
        self, uv_curve: Curve2D, tol: Optional[_ToleranceLike] = None
    ) -> Curve:
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
        self, uv_curve: TrimmedCurve2D, tol: Optional[_ToleranceLike] = None
    ) -> TrimmedCurve: ...
    @overload
    def curve_3d_to_uv(
        self, curve: Curve, tol: Optional[_ToleranceLike] = None
    ) -> Curve2D:
        """The exact UV curve of `curve`, a 3D curve lying on this surface: the
        curve in the surface's (u, v) parameter space, in the same parameter,
        whose point at `t` is the (u, v) of the curve's point at `t`.

        Closed forms: any line, circle, ellipse, parabola, hyperbola, B-spline
        or Bezier curve on a plane, as the 2D curve of the same kind; a line
        along the axis of a cylinder or a generator of a cone, as a line along
        v; a circle about the axis of a cylinder or cone, a parallel of a
        sphere, or a poloidal or toroidal circle of a torus, as a line along u
        or v; and on an extrusion or revolution, a generator line, a parallel
        circle, or the profile translated or turned in its own parameter, as a
        line along u or v. A `TrimmedCurve` maps to the `TrimmedCurve2D` of its
        basis's image over the same range. A sphere's meridian runs through
        both poles, where its u jumps, so only a `TrimmedCurve` of it from
        pole to pole at most maps, to a line along v. Raises
        `Curve3DToUvError` for any other pair, or when the curve does not
        lie on the surface within `tol` (`Tolerance.default()` when
        omitted).
        """
    @overload
    def curve_3d_to_uv(
        self, curve: TrimmedCurve, tol: Optional[_ToleranceLike] = None
    ) -> TrimmedCurve2D: ...
    def curve_3d_to_uv_within(
        self,
        curve: _CurveMethods,
        max_deviation: float,
        tol: Optional[_ToleranceLike] = None,
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
        uv_curve: _Curve2DMethods,
        max_deviation: float,
        tol: Optional[_ToleranceLike] = None,
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
    def trimmed(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> TrimmedSurface:
        """The patch of the surface over the parameter ranges
        `u_range` and `v_range`, each `(first, last)`, as a
        `TrimmedSurface`. A periodic direction's range may cross the
        seam, up to one period. Raises
        `TrimmedSurfaceConstructionError` for non-finite ends, an
        empty range, a range outside a non-periodic direction's
        bounds, or one longer than a period.
        """
    @overload
    def intersect(
        self, other: _CurveMethods, tol: Optional[_ToleranceLike] = None
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
        self, other: _SurfaceMethods[Any, Any] | TrimmedSurface, tol: Optional[_ToleranceLike] = None
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
        self, point: _PointLike, tol: Optional[_ToleranceLike] = None
    ) -> List[SurfaceProjection]:
        """The local minima of the distance from `point` to the
        surface, as `SurfaceProjection` objects ordered by ascending
        distance.
        """


@final
class Plane(_SurfaceMethods[Line, Line]):
    """An infinite plane, parametrized by in-plane coordinates (u, v). Its
    normal is the `normal` it was built with.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, point: _PointLike, normal: Vector) -> Plane: ...
    @staticmethod
    def from_frame(frame: Frame) -> Plane:
        """Build a plane from a placement frame (normal = frame z direction)."""
    @staticmethod
    def from_three_points(p1: _PointLike, p2: _PointLike, p3: _PointLike) -> Plane:
        """Build a plane through three non-collinear points."""
    @staticmethod
    def from_coefficients(a: float, b: float, c: float, d: float) -> Plane:
        """Build a plane from the equation `a*x + b*y + c*z + d = 0`."""
    def axis(self) -> Axis:
        """The normal axis of the plane: through the frame's origin, along
        the unit normal, which reverses with `reversed_u()` or
        `reversed_v()`. `normal(u, v)` gives the same direction at
        every point.
        """
    def signed_distance(self, point: _PointLike) -> float:
        """The signed distance from the plane to `point`: positive on the side
        the normal points to, negative on the other.
        """
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface:
        """The patch `u_range × v_range` as an exact bilinear B-spline
        surface (raises `ValueError` for a bad range). It keeps the
        parameters: the same `(u, v)` gives the same point on both.
        """
    def frame(self) -> Frame:
        """The placement frame of the plane."""

@final
class Cylinder(_SurfaceMethods[Line, Circle]):
    """An infinite circular cylinder; u is the angle around the axis, v the
    height along it. Its normal points away from the axis.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, center: _PointLike, axis_direction: Vector, radius: float) -> Cylinder: ...
    @staticmethod
    def from_frame(frame: Frame, radius: float) -> Cylinder:
        """Build a cylinder from a placement frame (axis = frame z direction)."""
    @staticmethod
    def from_axis(axis: Axis, radius: float) -> Cylinder:
        """Build a cylinder from an axis and a radius."""
    @staticmethod
    def from_circle(circle: Circle) -> Cylinder:
        """Build the cylinder containing a circle (axis through the circle
        center, along its normal).
        """
    def radius(self) -> float:
        """The cylinder radius."""
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface:
        """The patch `u_range × v_range` as an exact rational B-spline
        surface (raises `ValueError` for a bad range).
        Its rational arcs do not keep the angles as parameters: the same
        `(u, v)` gives a different point on each surface between the
        knots, so map a point across with `bspline.parameters_of(point)`.
        """
    def frame(self) -> Frame:
        """The placement frame of the cylinder."""
    def axis(self) -> Axis:
        """The axis of the cylinder: its origin and unit direction."""

@final
class Cone(_SurfaceMethods[Line, Circle]):
    """An infinite cone; u is the angle around the axis, v the distance along a
    generator from the reference circle. Its normal points away from the
    axis.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        center: _PointLike,
        axis_direction: Vector,
        half_angle: float,
        reference_radius: float,
    ) -> Cone: ...
    @staticmethod
    def from_frame(frame: Frame, half_angle: float, reference_radius: float) -> Cone:
        """Build a cone from a placement frame, half-angle and reference radius."""
    @staticmethod
    def from_two_points_and_radii(
        p1: _PointLike, p2: _PointLike, r1: float, r2: float
    ) -> Cone:
        """Build a cone through two circular sections given by axis points and radii."""
    @staticmethod
    def from_apex(apex: _PointLike, axis_direction: Vector, half_angle: float) -> Cone:
        """Build the cone of apex `apex` opening along `axis_direction` at
        half-angle `half_angle`; its reference circle shrinks to the apex.
        """
    def axis(self) -> Axis:
        """The axis of the cone: through the center of the reference section
        along the frame's z direction.
        """
    def apex(self) -> Point:
        """The cone apex."""
    def half_angle(self) -> float:
        """The angle between the axis and a generator, in radians, with
        `0 < |half_angle| < pi/2`; a negative half-angle makes the radius
        shrink, not grow, along `axis_direction`.
        """
    def reference_radius(self) -> float:
        """The radius of the reference section at v = 0."""
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface:
        """The patch `u_range × v_range` as an exact rational B-spline
        surface (raises `ValueError` for a bad range).
        Its rational arcs do not keep the angles as parameters: the same
        `(u, v)` gives a different point on each surface between the
        knots, so map a point across with `bspline.parameters_of(point)`.
        """
    def frame(self) -> Frame:
        """The placement frame of the cone."""

@final
class Sphere(_SurfaceMethods[Circle, Circle]):
    """A sphere; u is the longitude in [0, 2*pi), v the latitude in [-pi/2, pi/2].
    Its normal points out of the sphere.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, center: _PointLike, radius: float) -> Sphere: ...
    @staticmethod
    def from_frame(frame: Frame, radius: float) -> Sphere:
        """Build a sphere from a placement frame and radius."""
    def center(self) -> Point:
        """The sphere center."""
    def radius(self) -> float:
        """The sphere radius."""
    def area(self) -> float:
        """The surface area (`4*PI*radius^2`)."""
    def volume(self) -> float:
        """The enclosed volume (`4/3*PI*radius^3`)."""
    def to_bspline(self) -> BSplineSurface:
        """The whole sphere as an exact rational B-spline surface.
        Its rational arcs do not keep the angles as parameters: the same
        `(u, v)` gives a different point on each surface between the
        knots, so map a point across with `bspline.parameters_of(point)`.
        """
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface:
        """The patch `u_range × v_range` as an exact rational B-spline
        surface (raises `ValueError` for a bad range).
        Its rational arcs do not keep the angles as parameters: the same
        `(u, v)` gives a different point on each surface between the
        knots, so map a point across with `bspline.parameters_of(point)`.
        """
    def frame(self) -> Frame:
        """The placement frame of the sphere."""

@final
class Torus(_SurfaceMethods[Circle, Circle]):
    """A torus; u is the angle around the main axis, v the angle around the tube.
    Its normal points out of the tube.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, center: _PointLike, axis_direction: Vector, major_radius: float, minor_radius: float
    ) -> Torus: ...
    @staticmethod
    def from_frame(frame: Frame, major_radius: float, minor_radius: float) -> Torus:
        """Build a torus from a placement frame and the two radii."""
    def major_radius(self) -> float:
        """The distance from the center to the tube center line."""
    def minor_radius(self) -> float:
        """The tube radius."""
    def to_bspline(self) -> BSplineSurface:
        """The whole torus as an exact rational B-spline surface.
        Its rational arcs do not keep the angles as parameters: the same
        `(u, v)` gives a different point on each surface between the
        knots, so map a point across with `bspline.parameters_of(point)`.
        """
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface:
        """The patch `u_range × v_range` as an exact rational B-spline
        surface (raises `ValueError` for a bad range).
        Its rational arcs do not keep the angles as parameters: the same
        `(u, v)` gives a different point on each surface between the
        knots, so map a point across with `bspline.parameters_of(point)`.
        """
    def axis(self) -> Axis:
        """The axis the tube turns about: through the center along the
        frame's z direction.
        """
    def frame(self) -> Frame:
        """The placement frame of the torus."""
    def center(self) -> Point:
        """The center of the torus."""

@final
class ExtrusionSurface(_SurfaceMethods[Line, Curve]):
    """A profile curve swept along a direction: `profile(u) + v * direction`,
    `v` the signed distance travelled.

    `ExtrusionSurface(profile, direction)` normalizes `direction`, so its
    length is dropped: `Vector(0, 0, 2)` sweeps the same surface as
    `Vector(0, 0, 1)`, and `v` stays the distance along it. Raises
    `ValueError` for a zero direction and for a profile lying on one line
    along it.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, profile: Curve, direction: Vector) -> ExtrusionSurface: ...
    def profile(self) -> Curve:
        """The profile, the curve at `v = 0`."""
    def direction(self) -> Vector:
        """The unit direction of the sweep."""
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface:
        """The patch `u_range × v_range` as an exact B-spline surface
        (raises `ValueError` for a bad range). It keeps `v`, and keeps `u`
        when the swept curve's `to_bspline` keeps its parameter; otherwise
        map a point across with `bspline.parameters_of(point)`.
        """

@final
class RevolutionSurface(_SurfaceMethods[Curve, Circle]):
    """A profile curve turned about an axis: the profile point at `v` rotated
    by the angle `u`.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(cls, profile: Curve, axis: Axis) -> RevolutionSurface: ...
    def profile(self) -> Curve:
        """The profile, the curve at `u = 0`."""
    def axis(self) -> Axis:
        """The axis the profile turns about."""
    def to_bspline(self) -> BSplineSurface:
        """The whole surface as an exact rational B-spline surface (raises
        `ValueError` for an unbounded profile).
        Its rational arcs do not keep the angles as parameters: the same
        `(u, v)` gives a different point on each surface between the
        knots, so map a point across with `bspline.parameters_of(point)`.
        """
    def to_bspline_patch(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface:
        """The patch `u_range × v_range` as an exact rational B-spline
        surface (raises `ValueError` for a bad range).
        Its rational arcs do not keep the angles as parameters: the same
        `(u, v)` gives a different point on each surface between the
        knots, so map a point across with `bspline.parameters_of(point)`.
        """

@final
class BSplineSurface(
    _SurfaceMethods[BSplineCurve, BSplineCurve],
    _SidedSurfaceEvaluationMethods,
    _SurfaceExtremaMethods,
):
    """A tensor-product B-spline (optionally rational, optionally periodic) surface."""
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls,
        u_degree: int,
        v_degree: int,
        poles: _PointGrid,
        u_knots: _Floats,
        u_multiplicities: _Ints,
        v_knots: _Floats,
        v_multiplicities: _Ints,
        u_periodic: bool,
        v_periodic: bool,
    ) -> BSplineSurface: ...
    @staticmethod
    def new_rational(
        u_degree: int,
        v_degree: int,
        poles: _PointGrid,
        weights: Sequence[_Floats],
        u_knots: _Floats,
        u_multiplicities: _Ints,
        v_knots: _Floats,
        v_multiplicities: _Ints,
        u_periodic: bool,
        v_periodic: bool,
    ) -> BSplineSurface:
        """Build a rational B-spline (NURBS) surface with a weight per pole."""
    def u_degree(self) -> int:
        """Polynomial degree in the u direction."""
    def v_degree(self) -> int:
        """Polynomial degree in the v direction."""
    def is_rational(self) -> bool:
        """Whether the surface carries rational weights."""
    def u_pole_count(self) -> int:
        """The number of poles in the `u` direction (rows of the net)."""
    def v_pole_count(self) -> int:
        """The number of poles in the `v` direction (columns of the net)."""
    def poles(self) -> List[List[Point]]:
        """The control poles as `poles[i_u][i_v]`."""
    def weights(self) -> Optional[List[List[float]]]:
        """The per-pole weights laid out as `poles`, or `None` when the
        surface is not rational.
        """
    def u_knots(self) -> List[float]:
        """The distinct knots in the `u` direction."""
    def u_multiplicities(self) -> List[int]:
        """The multiplicity of each `u` knot."""
    def v_knots(self) -> List[float]:
        """The distinct knots in the `v` direction."""
    def v_multiplicities(self) -> List[int]:
        """The multiplicity of each `v` knot."""
    def u_flat_knots(self) -> List[float]:
        """The flat knot vector in the `u` direction: each distinct knot
        repeated by its multiplicity.
        """
    def v_flat_knots(self) -> List[float]:
        """The flat knot vector in the `v` direction: each distinct knot
        repeated by its multiplicity.
        """
    @staticmethod
    def from_flat_knots(
        u_degree: int,
        v_degree: int,
        poles: _PointGrid,
        u_flat_knots: _Floats,
        v_flat_knots: _Floats,
        weights: Optional[Sequence[_Floats]] = None,
        u_periodic: bool = False,
        v_periodic: bool = False,
    ) -> BSplineSurface:
        """Build a surface from its flat knot vectors, the form textbooks
        use: each distinct knot repeated by its multiplicity, so a
        non-periodic direction of `n` poles takes `n + degree + 1` knots.
        `weights`, a grid laid out as `poles`, make it rational. The flat
        knots are as `u_flat_knots()` and `v_flat_knots()` give them.
        """
    @staticmethod
    def approximate(
        points: _PointGrid,
        u_degree: int,
        v_degree: int,
        u_pole_count: int,
        v_pole_count: int,
        parametrization: _Parametrization = "centripetal",
    ) -> BSplineSurface:
        """Least-squares B-spline fit through the grid `points[i][j]`.

        `points` is rows of `Point` objects or an `(n, m, 3)` array-like.

        `parametrization` is a `PointParametrization` or its string value,
        `"centripetal"` (the default) or `"chordal"`.
        """
    @staticmethod
    def approximate_within(
        points: _PointGrid,
        u_degree: int,
        v_degree: int,
        max_deviation: float,
        parametrization: _Parametrization = "centripetal",
    ) -> SurfaceApproximation:
        """Least-squares B-spline fit through the grid `points[i][j]` with as
        few poles as keep every point within `max_deviation` of the
        surface.

        Returns the surface with its `deviation`, the largest distance
        between a point and the surface at its parameters, proven.

        `points` is rows of `Point` objects or an `(n, m, 3)` array-like.

        `parametrization` is a `PointParametrization` or its string value,
        `"centripetal"` (the default) or `"chordal"`.
        """
    def insert_u_knot(self, u: float, times: int) -> BSplineSurface:
        """The surface with the knot `u` inserted `times` times in the `u`
        direction (same shape, more poles).
        """
    def insert_v_knot(self, v: float, times: int) -> BSplineSurface:
        """The surface with the knot `v` inserted `times` times in the `v`
        direction (same shape, more poles).
        """
    def refine_u_knots(self, knots: _Floats) -> BSplineSurface:
        """The surface with each of `knots` inserted once in the `u`
        direction (same shape, more poles).
        """
    def refine_v_knots(self, knots: _Floats) -> BSplineSurface:
        """The surface with each of `knots` inserted once in the `v`
        direction (same shape, more poles).
        """
    def elevate_u_degree(self, amount: int) -> BSplineSurface:
        """The surface of `u` degree raised by `amount` (same shape)."""
    def elevate_v_degree(self, amount: int) -> BSplineSurface:
        """The surface of `v` degree raised by `amount` (same shape)."""
    def remove_u_knot(
        self, u: float, times: int, max_deviation: float
    ) -> Tuple[BSplineSurface, int]:
        """The surface with the knot `u` removed up to `times` times in the
        `u` direction, and how many copies were removed: a removal is kept
        only while every point stays within `max_deviation` of this one.
        """
    def remove_v_knot(
        self, v: float, times: int, max_deviation: float
    ) -> Tuple[BSplineSurface, int]:
        """The surface with the knot `v` removed up to `times` times in the
        `v` direction, and how many copies were removed: a removal is kept
        only while every point stays within `max_deviation` of this one.
        """
    def remove_u_knots(
        self, knots: _Floats, max_deviation: float
    ) -> Tuple[BSplineSurface, int]:
        """The surface with each of `knots` removed once in the `u`
        direction, and how many copies were removed in all.
        """
    def remove_v_knots(
        self, knots: _Floats, max_deviation: float
    ) -> Tuple[BSplineSurface, int]:
        """The surface with each of `knots` removed once in the `v`
        direction, and how many copies were removed in all.
        """
    def to_bezier_patches(self) -> List[List[BezierSurface]]:
        """The polynomial (or rational) patch of each pair of knot spans as
        a `BezierSurface`: `patches[i][j]` covers the `i`-th `u` span and
        the `j`-th `v` span, which it keeps as its bounds.
        """
    def to_bspline(self) -> BSplineSurface:
        """The surface itself, as the other surfaces give their exact
        B-spline with `to_bspline()`.
        """
    def extract(
        self, u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> BSplineSurface:
        """The piece over `u_range × v_range` as a clamped, non-periodic
        B-spline (a range may cross the seam of a periodic direction, up
        to one period).
        """
    def extract_u(self, start: float, end: float) -> BSplineSurface:
        """The piece over `[start, end]` in `u`, the `v` direction
        unchanged.
        """
    def extract_v(self, start: float, end: float) -> BSplineSurface:
        """The piece over `[start, end]` in `v`, the `u` direction
        unchanged.
        """
    def split_u_at(self, u: float) -> Tuple[BSplineSurface, BSplineSurface]:
        """The two halves split at `u`, each a clamped, non-periodic B-spline."""
    def split_v_at(self, v: float) -> Tuple[BSplineSurface, BSplineSurface]:
        """The two halves split at `v`, each a clamped, non-periodic B-spline."""

@final
class BezierSurface(
    _SurfaceMethods[BezierCurve, BezierCurve],
    _SidedSurfaceEvaluationMethods,
    _SurfaceExtremaMethods,
):
    """A Bézier surface (optionally rational): a grid of
    `(u_degree + 1) x (v_degree + 1)` control poles over one parameter
    rectangle, `[0, 1] x [0, 1]` when built from poles. It has no knot
    vectors: `to_bspline` gives the `BSplineSurface`.
    """
    def __repr__(self) -> str: ...
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, u_degree: int, v_degree: int, poles: _PointGrid
    ) -> BezierSurface: ...
    @staticmethod
    def new_rational(
        u_degree: int,
        v_degree: int,
        poles: _PointGrid,
        weights: Sequence[_Floats],
    ) -> BezierSurface:
        """Build a rational Bézier surface with a positive weight per pole,
        laid out as `poles`.
        """
    @staticmethod
    def from_bspline(surface: BSplineSurface) -> BezierSurface:
        """The Bézier surface a `BSplineSurface` is: one clamped,
        non-periodic span in each direction. Raises `ConstructionError`
        otherwise.
        """
    def to_bspline(self) -> BSplineSurface:
        """The surface as a `BSplineSurface` over the same bounds, with the
        same point at each `(u, v)`.
        """
    def u_degree(self) -> int:
        """The surface's degree in the `u` direction."""
    def v_degree(self) -> int:
        """The surface's degree in the `v` direction."""
    def is_rational(self) -> bool:
        """Whether the surface carries rational weights."""
    def u_pole_count(self) -> int:
        """The number of pole rows, `u_degree + 1`."""
    def v_pole_count(self) -> int:
        """The number of poles per row, `v_degree + 1`."""
    def poles(self) -> List[List[Point]]:
        """The control poles as `poles[i_u][i_v]`."""
    def weights(self) -> Optional[List[List[float]]]:
        """The per-pole weights laid out as `poles`, or `None` when the
        surface is not rational.
        """
    def split_u_at(self, u: float) -> Tuple[BezierSurface, BezierSurface]:
        """The two pieces either side of `u`, over `[first, u]` and
        `[u, last]` in `u`, each a `BezierSurface`.
        """
    def split_v_at(self, v: float) -> Tuple[BezierSurface, BezierSurface]:
        """The two pieces either side of `v`, as `split_u_at` for `v`."""
    def elevate_u_degree(self, amount: int) -> BezierSurface:
        """The surface with its `u` degree raised by `amount` (same shape)."""
    def elevate_v_degree(self, amount: int) -> BezierSurface:
        """The surface with its `v` degree raised by `amount` (same shape)."""

@final
class TrimmedSurface:
    """A surface restricted to a parameter rectangle `u_range` by `v_range`
    of a basis surface: the patch a face of a model lies on. Geometry only,
    with no boundary curves or orientation.

    The basis keeps its parametrisation, and `u_bounds()` and `v_bounds()`
    are the ranges; evaluation does not check a parameter against them. On
    a periodic basis direction (a cylinder's angle) the range may cross the
    seam, up to one period. Pass it to `intersect` to get only the curves
    and points inside the ranges of the trimmed sides.
    """
    __hash__: ClassVar[None]  # type: ignore[assignment]
    def __new__(
        cls, basis: _SurfaceMethods[Any, Any], u_range: Tuple[float, float], v_range: Tuple[float, float]
    ) -> TrimmedSurface: ...
    def basis(self) -> Surface:
        """The basis surface."""
    def u_range(self) -> Tuple[float, float]:
        """The `(first, last)` range of `u`."""
    def v_range(self) -> Tuple[float, float]:
        """The `(first, last)` range of `v`."""
    def u_bounds(self) -> Tuple[float, float]:
        """The `(first, last)` `u` parameter bounds: the `u` range."""
    def v_bounds(self) -> Tuple[float, float]:
        """The `(first, last)` `v` parameter bounds: the `v` range."""
    def point_at(self, u: float, v: float) -> Point:
        """The basis's point at surface parameters (u, v), inside the ranges
        or not.
        """
    def points_at(self, uvs: _Pairs) -> List[Point]:
        """Evaluate many (u, v) pairs at once."""
    def first_derivatives(self, u: float, v: float) -> Tuple[Point, Tuple[Vector, Vector]]:
        """The point and the first partial derivatives `(du, dv)` at `(u, v)`."""
    def normal(
        self, u: float, v: float, tol: Optional[_ToleranceLike] = None
    ) -> Optional[Vector]:
        """The basis's unit normal at `(u, v)` (`Tolerance.default()` when
        omitted), as `Surface.normal`.
        """
    def project_point(
        self, point: _PointLike, tol: Optional[_ToleranceLike] = None
    ) -> SurfaceProjection:
        """Projects `point` onto the patch: the basis's closest point when it
        lies in the ranges, otherwise the closest point on the boundary of
        the ranges, as a `SurfaceProjection`. A local minimum of the
        distance inside the ranges other than the basis's closest point is
        not searched.
        """
    def contains(self, point: _PointLike, tol: Optional[_ToleranceLike] = None) -> bool:
        """Returns whether `point` lies on the patch within `tol`
        (`Tolerance.default()` when omitted).
        """
    def parameters_of(
        self, point: _PointLike, tol: Optional[_ToleranceLike] = None
    ) -> Tuple[float, float]:
        """The basis's `(u, v)` parameters of `point`, a periodic parameter
        shifted by whole periods into its range where it can be
        (`Tolerance.default()` when omitted).
        """
    def bounding_box(self) -> Optional[BoundingBox]:
        """The axis-aligned box around the patch: that of its exact rational
        B-spline patch, or the basis's box for a B-spline or Bézier basis.
        """
    def intersect(
        self, other: _SurfaceMethods[Any, Any] | TrimmedSurface, tol: Optional[_ToleranceLike] = None
    ) -> SurfaceSurfaceIntersection:
        """Intersects with any surface or trimmed surface, in either order.

        Returns a `SurfaceSurfaceIntersection` holding this object as
        `first` and `other` as `second`: the intersection of the bases,
        its curves cut to the pieces inside the ranges of each trimmed side
        (each curve's `parameter_range`), and its contact points kept when
        inside them. Raises `TypeError` for a curve, and for an object that
        is not a geomcore surface. Other Python threads run meanwhile.
        """
    def __repr__(self) -> str: ...


__all__ = [
    "GeomcoreError",
    "ConstructionError",
    "FrameConstructionError",
    "TrimmedCurveConstructionError",
    "TrimmedSurfaceConstructionError",
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
    "QuarticBranch",
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
    "TrimmedSurface",
    "curve_3d_to_uv_within",
    "uv_curve_to_3d_within",
    "curve_curve_extrema",
    "Curve",
    "Curve2D",
    "Surface",
]
