//! The enumerations Python sees as real `enum` classes.
//!
//! A PyO3 enum class cannot list its members, so each enumeration here is a
//! plain Rust enum paired with a Python class built with the functional API
//! of the `enum` module: `IntEnum` for the ones that were integers, and an
//! `Enum` mixed with `str` for `PointParametrization`, whose members compare
//! equal to the strings it also accepts. A Rust value converts to its member
//! by name; a Python argument converts by calling the class on it, which
//! takes a member or a member's value, as `enum` does.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyDict, PyModule};

use geomcore::{IncompleteReason, IntersectionMethod, PointParametrization, Side};

/// Declares a Rust enum and its Python `enum` class: the class's name, its
/// base in the `enum` module, an optional builtin mixed in, and each
/// variant's member name and value.
macro_rules! python_enum {
    (
        $(#[doc = $doc:literal])*
        $rust:ident => $class:literal ($base:literal, mixin $mixin:expr) {
            $($variant:ident => $member:literal = $value:expr),+ $(,)?
        }
    ) => {
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
        pub(crate) enum $rust {
            $($variant),+
        }

        impl $rust {
            /// The member's name, as written after the dot.
            fn member_name(self) -> &'static str {
                match self {
                    $($rust::$variant => $member),+
                }
            }

            /// The member as `Class.NAME`, as the result reprs write it.
            #[allow(dead_code, reason = "only the result enumerations are in a repr")]
            pub(crate) fn repr(self) -> &'static str {
                match self {
                    $($rust::$variant => concat!($class, ".", $member)),+
                }
            }

            /// The Python class, built on first use.
            fn class(py: Python<'_>) -> PyResult<&Bound<'_, PyAny>> {
                static CLASS: PyOnceLock<Py<PyAny>> = PyOnceLock::new();
                CLASS
                    .get_or_try_init(py, || {
                        let doc: &[&str] = &[$($doc),*];
                        let mixin: Option<&str> = $mixin;
                        build_class(py, $class, $base, mixin, doc, vec![
                            $(($member, $value.into_pyobject(py)?.into_any().unbind())),+
                        ])
                    })
                    .map(|class| class.bind(py))
            }

            /// Adds the class to `module`.
            fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
                module.add($class, $rust::class(module.py())?)
            }
        }

        impl<'py> IntoPyObject<'py> for $rust {
            type Target = PyAny;
            type Output = Bound<'py, PyAny>;
            type Error = PyErr;

            fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
                $rust::class(py)?.getattr(self.member_name())
            }
        }

        impl<'py> FromPyObject<'py> for $rust {
            fn extract_bound(object: &Bound<'py, PyAny>) -> PyResult<Self> {
                let member = $rust::class(object.py())?.call1((object,))?;
                let name: String = member.getattr("name")?.extract()?;
                match name.as_str() {
                    $($member => Ok($rust::$variant),)+
                    _ => Err(PyValueError::new_err(format!(
                        "{name} is not a member of {}",
                        $class
                    ))),
                }
            }
        }
    };
}

/// `__str__` and `__format__` for every enumeration: the member as
/// `Class.NAME`. Since Python 3.11 an `IntEnum` prints its number from
/// `str()`, and before it from an f-string, so a printed status read `2`.
const MEMBER_TEXT: &std::ffi::CStr = c"
def __str__(self):
    return f'{type(self).__name__}.{self.name}'

def __format__(self, format_spec):
    return format(str(self), format_spec)

def _missing_(cls, value):
    valid = ', '.join(repr(member.value) for member in cls)
    raise ValueError(f'{value!r} is not a valid {cls.__name__}; use a member or one of {valid}')
";

/// The functions of [`MEMBER_TEXT`], compiled once.
fn member_text(py: Python<'_>) -> PyResult<&Bound<'_, PyModule>> {
    static MODULE: PyOnceLock<Py<PyModule>> = PyOnceLock::new();
    MODULE
        .get_or_try_init(py, || {
            Ok::<_, PyErr>(
                PyModule::from_code(py, MEMBER_TEXT, c"geomcore_enums.py", c"geomcore_enums")?
                    .unbind(),
            )
        })
        .map(|module| module.bind(py))
}

/// Builds `enum.<base>(class, members, module="geomcore", type=mixin)`,
/// gives it the docstring `doc`, and makes its members print as
/// `Class.NAME` (see [`MEMBER_TEXT`]).
fn build_class(
    py: Python<'_>,
    class: &str,
    base: &str,
    mixin: Option<&str>,
    doc: &[&str],
    members: Vec<(&str, Py<PyAny>)>,
) -> PyResult<Py<PyAny>> {
    let options = PyDict::new(py);
    options.set_item("module", "geomcore")?;
    options.set_item("qualname", class)?;
    if let Some(mixin) = mixin {
        options.set_item("type", py.import("builtins")?.getattr(mixin)?)?;
    }
    let built = py
        .import("enum")?
        .getattr(base)?
        .call((class, members), Some(&options))?;
    let lines: Vec<&str> = doc
        .iter()
        .map(|line| line.strip_prefix(' ').unwrap_or(line))
        .collect();
    built.setattr("__doc__", lines.join("\n"))?;
    let text = member_text(py)?;
    for method in ["__str__", "__format__"] {
        built.setattr(method, text.getattr(method)?)?;
    }
    let classmethod = py.import("builtins")?.getattr("classmethod")?;
    built.setattr(
        "_missing_",
        classmethod.call1((text.getattr("_missing_")?,))?,
    )?;
    Ok(built.unbind())
}

python_enum! {
    /// Which side of a parameter to evaluate from, for one-sided limits at
    /// B-spline knots and periodic seams.
    ///
    /// `BELOW` approaches from the span ending at the parameter, `ABOVE`
    /// from the span starting there (the side plain evaluation takes,
    /// except at the upper bound). Away from knots and seams both sides
    /// agree.
    PySide => "Side" ("IntEnum", mixin None) {
        Below => "BELOW" = 0u8,
        Above => "ABOVE" = 1u8,
    }
}

impl From<PySide> for Side {
    fn from(side: PySide) -> Self {
        match side {
            PySide::Below => Side::Below,
            PySide::Above => Side::Above,
        }
    }
}

python_enum! {
    /// How the parameters of fitted points are spaced: `CHORDAL` by chord
    /// length, `CENTRIPETAL` by its square root (tamer at sharp turns, the
    /// default). A member equals its lowercase name, which is accepted in
    /// its place.
    PyPointParametrization => "PointParametrization" ("Enum", mixin Some("str")) {
        Chordal => "CHORDAL" = "chordal",
        Centripetal => "CENTRIPETAL" = "centripetal",
    }
}

impl From<PyPointParametrization> for PointParametrization {
    fn from(parametrization: PyPointParametrization) -> Self {
        match parametrization {
            PyPointParametrization::Chordal => PointParametrization::Chordal,
            PyPointParametrization::Centripetal => PointParametrization::Centripetal,
        }
    }
}

python_enum! {
    /// How an intersection was computed.
    ///
    /// `EXACT` is a closed form, always complete. `NUMERIC` is a search,
    /// whose `completeness` tells whether it proved that nothing was
    /// missed.
    PyIntersectionMethod => "IntersectionMethod" ("IntEnum", mixin None) {
        Exact => "EXACT" = 0u8,
        Numeric => "NUMERIC" = 1u8,
    }
}

impl From<IntersectionMethod> for PyIntersectionMethod {
    fn from(method: IntersectionMethod) -> Self {
        match method {
            IntersectionMethod::Exact => PyIntersectionMethod::Exact,
            IntersectionMethod::Numeric => PyIntersectionMethod::Numeric,
        }
    }
}

python_enum! {
    /// Whether a numeric search proved that it found everything.
    ///
    /// - `COMPLETE`: nothing was missed (always the case for an exact
    ///   result).
    /// - `INCOMPLETE`: the search may have missed something inside its
    ///   `regions`, and found everything outside them.
    /// - `SAMPLED`: the result comes from samples; nothing is proved.
    PyCompletenessStatus => "CompletenessStatus" ("IntEnum", mixin None) {
        Complete => "COMPLETE" = 0u8,
        Incomplete => "INCOMPLETE" = 1u8,
        Sampled => "SAMPLED" = 2u8,
    }
}

python_enum! {
    /// Why an incomplete search stopped.
    ///
    /// - `BUDGET_EXHAUSTED`: it ran out of subdivisions.
    /// - `UNSEPARATED`: pieces of the operands stayed within tolerance of
    ///   each other down to the smallest piece, as along an overlap or at a
    ///   tangent contact of a curve with a surface.
    /// - `UNBOUNDED`: a branch of two unbounded surfaces runs beyond the
    ///   box it was traced in, and continues past the regions reported.
    PyIncompleteReason => "IncompleteReason" ("IntEnum", mixin None) {
        BudgetExhausted => "BUDGET_EXHAUSTED" = 0u8,
        Unseparated => "UNSEPARATED" = 1u8,
        Unbounded => "UNBOUNDED" = 2u8,
    }
}

impl From<IncompleteReason> for PyIncompleteReason {
    fn from(reason: IncompleteReason) -> Self {
        match reason {
            IncompleteReason::BudgetExhausted => PyIncompleteReason::BudgetExhausted,
            IncompleteReason::Unseparated => PyIncompleteReason::Unseparated,
            IncompleteReason::Unbounded => PyIncompleteReason::Unbounded,
        }
    }
}

/// Adds every enumeration class to the root module.
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    PySide::register(module)?;
    PyPointParametrization::register(module)?;
    PyIntersectionMethod::register(module)?;
    PyCompletenessStatus::register(module)?;
    PyIncompleteReason::register(module)
}
