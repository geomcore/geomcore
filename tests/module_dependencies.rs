//! Dependency rules between the crate's modules, checked on the sources.
//!
//! The internal math sits under the shapes: no production file under
//! `src/math/` names the `geometry`, `operations` or `bounding_box`
//! modules, nor any name the crate root re-exports from them. It works on
//! the values of `base/` and on its own types.
//!
//! The values sit under everything else: no production file under
//! `src/base/` names the `math`, `geometry`, `operations` or `bounding_box`
//! modules, nor any name the crate root re-exports from them. The pole
//! algebra the B-spline engine runs on lives in `base/pole.rs`, and the
//! box types in `base/bounding_box.rs` and `base/bounding_box_2d.rs`, so
//! the tier trait requires them without reaching up.
//!
//! Curves do not depend on surfaces: a surface owns the mapping between its
//! UV curves and 3D curves both ways (ADRs 0030, 0031), so no production
//! file under `src/geometry/curves/` names `crate::geometry::surfaces`, nor
//! any name the crate root re-exports from it, such as `Plane` reached
//! through `crate::Plane`.
//!
//! The shapes do not depend on the intersection router: no production
//! file under `src/geometry/` names `intersect::routing`. The rules of a
//! parameter range they share with it live in `src/math/` (#291).
//!
//! Only the files under `src/geometry/curves/` are read for the curves.
//! `QuarticBranch` lives beside the quartic it is a piece of, under
//! `src/operations/intersect/exact/quartic/`, and holds the cones or the
//! cone and cylinder it lies on (ADR 0042): through that variant the
//! `Curve` enum reaches surfaces, and this check does not see it.

use std::fs;
use std::path::{Path, PathBuf};

/// Every production `.rs` file under `directory`, recursively. Test code, a
/// `tests.rs` file or anything under a `tests/` directory, may build
/// surfaces.
fn production_files(directory: &Path) -> Vec<PathBuf> {
    let mut files = Vec::new();
    for entry in fs::read_dir(directory).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_str().unwrap();
        if path.is_dir() {
            if name != "tests" {
                files.extend(production_files(&path));
            }
        } else if name.ends_with(".rs") && name != "tests.rs" {
            files.push(path);
        }
    }
    files
}

/// The names `src/lib.rs` re-exports from the module at `path` (such as
/// `geometry::surfaces` or `operations`), read from its `pub use` statements,
/// braced or not.
fn reexported_names(root: &Path, path: &str) -> Vec<String> {
    let lib = fs::read_to_string(root.join("src/lib.rs")).unwrap();
    let marker = format!("pub use {path}::");
    let mut names = Vec::new();
    let mut remaining = lib.as_str();
    while let Some(start) = remaining.find(&marker) {
        let statement = &remaining[start + marker.len()..];
        let end = statement.find(';').unwrap();
        let body = &statement[..end];
        let list = match body.find('{') {
            Some(open) => body[open + 1..].trim_end_matches('}'),
            None => body.rsplit("::").next().unwrap(),
        };
        names.extend(
            list.split(',')
                .map(str::trim)
                .filter(|name| !name.is_empty())
                .map(String::from),
        );
        remaining = &statement[end..];
    }
    names
}

/// The code of `source` without its comments: every line from its first
/// `//` on is dropped, so a doc comment may name the shape a function
/// serves.
fn code_of(source: &str) -> String {
    source
        .lines()
        .map(|line| line.find("//").map_or(line, |at| &line[..at]))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Whether `word` appears in `text` as a whole identifier.
fn names_word(text: &str, word: &str) -> bool {
    let is_identifier = |c: char| c.is_alphanumeric() || c == '_';
    text.match_indices(word).any(|(at, _)| {
        let before = text[..at].chars().next_back();
        let after = text[at + word.len()..].chars().next();
        !before.is_some_and(is_identifier) && !after.is_some_and(is_identifier)
    })
}

#[test]
fn no_curve_source_names_the_surfaces_module_or_a_surface_type() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let surface_names = reexported_names(root, "geometry::surfaces");
    // A parse that found nothing would pass every file, so pin the list.
    for expected in [
        "Surface",
        "ElementarySurface",
        "ParametricSurface",
        "Plane",
        "Sphere",
        "ExtrusionSurface",
        "RevolutionSurface",
        "BSplineSurface",
    ] {
        assert!(
            surface_names.iter().any(|name| name == expected),
            "{expected} missing from {surface_names:?}"
        );
    }

    let files = production_files(&root.join("src/geometry/curves"));
    assert!(files.len() > 10, "only {} curve files", files.len());
    let mut offenders = Vec::new();
    for path in &files {
        let source = fs::read_to_string(path).unwrap();
        if source.contains("crate::geometry::surfaces") {
            offenders.push(format!("{}: crate::geometry::surfaces", path.display()));
        }
        for name in &surface_names {
            if names_word(&source, name) {
                offenders.push(format!("{}: {name}", path.display()));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "curves must not depend on surfaces: {offenders:#?}"
    );
}

#[test]
fn no_math_source_names_a_shape_or_an_operation() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let modules = ["geometry", "operations", "bounding_box"];
    let names: Vec<String> = modules
        .iter()
        .flat_map(|module| reexported_names(root, module))
        .collect();
    // A parse that found nothing would pass every file, so pin the list.
    for expected in [
        "BSplineCurve",
        "Curve",
        "Plane",
        "CurveApproximation",
        "HasBoundingBox",
        "Discretize",
        "Intersect",
        "Completeness",
        "CurveProjection",
        "curve_curve_extrema",
    ] {
        assert!(
            names.iter().any(|name| name == expected),
            "{expected} missing from {names:?}"
        );
    }

    let files = production_files(&root.join("src/math"));
    assert!(files.len() > 10, "only {} math files", files.len());
    let mut offenders = Vec::new();
    for path in &files {
        let code = code_of(&fs::read_to_string(path).unwrap());
        for module in modules {
            if code.contains(&format!("crate::{module}")) {
                offenders.push(format!("{}: crate::{module}", path.display()));
            }
        }
        for name in &names {
            if names_word(&code, name) {
                offenders.push(format!("{}: {name}", path.display()));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "math must not depend on the shapes or the operations: {offenders:#?}"
    );
}

#[test]
fn no_base_source_names_a_module_above_it() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let modules = ["math", "geometry", "operations", "bounding_box"];
    let names: Vec<String> = modules
        .iter()
        .flat_map(|module| reexported_names(root, module))
        .collect();
    // A parse that found nothing would pass every file, so pin the list.
    for expected in [
        "BSplineConstructionError",
        "KnotRemovalError",
        "BSplineCurve",
        "Curve",
        "Plane",
        "HasBoundingBox",
        "HasBoundingBox2D",
        "Discretize",
        "Intersect",
        "CurveProjection",
    ] {
        assert!(
            names.iter().any(|name| name == expected),
            "{expected} missing from {names:?}"
        );
    }

    let files = production_files(&root.join("src/base"));
    assert!(files.len() > 5, "only {} base files", files.len());
    let mut offenders = Vec::new();
    for path in &files {
        let code = code_of(&fs::read_to_string(path).unwrap());
        for module in modules {
            if code.contains(&format!("crate::{module}")) {
                offenders.push(format!("{}: crate::{module}", path.display()));
            }
        }
        for name in &names {
            if names_word(&code, name) {
                offenders.push(format!("{}: {name}", path.display()));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "base must not depend on a module above it: {offenders:#?}"
    );
}

#[test]
fn no_shape_source_names_the_intersection_router() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let files = production_files(&root.join("src/geometry"));
    assert!(files.len() > 20, "only {} shape files", files.len());
    let offenders: Vec<String> = files
        .iter()
        .filter(|path| code_of(&fs::read_to_string(path).unwrap()).contains("intersect::routing"))
        .map(|path| path.display().to_string())
        .collect();
    assert!(
        offenders.is_empty(),
        "the shapes must not depend on the intersection router: {offenders:#?}"
    );
}
