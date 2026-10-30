//! Dependency rules between the crate's modules, checked on the sources.
//! Curves do not depend on surfaces: a surface owns the mapping between its
//! UV curves and 3D curves both ways (ADRs 0030, 0031), so no production
//! file under `src/curves/` names `crate::surfaces`, nor any name the crate
//! root re-exports from it, such as `Plane` reached through `crate::Plane`.

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

/// The names `src/lib.rs` re-exports from `surfaces`, read from its
/// `pub use surfaces::...;` statements, braced or not.
fn surface_names(root: &Path) -> Vec<String> {
    let lib = fs::read_to_string(root.join("src/lib.rs")).unwrap();
    let marker = "pub use surfaces::";
    let mut names = Vec::new();
    let mut remaining = lib.as_str();
    while let Some(start) = remaining.find(marker) {
        let statement = &remaining[start + marker.len()..];
        let end = statement.find(';').unwrap();
        let list = statement[..end].trim_matches(|c| c == '{' || c == '}');
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
    let surface_names = surface_names(root);
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

    let files = production_files(&root.join("src/curves"));
    assert!(files.len() > 10, "only {} curve files", files.len());
    let mut offenders = Vec::new();
    for path in &files {
        let source = fs::read_to_string(path).unwrap();
        if source.contains("crate::surfaces") {
            offenders.push(format!("{}: crate::surfaces", path.display()));
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
