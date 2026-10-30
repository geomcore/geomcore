//! Dependency rules between the crate's modules, checked on the sources.
//! Curves do not depend on surfaces: a surface owns the mapping between its
//! UV curves and 3D curves both ways (ADR 0030), so no file under
//! `src/curves/` names `crate::surfaces`.

use std::fs;
use std::path::{Path, PathBuf};

/// Every `.rs` file under `directory`, recursively.
fn rust_files(directory: &Path) -> Vec<PathBuf> {
    let mut files = Vec::new();
    for entry in fs::read_dir(directory).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            files.extend(rust_files(&path));
        } else if path.extension().is_some_and(|extension| extension == "rs") {
            files.push(path);
        }
    }
    files
}

#[test]
fn no_curve_source_names_the_surfaces_module() {
    let curves = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/curves");
    let files = rust_files(&curves);
    assert!(
        files.len() > 10,
        "only {} files under {curves:?}",
        files.len()
    );
    let offenders: Vec<_> = files
        .iter()
        .filter(|path| {
            fs::read_to_string(path)
                .unwrap()
                .contains("crate::surfaces")
        })
        .collect();
    assert!(
        offenders.is_empty(),
        "curves must not depend on surfaces: {offenders:?}"
    );
}
