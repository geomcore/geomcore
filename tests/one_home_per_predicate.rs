//! Each geometric question has one home (ARCHITECTURE.md, rule 7), checked
//! on the sources: no production file outside their homes, [`HOMES`],
//! writes out, on one line,
//!
//! - a parallel test, a cross product's length against an angular bound
//!   (`Direction::is_parallel`, `Vector::is_parallel`);
//! - a perpendicular test, a dot product's absolute value against an
//!   angular bound (`Direction::is_perpendicular`,
//!   `Vector::is_perpendicular`);
//! - the part of `v` across `d` as `v - d * v.dot(d)`
//!   (`Vector::component_perpendicular_to`).
//!
//! The check reads one line at a time, so it catches the copies as they are
//! usually written, not every way of spelling them.

use std::fs;
use std::path::{Path, PathBuf};

/// The files the predicates live in: the vector forms, the direction forms
/// and the tier-generic form the curves of both tiers share.
const HOMES: [&str; 3] = [
    "src/base/vector.rs",
    "src/base/direction.rs",
    "src/base/tier.rs",
];

/// Every production `.rs` file under `directory`, recursively, test code
/// (a `tests.rs` file, a `tests_*.rs` file or anything under a `tests/`
/// directory) left out.
fn production_files(directory: &Path) -> Vec<PathBuf> {
    let mut files = Vec::new();
    for entry in fs::read_dir(directory).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_str().unwrap();
        if path.is_dir() {
            if name != "tests" {
                files.extend(production_files(&path));
            }
        } else if name.ends_with(".rs") && name != "tests.rs" && !name.starts_with("tests_") {
            files.push(path);
        }
    }
    files
}

/// The identifier that ends just before byte `end` of `line`.
fn identifier_before(line: &str, end: usize) -> &str {
    let start = line[..end]
        .rfind(|c: char| !(c.is_alphanumeric() || c == '_'))
        .map_or(0, |at| at + 1);
    &line[start..end]
}

/// The copy of a predicate `line` writes out by hand, if any.
fn hand_written_predicate(line: &str) -> Option<&'static str> {
    if line.contains(".cross(") && line.contains("magnitude()") && line.contains("angular") {
        return Some("a parallel test: use Direction::is_parallel or Vector::is_parallel");
    }
    if line.contains(".dot(") && line.contains(".abs()") && line.contains("angular") {
        return Some(
            "a perpendicular test: use Direction::is_perpendicular or Vector::is_perpendicular",
        );
    }
    for (at, _) in line.match_indices(".dot(") {
        let vector = identifier_before(line, at);
        let argument = &line[at + ".dot(".len()..];
        let Some(close) = argument.find(')') else {
            continue;
        };
        let direction = &argument[..close];
        if !vector.is_empty()
            && line.contains(&format!(
                "{vector} - {direction} * {vector}.dot({direction})"
            ))
        {
            return Some("the part across a direction: use Vector::component_perpendicular_to");
        }
    }
    None
}

#[test]
fn the_check_recognises_each_copy() {
    assert!(hand_written_predicate("if a.cross(b).magnitude() > tol.angular {").is_some());
    assert!(hand_written_predicate("x.dot(y).abs() <= tol.angular").is_some());
    assert!(hand_written_predicate("let p = w - a1 * w.dot(a1);").is_some());
    assert!(hand_written_predicate("let p = w.component_perpendicular_to(a1);").is_none());
    assert!(hand_written_predicate("if !a.is_parallel(b, tol.angular) {").is_none());
}

#[test]
fn no_solver_writes_a_predicate_by_hand() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let files = production_files(&root.join("src"));
    assert!(files.len() > 50, "only {} source files", files.len());
    let mut offenders = Vec::new();
    for path in &files {
        if HOMES.iter().any(|home| path.ends_with(home)) {
            continue;
        }
        let source = fs::read_to_string(path).unwrap();
        for (number, line) in source.lines().enumerate() {
            let code = line.find("//").map_or(line, |at| &line[..at]);
            if let Some(what) = hand_written_predicate(code) {
                offenders.push(format!("{}:{}: {what}", path.display(), number + 1));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "geometric predicates written by hand: {offenders:#?}"
    );
}
