#!/usr/bin/env python3
"""Instruction budget for the criterion benchmarks.

Counts the instructions each benchmark spends inside criterion's
`Bencher::iter` with callgrind, and compares them with `benches/budget.tsv`.
Only the routine is counted: the harness's fixed cost is outside the toggle.

  measure [ids ...]  print the instruction count of each id
                     (default: the budgeted ids; --all for every id)
  check              fail if a budgeted id is more than 1% above its budget
                     (--base REF also prints a per-function diff on failure)
  update [ids ...]   re-measure the budgeted ids (or the given ids, added to
                     the budget) and rewrite benches/budget.tsv
                     (--all budgets every id)

Requires valgrind and callgrind_annotate on the PATH.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BUDGET_PATH = REPO_ROOT / "benches" / "budget.tsv"
# Collect only while criterion's Bencher::iter (and its callees) run.
COLLECT_TOGGLE = "criterion::bencher::Bencher<M>::iter*"
TOLERANCE_PERCENT = 1
DIFF_ROWS = 15
# The count depends on the binary's path, so every run uses this one (see
# stage_target). The environment is fixed too.
STAGING_DIRECTORY = Path("/tmp/geomcore-benchmark-budget")
CLEAN_ENVIRONMENT = {"PATH": "/usr/bin:/bin"}
# The bench binaries are linked without RELRO. LLD, rust's default linker,
# can make the writable segment all RELRO, and valgrind 3.22 then reads no
# symbols of the binary, so the toggle collects nothing. They build in their
# own target directory, so the flag does not rebuild the normal one.
BUILD_FLAGS = "-C link-arg=-Wl,-z,norelro"
BUILD_DIRECTORY = Path("target") / "benchmark-budget"


def run(command, cwd, environment=None):
    """Runs a command and returns its stdout, or exits with its stderr."""
    completed = subprocess.run(
        command, cwd=cwd, env=environment, capture_output=True, text=True
    )
    if completed.returncode != 0:
        command_text = " ".join(map(str, command))
        sys.exit(f"command failed: {command_text}\n{completed.stderr}")
    return completed.stdout


def all_ids(tree):
    """Maps every benchmark id of `tree` to (bench target, binary path)."""
    environment = dict(os.environ)
    environment["RUSTFLAGS"] = f"{environment.get('RUSTFLAGS', '')} {BUILD_FLAGS}".strip()
    environment["CARGO_TARGET_DIR"] = str(Path(tree) / BUILD_DIRECTORY)
    output = run(
        ["cargo", "bench", "--no-run", "--message-format=json"], tree, environment
    )
    ids = {}
    for line in output.splitlines():
        message = json.loads(line)
        if message.get("reason") != "compiler-artifact" or not message.get("executable"):
            continue
        if message["target"]["kind"] != ["bench"]:
            continue
        target, executable = message["target"]["name"], message["executable"]
        listing = run([executable, "--bench", "--list"], REPO_ROOT)
        for listed in listing.splitlines():
            if listed.endswith(": benchmark"):
                ids[listed[: -len(": benchmark")]] = (target, executable)
    return ids


def stage_target(target, executable):
    """Copies a bench binary to the fixed staging path and returns that path.

    The count depends on where the binary lives (the heap layout shifts the
    malloc work inside the routine), so every run uses the same path.
    """
    STAGING_DIRECTORY.mkdir(parents=True, exist_ok=True)
    staged = STAGING_DIRECTORY / target
    shutil.copy2(executable, staged)
    return staged


def profile(executable, bench_id, out_file):
    """Runs one benchmark under callgrind and returns its instruction count."""
    run(
        [
            "valgrind",
            "--tool=callgrind",
            "--collect-atstart=no",
            f"--toggle-collect={COLLECT_TOGGLE}",
            f"--callgrind-out-file={out_file}",
            executable,
            "--bench",
            "--test",
            "--exact",
            bench_id,
        ],
        REPO_ROOT,
        environment=CLEAN_ENVIRONMENT,
    )
    for line in Path(out_file).read_text().splitlines():
        if line.startswith("summary:"):
            count = int(line.split(":", 1)[1])
            if count == 0:
                sys.exit(f"{bench_id}: the toggle collected no instructions")
            return count
    sys.exit(f"{bench_id}: callgrind wrote no summary")


def measure(ids, known, jobs):
    """Counts each id in parallel: {bench id: count}."""
    staged = {}
    for bench_id in ids:
        if bench_id not in known:
            sys.exit(f"unknown benchmark id: {bench_id}")
        target, executable = known[bench_id]
        if target not in staged:
            staged[target] = stage_target(target, executable)
    with tempfile.TemporaryDirectory() as directory:

        def task(index_and_id):
            index, bench_id = index_and_id
            out_file = Path(directory) / f"{index}.callgrind"
            return bench_id, profile(staged[known[bench_id][0]], bench_id, out_file)

        with ThreadPoolExecutor(max_workers=jobs) as pool:
            return dict(pool.map(task, enumerate(ids)))


def read_budget():
    """Returns {bench id: (bench target, budgeted instructions)}."""
    if not BUDGET_PATH.exists():
        return {}
    rows = {}
    for line in BUDGET_PATH.read_text().splitlines()[1:]:
        target, bench_id, instructions = line.split("\t")
        rows[bench_id] = (target, int(instructions))
    return rows


def write_budget(rows):
    lines = ["binary\tid\tinstructions"]
    for bench_id in sorted(rows, key=lambda name: (rows[name][0], name)):
        target, instructions = rows[bench_id]
        lines.append(f"{target}\t{bench_id}\t{instructions}")
    BUDGET_PATH.write_text("\n".join(lines) + "\n")


def command_measure(arguments, jobs):
    known = all_ids(REPO_ROOT)
    if arguments.all:
        ids = sorted(known)
    else:
        ids = arguments.ids or sorted(read_budget())
    counts = measure(ids, known, jobs)
    for bench_id in ids:
        print(f"{counts[bench_id]:>14,}  {bench_id}")
    return 0


def command_check(arguments, jobs):
    budget = read_budget()
    if not budget:
        sys.exit(f"no budget in {BUDGET_PATH}: run `update` first")
    known = all_ids(REPO_ROOT)
    ids = sorted(budget)
    for bench_id in ids:
        if bench_id not in known:
            sys.exit(f"budgeted id no longer exists: {bench_id} (run `update`)")
        if known[bench_id][0] != budget[bench_id][0]:
            sys.exit(f"{bench_id} moved to another bench (run `update`)")
    counts = measure(ids, known, jobs)

    failures, lowered = [], []
    print(f"{'id':<58}{'budget':>15}{'now':>15}{'change':>10}")
    for bench_id in ids:
        budgeted, count = budget[bench_id][1], counts[bench_id]
        change = 100.0 * (count - budgeted) / budgeted
        print(f"{bench_id:<58}{budgeted:>15,}{count:>15,}{change:>+9.2f}%")
        if count * 100 > budgeted * (100 + TOLERANCE_PERCENT):
            failures.append(bench_id)
        elif count < budgeted:
            lowered.append(bench_id)

    if lowered:
        print(f"\n{len(lowered)} id(s) are below budget: lower the budget in this PR "
              "with `python3 scripts/benchmark_budget.py update` (the ratchet).")
    if not failures:
        print(f"\nall {len(ids)} budgeted ids are within budget (+{TOLERANCE_PERCENT}%)")
        return 0
    print(f"\n{len(failures)} id(s) are more than {TOLERANCE_PERCENT}% above budget:")
    print("\n".join(f"  {bench_id}" for bench_id in failures))
    if arguments.base:
        print_base_diff(arguments.base, failures, counts)
    else:
        print("\nRe-run with --base origin/<base branch> for a per-function diff.")
    return 1


def function_costs(out_file):
    """Reads a callgrind profile as {function: instructions} via callgrind_annotate."""
    text = run(
        ["callgrind_annotate", "--auto=no", "--threshold=100", "--show=Ir", out_file],
        REPO_ROOT,
    )
    costs = {}
    row = re.compile(r"^\s*([\d,]+) \(\s*[\d.]+%\)\s+(.+?)\s*$")
    for line in text.splitlines():
        match = row.match(line)
        if match:
            name = re.sub(r"\s+\[[^\]]*\]$", "", match.group(2))
            costs[name] = costs.get(name, 0) + int(match.group(1).replace(",", ""))
    return costs


def print_base_diff(base_ref, failures, pr_counts):
    """Builds `base_ref` in a temporary worktree and prints per-function diffs."""
    with tempfile.TemporaryDirectory() as directory:
        tree = Path(directory) / "base"
        run(["git", "worktree", "add", "--detach", str(tree), base_ref], REPO_ROOT)
        try:
            # Cargo.lock is not tracked; the base gets the PR's dependency versions.
            if (REPO_ROOT / "Cargo.lock").exists():
                shutil.copy2(REPO_ROOT / "Cargo.lock", tree / "Cargo.lock")
            base_known = all_ids(tree)
            pr_known = all_ids(REPO_ROOT)
            for bench_id in failures:
                print(f"\n=== {bench_id}")
                if bench_id not in base_known:
                    print(f"not in {base_ref}; PR count {pr_counts[bench_id]:,}")
                    continue
                with tempfile.TemporaryDirectory() as profiles:
                    base_profile, pr_profile = Path(profiles) / "base", Path(profiles) / "pr"
                    base_count = profile(
                        stage_target(*base_known[bench_id]), bench_id, base_profile
                    )
                    pr_count = profile(
                        stage_target(*pr_known[bench_id]), bench_id, pr_profile
                    )
                    print(f"base {base_count:,}  PR {pr_count:,}  "
                          f"({100.0 * (pr_count - base_count) / base_count:+.2f}%)")
                    print_function_table(function_costs(base_profile), function_costs(pr_profile))
        finally:
            run(["git", "worktree", "remove", "--force", str(tree)], REPO_ROOT)


def print_function_table(base, pr):
    """Prints the functions whose instruction counts changed most."""
    changes = sorted(
        ((pr.get(name, 0) - base.get(name, 0), name) for name in set(base) | set(pr)),
        key=lambda change: -abs(change[0]),
    )
    print(f"{'base':>14}{'PR':>14}{'change':>14}  function")
    for delta, name in [change for change in changes if change[0] != 0][:DIFF_ROWS]:
        print(f"{base.get(name, 0):>14,}{pr.get(name, 0):>14,}{delta:>+14,}  {name}")


def command_update(arguments, jobs):
    rows = read_budget()
    known = all_ids(REPO_ROOT)
    if arguments.all:
        ids = sorted(known)
    else:
        ids = arguments.ids or sorted(rows)
    if not ids:
        sys.exit("no ids to budget: pass them as arguments or use --all")
    counts = measure(ids, known, jobs)
    for bench_id in ids:
        rows[bench_id] = (known[bench_id][0], counts[bench_id])
    write_budget(rows)
    print(f"wrote {len(rows)} budgeted ids to {BUDGET_PATH.relative_to(REPO_ROOT)}")
    return 0


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--jobs", type=int, default=os.cpu_count(), help="parallel valgrind runs")
    commands = parser.add_subparsers(dest="command", required=True)

    measure_parser = commands.add_parser("measure", help="print instruction counts")
    measure_parser.add_argument("ids", nargs="*")
    measure_parser.add_argument("--all", action="store_true", help="every benchmark id")
    measure_parser.set_defaults(run=command_measure)

    check_parser = commands.add_parser("check", help="compare with benches/budget.tsv")
    check_parser.add_argument("--base", help="git ref to diff against on failure")
    check_parser.set_defaults(run=command_check)

    update_parser = commands.add_parser("update", help="rewrite benches/budget.tsv")
    update_parser.add_argument("ids", nargs="*", help="ids to add or re-measure")
    update_parser.add_argument("--all", action="store_true", help="budget every benchmark id")
    update_parser.set_defaults(run=command_update)

    arguments = parser.parse_args()
    if shutil.which("valgrind") is None or shutil.which("callgrind_annotate") is None:
        sys.exit("valgrind and callgrind_annotate must be on the PATH")
    sys.exit(arguments.run(arguments, arguments.jobs))


if __name__ == "__main__":
    main()
