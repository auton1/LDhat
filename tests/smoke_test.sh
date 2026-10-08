#!/usr/bin/env bash
#
# Smoke test for LDhat. Not a correctness/regression test suite (LDhat has
# none) -- this just builds every program and runs each one through a
# minimal, non-interactive, real invocation to confirm it starts, reads its
# inputs, runs to completion and writes output, i.e. that nothing crashes,
# hangs waiting on stdin, or silently fails after a change to the code.
#
# Usage:
#   tests/smoke_test.sh            # build + run all smoke tests
#   tests/smoke_test.sh --no-build # skip `make`, test whatever is already built
#   tests/smoke_test.sh --keep     # keep the temporary output directory around
#
# Exit status is 0 if every test passed, 1 otherwise.

set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
EXAMPLE_DIR="$REPO_DIR/Example"

DO_BUILD=1
KEEP_OUTPUT=0
for arg in "$@"; do
	case "$arg" in
		--no-build) DO_BUILD=0 ;;
		--keep) KEEP_OUTPUT=1 ;;
		-h|--help)
			sed -n '2,15p' "$0"
			exit 0
			;;
		*)
			echo "Unknown option: $arg" >&2
			exit 2
			;;
	esac
done

OUT_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ldhat_smoke.XXXXXX")"
LOG_DIR="$OUT_DIR/logs"
mkdir -p "$LOG_DIR"

cleanup() {
	if [ "$KEEP_OUTPUT" -eq 1 ]; then
		echo "Output kept at: $OUT_DIR"
	else
		rm -rf "$OUT_DIR"
	fi
}
trap cleanup EXIT

PASS=0
FAIL=0
FAILED_NAMES=()

# run_test <name> <timeout_seconds> <log_basename> <command...>
# Runs a command with a timeout, records pass/fail, always continues.
run_test() {
	local name="$1"; shift
	local secs="$1"; shift
	local logfile="$LOG_DIR/$1"; shift
	if timeout "$secs" "$@" >"$logfile" 2>&1; then
		echo "PASS: $name"
		PASS=$((PASS + 1))
		return 0
	else
		local rc=$?
		if [ "$rc" -eq 124 ]; then
			echo "FAIL: $name (timed out after ${secs}s -- see $logfile)"
		else
			echo "FAIL: $name (exit $rc -- see $logfile)"
		fi
		FAIL=$((FAIL + 1))
		FAILED_NAMES+=("$name")
		return 1
	fi
}

# check_file <name> <path>
# Fails the most recently run test's name if the file is missing/empty.
require_file() {
	local name="$1" path="$2"
	if [ ! -s "$path" ]; then
		echo "FAIL: $name (expected output file missing or empty: $path)"
		FAIL=$((FAIL + 1))
		FAILED_NAMES+=("$name")
		PASS=$((PASS - 1))
	fi
}

require_stdout_contains() {
	local name="$1" logfile="$2" needle="$3"
	if ! grep -q "$needle" "$logfile"; then
		echo "FAIL: $name (expected to see '$needle' in output -- see $logfile)"
		FAIL=$((FAIL + 1))
		FAILED_NAMES+=("$name")
		PASS=$((PASS - 1))
	fi
}

cd "$REPO_DIR" || exit 1

if [ "$DO_BUILD" -eq 1 ]; then
	echo "== Building =="
	rm -f ./*.o convert pairwise interval lkgen complete stat fin rhomap
	if ! make -k >"$LOG_DIR/build.log" 2>&1; then
		echo "FAIL: build (see $LOG_DIR/build.log)"
		echo "Cannot continue without the binaries."
		exit 1
	fi
	WARN_COUNT=$(grep -c "warning:" "$LOG_DIR/build.log" || true)
	echo "Build finished with $WARN_COUNT compiler warning(s)."
fi

BINARIES=(convert pairwise interval lkgen complete stat fin rhomap)
MISSING=0
for b in "${BINARIES[@]}"; do
	if [ ! -x "$REPO_DIR/$b" ]; then
		echo "FAIL: missing binary '$b'"
		MISSING=1
	fi
done
if [ "$MISSING" -eq 1 ]; then
	echo "Cannot continue: not all binaries were built."
	exit 1
fi

echo
echo "== Usage (-h) smoke tests =="
for b in "${BINARIES[@]}"; do
	run_test "$b -h" 10 "${b}_help.log" "./$b" -h
done

echo
echo "== Functional smoke tests (Example/lpl_fn dataset + synthetic inputs) =="

# --- convert: synthetic tiny FASTA-style input (own header + sequence pairs,
# same format read_fasta expects, matching Example/lpl_fn.sites' layout) ---
CONVERT_IN="$OUT_DIR/tiny.sites"
cat >"$CONVERT_IN" <<'EOF'
6 20 1
>Seq1
ACGTACGTACGTACGTACGT
>Seq2
ACGTACGTACGTACGTACGA
>Seq3
ACGTACGTACGTACGTACGT
>Seq4
ACGAACGTACGTACGTACGT
>Seq5
ACGTACGTACGTACGTACGT
>Seq6
ACGTACGTACGTACGTACGA
EOF
if run_test "convert" 15 "convert.log" \
	./convert -seq "$CONVERT_IN" -prefix "$OUT_DIR/conv_"; then
	require_file "convert" "${OUT_DIR}/conv_sites.txt"
	require_file "convert" "${OUT_DIR}/conv_locs.txt"
fi

# --- lkgen: Example/lpl_fn.lk is NOT an exhaustive table (it only holds
# pair-types actually observed in one prior run: header says "48 400" but an
# exhaustive table for n=48 needs 2900), so it can't be used as -lk with real
# data that has missing bases, and it can't be shrunk further either. Instead
# derive a genuinely exhaustive n=48 table by shrinking the exhaustive n=50
# table shipped in lk_files/ -- this doubles as the lkgen functional test and
# feeds pairwise/interval/rhomap below. ---
LK50_GZ="$REPO_DIR/lk_files/lk_n50_t0.001.gz"
LK50="$OUT_DIR/lk_n50.txt"
LK48="$OUT_DIR/lk48_new_lk.txt"
if [ ! -f "$LK50_GZ" ]; then
	echo "FAIL: lkgen (missing $LK50_GZ)"
	FAIL=$((FAIL + 1))
	FAILED_NAMES+=("lkgen")
elif ! gunzip -c "$LK50_GZ" >"$LK50" 2>"$LOG_DIR/gunzip.log"; then
	echo "FAIL: lkgen (could not decompress $LK50_GZ -- see $LOG_DIR/gunzip.log)"
	FAIL=$((FAIL + 1))
	FAILED_NAMES+=("lkgen")
elif run_test "lkgen" 60 "lkgen.log" \
	./lkgen -lk "$LK50" -nseq 48 -prefix "$OUT_DIR/lk48_"; then
	require_file "lkgen" "$LK48"
fi

# --- pairwise: needs six interactive yes/no answers regardless of flags
# (hardcoded prompts in pairdip.c) -- answer "no" to all of them. Do NOT pass
# -exact: the LPL sequences contain missing bases ('-'), so even an
# exhaustive table will show pnew/miss > 0, and -exact asserts there should
# be none of those (a different, legitimate code path handles resolving
# missing data against an exhaustive table when -exact is left off). ---
if [ -s "$LK48" ]; then
	if printf '0\n0\n0\n0\n0\n0\n' | timeout 30 ./pairwise \
		-seq "$EXAMPLE_DIR/lpl_fn.sites" -loc "$EXAMPLE_DIR/lpl_fn.locs" \
		-lk "$LK48" -prefix "$OUT_DIR/pw_" \
		>"$LOG_DIR/pairwise.log" 2>&1
	then
		echo "PASS: pairwise"
		PASS=$((PASS + 1))
		require_file "pairwise" "${OUT_DIR}/pw_type_table.txt"
		require_file "pairwise" "${OUT_DIR}/pw_fit.txt"
	else
		rc=$?
		echo "FAIL: pairwise (exit $rc -- see $LOG_DIR/pairwise.log)"
		FAIL=$((FAIL + 1))
		FAILED_NAMES+=("pairwise")
	fi
else
	echo "SKIP: pairwise (no exhaustive n=48 lookup table from lkgen step)"
fi

# --- interval: MCMC over the LPL example; -its is kept at the minimum
# interval accepts (BURNIN*2 in ldhat.h) but it's still well under a second
# for this dataset ---
if [ -s "$LK48" ]; then
	if run_test "interval" 60 "interval.log" \
		./interval -seq "$EXAMPLE_DIR/lpl_fn.sites" -loc "$EXAMPLE_DIR/lpl_fn.locs" \
		-lk "$LK48" -bpen 5 -its 200000 -samp 2000 \
		-prefix "$OUT_DIR/int_"; then
		require_file "interval" "${OUT_DIR}/int_rates.txt"
	fi
else
	echo "SKIP: interval (no exhaustive n=48 lookup table from lkgen step)"
fi

# --- stat: summarise the rates.txt that interval just produced ---
if [ -s "$OUT_DIR/int_rates.txt" ]; then
	run_test "stat" 15 "stat.log" \
		./stat -input "$OUT_DIR/int_rates.txt" -loc "$EXAMPLE_DIR/lpl_fn.locs" \
		-burn 0 -prefix "$OUT_DIR/st_"
else
	echo "SKIP: stat (no rates.txt from interval to summarise)"
fi

# --- complete: generate a small complete-enumeration lookup table ---
if run_test "complete" 60 "complete.log" \
	./complete -n 10 -rhomax 10 -n_pts 11 -theta 0.01 -prefix "$OUT_DIR/cp_"; then
	require_file "complete" "${OUT_DIR}/cp_new_lk.txt"
fi

# --- fin: simulate a small dataset (writes results to stdout by default) ---
if run_test "fin" 30 "fin.log" \
	./fin -nsamp 10 -len 20 -theta 0.01 -R 5 -prefix "$OUT_DIR/fs_"; then
	require_stdout_contains "fin" "$LOG_DIR/fin.log" "Results of"
fi

# --- rhomap: hotspot MCMC over the LPL example, kept short for a smoke test ---
if [ -s "$LK48" ]; then
	if run_test "rhomap" 60 "rhomap.log" \
		./rhomap -seq "$EXAMPLE_DIR/lpl_fn.sites" -loc "$EXAMPLE_DIR/lpl_fn.locs" \
		-lk "$LK48" -bpen 5 -hpen 5 -its 2000 -samp 200 -burn 0 \
		-prefix "$OUT_DIR/rh_" -nonverbose; then
		require_file "rhomap" "${OUT_DIR}/rh_rates.txt"
	fi
else
	echo "SKIP: rhomap (no exhaustive n=48 lookup table from lkgen step)"
fi

echo
echo "== Summary =="
echo "$PASS passed, $FAIL failed"
if [ "$FAIL" -ne 0 ]; then
	printf 'Failed: %s\n' "${FAILED_NAMES[@]}"
	exit 1
fi
exit 0
