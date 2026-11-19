#!/usr/bin/env bash

# End-to-end test of phycoeval's checkpointing and resuming.
#
# Checks that
#
#   A. an uninterrupted run (no checkpointing options) with a fixed seed,
#   B. the same run stopped part-way with SIGTERM (and, separately, killed
#      with SIGKILL after a periodic checkpoint) and then resumed with
#      --resume,
#
# produce byte-identical state, tree, and operator logs, and (optionally)
#
#   C. a run without any of the new options is byte-identical to a run of an
#      unmodified build (e.g., upstream v1.1.1).
#
# It also checks the exit statuses, that a mismatched config is refused, and
# that a SIGINT the run was started ignoring stays ignored.
#
# usage: dev-tools/test-checkpoint-resume.sh PHYCOEVAL [REFERENCE-PHYCOEVAL]
#
# PHYCOEVAL is the phycoeval executable to test. REFERENCE-PHYCOEVAL is an
# executable built from the unmodified upstream code (e.g., from a separate
# build directory of the v1.1.1 tag); if it is not given, test C is skipped.
# Only synthetic test data from this repository is used.

set -u

usage () {
    echo ""
    echo "usage: $0 PHYCOEVAL [REFERENCE-PHYCOEVAL]"
    echo ""
}

if [ $# -lt 1 ] || [ $# -gt 2 ]
then
    usage
    exit 2
fi

abspath () {
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

phycoeval="$(abspath "$1")"
reference=""
if [ $# -eq 2 ]
then
    reference="$(abspath "$2")"
fi

this_dir="$(cd "$(dirname "$0")" && pwd)"
data_dir="${this_dir}/../test/data"
config_name="config-checkpoint-resume.yml"
data_name="species-6-genomes-2-chars-300-simulated.yml"
stem="config-checkpoint-resume"
seed=4242

work_dir="$(mktemp -d "${TMPDIR:-/tmp}/phycoeval-checkpoint-test.XXXXXX")"
echo "Working directory: ${work_dir}"

failures=0

pass () { echo "PASS: $1"; }
fail () { echo "FAIL: $1"; failures=$((failures + 1)); }

new_run_dir () {
    mkdir -p "${work_dir}/$1"
    cp "${data_dir}/${config_name}" "${data_dir}/${data_name}" "${work_dir}/$1/"
}

compare_logs () {
    # $1 = description, $2 = expected dir, $3 = actual dir
    local ok=1
    for suffix in state-run-1.log trees-run-1.nex operator-run-1.log
    do
        if cmp "${2}/${stem}-${suffix}" "${3}/${stem}-${suffix}"
        then
            echo "    identical: ${suffix}"
        else
            echo "    DIFFERENT: ${suffix}"
            ok=0
        fi
    done
    if [ $ok -eq 1 ]
    then
        pass "$1"
    else
        fail "$1"
    fi
}

# Wait (up to ~10 minutes) until a file contains a line matching a pattern.
# Gives up at once if the process PID (third argument) has exited, e.g. when
# phycoeval rejects its options.
wait_for_line () {
    local file="$1"
    local pattern="$2"
    local pid="$3"
    local tries=0
    while [ $tries -lt 6000 ]
    do
        if [ -f "$file" ] && grep -q -E "$pattern" "$file"
        then
            return 0
        fi
        if ! kill -0 "$pid" 2>/dev/null
        then
            return 1
        fi
        sleep 0.1
        tries=$((tries + 1))
    done
    return 1
}

# Wait for a process to exit; sets exit_status
wait_for_exit () {
    wait "$1"
    exit_status=$?
}

# Run a command with SIGINT handled the default way. A script without job
# control starts background commands with SIGINT ignored, and phycoeval
# leaves an ignored SIGINT ignored. Only for background commands ("... &"):
# the exec makes $! the PID of the command itself.
with_sigint () {
    exec perl -e '$SIG{INT} = "DEFAULT"; exec @ARGV or die "exec: $!\n"' "$@"
}

#############################################################################
# A. Uninterrupted run (no checkpointing)
#############################################################################
echo ""
echo "=== A: uninterrupted run ==="
new_run_dir A
( cd "${work_dir}/A" && "${phycoeval}" --seed ${seed} ${config_name} > stdout.txt 2>&1 )
if [ $? -eq 0 ] && tail -c 4 "${work_dir}/A/${stem}-trees-run-1.nex" | grep -q "END;"
then
    pass "A: uninterrupted run finished, trees file ends with END;"
else
    fail "A: uninterrupted run"
fi
if [ -e "${work_dir}/A/${stem}-checkpoint-run-1.txt" ]
then
    fail "A: no checkpoint file should be written without checkpointing options"
else
    pass "A: no checkpoint file without checkpointing options"
fi

#############################################################################
# B1. Stopped with SIGTERM, then resumed
#############################################################################
echo ""
echo "=== B1: SIGTERM, then resume ==="
new_run_dir B1
cd "${work_dir}/B1"
"${phycoeval}" --seed ${seed} --checkpoint-every=500 ${config_name} > stdout1.txt 2>&1 &
pid=$!
# The checkpoint at generation 500 is written right after its sample; wait
# until the run is a bit past that
if ! wait_for_line "${stem}-state-run-1.log" "^600	" $pid
then
    kill -9 $pid
    fail "B1: run did not reach generation 600"
fi
kill -TERM $pid
wait_for_exit $pid
if [ $exit_status -eq 75 ]
then
    pass "B1: stopped by SIGTERM with exit status 75"
else
    fail "B1: expected exit status 75 after SIGTERM, got ${exit_status}"
fi
if grep -q -- "--resume=${stem}-checkpoint-run-1.txt" stdout1.txt
then
    pass "B1: message says how to resume"
else
    fail "B1: message about how to resume not found"
fi
if tail -c 4 "${stem}-trees-run-1.nex" | grep -q "END;"
then
    fail "B1: run was not stopped part-way"
fi

# A resume with a different config must be refused, and must not touch the logs
cp "${stem}-state-run-1.log" state-before.log
sed 's/chain_length: 2000/chain_length: 2500/' ${config_name} > changed.yml
"${phycoeval}" --seed ${seed} --resume=${stem}-checkpoint-run-1.txt changed.yml > stdout-changed.txt 2>&1
changed_status=$?
if [ $changed_status -ne 0 ] && cmp -s state-before.log "${stem}-state-run-1.log"
then
    pass "B1: resuming with a different config is refused"
else
    fail "B1: resuming with a different config was not refused (status ${changed_status})"
fi
"${phycoeval}" --seed 1 --resume=${stem}-checkpoint-run-1.txt ${config_name} > stdout-seed.txt 2>&1
if [ $? -ne 0 ] && cmp -s state-before.log "${stem}-state-run-1.log"
then
    pass "B1: resuming with a different seed is refused"
else
    fail "B1: resuming with a different seed was not refused"
fi

# Resume, and stop again (second signal after resuming), then resume to the end
with_sigint "${phycoeval}" --resume=${stem}-checkpoint-run-1.txt ${config_name} > stdout2.txt 2>&1 &
pid=$!
if ! wait_for_line "${stem}-state-run-1.log" "^1[0-9]00	" $pid
then
    kill -9 $pid
    fail "B1: resumed run did not reach generation 1000"
fi
kill -INT $pid
wait_for_exit $pid
if [ $exit_status -eq 75 ]
then
    pass "B1: resumed run stopped by SIGINT with exit status 75"
else
    fail "B1: expected exit status 75 after SIGINT, got ${exit_status}"
fi
"${phycoeval}" --resume=${stem}-checkpoint-run-1.txt ${config_name} > stdout3.txt 2>&1
if [ $? -eq 0 ]
then
    pass "B1: run finished after being resumed twice"
else
    fail "B1: final resume failed"
fi
if [ -e "${stem}-checkpoint-run-1.txt" ]
then
    fail "B1: checkpoint file should be removed when the run finishes"
else
    pass "B1: checkpoint file removed when the run finished"
fi
cd "${work_dir}"
compare_logs "B1: SIGTERM/SIGINT + resume gives the same logs as an uninterrupted run" A B1

#############################################################################
# B2. Killed with SIGKILL after a periodic checkpoint, then resumed
#############################################################################
echo ""
echo "=== B2: SIGKILL, then resume ==="
new_run_dir B2
cd "${work_dir}/B2"
"${phycoeval}" --seed ${seed} --checkpoint-every=500 ${config_name} > stdout1.txt 2>&1 &
pid=$!
# Past the checkpoint of generation 1000 and part-way to the next samples, so
# there is output after the checkpoint that has to be cut off
if ! wait_for_line "${stem}-state-run-1.log" "^1100	" $pid
then
    kill -9 $pid
    fail "B2: run did not reach generation 1100"
fi
kill -KILL $pid
wait_for_exit $pid
echo "    killed, state log has $(wc -l < ${stem}-state-run-1.log) lines"
if [ -e "${stem}-checkpoint-run-1.txt" ]
then
    pass "B2: checkpoint exists after kill -9"
else
    fail "B2: no checkpoint after kill -9"
fi
"${phycoeval}" --resume=${stem}-checkpoint-run-1.txt ${config_name} > stdout2.txt 2>&1
if [ $? -eq 0 ]
then
    pass "B2: run finished after being resumed"
else
    fail "B2: resume failed"
fi
cd "${work_dir}"
compare_logs "B2: SIGKILL + resume gives the same logs as an uninterrupted run" A B2

#############################################################################
# B3. Checkpointing alone does not change the output, and a SIGINT that was
#     ignored when the run started stays ignored
#############################################################################
echo ""
echo "=== B3: checkpointing without interruption ==="
new_run_dir B3
cd "${work_dir}/B3"
# Started in the background without job control, so SIGINT is ignored
"${phycoeval}" --seed ${seed} --checkpoint-every=500 ${config_name} > stdout.txt 2>&1 &
pid=$!
if wait_for_line "${stem}-state-run-1.log" "^600	" $pid
then
    kill -INT $pid
fi
wait_for_exit $pid
if [ $exit_status -eq 0 ]
then
    pass "B3: run with checkpointing finished, ignoring SIGINT"
else
    fail "B3: run with checkpointing failed (exit status ${exit_status})"
fi
cd "${work_dir}"
compare_logs "B3: run with checkpointing is identical to a run without" A B3

#############################################################################
# C. Without the new options, identical to an unmodified build
#############################################################################
echo ""
echo "=== C: comparison with unmodified build ==="
if [ -z "${reference}" ]
then
    echo "SKIPPED: no reference executable given"
else
    new_run_dir C
    ( cd "${work_dir}/C" && "${reference}" --seed ${seed} ${config_name} > stdout.txt 2>&1 )
    if [ $? -eq 0 ]
    then
        pass "C: reference run finished"
    else
        fail "C: reference run failed"
    fi
    compare_logs "C: output without the new options is identical to the unmodified build" C A
fi

echo ""
if [ $failures -eq 0 ]
then
    echo "All checks passed."
    rm -rf "${work_dir}"
    exit 0
fi
echo "${failures} check(s) FAILED. Output kept in ${work_dir}"
exit 1
