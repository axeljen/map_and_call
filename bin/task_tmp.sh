#!/usr/bin/env bash

# Source from a task shell; never change its working directory or stage outputs here.
set -euo pipefail

if [[ ${1:?Missing scratch enabled flag} == true ]]; then
    task_tmp_root=${2:?--use_scratch requires --scratch_dir or a profile defining it}
else
    task_tmp_root=$PWD
fi

if [[ ! -d "$task_tmp_root" || ! -w "$task_tmp_root" ]]; then
    printf 'Error: temporary directory root is missing or not writable: %s\n' "$task_tmp_root" >&2
    exit 1
fi

TASK_TMP_DIR=$(mktemp -d "${task_tmp_root%/}/map-and-call.XXXXXXXXXX") || exit 1
TASK_TMP_DIR=$(cd "$TASK_TMP_DIR" && pwd -P) || exit 1
readonly TASK_TMP_DIR
export TMPDIR="$TASK_TMP_DIR"
printf 'Temporary workspace: %s\n' "$TMPDIR" >&2

task_tmp_cleanup() {
    local status=$?
    trap - EXIT
    if ! rm -rf -- "$TASK_TMP_DIR"; then
        printf 'Error: could not remove task temporary directory: %s\n' "$TASK_TMP_DIR" >&2
        [[ $status -ne 0 ]] || status=1
    fi
    exit "$status"
}

trap task_tmp_cleanup EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
