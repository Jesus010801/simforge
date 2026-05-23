#!/bin/bash
set -e
bash "$(dirname "$0")/run_nvt.sh"
bash "$(dirname "$0")/run_npt.sh"
