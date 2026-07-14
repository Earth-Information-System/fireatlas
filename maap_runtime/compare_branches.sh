#!/bin/bash
# End-to-end cross-branch comparison test.
# Runs the FEDS pipeline on main, staging, and a dev branch (default: current branch)
# under both NRT and archive modes, then compares perimeter outputs.
#
# Usage:
#   bash maap_runtime/compare_branches.sh [DEV_BRANCH]
#
# If DEV_BRANCH is not provided, defaults to the currently checked-out branch.
# 
# Requires optional dev dependency group to be installed in working python environment 
# (cd fireatlas; pip install -e .) 

set -eo pipefail
export TZ="Etc/UTC"

SCRIPT_DIR="$(cd "$(dirname "$0")"; pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.."; pwd -P)"
TIMESTAMP=$(date -u +"%Y%m%d%H%M%S")

# ── Branch setup ──────────────────────────────────────────────────────────────
ORIGINAL_BRANCH=$(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD)
DEV_BRANCH="${1:-$ORIGINAL_BRANCH}"

echo "=== Cross-Branch E2E Comparison Test ==="
echo "Timestamp:   $TIMESTAMP"
echo "Prod branch: main"
echo "Staging:     staging"
echo "Dev branch:  $DEV_BRANCH"
echo ""

# Stash any local changes and restore on exit
STASH_CREATED=false
if ! git -C "$REPO_ROOT" diff --quiet || ! git -C "$REPO_ROOT" diff --cached --quiet; then
    echo "Stashing local changes..."
    git -C "$REPO_ROOT" stash push -m "compare_branches.sh auto-stash $TIMESTAMP"
    STASH_CREATED=true
fi

restore_branch() {
    echo ""
    echo "Restoring branch: $ORIGINAL_BRANCH"
    git -C "$REPO_ROOT" checkout "$ORIGINAL_BRANCH" 2>/dev/null || true
    if [ "$STASH_CREATED" = true ]; then
        echo "Popping stash..."
        git -C "$REPO_ROOT" stash pop || true
    fi
}
trap restore_branch EXIT

# ── Date derivation ───────────────────────────────────────────────────────────
# NRT: yesterday-10d → yesterday
YESTERDAY=$(date -u -d "yesterday" +"%Y-%m-%d" 2>/dev/null || date -u -v-1d +"%Y-%m-%d")
NRT_TED_DATE="$YESTERDAY"
NRT_TST_DATE=$(date -u -d "$YESTERDAY - 10 days" +"%Y-%m-%d" 2>/dev/null || date -u -j -f "%Y-%m-%d" -v-10d "$YESTERDAY" +"%Y-%m-%d")

NRT_TST_YEAR=$(echo $NRT_TST_DATE | cut -d'-' -f1)
NRT_TST_MONTH=$(echo $NRT_TST_DATE | cut -d'-' -f2 | sed 's/^0//')
NRT_TST_DAY=$(echo $NRT_TST_DATE | cut -d'-' -f3 | sed 's/^0//')
NRT_TED_YEAR=$(echo $NRT_TED_DATE | cut -d'-' -f1)
NRT_TED_MONTH=$(echo $NRT_TED_DATE | cut -d'-' -f2 | sed 's/^0//')
NRT_TED_DAY=$(echo $NRT_TED_DATE | cut -d'-' -f3 | sed 's/^0//')

# Archive: same window anchored in 2024
ARCH_TST_YEAR=2024
ARCH_TED_YEAR=2024
ARCH_TST_MONTH=$NRT_TST_MONTH
ARCH_TST_DAY=$NRT_TST_DAY
ARCH_TED_MONTH=$NRT_TED_MONTH
ARCH_TED_DAY=$NRT_TED_DAY

NRT_TST="[$NRT_TST_YEAR,$NRT_TST_MONTH,$NRT_TST_DAY,\"AM\"]"
NRT_TED="[$NRT_TED_YEAR,$NRT_TED_MONTH,$NRT_TED_DAY,\"AM\"]"
ARCH_TST="[$ARCH_TST_YEAR,$ARCH_TST_MONTH,$ARCH_TST_DAY,\"AM\"]"
ARCH_TED="[$ARCH_TED_YEAR,$ARCH_TED_MONTH,$ARCH_TED_DAY,\"AM\"]"

# Date strings for comparison script (YYYYMMDDAMPM)
NRT_DATE_STRING=$(printf "%04d%02d%02d" $NRT_TED_YEAR $NRT_TED_MONTH $NRT_TED_DAY)"AM"
ARCH_DATE_STRING=$(printf "%04d%02d%02d" $ARCH_TED_YEAR $ARCH_TED_MONTH $ARCH_TED_DAY)"AM"

echo "NRT window:     $NRT_TST → $NRT_TED  (date string: $NRT_DATE_STRING)"
echo "Archive window: $ARCH_TST → $ARCH_TED  (date string: $ARCH_DATE_STRING)"
echo ""

# ── Read prod FIRE_SOURCE default ─────────────────────────────────────────────
# Extract the default FIRE_SOURCE value from main's FireConsts.py so the
# NRT runs automatically track whatever production defaults to.
NRT_FIRE_SOURCE=$(git -C "$REPO_ROOT" show main:fireatlas/FireConsts.py \
    | python3 -c $'
import sys, re
content = sys.stdin.read()
m = re.search(r\'FIRE_SOURCE\\s*:.*?=\\s*Field\\(\\s*["\\\']+(SNPP|NOAA20|NOAA21|VIIRS|BAMOD)\', content)
print(m.group(1) if m else "NOAA20")
')
echo "NRT FIRE_SOURCE (from main): $NRT_FIRE_SOURCE"
echo ""

# ── NRT pre-flight: ensure input data exists, download if missing ──────────────
echo "=== Pre-flight: Checking NRT input data ==="
python3 -c "
from datetime import date, timedelta
import s3fs, sys

s3 = s3fs.S3FileSystem()
sat = '$NRT_FIRE_SOURCE'
tst_date = date($NRT_TST_YEAR, $NRT_TST_MONTH, $NRT_TST_DAY)
ted_date = date($NRT_TED_YEAR, $NRT_TED_MONTH, $NRT_TED_DAY)

from fireatlas import settings
base = settings.dirextdata.replace('s3://', '') if settings.dirextdata.startswith('s3://') else None

missing = []
d = tst_date
while d <= ted_date:
    fname = f'FIRMS_VIIRS_{sat}_NRT_{d.strftime(\"%Y%m%d\")}.csv'
    if base:
        path = f'{base}/VIIRS/FIRMS_VIIRS_{sat}_NRT/{fname}'
        if not s3.exists(path):
            missing.append(d)
    d += timedelta(days=1)

if missing:
    print(f'Missing NRT files for {len(missing)} day(s), attempting download...')
    import sys
    sys.path.insert(0, '.')
    from fireatlas.DataCheckUpdate import update_FIRMS
    failed = []
    for d in missing:
        result = update_FIRMS(d, sat, 'NRT')
        if result is None:
            failed.append(d)
    if failed:
        print(f'ERROR: Could not download NRT data for: {failed}')
        sys.exit(1)
    print('All missing files downloaded successfully.')
else:
    print('All NRT input files present.')
"
echo ""

# ── Run function ──────────────────────────────────────────────────────────────
run_branch() {
    local branch="$1"
    local suffix="$2"
    local regnm="$3"
    local tst="$4"
    local ted="$5"
    local fire_nrt="$6"
    local fire_source="$7"

    echo "--- Checking out $branch ---"
    git -C "$REPO_ROOT" checkout "$branch"

    echo "Running: regnm=$regnm  tst=$tst  ted=$ted  FIRE_NRT=$fire_nrt  FIRE_SOURCE=$fire_source"
    FEDS_FIRE_NRT="$fire_nrt" \
    FEDS_FIRE_SOURCE="$fire_source" \
    python3 "$REPO_ROOT/fireatlas/FireRunDaskCoordinator.py" \
        --regnm="$regnm" \
        --bbox="[-126,24,-61,49]" \
        --tst="$tst" \
        --ted="$ted" \
        --no-veda-copy
    echo "Completed: $regnm"
}

# ── NRT runs ──────────────────────────────────────────────────────────────────
echo "=== NRT Runs ==="
NRT_PROD_REGNM="CONUS_TEST_${TIMESTAMP}_PROD_NRT"
NRT_STAGING_REGNM="CONUS_TEST_${TIMESTAMP}_STAGING_NRT"
NRT_DEV_REGNM="CONUS_TEST_${TIMESTAMP}_DEV_NRT"

run_branch "main" "PROD"    "$NRT_PROD_REGNM"    "$NRT_TST" "$NRT_TED" "true" "$NRT_FIRE_SOURCE"
run_branch "staging"   "STAGING" "$NRT_STAGING_REGNM" "$NRT_TST" "$NRT_TED" "true" "$NRT_FIRE_SOURCE"
run_branch "$DEV_BRANCH" "DEV"   "$NRT_DEV_REGNM"     "$NRT_TST" "$NRT_TED" "true" "$NRT_FIRE_SOURCE"

# ── Archive runs ──────────────────────────────────────────────────────────────
echo ""
echo "=== Archive Runs ==="
ARCH_PROD_REGNM="CONUS_TEST_${TIMESTAMP}_PROD_ARCHIVE"
ARCH_STAGING_REGNM="CONUS_TEST_${TIMESTAMP}_STAGING_ARCHIVE"
ARCH_DEV_REGNM="CONUS_TEST_${TIMESTAMP}_DEV_ARCHIVE"

run_branch "main" "PROD"    "$ARCH_PROD_REGNM"    "$ARCH_TST" "$ARCH_TED" "false" "SNPP"
run_branch "staging"   "STAGING" "$ARCH_STAGING_REGNM" "$ARCH_TST" "$ARCH_TED" "false" "SNPP"
run_branch "$DEV_BRANCH" "DEV"   "$ARCH_DEV_REGNM"     "$ARCH_TST" "$ARCH_TED" "false" "SNPP"

# Restore original branch before verification steps
git -C "$REPO_ROOT" checkout "$ORIGINAL_BRANCH"

# ── S3 output verification ────────────────────────────────────────────────────
echo ""
echo "=== S3 Output Verification ==="
python3 -c "
import s3fs
s3 = s3fs.S3FileSystem()

prod_base    = 'maap-ops-workspace/shared/gsfc_landslides/FEDSoutput-v3'
staging_base = 'maap-ops-workspace/shared/zbecker/FEDSstaging/FEDSoutput-v3'
dev_base     = 'maap-ops-workspace/shared/gsfc_landslides/FEDS-staging/FEDSoutput-v3'

runs = [
    ('$NRT_PROD_REGNM',    prod_base),
    ('$NRT_STAGING_REGNM', staging_base),
    ('$NRT_DEV_REGNM',     dev_base),
    ('$ARCH_PROD_REGNM',   prod_base),
    ('$ARCH_STAGING_REGNM',staging_base),
    ('$ARCH_DEV_REGNM',    dev_base),
]

all_ok = True
for regnm, base in runs:
    path = f'{base}/{regnm}'
    files = s3.glob(f'{path}/**')
    status = 'OK' if files else 'EMPTY/MISSING'
    if not files:
        all_ok = False
    print(f'  s3://{path}: {status}')

if not all_ok:
    print('WARNING: Some output folders are empty or missing.')
else:
    print('All output folders verified.')
"

# ── Comparison ────────────────────────────────────────────────────────────────
echo ""
echo "=== Running Output Comparison ==="
python3 "$SCRIPT_DIR/compare_outputs.py" \
    --prod-nrt-regnm="$NRT_PROD_REGNM" \
    --staging-nrt-regnm="$NRT_STAGING_REGNM" \
    --dev-nrt-regnm="$NRT_DEV_REGNM" \
    --prod-archive-regnm="$ARCH_PROD_REGNM" \
    --staging-archive-regnm="$ARCH_STAGING_REGNM" \
    --dev-archive-regnm="$ARCH_DEV_REGNM" \
    --nrt-date-string="$NRT_DATE_STRING" \
    --archive-date-string="$ARCH_DATE_STRING" \
    --timestamp="$TIMESTAMP"

echo ""
echo "=== Done ==="
