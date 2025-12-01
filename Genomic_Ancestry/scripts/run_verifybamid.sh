#!/usr/bin/env bash
# run_verifybamid.sh
#
# Wrapper script for running VerifyBamID2 on a single CRAM file.
# Invoked automatically by Snakemake.
#
# Usage:
#   run_verifybamid.sh <cram> <reference_fa> <sample_id>
#
# Produces:
#   output/verifybamid/<sample>.selfSM
#   output/verifybamid/<sample>.Ancestry

set -euo pipefail

CRAM="$1"
REF="$2"
SAMPLE="$3"

mkdir -p output/verifybamid

# Path prefix for the 1000 Genomes GRCh38 SVD dataset used by VerifyBamID2.
# The pipeline provides at least:
#   data/1000G/1000g.phase3.100k.b38.vcf.gz.dat.V
SVD_PREFIX="data/1000G/1000g.phase3.100k.b38.vcf.gz.dat"

verifybamid2 \
    --BamFile "${CRAM}" \
    --Reference "${REF}" \
    --SVDPrefix "${SVD_PREFIX}" \
    --NumPC 4 \
    --Output "output/verifybamid/${SAMPLE}"

echo "VerifyBamID2 completed for ${SAMPLE}"
