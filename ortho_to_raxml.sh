#!/usr/bin/env bash

###################################################################################
# Description:
#   A pipeline to reconstruct a species tree from the Single-Copy Orthologues (SCOs)
#   identified by OrthoFinder. The pipeline performs the necessary steps, starting from
#   the Orthofinder results and produces the input needed for generax.
#   In sequense:
#     1. Extracts SCO sequences from OrthoFinder results
#     2. Replaces non-standard amino acid character U with X, for trimAl errors
#     3. Aligns each SCO with MAFFT
#     4. Converts multi-line FASTA alignments to one-line format
#     5. Concatenates all alignments into a supermatrix
#     6. Trims the supermatrix with trimAl
#     7. Selects the best substitution model with ModelTest-NG
#
# Dependencies:
#   - MAFFT    (https://mafft.cbrc.jp/alignment/software/)
#   - trimAl   (http://trimal.cgenomics.org/)
#   - ModelTest-NG (https://github.com/ddarriba/modeltest)
#   - RAxML-NG (https://github.com/amkozlov/raxml-ng)  [downstream use]
#
# Usage:
#   bash species_tree_pipeline.sh
# Arguments:
#   -p  Path to the OrthoFinder Results directory (like: /..../OrthoFinder/Results)
#   -n  Number of CPU cores to use for parallel steps
#
# Notes:
#   - Edit the Parameters section below before running.
#   - All output is written under ./sp_tree/.
#   - The script assumes OrthoFinder was run beforehand and its results are
#     available at the paths specified in the Parameters section.
###################################################################################

set -euo pipefail

####################################################################################
# Arguments

usage() {
    echo "Usage: $0 -p <orthofinder_results_path> -n <ncores>"
    echo ""
    echo "  -p  Path to the OrthoFinder Results directory"
    echo "      (sth. like /../../OrthoFinder/Results)"
    echo "  -n  Number of CPU cores to use (e.g. 40)"
    exit 1
}

ncores=""
orthofinder_results=""

while getopts ":p:n:" opt; do
    case ${opt} in
        p) orthofinder_results="${OPTARG}" ;;
        n) ncores="${OPTARG}" ;;
        :) echo "Error: Option -${OPTARG} requires an argument." >&2; usage ;;
        \?) echo "Error: Unknown option -${OPTARG}." >&2; usage ;;
    esac
done

if [[ -z "${orthofinder_results}" || -z "${ncores}" ]]; then
    echo "Error: Both -p and -n are required." >&2
    usage
fi

if [[ ! -d "${orthofinder_results}" ]]; then
    echo "Error: OrthoFinder Results directory not found: ${orthofinder_results}" >&2
    exit 1
fi

if ! [[ "${ncores}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: -n must be a positive integer, got: ${ncores}" >&2
    exit 1
fi

####################################################################################
# Parameters

# Path to the file of all HOGs that are Single-Copy Orthologues (SCOs)
SCOs_HOGs_file="${orthofinder_results}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"

# Directory containing the per-HOG FASTA sequence files from OrthoFinder
SCOs_seqs_path="${orthofinder_results}/Single_Copy_Orthologue_Sequences"

if [[ ! -f "${SCOs_HOGs_file}" ]]; then
    echo "Error: SCOs HOGs file not found: ${SCOs_HOGs_file}" >&2
    exit 1
fi

if [[ ! -d "${SCOs_seqs_path}" ]]; then
    echo "Error: SCO sequences directory not found: ${SCOs_seqs_path}" >&2
    exit 1
fi

# echo "OrthoFinder Results : ${orthofinder_results}"
# echo "CPU cores           : ${ncores}"
# echo "SCOs HOGs file      : ${SCOs_HOGs_file}"
# echo "SCOs sequences dir  : ${SCOs_seqs_path}"

####################################################################################
# 1. Copy SCO sequence files into local working directory

echo "1. Copying SCO sequence files"

mkdir -p ./sp_tree/SCseqs/

while read -r hogs; do
    echo "HOG: ${hogs}"
    filepath=$(find "${SCOs_seqs_path}" -type f -name "${hogs}.fa")
    if [[ -z "${filepath}" ]]; then
        echo "WARNING: FASTA file not found for ${hogs}, skipping." >&2
        continue
    fi
    echo "Found: ${filepath}"
    cp "${filepath}" ./sp_tree/SCseqs/
done < "${SCOs_HOGs_file}"

echo "Copying complete"

####################################################################################
# 2. Replace non-standard residue U with X, avoids downstream errors, esp. with trimAl

echo "2. Cleaning sequences"

mkdir -p ./sp_tree/seqs_moded/

for f in ./sp_tree/SCseqs/*; do
    base=$(basename "${f}")
    echo "Processing: ${base}"
    sed '/^>/! s/U/X/g' "${f}" > "./sp_tree/seqs_moded/${base}"
done

echo "Cleaning complete"

####################################################################################
# 3. Multiple sequence alignment of each SCO fasta with MAFFT (--auto)

echo "3. Running MAFFT alignments"

mod_dir="./sp_tree/seqs_moded"
msas_dir="./sp_tree/msas"
mkdir -p "${msas_dir}"

i=0
ending=$(ls "${mod_dir}" | wc -l)

for fasta in "${mod_dir}"/*.fa; do
    ((i++))
    echo "  Starting ${i} / ${ending}: $(basename "${fasta}")"
    name=$(basename "${fasta}")
    mafft --auto --thread "${ncores}" "${fasta}" > "${msas_dir}/${name}"
    echo "  Done"
done

echo "Alignments complete"

#####################################################################################
# 4. Convert multi-line FASTA to one-line FASTA

echo "=== Step 4: Converting alignments to one-line FASTA ==="

oneline_msas_dir="./sp_tree/oneline_msas"
mkdir -p "${oneline_msas_dir}"

i=0
ending=$(ls "${msas_dir}" | wc -l)

for align in "${msas_dir}"/*.fa; do
    ((i++))
    name=$(basename "${align}" .fa)
    echo "  NO.${name} — Starting ${i} / ${ending}"
    perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' "${align}" \
        > "${oneline_msas_dir}/NO.${name}.fasta"
done

echo "=== Step 4 complete ==="

#####################################################################################
# 5. Concatenate alignments into a supermatrix

echo "5. Concatenating alignments"
 
ulimit -n 2400
concat_msa="./sp_tree/concat_msas.fasta"
tmp_dir=$(mktemp -d)
 
# Save headers from the first alignment
first_file=$(ls "${oneline_msas_dir}"/*.fasta | head -1)
grep "^>" "${first_file}" > "${tmp_dir}/headers.txt"
 
# Extract only sequence lines from each alignment
for f in "${oneline_msas_dir}"/*.fasta; do
    name=$(basename "${f}")
    grep -v "^>" "${f}" > "${tmp_dir}/${name}.seqs"
done
 
# Paste sequence lines side-by-side across all alignments
paste -d "" "${tmp_dir}"/*.seqs > "${tmp_dir}/concat_seqs.txt"
 
# Interleave headers and concatenated sequences into final FASTA
paste -d '\n' "${tmp_dir}/headers.txt" "${tmp_dir}/concat_seqs.txt" > "${concat_msa}"
 
rm -rf "${tmp_dir}"
 
echo "Supermatrix written to: ${concat_msa}"

#####################################################################################
# 6. Trim the supermatrix with trimAl

echo "6. Trimming supermatrix with trimAl"

trimmed_msa="./sp_tree/concat_trimmed.fasta"
trimal -in "${concat_msa}" -out "${trimmed_msa}" -gt 0.9 -cons 65

echo "Trimmed supermatrix written to: ${trimmed_msa}"

#####################################################################################
#7. Model selection with ModelTest-NG
# Identify the best-fit amino acid substitution model for use with RAxML-NG.

echo "7. Running ModelTest-NG"

model_path="./sp_tree/models"
mkdir -p "${model_path}"

modeltest-ng \
    -i "${trimmed_msa}" \
    -d "aa" \
    -p "${ncores}" \
    -o "${model_path}/concat_trimmed_mod.fasta" \
    -T "raxml"

echo "  Model selection results written to: ${model_path}/"

echo "Pipeline finished. Use the ModelTest-NG output to configure RAxML-NG for the final tree inference step."