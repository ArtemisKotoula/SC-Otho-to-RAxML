# SC-Otho-to-RAxML
A Bash pipeline for reconstructing a species tree from **Single-Copy Orthologues (SCOs)** identified by [OrthoFinder](https://github.com/davidemms/OrthoFinder). 

## Description
Starting from OrthoFinder output, the pipeline aligns, concatenates, trims, prepares a supermatrix, and selects the best model for phylogenetic inference with [RAxML-NG](https://github.com/amkozlov/raxml-ng).

## Usage

   ```bash
   bash ortho_to_raxml.sh -p <orthofinder_results_path> -n <ncores>
   ```
 
   | Argument | Description |
   |----------|-------------|
   | `-p` | Path to the OrthoFinder `Results/` directory |
   | `-n` | Number of CPU cores for MAFFT and ModelTest-NG |

### Input

The script uses as input the `Orthogroups_SingleCopyOrthologues.txt` file and `Single_Copy_Orthologue_Sequences/` directory, created by OrthoFinder. The path to the OrthoFinder `Results/` directory is needed to locate the required inputs.

### Tool parameters

- MAFFT is called with `--auto`.
- trimAl is run with: `-gt 0.9` and `-cons 65`.

Tool paraneters can be editied in their respecitve steps within the script.

## Output

All output is written under `./sp_tree/`:

```
sp_tree/
├── SCseqs/                   # Raw SCO FASTA files copied from OrthoFinder
├── seqs_moded/               # U → X substituted FASTA files
├── msas/                     # Per-gene MAFFT alignments
├── oneline_msas/             # One-line FASTA versions of alignments
├── concat_msas.fasta         # Concatenated supermatrix
├── concat_trimmed.fasta      # trimAl-trimmed supermatrix
└── models/                   # ModelTest-NG output (model selection results)
```

The trimmed supermatrix (`concat_trimmed.fasta`) and ModelTest-NG results (`models/`) are the direct inputs needed for a downstream RAxML-NG tree inference analysis.


## Dependencies

| Tool | Version | Purpose |
|------|---------------|---------|
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | ≥ 7.0 | Multiple sequence alignment |
| [trimAl](http://trimal.cgenomics.org/) | ≥ 1.4 | Alignment trimming |
| [ModelTest-NG](https://github.com/ddarriba/modeltest) | ≥ 0.1 | Substitution model selection |
| Perl | ≥ 5 | One-line FASTA conversion |

All tools must be available in your `$PATH`.



