<img width="2578" height="1334" alt="image" src="https://github.com/user-attachments/assets/39deeb99-88c6-4c93-a45b-95d41f8f72ca" />

# avaclust

All-vs-all comparison and clustering of molecular structures, complexes and interfaces.


## Installation

```bash
pip install git+https://github.com/sahakyanhk/avaclust.git

conda install -c bioconda usalign     
```

This installs the `avaclust` command (plus `avaclust-cluster`, `avaclust-plot`,
and `avaclust-interface` for running individual steps).

## Usage

```bash
avaclust -i <pdb_dir> [-o <output_dir>] [--cutoff <TM>] [options]
```



### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Directory of `.pdb`/`.ent` files, or a tar archive (`.tar.gz`, `.tgz`, `.tar`, `.tar.bz2`, `.tar.xz`) | *required* |
| `-o, --output` | Output directory | `avaclust_out` |
| `--cutoff` | TM-score cutoff (cut tree at `1 - cutoff`); omit for auto-detection | `auto` |
| `--k` | Force a specific number of clusters | — |
| `--linkage` | `average` (UPGMA), `complete`, or `single` | `average` |
| `--chains` | Chains to align, e.g. `A,B` or `A`; if omitted, USalign uses the **first chain** of each structure | first chain |
| `--threads` | Threads for USalign | all cores |
| `--interface-cutoff` | Extract inter-chain interface residues (Å) before clustering | off |
| `--min-cluster-size` | Min members for a cluster to be structurally aligned | `3` |
| `--labels` | Label the dendrogram leaves with structure names | off |
| `--no-align` / `--no-plot` | Skip the alignment / figure step | — |

### Examples

```bash
# Cut the tree at a fixed TM-score cutoff and label the leaves
avaclust -i examples/c2_folds.tgz -o examples/clustering_results/ --cutoff 0.75 

# Cluster interfaces between chains A and B (residues within 15 Å) 
avaclust -i complexes/ -o results --chains A,B --cutoff 0.7 --interface-cutoff 15
```


