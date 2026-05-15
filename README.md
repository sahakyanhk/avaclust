# RNPclust

Structural clustering of ribonucleoproteins and RNA-protein interfaces using USalign.


## Dependencies 

- [USalign](https://github.com/pylelab/USalign)
- [GNU parallel](https://www.gnu.org/software/parallel/)
- Python 3.6+ with NumPy, SciPy, BioPython

## Installation

```bash
wget https://github.com/sahakyanhk/rnpclust/archive/refs/heads/main.zip
unzip main.zip && mv rnpclust-main rnpclust
chmod +x rnpclust/bin/*

export PATH="$PWD/rnpclust/bin:$PATH" 


# Install Python dependencies and USalign
pip install numpy scipy biopython
conda install -c bioconda usalign
```

## Usage

```bash
rnpclust -i <pdb_dir> [-o <output_dir>] [-c <cutoff>] [-m <method>]
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i` | Input directory with `.pdb` files | *required* |
| `-o` | Output directory | `rnpclust_out` |
| `-c` | TM-score cutoff | `auto` |
| `-m` | Clustering method: `0` = hierarchical, `1` = set-cover | `0` |

### Examples

```bash
# Hierarchical clustering with TM-score cutoff
rnpclust -i examples/pdb_interface/ -o examples/clustering_results -c 0.4  -m 0

# extract interfaces intercating with chain A within 15A, then cluster with cutoff 0.4
python bin/extract_interface.py -i examples/pdb -c A --cut 15 -o examples/pdb_interface 
rnpclust -i examples/pdb_interface/ -o examples/clustering_results -c 0.4  -m 0

```

### Clustering methods

**Hierarchical (`-m 0`, default)** — Average linkage clustering on a TM-score distance matrix (1 - TM)\
**Set-cover (`-m 1`)** — Greedy clustering. (Faster, but less accurate)

### Output

```
results/
  all_vs_all.tsv           # pairwise USalign scores (one row per PDB pair, best chain combo)
  clusters.dat             # cluster assignments (wide format, centroid first)
  aligned_clusters/        # superimposed PDBs per cluster
    clust0001/
    clust0002/
    ...
```

