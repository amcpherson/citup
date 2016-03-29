# Citup 

The following package implements the method described in
[Clonality inference in multiple tumor samples using phylogeny](http://dx.doi.org/10.1093/bioinformatics/btv003)

Citup estimates the clone phylogeny and clonal genotypes for deep sequencing of SNVs in multiple tumour biopsies.
The input are cellular frequencies of mutations as estimated from the deep sequencing data.  The method infers an
evolutionary tree, and an assignment of mutations to nodes in the tree.

## Installation

The supported installation method for citup is using [conda](http://conda.pydata.org/docs/).  If you do not already
use conda, we recommend installing [anaconda](https://www.continuum.io/downloads) or
[miniconda](http://conda.pydata.org/miniconda.html).

You will need to add my conda channel using:

```
conda config --add channels http://conda.anaconda.org/dranew
```

Cplex is required to run citup.  For licensing reasons this could not be included in my public conda repo.
A license is required and is free for academic use.  The recommended installation procedure involves building
and installing cplex using conda.  A [conda build recipe for cplex](https://bitbucket.org/dranew/conda_recipes/src)
is provided.

With cplex installed you should be able to type:

```
conda install citup
```

## Usage

### Citup Iter

To run the iterative version of citup, provide a table of frequencies.  The input format is tab/whitespace separated, with each row
representing the frequency of a mutation and each column is a tumour sample.  No header is required.

For example, the following would be input for 3 mutations in 2 samples.

```
0.2 0.1
0.4 0.3
0.5 0.1
```

Given mutation frequences `freq.txt`, run citup iter using the following command.

```
python run_citup_iter.py freq.txt results.h5
```

The above command will run citup iter and produce results in pandas hdf5 format in `results.h5`.

For additional options run `python citup_iter.py -h`.

### Citup QIP

To run the QIP version of citup, provide a table of mutation frequencies, and a table of mutation clusters.
The input format for mutation frequencies is described above.  The mutation clusters is a single line per mutation,
containing a 0 based integer cluster index for that mutation.  For example, the following specifies 3 mutations, the
first two in the same cluster the last in a different cluster.

```
0
0
1
```

Given mutation frequences `freq.txt` and mutation clusters `clusters.txt` run citup QIP using the following command.

```
python run_citup_qip.py freq.txt clusters.txt results.h5
```

The above command will run citup iter and produce results in pandas hdf5 format in `results.h5`.

For additional options run `python citup.py -h`.

## Output Format

The output of citup is [pandas hdf5](http://pandas.pydata.org/pandas-docs/version/0.18.0/generated/pandas.read_hdf.html) format.

The results store contains the following pandas series, with an entry per tree solution:

 - `/results/bic`
 - `/results/error_rate`
 - `/results/likelihood`
 - `/results/num_mutations`
 - `/results/num_nodes`
 - `/results/num_samples`
 - `/results/objective_value`
 - `/results/optimal`
 - `/results/tree_id`
 - `/results/tree_index`
 - `/results/tree_string`

The store also contains the following pandas series, with one series per tree solution:

 - `/trees/{tree_solution}/cluster_assignment`
 - `/trees/{tree_solution}/objective_value`

Finally, the store contains the following pandas dataframes, with one frame per tree solution:

 - `/trees/{tree_solution}/adjacency_list`
 - `/trees/{tree_solution}/clade_freq`
 - `/trees/{tree_solution}/clone_freq`
 - `/trees/{tree_solution}/gamma_matrix`

