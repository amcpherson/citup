# Citup 

The following package implements the method described in
[Clonality inference in multiple tumor samples using phylogeny](http://dx.doi.org/10.1093/bioinformatics/btv003)

Citup estimates the clone phylogeny and clonal genotypes for deep sequencing of SNVs in multiple tumour biopsies.
The input are cellular frequencies of mutations as estimated from the deep sequencing data.  The method infers an
evolutionary tree, and an assignment of mutations to nodes in the tree.

## Prerequisites

### CPLEX

CPLEX is required to run citup.  A license is required and is free for academic use.  Install CPLEX and set the `CPLEX_DIRECTORY`
and `CPLEX_BUILD` environment variables.

For example:

```
CPLEX_DIRECTORY=/Applications/IBM/ILOG/CPLEX_Studio125
CPLEX_BUILD=x86-64_darwin
```

Specifies that you have installed to `/Applications/IBM/ILOG/` and are using the 64 bit mac binaries.  The available builds will
be in subdirectories of `/Users/amcphers/Applications/IBM/ILOG/CPLEX_Studio125/cplex/lib/`

### Boost C++

A boost c++ header only installation is required.  If boost is not installed on your system, [download](http://www.boost.org/users/download/)
and unpack, and set the following environment variable to the resulting directory.

```
BOOST_DIRECTORY=
```

## Build

Run the following commands to build citup.

```
cd src
make citupiter citupqip
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
python citup_iter.py freq.txt output_solution.txt output_alltrees.tsv
```

The above command will run citup iter and produce a solution file in `output_solution.txt` and a table
of results for all trees in `output_alltrees.tsv` in tab separated format.

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
python citup.py freq.txt clusters.txt output_solution.txt output_alltrees.tsv
```

The above command will run citup iter and produce a solution file in `output_solution.txt` and a table
of results for all trees in `output_alltrees.tsv` in tab separated format.

For additional options run `python citup.py -h`.




