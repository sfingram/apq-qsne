## apq-qsne

A collection of C/C++ code implementing the algorithms from the paper "[Dimensionality Reduction for Documents with Nearest Neighbor Queries](http://www.cs.ubc.ca/labs/imager/tr/2014/QSNE/)" by Ingram and Munzner, 2014.

This includes the All-Pairs Query, or *APQ*, nearest-neighbor algorithm and the *Q-SNE*, dimensionality reduction algorithm.

## Building

The code requires a C/C++ compiler.  I have successfully used both [clang](http://clang.llvm.org/) and [gcc](https://gcc.gnu.org/).

### Prerequisites

You will first need to have [GLib](https://developer.gnome.org/glib/) installed on your system.

For Mac (with [homebrew](http://brew.sh/) installed):

`brew install glib`

For Ubuntu Linux:

`sudo apt-get install libglib2.0-dev`

For other Linux [build from source](http://www.linuxfromscratch.org/blfs/view/svn/general/glib2.html).

## Compiling

For building the APQ algorithm just type

`make`

which will make the executable file `testapq.` For building BH-SNE with arbitrary input (use input from apq for Q-SNE) type

`make tsne`

As a fun bonus, there are some other executables you can build if you dig around in the makefile.

## Input File Format

APQ processes document _term-vectors_ into _nearest-neighbor_ files (both in `vec` format).

QSNE processes nearest-neighbor files into _coordinate_ files (in `csv` file).

### Term-Vector File Description

Term vectors represent documents as a vector of term counts in a term vectors space.  Each dimension represents a different term and the value represents the [TFIDF](http://en.wikipedia.org/wiki/Tf%E2%80%93idf) weight of that term.  Because most terms don't appear in a given document, APQ expects _sparse_ vectors as input.

### Nearest-Neighbor File Description

APQ computes the approximate k most similar documents for each document in the input as well as the approximate [cosine distance](http://en.wikipedia.org/wiki/Cosine_distance) between the specified document pairs.

### Vec File Format

Both term term vectors and nearest neighbors are stored in `vec` format.  Each line describes a document.

```
(index1,value1) (index2,value2)
```
For term-vectors, indices represent term numbers and values.  For nearest neighbors, indices represent document numbers.

## Running

(soon)

### Command Line Options

for help running the APQ algorithm

`testapq -?`

for help running BH-SNE (and Q-SNE)

`tsne -?`

### Some Examples

(soon)

## About the Algorithms

(soon)

### APQ

### Q-SNE

