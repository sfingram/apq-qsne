apq-qsne
========

A collection of C/C++ code implementing the algorithms from the paper "Dimensionality Reduction for Documents with Nearest Neighbor Queries" by Ingram and Munzner, 2014.

requires:

[GLib](https://developer.gnome.org/glib/)

build:

for building the APQ algorithm 

`make`

for building BH-SNE with arbitrary input (use input from apq for Q-SNE)

`make tsne`

run:

for help running the APQ algorithm

`testapq -?`

for help running BH-SNE (and Q-SNE)

`tsne -?`