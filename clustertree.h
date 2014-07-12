#include <glib.h>
#include "vec.h"
#include "disjointset.h"

#ifndef CLUSTERTREE_H
#define CLUSTERTREE_H

typedef struct _EDGE {
	int src;
	int dst;
	float weight;
} EDGE;

typedef struct _VERTEX {
	int id;
	DisjointSet* set;
	struct _VERTEX* a;
	struct _VERTEX* b;
	struct _VERTEX* parent;
	int size;
} VERTEX;

GArray* make_leaves(int N);
GArray* make_edges( DMAT* distance_matrix );
VERTEX* build_tree( GArray* leaves, GArray *edges );
GHashTable* build_k_tree( GArray* leaves, GArray *edges, int k );
GHashTable* build_thresh_tree( GArray* leaves, GArray *edges, float thresh );
int* get_cluster_membership( GHashTable* roots, GArray* leaves );

#endif
