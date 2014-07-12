#ifndef DISJOINTSET_H
#define DISJOINTSET_H

typedef struct _DisjointSet {
	struct _DisjointSet* parent;
	int rank;
} DisjointSet;
DisjointSet* ds_find( DisjointSet* set );
void ds_union( DisjointSet* x, DisjointSet* y );

#endif
