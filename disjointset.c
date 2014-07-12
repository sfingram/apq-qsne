#include <stdlib.h>
#include "disjointset.h"

DisjointSet* ds_find( DisjointSet* set ) {
	if( set->parent != set ) {
		set->parent = ds_find( set->parent );
	}
	return set->parent;
}

void ds_union( DisjointSet* x, DisjointSet* y ) {
	DisjointSet *xRoot = NULL;
	DisjointSet *yRoot = NULL;

	yRoot = ds_find(y);
	xRoot = ds_find(x);
	
		if( xRoot != yRoot ) {
		if(xRoot->rank < yRoot->rank ) {
			xRoot->parent = yRoot;
		}
		else if( xRoot->rank > yRoot->rank ) {
			yRoot->parent = xRoot;
		}
		else {
			yRoot->parent = xRoot;
			xRoot->rank = xRoot->rank + 1;
		}
	}
}
