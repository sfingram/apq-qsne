#include <stdlib.h>
#include <stdio.h>
#include "clustertree.h"

VERTEX* make_leaf( int id ) {
	VERTEX* v = (VERTEX*)malloc(sizeof(VERTEX));
	v->id = id;
	v->parent = NULL;
	v->size = 1;
	v->a = NULL;
	v->b = NULL;
	v->set = (DisjointSet*)malloc(sizeof(DisjointSet));
	v->set->rank=0;
	v->set->parent=v->set;
	return v;
}

VERTEX* make_branch( VERTEX* a, VERTEX* b ) {	
	VERTEX* v = (VERTEX*)malloc(sizeof(VERTEX));
	v->parent = NULL;
	v->id = -1;
	v->a = a;
	v->b = b;
	a->parent = v;
	b->parent = v;
	v->set = (DisjointSet*)malloc(sizeof(DisjointSet));
	v->set->rank=0;
	v->set->parent = v->set;
	ds_union(v->set,a->set);
	ds_union(v->set,b->set);
	v->size = a->size+b->size;
	return v;
}

VERTEX* get_root( VERTEX* v) {
	if(v->parent == NULL)
		return v;
	else
		return get_root(v->parent);
	
	return NULL;
}

int compare_edges(gpointer a, gpointer b) {
	EDGE** x = (EDGE**)a;
	EDGE** y = (EDGE**)b;
	return ((*x)->weight < (*y)->weight)?-1:1;
}

GArray* make_edges( DMAT* distance_matrix ) {
	
	int i = 0;
	int j = 0;
	GHashTable* row;
	GArray* edges = g_array_new(FALSE,FALSE,sizeof(EDGE*));
	EDGE* e = NULL;
	
	for( i = 0; i < distance_matrix->N; i++ ) {
		row = g_array_index(distance_matrix->rows,GHashTable*,i);
		GList* indices = g_hash_table_get_keys(row);
		while(indices != NULL) {
			j = *((int*)(indices->data));
			if( i != j) {
				e = (EDGE*)malloc(sizeof(EDGE));
				e->src = i;
				e->dst = j;
				e->weight = *((float*)g_hash_table_lookup(row,indices->data));
				g_array_append_val(edges,e);
			}
			indices = indices->next;
		}
	}
	
	g_array_sort(edges,(GCompareFunc)compare_edges);
	
	return edges;
}

GArray* make_leaves(int N) {
	int i;
	VERTEX* v;
	GArray *leaves = g_array_new(FALSE,FALSE,sizeof(VERTEX*));
	for( i = 0; i < N; i++ ) {
		v = make_leaf(i);
		g_array_append_val(leaves,v);
	}
	return leaves;
}

int* get_cluster_membership( GHashTable* roots, GArray* leaves ) {
	
	GHashTable* root_lookup = g_hash_table_new(g_direct_hash,g_direct_equal);
	
	GList* root = g_hash_table_get_keys(roots);
	int* root_counts = (int*)malloc(g_hash_table_size(roots)*sizeof(int));
	int i = 0;
	int* cluster_numbers = (int*)malloc(sizeof(int)*(leaves->len));
	
	for( i = 0; i < (int)g_hash_table_size(roots); i++ )
		root_counts[i]=i;
		
	int root_count=0;
	while( root != NULL ) {
		g_hash_table_insert(root_lookup,root->data,root_counts + root_count);
		root = root->next;
		root_count++;
	}
	
	VERTEX*v=NULL;
	for( i = 0; i < (int)(leaves->len); i++ ) {
		
		v=get_root(g_array_index(leaves,VERTEX*,i));
		cluster_numbers[i] = *((int*)g_hash_table_lookup(root_lookup,v));
	}
	
	return cluster_numbers;
}

GHashTable* build_thresh_tree( GArray* leaves, GArray *edges, float thresh ) {

	unsigned int i = 0;
	VERTEX* v = NULL;
	VERTEX* src = NULL;
	VERTEX* dst = NULL;
	VERTEX* src_root = NULL;
	VERTEX* dst_root = NULL;
	VERTEX* new_branch = NULL;

	// make root list
	GHashTable* root_list = g_hash_table_new(g_direct_hash,g_direct_equal);
	for( i = 0; i < leaves->len; i++ ) {
		v = g_array_index(leaves,VERTEX*,i);
		g_hash_table_insert(root_list,v,v);
	}

	// process edges
	for( i = 0; i < edges->len; i++ ) {
		
		EDGE* e = g_array_index(edges,EDGE*,i);
		if( e->weight > thresh ) {
			return root_list;
		}
		
		src = g_array_index(leaves,VERTEX*,e->src);
		dst = g_array_index(leaves,VERTEX*,e->dst);
		if( ds_find(src->set) != ds_find(dst->set) ) {
			src_root = get_root(src);
			dst_root = get_root(dst);
			if( g_hash_table_contains(root_list,src_root) ) {
				g_hash_table_remove(root_list,src_root);
			}
			if( g_hash_table_contains(root_list,dst_root) ) {
				g_hash_table_remove(root_list,dst_root);
			}
			new_branch = make_branch(src_root,dst_root);
			g_hash_table_insert(root_list,new_branch,new_branch);
		}
	}
	
	return root_list;
}

GHashTable* build_k_tree( GArray* leaves, GArray *edges, int k ) {

	unsigned int i = 0;
	VERTEX* v = NULL;
	VERTEX* src = NULL;
	VERTEX* dst = NULL;
	VERTEX* src_root = NULL;
	VERTEX* dst_root = NULL;
	VERTEX* new_branch = NULL;

	// make root list
	GHashTable* root_list = g_hash_table_new(g_direct_hash,g_direct_equal);
	for( i = 0; i < leaves->len; i++ ) {
		v = g_array_index(leaves,VERTEX*,i);
		g_hash_table_insert(root_list,v,v);
	}

	// process edges
	for( i = 0; i < edges->len; i++ ) {

		EDGE* e = g_array_index(edges,EDGE*,i);
		src = g_array_index(leaves,VERTEX*,e->src);
		dst = g_array_index(leaves,VERTEX*,e->dst);
		if( ds_find(src->set) != ds_find(dst->set) ) {
			src_root = get_root(src);
			dst_root = get_root(dst);
			if( g_hash_table_contains(root_list,src_root) ) {
				g_hash_table_remove(root_list,src_root);
			}
			if( g_hash_table_contains(root_list,dst_root) ) {
				g_hash_table_remove(root_list,dst_root);
			}
			new_branch = make_branch(src_root,dst_root);
			g_hash_table_insert(root_list,new_branch,new_branch);
			if( k > 0 && g_hash_table_size(root_list)==((unsigned int)k) )
				return root_list;
		}
	}
	
	return root_list;
}


VERTEX* build_tree( GArray* leaves, GArray *edges ) {
	
	GHashTable* root_list =  build_k_tree(leaves,edges,-1);
	
	// process remaining roots
	GList* root = g_hash_table_get_keys(root_list);
	VERTEX* root_vertex = (VERTEX*)g_hash_table_lookup(root_list,root->data);
	GList* next_root = root->next;
	int root_count = 1;
	while( next_root != NULL ) {
		root_count++;
		VERTEX* next_root_vertex = (VERTEX*)g_hash_table_lookup(root_list,next_root->data);
		root_vertex = make_branch(root_vertex,next_root_vertex);
		next_root = next_root->next;
	}
	
	// printf("processed %d roots\n",root_count);
	return root_vertex;
}
