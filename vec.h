#include <glib.h>

#ifndef VEC_H
#define VEC_H

typedef struct _TERM {
	int num;
	float val;
} TERM;

typedef struct _VECFILE {
	int max_dim;
	int max_nnz;
	TERM* terms;
	GArray* points;
	int term_count;
} VECFILE;

typedef struct _DMAT{
	int N;
	GArray* rows;
	float* vals;
	int* idxs;
} DMAT;

VECFILE* load_vec( const char* filename );
DMAT* load_dmat( const char* filename );
#ifdef __cplusplus 
extern "C" {
#endif
void output_dmat( DMAT* dmat, const char* filename );
#ifdef __cplusplus
}
#endif

#endif
