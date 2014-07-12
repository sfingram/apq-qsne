#include <glib.h>
#include <glib/gprintf.h>
#include "pqueue.h"
#include "vec.h"

#ifndef APQ_H
#define APQ_H

typedef struct _ACCUM {
	int pt;		// point index
	float acc;  // value to accumulate
} ACCUM;

typedef struct _IFQUAD {
	int cv;		// counter voule
	int iv;		// index into terms
	float pv;	// priority value
	float vv;	// vector value
	size_t pos;
} IFQUAD;

typedef struct _IFPAIR {
	int iv;		// index into points 
	float fv;	// weight balue
} IFPAIR;

typedef struct _APQ {
	int N;
	GHashTable** self_lengths;
	GHashTable* self_lengths_pt;
	GArray** inverted_file;
	pqueue_t** impact_file;
	pqueue_t* impact_file_pt;
	GHashTable** accumulator_file;
	GHashTable* accumulator_file_pt;
	ACCUM* accumulatorstore;
	IFQUAD* quadstore;
	float* false_lengths;
	float false_length;
} APQ;

void accumulate(APQ* apq);
APQ* build_apq(VECFILE* vecfile);
DMAT* apq_2_dmat(APQ* apq,int k);
void output_nn( APQ* apq, int k, const char* filename );
void full_apq( APQ* apq );
void build_apq_pt( VECFILE* vecfile, int accumulator_runs, int k, char* dmat_filename, char* nn_filename );

#endif
