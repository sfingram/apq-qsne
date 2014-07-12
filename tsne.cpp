/*
 *  tsne.cpp
 *  Implementation of both standard and Barnes-Hut-SNE.
 *
 *  Created by Laurens van der Maaten.
 *  Copyright 2012, Delft University of Technology. All rights reserved.
 *
 *  Modified by Stephen ingram in 2013 to take truncated distance matrix input
 *  truncated distance matrices are arrays of hashtables, where distances are keyed by index
 *  arrays and hashtable data structures are built using GLib
 *
 */

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <time.h>
#include "quadtree.h"
#include "tsne.h"

using namespace std;

const char *g_inputfile = NULL;
const char *g_outputfile = NULL;
const char *g_initstartfile = NULL;

TSNE::TSNE() {
	dmat		= NULL;
	theta		= 0.5;
	perplexity	= 3.0;
	Y			= NULL;
	costs		= NULL;
	random_start = true;
	max_iter = 3000;
	stop_lying_iter=1000;
	mom_switch_iter = 250;
}

// Perform t-SNE
void TSNE::run( ) {
    
	int N = dmat->N;
	Y = (double*) malloc(dmat->N * 2 * sizeof(double));
	costs = (double*) calloc(dmat->N, sizeof(double));
    if(Y == NULL || costs == NULL) { printf("A: Memory allocation failed!\n"); exit(1); }

    printf("perplexity = %f, and theta = %f\n", perplexity, theta);
    
    // Set learning parameters
    float total_time = .0;
    clock_t start, end;
	//int max_iter = 1000, stop_lying_iter = 250, mom_switch_iter = 250;
	//int max_iter = 3000, stop_lying_iter = 1000, mom_switch_iter = 250;
	double momentum = .5, final_momentum = .8;
	double eta = 200.0;
    
    // Allocate some memory
	int no_dims = 2;
	
	// printf("N=%d no_dims=%d\n",N,no_dims);
	
    double* dY    = (double*) malloc(N * no_dims * sizeof(double));
    double* uY    = (double*) malloc(N * no_dims * sizeof(double));
    double* gains = (double*) malloc(N * no_dims * sizeof(double));
    if(dY == NULL || uY == NULL || gains == NULL) { printf("B: Memory allocation failed!\n"); exit(1); }
    for(int i = 0; i < N * no_dims; i++)    uY[i] =  .0;
    for(int i = 0; i < N * no_dims; i++) gains[i] = 1.0;
    
    // Normalize input data (to prevent numerical problems)
    printf("Computing input similarities...\n");
    start = clock();
    
    // Compute asymmetric pairwise input similarities

    int* row_P; int* col_P; double* val_P;
    computeGaussianPerplexity(&row_P, &col_P, &val_P);
    
    // Symmetrize input similarities
    symmetrizeMatrix(&row_P, &col_P, &val_P);
    double sum_P = .0;
    for(int i = 0; i < row_P[N]; i++) sum_P += val_P[i];
    for(int i = 0; i < row_P[N]; i++) val_P[i] /= sum_P;

    end = clock();
    
	// double lie_value = 4.00;
	double lie_value = 12.00;

    // Lie about the P-values
	for(int i = 0; i < row_P[N]; i++) val_P[i] *= lie_value;

	// Initialize solution (randomly)
	if( init_startfile ) {
	  // load the file into Y
	  //load_init_startfile(Y);
	}
	else if( random_start ) {
		for(int i = 0; i < N * no_dims; i++) Y[i] = randn() * .0001;
	}
	else{ 
		VERTEX* tree = build_tree(make_leaves(dmat->N),make_edges(dmat));
		dendrogram_traverse( tree, true, -1.0, -1.0, 2.0, 2.0 );
	}
	
	// Perform main training loop
    printf("Done in %4.2f seconds (sparsity = %f)!\nLearning embedding...\n", (float) (end - start) / CLOCKS_PER_SEC, (double) row_P[N] / ((double) N * (double) N));
    start = clock();
	for(int iter = 0; iter < max_iter; iter++) {
        
        // Compute (approximate) gradient
        computeGradient(row_P, col_P, val_P, dY );
        
        // Update gains
        for(int i = 0; i < N * no_dims; i++) gains[i] = (sign(dY[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
        for(int i = 0; i < N * no_dims; i++) if(gains[i] < .01) gains[i] = .01;
            
        // Perform gradient update (with momentum and gains)
        for(int i = 0; i < N * no_dims; i++) uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
		for(int i = 0; i < N * no_dims; i++)  Y[i] = Y[i] + uY[i];
        
        // Make solution zero-mean
		zeroMean(Y, N, 2);
        
        // Stop lying about the P-values after a while, and switch momentum
        if(iter == stop_lying_iter) {
            for(int i = 0; i < row_P[N]; i++) val_P[i] /= lie_value;
        }
        if(iter == mom_switch_iter) momentum = final_momentum;
        
        // Print out progress
        if((iter > 0 && iter % 50 == 0) || iter == max_iter - 1) {
            end = clock();
            double C = .0;
			C = evaluateError(row_P, col_P, val_P);  // doing approximate computation here!
            if(iter == 0)
                printf("Iteration %d: error is %f\n", iter + 1, C);
            else {
                total_time += (float) (end - start) / CLOCKS_PER_SEC;
                printf("Iteration %d: error is %f (50 iterations in %4.2f seconds)\n", iter, C, (float) (end - start) / CLOCKS_PER_SEC);
            }
			start = clock();
        }
    }
    end = clock(); total_time += (float) (end - start) / CLOCKS_PER_SEC;
    
    // Clean up memory
    free(dY);
    free(uY);
    free(gains);
    free(row_P); row_P = NULL;
    free(col_P); col_P = NULL;
    free(val_P); val_P = NULL;
    printf("Fitting performed in %4.2f seconds.\n", total_time);
}


void TSNE::zeroMean(double* X, int N, int D) {
	
	// Compute data mean
	double* mean = (double*) calloc(D, sizeof(double));
    if(mean == NULL) { printf("C:Memory allocation failed!\n"); exit(1); }
	for(int n = 0; n < N; n++) {
		for(int d = 0; d < D; d++) {
			mean[d] += X[n * D + d];
		}
	}
	for(int d = 0; d < D; d++) {
		mean[d] /= (double) N;
	}
	
	// Subtract data mean
	for(int n = 0; n < N; n++) {
		for(int d = 0; d < D; d++) {
			X[n * D + d] -= mean[d];
		}
	}
    free(mean); mean = NULL;
}

// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
void TSNE::computeGradient(int* inp_row_P, int* inp_col_P, double* inp_val_P, double* dC)
{
	int D = 2;
	int N = dmat->N;
    
    // Construct quadtree on current map
    QuadTree* tree = new QuadTree(Y, N);
    
    // Compute all terms required for t-SNE gradient
    double sum_Q = .0;
    double* pos_f = (double*) calloc(N * D, sizeof(double));
    double* neg_f = (double*) calloc(N * D, sizeof(double));
    if(pos_f == NULL || neg_f == NULL) { printf("D:Memory allocation failed!\n"); exit(1); }
    tree->computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, N, pos_f);
    for(int n = 0; n < N; n++) tree->computeNonEdgeForces(n, theta, neg_f + n * D, &sum_Q);
    
    // Compute final t-SNE gradient
    for(int i = 0; i < N * D; i++) {
        dC[i] = pos_f[i] - (neg_f[i] / sum_Q);
    }
    free(pos_f);
    free(neg_f);
    delete tree;
}


// Evaluate t-SNE cost function (approximately)
double TSNE::evaluateError(int* row_P, int* col_P, double* val_P)
{
	int N = dmat->N;
    
    // Get estimate of normalization term
    const int QT_NO_DIMS = 2;
    QuadTree* tree = new QuadTree(Y, N);
    double buff[QT_NO_DIMS] = {.0, .0};
    double sum_Q = .0;
    for(int n = 0; n < N; n++) tree->computeNonEdgeForces(n, theta, buff, &sum_Q);
    
    // Loop over all edges to compute t-SNE error
    int ind1, ind2;
    double C = .0, Q;
    for(int n = 0; n < N; n++) {
        ind1 = n * QT_NO_DIMS;
        for(int i = row_P[n]; i < row_P[n + 1]; i++) {
            Q = .0;
            ind2 = col_P[i] * QT_NO_DIMS;
            for(int d = 0; d < QT_NO_DIMS; d++) buff[d]  = Y[ind1 + d];
            for(int d = 0; d < QT_NO_DIMS; d++) buff[d] -= Y[ind2 + d];
            for(int d = 0; d < QT_NO_DIMS; d++) Q += buff[d] * buff[d];
            Q = (1.0 / (1.0 + Q)) / sum_Q;
            C += val_P[i] * log((val_P[i] + FLT_MIN) / (Q + FLT_MIN));
        }
    }

    delete tree;

    return C;
}


// Compute input similarities with a fixed perplexity using ball trees (this function allocates memory another function should free)
void TSNE::computeGaussianPerplexity(int** _row_P, int** _col_P, double** _val_P) {
    
	int N = dmat->N;
	
	// compute memory requirements
	int i,j,m;
	double w;
	GList* listhead=NULL;
	GList* indices=NULL;
	int num_vals=0;
	for( i = 0; i < dmat->N; i++ ) {
		num_vals += g_hash_table_size(g_array_index(dmat->rows,GHashTable*,i));
	}

    // Allocate the memory we need
    *_row_P = (int*)    malloc((N + 1) * sizeof(int));
    *_col_P = (int*)    calloc(num_vals, sizeof(int));
    *_val_P = (double*) calloc(num_vals, sizeof(double));
    if(*_row_P == NULL || *_col_P == NULL || *_val_P == NULL) { printf("E:Memory allocation failed!%d,%d\n",N+1,num_vals); exit(1); }
    int* row_P = *_row_P;
    int* col_P = *_col_P;
    double* val_P = *_val_P;
    double* cur_P = (double*) malloc((N - 1) * sizeof(double));
    if(cur_P == NULL) { printf("F:Memory allocation failed!\n"); exit(1); }
    row_P[0] = 0;
    for(int n = 0; n < N; n++) row_P[n + 1] = row_P[n] + g_hash_table_size(g_array_index(dmat->rows,GHashTable*,n));
        
    for(int n = 0; n < N; n++) {
        
        
        // Initialize some variables for binary search
		bool found = false;
		double beta = 1.0;
		double min_beta = -DBL_MAX;
		double max_beta =  DBL_MAX;
		double tol = 1e-5;
		int K;
		GHashTable* row = g_array_index(dmat->rows,GHashTable*,n);
		
		// Iterate until we found a good perplexity
		int iter = 0; double sum_P;
		while(!found && iter < 200) {
			
			// Compute Gaussian kernel row
			m = 0;
			indices = g_hash_table_get_keys(row);
			listhead = indices;
			while(indices != NULL) {
				w = *((float*)g_hash_table_lookup(row,indices->data));
				cur_P[m] = exp(-beta * w);
				m++;
				indices = indices->next;
			}
			K = m;
			g_list_free(listhead);
			
			// Compute entropy of current row
			sum_P = DBL_MIN;
			for(m = 0; m < K; m++) sum_P += cur_P[m];
			double H = .0;
			m=0;
			indices = g_hash_table_get_keys(row);
			listhead = indices;
			while(indices != NULL) {
				w = *((float*)g_hash_table_lookup(row,indices->data));
				H += beta * (w * cur_P[m]);
				m++;
				indices = indices->next;
			}
			H = (H / sum_P) + log(sum_P);
			g_list_free(listhead);
			
			// Evaluate whether the entropy is within the tolerance level
			double Hdiff = H - log(perplexity);
			if(Hdiff < tol && -Hdiff < tol) {
				found = true;
			}
			else {
				if(Hdiff > 0) {
					min_beta = beta;
					if(max_beta == DBL_MAX || max_beta == -DBL_MAX)
						beta *= 2.0;
					else
						beta = (beta + max_beta) / 2.0;
				}
				else {
					max_beta = beta;
					if(min_beta == -DBL_MAX || min_beta == DBL_MAX)
						beta /= 2.0;
					else
						beta = (beta + min_beta) / 2.0;
				}
			}
			
			// Update iteration counter
			iter++;
		}
		
		// Row-normalize current row of P and store in matrix
        for(m = 0; m < K; m++) cur_P[m] /= sum_P;
		m=0;
		indices = g_hash_table_get_keys(row);
		listhead = indices;
		while(indices != NULL) {
			j = *((int*)(indices->data));
			col_P[row_P[n] + m] = j;
			val_P[row_P[n] + m] = cur_P[m];
			m++;
			indices = indices->next;
		}
		g_list_free(listhead);
    }
    
    // Clean up memory
    free(cur_P);
}


void TSNE::symmetrizeMatrix(int** _row_P, int** _col_P, double** _val_P) {
    
	int N = dmat->N;
	
    // Get sparse matrix
    int* row_P = *_row_P;
    int* col_P = *_col_P;
    double* val_P = *_val_P;

    // Count number of elements and row counts of symmetric matrix
    int* row_counts = (int*) calloc(N, sizeof(int));
    if(row_counts == NULL) { printf("G:Memory allocation failed!\n"); exit(1); }
    for(int n = 0; n < N; n++) {
        for(int i = row_P[n]; i < row_P[n + 1]; i++) {
            
            // Check whether element (col_P[i], n) is present
            bool present = false;
            for(int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
                if(col_P[m] == n) present = true;
            }
            if(present) row_counts[n]++;
            else {
                row_counts[n]++;
                row_counts[col_P[i]]++;
            }
        }
    }
    int no_elem = 0;
    for(int n = 0; n < N; n++) no_elem += row_counts[n];
    
    // Allocate memory for symmetrized matrix
    int*    sym_row_P = (int*)    malloc((N + 1) * sizeof(int));
    int*    sym_col_P = (int*)    malloc(no_elem * sizeof(int));
    double* sym_val_P = (double*) malloc(no_elem * sizeof(double));
    if(sym_row_P == NULL || sym_col_P == NULL || sym_val_P == NULL) { printf("H:Memory allocation failed!\n"); exit(1); }
    
    // Construct new row indices for symmetric matrix
    sym_row_P[0] = 0;
    for(int n = 0; n < N; n++) sym_row_P[n + 1] = sym_row_P[n] + row_counts[n];
    
    // Fill the result matrix
    int* offset = (int*) calloc(N, sizeof(int));
    if(offset == NULL) { printf("I:Memory allocation failed!\n"); exit(1); }
    for(int n = 0; n < N; n++) {
        for(int i = row_P[n]; i < row_P[n + 1]; i++) {                                  // considering element(n, col_P[i])
            
            // Check whether element (col_P[i], n) is present
            bool present = false;
            for(int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
                if(col_P[m] == n) {
                    present = true;
                    if(n <= col_P[i]) {                                                 // make sure we do not add elements twice
                        sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
                        sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
                        sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i] + val_P[m];
                        sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i] + val_P[m];
                    }
                }
            }
            
            // If (col_P[i], n) is not present, there is no addition involved
            if(!present) {
                sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
                sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
                sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i];
                sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i];
            }
            
            // Update offsets
            if(!present || (present && n <= col_P[i])) {
                offset[n]++;
                if(col_P[i] != n) offset[col_P[i]]++;               
            }
        }
    }
    
    // Divide the result by two
    for(int i = 0; i < no_elem; i++) sym_val_P[i] /= 2.0;
    
    // Return symmetrized matrices
    free(*_row_P); *_row_P = sym_row_P;
    free(*_col_P); *_col_P = sym_col_P;
    free(*_val_P); *_val_P = sym_val_P;
    
    // Free up some memery
    free(offset); offset = NULL;
    free(row_counts); row_counts  = NULL;
}

// Generates a Gaussian random number
double TSNE::randn() {
	double x, y, radius;
	do {
		x = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
		y = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
		radius = (x * x) + (y * y);
	} while((radius >= 1.0) || (radius == 0.0));
	radius = sqrt(-2 * log(radius) / radius);
	x *= radius;
	y *= radius;
	return x;
}

/* output data to a csv */
void TSNE::output_csv( const char* filename ) {
	
	FILE* fp = NULL;
	int i;
	
	if((fp = fopen(filename,"w"))==NULL) {
		fprintf(stderr,"Error: cannot open %s\n",filename);
		exit(0);
	}
	
	for( i = 0; i < dmat->N; i++ )
		fprintf(fp,"%f,%f\n",Y[i*2],Y[i*2+1]);
}

void TSNE::dendrogram_traverse( VERTEX* v, bool isVertical, double top, double left, double width, double height ) {
	
	if( v->a == NULL && v->b == NULL) {
		
		// pick a center spot in the rect
		Y[(v->id)*2+0] = left + width*(rand() / ((double) RAND_MAX + 1));
		Y[(v->id)*2+1] = top + height*(rand() / ((double) RAND_MAX + 1));;
	}
	else {
		
		double proportion = ((double)(v->a->size)) / ((double)(v->size));
		double a_top,b_top,a_left,b_left,a_width,b_width,a_height,b_height;
		if( isVertical ) {
			
			a_top = top;
			a_left = left;
			a_width = width;
			a_height = height*proportion;
			b_top = top+height*proportion;
			b_left = left;
			b_width = width;
			b_height = height*(1.0-proportion);
		}
		else {
			a_top = top;
			a_left = left;
			a_width = width*proportion;
			a_height = height;
			b_top = top;
			b_left = left+width*proportion;
			b_width = width*(1.0-proportion);
			b_height = height;
		}
		
		dendrogram_traverse(v->a, !isVertical, a_top, a_left, a_width, a_height);
		dendrogram_traverse(v->b, !isVertical, b_top, b_left, b_width, b_height);
	}
}


void usage() {
	printf("tsne [options]\n");
	printf("\ti - input truncated distance matrix (vec format)\n");
	printf("\to - input coordinates (csv format)\n");
	printf("\tt - Barnes hut approximation factor theta [0.5]\n");
	printf("\tp - perplexity (density) [3.0]\n");
	printf("\tN - normal iterations\n");
	printf("\tF - high energy iterations\n");
	printf("\th - hierarchical start config [NO ARGUMENT]\n");
}

/* catch buggy input params */
void check_args() {
	if( g_inputfile == NULL ) {
		fprintf(stderr,"Need input filename, use -i\n");
		exit(0);
	}
}

/* parse user input */
void proc_command_args( int argc, char **argv, TSNE* tsne ) {

	int i = 0; 
	char *argument = NULL;

	while( i < argc ) {
		if( ( argv[i][0] == '-' ) && (strlen( argv[i] ) > 1 ) ){
			
			if( argv[i][2] != '\0')
				argument = &(argv[i][2]);
			else
				argument = argv[i+1];
			
			if( argv[i][1] == '?' ) {
				usage(); exit(0);
			} else if( argv[i][1] == 'i' ) {
				g_inputfile = argument;
				tsne->dmat = load_dmat( g_inputfile );
			} else if( argv[i][1] == 't' ) {
				tsne->theta = atof(argument);
			} else if( argv[i][1] == 'p' ) {
				tsne->perplexity = atof(argument);
			} else if( argv[i][1] == 'N' ) {
				tsne->max_iter = atoi(argument);
			} else if( argv[i][1] == 'F' ) {
				tsne->stop_lying_iter = atoi(argument);
			} else if( argv[i][1] == 'o' ) {
				g_outputfile = argument;
			} else if( argv[i][1] == 'h' ) {
				tsne->random_start = false;
			} else if( argv[i][1] == 's') {
			  tsne->init_startfile = true;
			}
		}
		i++;
	}
	check_args();
}

// Function that runs the Barnes-Hut implementation of t-SNE
int main(int argc, char** argv) {
    
    TSNE* tsne = new TSNE();

	srand(time(NULL));
	proc_command_args(argc,argv,tsne);

    // Read the parameters and the dataset
	if(tsne->dmat != NULL) {
        
		// Now fire up the SNE implementation
		tsne->run( );
		
		// write to disk
		if( g_outputfile != NULL)
			tsne->output_csv(g_outputfile);
    }
}
