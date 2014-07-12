/*
 *  tsne.h
 *  Header file for t-SNE.
 *
 *  Created by Laurens van der Maaten.
 *  Copyright 2012, Delft University of Technology. All rights reserved.
 *
 *  Modified by Stephen Ingram in 2013 for sparse distance matrix input
 */

#include "vec.h"
#include "clustertree.h"

#ifndef TSNE_H
#define TSNE_H


static inline double sign(double x) { return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0)); }

class TSNE
{    
public:
	
	TSNE();
	
	DMAT* dmat;
	double theta;
	double perplexity;
	double* Y;
	double* costs;
	bool random_start;
	bool init_startfile;
	int max_iter;
	int stop_lying_iter;
	int mom_switch_iter;
	
    void run();
	void output_csv( const char* filename );

    
private:
	
	void dendrogram_traverse( VERTEX* v, bool isVertical, double top, double left, double width, double height );
    void symmetrizeMatrix(int** row_P, int** col_P, double** val_P);
	void zeroMean(double* X, int N, int D);
    void computeGradient(int* inp_row_P, int* inp_col_P, double* inp_val_P, double* dC );
    double evaluateError(int* row_P, int* col_P, double* val_P);
    void computeGaussianPerplexity(int** _row_P, int** _col_P, double** _val_P);
    double randn();
};

#endif

