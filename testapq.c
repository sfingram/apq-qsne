#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <glib.h>
#include "apq.h"
#include "vec.h"
#include "clustertree.h"

/* apq results generator */

char *input_vec_filename = NULL;	// input vec file
char *input_dmat_filename = NULL;	// input distance matrix file

char *output_dmat_filename = NULL;  // output distance matrix file
char *output_nn_filename = NULL;    // output nearest neighbor filename
char *output_pts_filename = NULL;	// output points filename
char *cluster_filename = NULL;		// cluster partition filename

int k_nearest = 5;					// how many nearest neighbors to compute
int perplexity = 3;					// tSNE perplexity (if computed)
int accumulator_runs = 100;			// accumulator runs
char b_compute_full = 0;			// compute the full distance matrix
int cluster_count = 0;				// computes clusters if > 0
float cluster_threshold=0.0;		// threshold for connected components clustering
char b_pointwise = 0;                   

/* command line help */
void usage() {
	printf("testapq [options]\n");
	printf("\tv - vec file input\n");
	printf("\td - dmat file input (vec format)\n");
	printf("\to - dmat file output (vec format)\n");
	printf("\tn - nn file output (octave table format)\n");
	printf("\tp - coordinate file output (octave table format)\n");
	printf("\tk - number of nearest neighbors to output [5]\n");
	printf("\ta - number of accumulator runs to produce [100]\n");
	printf("\tf - flags that we want to compute the full distance matrix\n");
	printf("\ts - compute a hierarchical clustering of at most s clusters\n");
	printf("\tw - compute a hierarchical clustering of of edgeweights greater than w\n");
	printf("\th - clustering output file only works with -s or -w\n");
	printf("\tz - calculate apq pointwise\n");
}

/* catch buggy input params */
void check_args() {
	if( input_vec_filename != NULL && input_dmat_filename != NULL ) {
		fprintf(stderr,"input dmat overrides input vec file\n");
		input_vec_filename = NULL;
	}
	if( input_vec_filename == NULL && input_dmat_filename == NULL) {
		fprintf(stderr,"Need an input filename. Use -v or -d\n");
		exit(0);
	}
	// if( input_dmat_filename != NULL && output_pts_filename == NULL ) {
	// 	fprintf(stderr,"Need coordinate output filename. Use -p\n");
	// 	exit(0);
	// }
}

void output_cluster_file(	const char* filename, 
							int N, 
							int* clusters) {
	int i;
	FILE* fp = NULL;
	fp = fopen(filename,"w");
	for( i = 0; i < N; i++ ) {
		fprintf(fp,"%d\n",clusters[i]);
	}
	fclose(fp);
}

/* parse user input */
void proc_command_args( int argc, char **argv ) {

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
			} else if( argv[i][1] == 'v' ) {
				input_vec_filename = argument;
			} else if( argv[i][1] == 'k' ) {
				k_nearest = atoi(argument);
			} else if( argv[i][1] == 'd' ) {
				input_dmat_filename = argument;
			} else if( argv[i][1] == 'o' ) {
				output_dmat_filename = argument;
			} else if( argv[i][1] == 'n' ) {
				output_nn_filename = argument;
			} else if( argv[i][1] == 'p' ) {
				output_pts_filename = argument;
			} else if( argv[i][1] == 'a' ) {
				accumulator_runs = atoi(argument);
			} else if( argv[i][1] == 'f' ) {
				b_compute_full = 1;
			} else if( argv[i][1] == 's' ) {
				cluster_count = atoi(argument);
			} else if( argv[i][1] == 'h' ) {
				cluster_filename = argument;
			} else if( argv[i][1] == 'w' ) {
				cluster_threshold = atof(argument);
			} else if( argv[i][1] == 'z') {
			  b_pointwise = 1;
			}
			
		}
		i++;
	}
	check_args();
}


int main(int argc, char** argv) {

	int i;
	VECFILE* vecfile = NULL;
	DMAT* dmat = NULL;
	APQ* apq;
	gint64 apq_start,apq_stop,apq_interval = 0;
	gint64 tree_start,tree_stop,tree_interval = 0;
	GHashTable* roots = NULL;
	
	proc_command_args(argc,argv);
	
	/* compute apq */
	
	if( input_vec_filename != NULL ) {
		
	  //printf("A\n");
	  //fflush(stdout);
		vecfile = load_vec(input_vec_filename);
		fprintf(stderr,"N = %d, M = %d, avg NNZ = %0.3f\n",vecfile->points->len,vecfile->max_dim,((double)(vecfile->term_count))/((double)(vecfile->points->len)));

		if( b_pointwise ) {

		  // builds the nearest neighbors point-by-point (much more memory efficient)
		  build_apq_pt( vecfile, 
				accumulator_runs, 
				k_nearest, 
				output_dmat_filename, 
				output_nn_filename );

		  return 0;
		}
		else {
		  apq_start = g_get_monotonic_time();
		  //printf("B\n");
		  fflush(stdout);
		  apq = build_apq(vecfile);
		  //printf("C\n");
		  fflush(stdout);
		  if( ! b_compute_full ) {
		    for( i = 0; i < accumulator_runs; i++ ) {
		      g_print("accumulator run:%d %ldms elapsed\n",i,g_get_monotonic_time()-apq_start);
		      accumulate(apq);
		    }
		  }
		  else {  
		    full_apq(apq);
		  }
		  apq_stop = g_get_monotonic_time();
		  apq_interval = apq_stop - apq_start;
		  fprintf(stdout,"%ld\n", apq_interval/1000);
		  if( ! b_compute_full )
		    dmat = apq_2_dmat(apq,k_nearest);
		  else
		    dmat = apq_2_dmat(apq,-1);

		  if( output_dmat_filename != NULL) {
		    output_dmat( dmat, output_dmat_filename );
		  }
		  if(output_nn_filename != NULL) {
		    output_nn( apq, k_nearest, output_nn_filename );
		  }
		}
	}
	
	/* compute cluster tree */

	if( ( cluster_count == 0 && cluster_threshold > 0.0 ) ||
	    ( cluster_count > 0 && cluster_threshold <= 0.0 )) {
		if( input_dmat_filename ) {
			dmat = load_dmat( input_dmat_filename );
		}
	
		tree_start = g_get_monotonic_time();
		GArray* leaves = make_leaves(dmat->N);
		if( cluster_threshold > 0.0 ) {
			roots = build_thresh_tree(leaves,make_edges(dmat),cluster_threshold);
		}
		else if ( cluster_count > 0 ) {
			roots = build_k_tree(leaves,make_edges(dmat),cluster_count);
		}
		tree_stop = g_get_monotonic_time();
		tree_interval = tree_stop - tree_start;
		
		fprintf(stderr,"Cluster tree time = %ld\n", tree_interval/1000);
		
		if( cluster_filename != NULL ) {
			output_cluster_file(	cluster_filename, 
									leaves->len, 
									get_cluster_membership( roots, leaves ) );
		}
	}
	
	/* compute h-bh-tsne */
	
	// results_list = (GList**)malloc(sizeof(GList*)*N);
	// for( i = 0; i < N; i++ ) {
	// 	results_list[i] = g_hash_table_get_values(accumulator_file[i]);
	// 	results_list[i] = g_list_sort(results_list[i],(GCompareFunc)compare_acc);
	// }
	// gint64 second_time =  g_get_monotonic_time();
	// 
	// fprintf(stderr,"Total time = %ld\n", (second_time - first_time)/1000);

	return 0;
}
