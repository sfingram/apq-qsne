#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <time.h>
#include <iostream>
#include <sys/time.h>
#include <ctime>
#include <vector>
#include "vptree.h"
#include "vec.h"

using namespace std;

char *g_dmat_filename = NULL; // for ground truth
char *g_csv_filename  = NULL; // for testing accuracy
char *g_true_nn_filename = NULL; // for ground truth
char *g_test_nn_filename = NULL;
int k_nearest = 5;

typedef struct _DMATPAIR {
	int key;
	float value;
} DMATPAIR;

void usage() {
	printf("layoutquality [options]\n");
	printf("\td - input truncated distance matrix (vec format)\n");
	printf("\ti - input coordinates (csv format)\n");
	printf("\tn - input true nn file\n");
	printf("\tt - input test nn file\n");
	printf("\tk - k nearest neighbors to check (if they exist in the distance matrix)\n");
}

/* catch buggy input params */
void check_args() {
	if( g_dmat_filename == NULL && g_true_nn_filename == NULL ) {
		fprintf(stderr,"Need input distance matrix or true nn file, use -d or -n\n");
		exit(0);
	}
	if( g_csv_filename == NULL  && g_test_nn_filename == NULL ) {
		fprintf(stderr,"Need input coordinates or test nn file , use -i or -t\n");
		exit(0);
	}
}

/* parse user input */
void proc_command_args( int argc, char **argv) {

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
			} else if( argv[i][1] == 'd' ) {
				g_dmat_filename = argument;
			} else if( argv[i][1] == 'i' ) {
				g_csv_filename = argument;
			} else if( argv[i][1] == 'n' ) {
				g_true_nn_filename = argument;
			} else if( argv[i][1] == 't' ) {
				g_test_nn_filename = argument;
			} else if( argv[i][1] == 'k' ) {
				k_nearest = atoi(argument);
			}
		}
		i++;
	}
	check_args();
}

int count_nns(const char* filename) {
	
	char line[65536];
	int line_num = 0;
	FILE* fp = NULL;
	fp = fopen( filename, "r" );
	
	while( fgets( line, 65535, fp ) != NULL ) {
		line_num++;
	}
	
	fclose(fp);
	
	return line_num;
}

int* load_nns( int N, int k, const char* filename ) {
	
	char line[65536];
	char item[256];
	int line_num = 0;
	int i,j,kk;
	FILE* fp = NULL;
	fp = fopen( filename, "r" );
	
	int* nns = (int*) calloc( N*k, sizeof(int) );
	
	while( fgets( line, 65535, fp ) != NULL ) {
		i = 0;
		j=0;
		kk = line_num*k;
		while( line[i] != '\0' && line[i] != '\n' ) {
			if( line[i] == ',') {
				item[j]='\0';
				nns[kk++]=atoi(item);
				j=0;
			}
			else {
				item[j++]=line[i];
			}
			i++;
		}
		item[j]='\0';
		nns[kk]=atoi(item);
		line_num++;
	}
	
	fclose(fp);
	
	return nns;
}

float* load_coordinates( int N, const char* filename ) {
	
	char line[65536];
	char item[256];
	int line_num = 0;
	int i,j;
	FILE* fp = NULL;
	fp = fopen( filename, "r" );
	
	float* coordinates = (float*) malloc( N*2*sizeof(float) );
	
	while( fgets( line, 65535, fp ) != NULL ) {
		i = 0;
		j=0;
		while( line[i] != '\0' && line[i] != '\n' ) {
			if( line[i] == ',') {
				item[j]='\0';
				coordinates[line_num*2]=atof(item);
				j=0;
			}
			item[j++]=line[i];
			i++;
		}
		item[j]='\0';
		coordinates[line_num*2+1]=atof(item);
		line_num++;
	}
	
	fclose(fp);
	
	return coordinates;
}

int compare_dp(gpointer a, gpointer b) {
	DMATPAIR* x = (DMATPAIR*)a;
	DMATPAIR* y = (DMATPAIR*)b;
	return (x->value < y->value)?-1:1;
}

inline int intmin( int a, int b) {return (a<b)?a:b;}
inline int intmax( int a, int b) {return (a>b)?a:b;}

int main (int argc, char* argv[]) {

	int N = 0;

	proc_command_args(argc,argv);
	
	bool b_use_dmat = (g_dmat_filename != NULL);
	bool b_use_coords = (g_csv_filename != NULL);
	
	int* ground_truth_nn = NULL;
	int* test_set_nn = NULL;
	
	// load ground truth nns
	
	if( b_use_dmat ) {
		DMAT* dmat = load_dmat( g_dmat_filename );
		N = dmat->N;
		ground_truth_nn = (int*)calloc(dmat->N*k_nearest,sizeof(int));
		for( int i = 0; i < dmat->N; i++ ) {

			GHashTable* row = g_array_index(dmat->rows,GHashTable*,i);
			GArray* pairs = g_array_sized_new( FALSE, FALSE, sizeof(DMATPAIR), g_hash_table_size(row));
			GList* indices = g_hash_table_get_keys(row);
			while(indices != NULL) {
				if( *((int*)(indices->data)) != i ) {
					DMATPAIR dmp = { *((int*)(indices->data)) ,*((float*)g_hash_table_lookup(row,(int*)(indices->data)))};
					g_array_append_val(pairs,dmp);
				}
				indices = indices->next;
			}
			g_array_sort(pairs,(GCompareFunc)compare_dp);
			for( int j = 0; j < k_nearest && j < pairs->len; j++ ) 
				ground_truth_nn[i*k_nearest+j] = g_array_index(pairs,DMATPAIR,j).key;
			g_array_free(pairs,TRUE);
		}
	}
	else {
		N = count_nns(g_true_nn_filename);
		ground_truth_nn = load_nns(N,k_nearest,g_true_nn_filename);
	}
	
	// load test set nns
	
	if( b_use_coords ) {
		
		float* coords = load_coordinates( N, g_csv_filename );
		
		// load the data into a structure
	
	    VpTree<DataPoint, euclidean_distance>* tree = new VpTree<DataPoint, euclidean_distance>();
	    vector<DataPoint> obj_X(N, DataPoint(2, -1, coords));
	    for(int n = 0; n < N; n++) obj_X[n] = DataPoint(2, n, coords + (n * 2));
		test_set_nn = new int[k_nearest*N];

	    // Loop over all points to find nearest neighbors in coords

	    // cerr << "Building tree...";
	    tree->create(obj_X);
		// cerr << "done." << endl;
		// cerr << "Computing coordinate knn...";
	    vector<DataPoint> indices;
	    vector<float> distances;
	    for(int n = 0; n < N; n++) {
        
	        // Find nearest neighbors
	        indices.clear();
	        distances.clear();
	        tree->search(obj_X[n], k_nearest+1, &indices, &distances);
			int found_self = 0;
			int ids_idx = 0;
			for( int kk = 0; kk < k_nearest+1; kk++) {
				if( indices[kk].index() == n )
					found_self++;
				else {
					test_set_nn[n*k_nearest+ids_idx] = indices[kk].index();
					ids_idx++;
				}
				if( !found_self && kk == k_nearest-1)
					break;
			}
		}	
	}
	else {
		
		test_set_nn = load_nns(N,k_nearest,g_test_nn_filename);
	}
	
	// cerr << "done." << endl;

	// calculate the size of the intersection of the coordinate nearest neighbors
	double* intersection_length = (double*) calloc( N, sizeof(double) );
	
	// output the result
	for( int i = 0; i < N; i++ ) {

		for( int j = 0; j < k_nearest; j++ ) {
			for( int k = 0; k < k_nearest; k++ ) {
				// if( i == 0) 
				// 	cout << "" << test_set_nn[i*k_nearest+j] << "==" << ground_truth_nn[i*k_nearest+k] << endl;
				if( test_set_nn[i*k_nearest+j] == ground_truth_nn[i*k_nearest+k] ) {
					intersection_length[i]+=1.0;
					break;
				}
			}
		}
		intersection_length[i] /= k_nearest;
	}
	
	double mean_intersection_length = 0.0;
	
	for( int i = 0; i < N; i++ ) {
		mean_intersection_length += (intersection_length[i]-mean_intersection_length)/((double)(i+1));
	}
	mean_intersection_length -= ((double)k_nearest)/((double)N - 1.0);
	
	cout << "" << k_nearest << "," << mean_intersection_length << endl;
	
	return 0;
}