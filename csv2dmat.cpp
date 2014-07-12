#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include "csv.h"
#include "vec.h"
#include "vptree.h"
#include <vector>

using namespace std;

gchar* csv_filename = NULL;
gchar* dmat_filename = NULL;
gint nn_count=6;

static GOptionEntry entries[] =
{
  { "input", 'i', 0, G_OPTION_ARG_FILENAME, &csv_filename, "Input CSV file name", "csvfilename" },
  { "dmat", 'd', 0, G_OPTION_ARG_FILENAME, &dmat_filename, "Output dmat filename", "dmatfilename" },
  { "nncount", 'n', 0, G_OPTION_ARG_INT, &nn_count, "Number of vptree neighbors", "N" },
  { NULL }
};

/* convert apq neighbors to distance matrix format */

DMAT* nn_2_dmat(DENSEFILE* dfile, int* ids, float* ds ) {
	
	int i,j;

	DMAT* dmat = (DMAT*)malloc(sizeof(DMAT));
	dmat->N = dfile->N;
	dmat->rows = g_array_new(FALSE,FALSE,sizeof(GHashTable*));
	dmat->vals = (float*)malloc(dfile->N*nn_count*sizeof(float));
	dmat->idxs = (int*)malloc(dfile->N*nn_count*sizeof(int));
	
	for( i = 0; i < dfile->N; i++ ) {
		GHashTable* row = g_hash_table_new(g_int_hash,g_int_equal);
		g_array_append_val(dmat->rows,row);
	}
	for( i = 0; i < dfile->N; i++ ) {
		GHashTable* row = g_array_index(dmat->rows,GHashTable*,i);
		for( j = 0; j < nn_count; j++ ) {
			if( i != ids[i*nn_count+j]) {
				dmat->idxs[i*nn_count+j] = ids[i*nn_count+j];
				dmat->vals[i*nn_count+j] = ds[i*nn_count+j];
				g_hash_table_insert(row,&(dmat->idxs[i*nn_count+j]),&(dmat->vals[i*nn_count+j]));
			}
		}
	}
	return dmat;
}


int main(int argc, char** argv) {
	
	GError *error = NULL;
	GOptionContext *context;

	context = g_option_context_new ("- compute nearest neighbors from csv file using a vptree");
	g_option_context_add_main_entries (context, entries, NULL);
	if (!g_option_context_parse (context, &argc, &argv, &error)) {
		g_print ("option parsing failed: %s\n", error->message);
		exit (1);
	}
	
	if( csv_filename == NULL ) {
		g_print("Input file missing.  Use -i\n");
		exit(0);
	}
	if( dmat_filename == NULL ) {
		g_print("Output file missing.  Use -d\n");
		exit(0);
	}
	
	DENSEFILE* dfile = load_csv( csv_filename );
    VpTree<DataPoint, euclidean_distance>* tree = new VpTree<DataPoint, euclidean_distance>();
    vector<DataPoint> obj_X(dfile->N, DataPoint(dfile->M, -1, dfile->data));
    for(int n = 0; n < dfile->N; n++) obj_X[n] = DataPoint( dfile->M, n, (dfile->data) + (n * (dfile->M)) );
	int* ids = new int[nn_count*(dfile->N)];
	float* ds = new float[nn_count*(dfile->N)];

    // Loop over all points to find nearest neighbors
    tree->create(obj_X);
    vector<DataPoint> indices;
    vector<float> distances;
    for(int n = 0; n < dfile->N; n++) {
        
        // Find nearest neighbors
        indices.clear();
        distances.clear();
        tree->search(obj_X[n], nn_count + 1, &indices, &distances);
		for( int kk = 0; kk < nn_count; kk++) {
			
			ids[n*nn_count+kk] = indices[kk].index();
			ds[n*nn_count+kk] = distances[kk];
		}
	}	
	
	DMAT* dmat = nn_2_dmat(dfile,ids,ds);
	output_dmat(dmat,dmat_filename);
}
