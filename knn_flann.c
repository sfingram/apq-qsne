#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include "flann/flann.h"

// TYPES

typedef struct COORD_FILE {
	float* dataset;
	int rows;
	int cols;
} COORD_FILE;

// GLOBALS

gchar* input_coordinates 	= NULL;
int nn = 5;
gdouble precision = 0.9;
gboolean output_indices = FALSE;
gboolean use_kdtree = TRUE;

// COMMAND LINE OPTS

static GOptionEntry entries[] =
{
  { "input_coordinates",	'i', 0, G_OPTION_ARG_FILENAME,	&input_coordinates, "Input coordinates", "input_coordinates" },
  { "nn", 					'n', 0, G_OPTION_ARG_INT, 	&nn, "Nearest Neighbors", "nn" },
  { "output_indices", 		'O', 0, G_OPTION_ARG_NONE, 	&output_indices, "Output indices to stdout", "output_indices" },
  { "precision", 			'p', 0, G_OPTION_ARG_DOUBLE, &precision, "Precision of nearest neighbors [0.9]", "output_indices" },
  { "use_kdtree", 			'k', 0, G_OPTION_ARG_NONE, &use_kdtree, "Use a kdtree for nn search", "use_kdtree" },
  { NULL }
};

/*

*/
COORD_FILE load_coords( const gchar *filename ) {
	
	COORD_FILE ret_file;
	
	// open the file
	
	GError *error = NULL;
	GIOChannel* input_channel = g_io_channel_new_file(filename,"r",&error);
	if( input_channel == NULL ) {
		g_error( "Error opening %s\n", filename );
		exit(0);
	}
	
	// parse lines
	
	int num_points = 0;
	int num_columns = 0;
	
	gchar* str_return = NULL;
	gsize length = 0;
	gsize terminator_pos = 0;
	GIOStatus g_status = g_io_channel_read_line(	input_channel,
													&str_return,
													&length,
													&terminator_pos,
													&error);
	while( str_return != NULL ) {
		
		if( num_columns == 0 ) {
			gchar** tokens = g_strsplit(	g_strstrip( str_return ),
											",",
											0);
			gint token_idx = 0;
			while( tokens[token_idx] != NULL ) {
				token_idx++;
			} 
			g_strfreev( tokens );
			num_columns = token_idx;
		}
		g_free( str_return );
		g_status = g_io_channel_read_line(	input_channel,
											&str_return,
											&length,
											&terminator_pos,
											&error);
		num_points++;
	}
	
	ret_file.dataset = (float*) malloc(sizeof(float)*num_points*num_columns);
	ret_file.rows = num_points;
	ret_file.cols = num_columns;
	
	// clean up
	
	g_status = g_io_channel_shutdown(	input_channel,
	                                    TRUE,
	                                    &error);
	
	// reopen the file
	
	input_channel = g_io_channel_new_file(filename,"r",&error);
	if( input_channel == NULL ) {
		g_error( "Error re-opening %s\n", filename );
		exit(0);
	}
	
	// parse lines
	
	num_points = 0;
	
	str_return = NULL;
	length = 0;
	terminator_pos = 0;
	g_status = g_io_channel_read_line(	input_channel,
													&str_return,
													&length,
													&terminator_pos,
													&error);
	while( str_return != NULL ) {
		
		num_columns = 0;
		gchar** tokens = g_strsplit(	g_strstrip( str_return ),
										",",
										0);
		gint token_idx = 0;
		while( tokens[token_idx] != NULL ) {
			ret_file.dataset[num_points*ret_file.cols+token_idx] = 
				(float)g_ascii_strtod( tokens[token_idx], NULL );
			token_idx++;
		} 
		g_strfreev( tokens );
		num_columns = token_idx;
		
		g_free( str_return );
		g_status = g_io_channel_read_line(	input_channel,
											&str_return,
											&length,
											&terminator_pos,
											&error);
		num_points++;
	}
	
	// clean up
	
	g_status = g_io_channel_shutdown(	input_channel,
	                                    TRUE,
	                                    &error);
	
	return ret_file;
}


int main (int argc, char* argv[]) {
	
	GError *error = NULL;
	
	// parse command line
	
	GOptionContext *context = g_option_context_new("- Compute mean precision recall values with optional confidence intervals (for sampled means).");
	g_option_context_add_main_entries (context, entries, NULL);
	if (!g_option_context_parse (context, &argc, &argv, &error)) {
		g_error ("option parsing failed: %s\n", error->message);
		exit ( 1 );
	}
	
	if( input_coordinates == NULL ) {
		g_error("Problem opening %s\n",input_coordinates);
		exit(1);
	}
	
	COORD_FILE coord_file = load_coords( input_coordinates );
	
	/* allocate memory for the nearest-neighbors indices */
	int* result = (int*) malloc(coord_file.rows*nn*sizeof(int));
	
	/* allocate memory for the distances */
	float* dists = (float*) malloc(coord_file.rows*nn*sizeof(float));
	
	/* index parameters are stored here */
	struct FLANNParameters p = DEFAULT_FLANN_PARAMETERS;
	p.algorithm = FLANN_INDEX_AUTOTUNED; 
	p.target_precision = precision; 
	if( use_kdtree ) {
		p.algorithm = FLANN_INDEX_KDTREE;
		p.target_precision = -1;
	}
	
	flann_find_nearest_neighbors(coord_file.dataset, coord_file.rows, coord_file.cols, coord_file.dataset, coord_file.rows,
		result, dists, nn, &p);
	
	if( output_indices ) {
		
		for( int i = 0; i < coord_file.rows; i++ ) {
			for( int j = 0; j < nn; j++ ) {
				printf("%d",result[i*nn+j]);
				if( j < nn-1 )
					printf(",");
			}
			printf("\n");
		}
	}
	
	free(coord_file.dataset);
	free(result);
	free(dists);
	
	return 0;
}