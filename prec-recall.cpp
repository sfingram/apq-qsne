/*
	Calculates mean precision/recall 
*/

#include <stdlib.h>
#include <glib.h>
#include <math.h>

// TYPES

typedef struct PREC_REC {
	gdouble mean_precision;
	gdouble mean_recall;
	gdouble precision_interval;
	gdouble recall_interval;
} PREC_REC;

// GLOBALS

gchar* input_neighbor_file 	= NULL;
gchar* output_neighbor_file	= NULL;
gint num_input_ngb	= 0;
gint num_output_ngb	= 0;
gint num_larger_set = 0;

// COMMAND LINE OPTS

static GOptionEntry entries[] =
{
  { "input_neighbor_file",		'P', 0, G_OPTION_ARG_FILENAME,	&input_neighbor_file, "Input nearest neighbor list", "input_neighbor_file" },
  { "output_neighbor_file", 	'Q', 0, G_OPTION_ARG_FILENAME,	&output_neighbor_file, "Output nearest neighbor list", "output_neighbor_file" },
  { "num_input_ngb", 			'r', 0, G_OPTION_ARG_INT, 	&num_input_ngb, "Input nearest neighbor count", "num_input_ngb" },
  { "num_output_ngb", 			'k', 0, G_OPTION_ARG_INT, 	&num_output_ngb, "Output nearest neighbor count", "num_output_ngb" },
  { "num_larger_set", 			'N', 0, G_OPTION_ARG_INT, 	&num_larger_set, "Size of larger set (for sampling)", "num_larger_set" },
  { NULL }
};

void hash_destroy(gpointer data) {
	GHashTable* hash_table = (GHashTable*)data;
	g_hash_table_destroy( hash_table );
}

/*
	Load an array of hashtables representing nearest neighbor sets from
	a csv file of nearest neighbor indices.
	
	Each line in a csv represents indices of a point whose number is
	indicated by the number of the line on which it appears. Neighbors
	are assumed to be entered in descending rank order.
	
	The num_read parameter only collects the first num_read indices 
	from the file
*/
GPtrArray *load_nns( const gchar *filename,  gint num_read ) {
	
	// create nearest neighbor array
	
	GPtrArray *ret_array = g_ptr_array_new_full( 1, hash_destroy );
	
	// open the file
	
	GError *error = NULL;
	GIOChannel* input_channel = g_io_channel_new_file(filename,"r",&error);
	if( input_channel == NULL ) {
		g_error( "Error opening %s\n", filename );
		exit(0);
	}
	
	// parse lines
	
	gchar* str_return = NULL;
	gsize length = 0;
	gsize terminator_pos = 0;
	GIOStatus g_status = g_io_channel_read_line(	input_channel,
													&str_return,
													&length,
													&terminator_pos,
													&error);
	while( str_return != NULL ) {
		
		GHashTable* table = g_hash_table_new( g_direct_hash, g_direct_equal );
		g_ptr_array_add(	ret_array, 
							table );
		gchar** tokens = g_strsplit(	g_strstrip( str_return ),
										",",
										0);
		gint token_idx = 0;
		while( tokens[token_idx] != NULL ) {
			gpointer ptr_data = GINT_TO_POINTER( (gint)g_ascii_strtoll( tokens[token_idx], NULL, 10 ) );
			if( token_idx < num_read )
				g_hash_table_insert( table, ptr_data, ptr_data );
			token_idx++;
		} 
		g_strfreev( tokens );
		g_free( str_return );
		g_status = g_io_channel_read_line(	input_channel,
											&str_return,
											&length,
											&terminator_pos,
											&error);
	}
	
	// clean up
	
	g_status = g_io_channel_shutdown(	input_channel,
	                                    TRUE,
	                                    &error);
	
	return ret_array;
}

/*
	Calculate mean precision and recall with confidence intervals if necessary
*/
PREC_REC calculate_prec_rec(	GPtrArray* P, 
								GPtrArray* Q, 
								gint r, 
								gint k,
								gint N ) {
	
	gint num_points = P->len;
	PREC_REC ret_prec_rec;
	
	// accumulate the mean and variance using Welford's
	// method (see http://www.johndcook.com/standard_deviation.html )
	
	gdouble M_prec = 0;		// init M and S vars
	gdouble M_rec = 0;
	gdouble S_prec = 0;
	gdouble S_rec = 0;
	
	for( gint i = 0; i < num_points; i++ ) {
		
		// calculate precision and recall for each point
		
		GHashTable* P_i = (GHashTable*)g_ptr_array_index( P, i );
		GHashTable* Q_i = (GHashTable*)g_ptr_array_index( Q, i );
		GList* Q_i_keys = g_hash_table_get_keys(Q_i);
		gint N_tp = 0;
		gint j = 0;
		while( Q_i_keys != NULL && j < k ) {

			// count the true positives

			if( g_hash_table_lookup( P_i, Q_i_keys->data ) != NULL )
				N_tp++;
			Q_i_keys = Q_i_keys->next;
			j++;
		}
		gdouble local_prec = ((gdouble)N_tp)/((gdouble)k);
		gdouble local_rec = ((gdouble)N_tp)/((gdouble)r);
		
		// accumulate prec/rec into rolling mean and variance measures
		
		if( i == 0 ) {
			M_prec = local_prec;
			M_rec = local_rec;
		}
		else {
			gdouble M_prec_old = M_prec;
			gdouble M_rec_old  = M_rec;
			gdouble S_prec_old = S_prec;
			gdouble S_rec_old  = S_rec;
			M_prec = M_prec_old + (local_prec - M_prec_old) / ((gdouble)(i+1));
			M_rec  = M_rec_old  + (local_rec  - M_rec_old)  / ((gdouble)(i+1));
			S_prec = S_prec_old + (local_prec - M_prec_old) * (local_prec - M_prec);
			S_rec  = S_rec_old  + (local_rec  - M_rec_old)  * (local_rec  - M_rec);
		}
	}
	
	ret_prec_rec.mean_precision	= M_prec;
	ret_prec_rec.mean_recall	= M_rec;

	if( N > 0 ) { // are we sampling?
		
		// correct final S accumulation
	
		S_prec = (num_points>1) ? (S_prec/(num_points - 1.0)) : 0.0;
		S_rec  = (num_points>1) ? (S_rec /(num_points - 1.0)) : 0.0;

		// compute 95 % confidence intervals

		ret_prec_rec.precision_interval = 1.96 * S_prec / sqrt( num_points );
		ret_prec_rec.recall_interval    = 1.96 * S_rec  / sqrt( num_points );
	}
	
	return ret_prec_rec;
}

/*
	Check the validity of the parameters
*/
gboolean check_args( ) {
	
	if( input_neighbor_file == NULL ) {
		g_error( "Need input_neighbor_file (use P)\n" );
		return FALSE;
	}
	if( output_neighbor_file == NULL ) {
		g_error( "Need output_neighbor_file (use Q)\n" );
		return FALSE;
	}
	if( num_input_ngb == 0) {
		g_error( "Need num_input_ngb > 0 (use r)\n" );
		return FALSE;
	}
	if( num_output_ngb == 0) {
		g_error( "Need num_output_ngb > 0 (use k)\n" );
		return FALSE;
	}
	
	return TRUE;
}

/*
	prec-recall main
*/
int main (int argc, char* argv[]) {

	GError *error = NULL;
	
	// parse command line
	
	GOptionContext *context = g_option_context_new("- Compute mean precision recall values with optional confidence intervals (for sampled means).");
	g_option_context_add_main_entries (context, entries, NULL);
	if (!g_option_context_parse (context, &argc, &argv, &error)) {
		g_error ("option parsing failed: %s\n", error->message);
		exit ( 1 );
	}
	if( check_args( ) == FALSE ) {
		exit( 1 );
	}

	// load the nearest neighbor sets from disk

	GPtrArray* P = load_nns( input_neighbor_file, num_input_ngb );
	GPtrArray* Q = load_nns( output_neighbor_file, num_output_ngb );

	// calculate precision and recall between P and Q

	PREC_REC prec_rec = calculate_prec_rec(	P,
											Q,
											num_input_ngb,
											num_output_ngb,
											num_larger_set);

	// output results to stdout

	if( num_larger_set > 0 )
		g_print("%lf,%lf,%lf,%lf\n",prec_rec.mean_precision,prec_rec.mean_recall,prec_rec.precision_interval,prec_rec.recall_interval);
	else
		g_print("%lf,%lf\n",prec_rec.mean_precision,prec_rec.mean_recall);

	// done

	return 0;
}