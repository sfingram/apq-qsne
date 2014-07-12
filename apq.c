#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "apq.h"

#define MAX_ACCUMULATORS 500

int accumptr	= 0;
int output_format = 0;

int compare_floats(gpointer a, gpointer b) {
	float* x = (float*)a;
	float* y = (float*)b;
	return (int)(*x - *y);
}
int compare_acc(gpointer a, gpointer b) {
	ACCUM* x = (ACCUM*)a;
	ACCUM* y = (ACCUM*)b;
	return (x->acc > y->acc)?-1:1;
}
int compare_if(gpointer a, gpointer b) {
	IFPAIR* x = (IFPAIR*)a;
	IFPAIR* y = (IFPAIR*)b;
	return (x->fv > y->fv)?-1:1;
}

pqueue_pri_t my_pqueue_get_pri_f( void* a ) {
	return ((IFQUAD*)a)->pv;
}
void my_pqueue_set_pri_f(void *a, pqueue_pri_t pri) {
	((IFQUAD*)a)->pv = (float)pri;
}
int my_pqueue_cmp_pri_f(pqueue_pri_t next, pqueue_pri_t curr) {
	return (int)(next < curr);
}
size_t my_pqueue_get_pos_f(void *a) {
	return ((IFQUAD*)a)->pos;
}
void my_pqueue_set_pos_f(void *a, size_t pos) {
	((IFQUAD*)a)->pos = pos;
}

void acc_destroyer(gpointer data) {

  g_slice_free(ACCUM,data);
}

/* traverse the whole inverted file without impact ordering and generate all accumulators */
void full_apq( APQ* apq ) {
	
	int i;
	unsigned int j;
	for( i = 0; i < apq->N; i++ ) {
		IFQUAD* val = NULL;
		while( (val = (IFQUAD*) pqueue_pop(apq->impact_file[i])) != NULL ) {
			for( j = 0; j < apq->inverted_file[val->iv]->len; j++ ) {
			
				IFPAIR a = g_array_index(apq->inverted_file[val->iv],IFPAIR,j);
				ACCUM temp = {a.iv,val->pv};
				if( a.iv != i ) {	// avoid computing self distances
			
					if( g_hash_table_contains(apq->accumulator_file[i],&(temp.pt))) {
						ACCUM* temp_accum = (ACCUM*) g_hash_table_lookup(apq->accumulator_file[i],&(temp.pt));
						temp_accum->acc += val->vv * a.fv;
					}
					else {
				
						if( g_hash_table_size(apq->accumulator_file[i]) < MAX_ACCUMULATORS ) {
							apq->accumulatorstore[accumptr].pt = a.iv;
							apq->accumulatorstore[accumptr].acc = val->vv * a.fv;
							g_hash_table_insert(apq->accumulator_file[i],&(apq->accumulatorstore[accumptr].pt),apq->accumulatorstore + accumptr);
							accumptr++;
						}
					}
				
					// maintain self lengths
					if( ! g_hash_table_contains(apq->self_lengths[i],&(val->iv))) {
						g_hash_table_insert(apq->self_lengths[i],&(val->iv),&(val->vv));
						apq->false_lengths[i] += (val->vv)*(val->vv);
					}
				}
			}
		}
	}	
}

/* 
compute the accumulators and self-lengths per point
 */
void pt_accumulate(APQ* apq, int num_runs, int pt_idx, FILE* raw_acc_file ) {

  int run_idx;
  for( run_idx = 0; run_idx < num_runs; run_idx++ ) {
    IFQUAD* val = (IFQUAD*)pqueue_peek(apq->impact_file_pt);
    if( val->pv <= 0.f )
      break;;
    IFPAIR a = g_array_index(apq->inverted_file[val->iv],IFPAIR,val->cv);
    ACCUM temp = {a.iv,val->pv};
    if( a.iv != pt_idx ) {	// avoid computing self distances
      
      if( g_hash_table_contains(apq->accumulator_file_pt,&(temp.pt))) {
	ACCUM* temp_accum = (ACCUM*) g_hash_table_lookup(apq->accumulator_file_pt,&(temp.pt));
	temp_accum->acc += temp.acc;
      }
      else {
	
	ACCUM* acc_temp = (ACCUM*) g_slice_new( ACCUM );
	acc_temp->pt = a.iv;
	acc_temp->acc = val->pv;
	g_hash_table_insert(apq->accumulator_file_pt,&(acc_temp->pt),acc_temp);
      }
      
      // maintain self lengths
      if( ! g_hash_table_contains(apq->self_lengths_pt,&(val->iv))) {
	g_hash_table_insert(apq->self_lengths_pt,&(val->iv),&(val->vv));
	apq->false_lengths[pt_idx] += (val->vv)*(val->vv);
      }
    }
    val->cv++;
    pqueue_pri_t tempprio = 0.0;
    if( val->cv < (int)(apq->inverted_file[val->iv]->len) )
      tempprio = (pqueue_pri_t)(val->vv * g_array_index(apq->inverted_file[val->iv],IFPAIR,val->cv).fv);
    pqueue_change_priority(apq->impact_file_pt,tempprio,val);
  }

  GList* results_list_head = g_hash_table_get_values(apq->accumulator_file_pt);
  results_list_head = g_list_sort(results_list_head,(GCompareFunc)compare_acc);
  GList* results_list = results_list_head;
  while( results_list != NULL ) {
    ACCUM* acc_temp = (ACCUM*) (results_list->data);
    fprintf(raw_acc_file,"%d,%f",acc_temp->pt,acc_temp->acc);
    results_list = results_list->next;
    if( results_list != NULL ) {
      fprintf(raw_acc_file,",");
    }
  }
  fprintf(raw_acc_file,"\n");
  g_list_free(results_list_head);

  // complete the length calculation
  apq->false_lengths[pt_idx] = sqrtf(apq->false_lengths[pt_idx]);
}

/* add an accumulator to the data */
void accumulate(APQ* apq) {
	int i;
	for( i = 0; i < apq->N; i++ ) {
		IFQUAD* val = (IFQUAD*)pqueue_peek(apq->impact_file[i]);
		if( val->pv <= 0.f )
			continue;
		IFPAIR a = g_array_index(apq->inverted_file[val->iv],IFPAIR,val->cv);
		ACCUM temp = {a.iv,val->pv};
		if( a.iv != i ) {	// avoid computing self distances
			
			if( g_hash_table_contains(apq->accumulator_file[i],&(temp.pt))) {
				ACCUM* temp_accum = (ACCUM*) g_hash_table_lookup(apq->accumulator_file[i],&(temp.pt));
				temp_accum->acc += temp.acc;
			}
			else {
				
			  ACCUM* acc_temp = (ACCUM*) g_slice_new( ACCUM );
			  acc_temp->pt = a.iv;
			  acc_temp->acc = val->pv;
			  g_hash_table_insert(apq->accumulator_file[i],&(acc_temp->pt),acc_temp);
			}
			
			// maintain self lengths
			if( ! g_hash_table_contains(apq->self_lengths[i],&(val->iv))) {
				g_hash_table_insert(apq->self_lengths[i],&(val->iv),&(val->vv));
				apq->false_lengths[i] += (val->vv)*(val->vv);
			}
		}
		val->cv++;
		pqueue_pri_t tempprio = 0.0;
		if( val->cv < (int)(apq->inverted_file[val->iv]->len) )
			tempprio = (pqueue_pri_t)(val->vv * g_array_index(apq->inverted_file[val->iv],IFPAIR,val->cv).fv);
		pqueue_change_priority(apq->impact_file[i],tempprio,val);
	}
}

void build_apq_pt( VECFILE* vecfile, int accumulator_runs, int k, char* dmat_filename, char* nn_filename ) {

  int i,j;
  int N = vecfile->points->len;
  int maxDim = vecfile->max_dim;
	
  gint64 apq_start = g_get_monotonic_time();
  
  APQ* apq = (APQ*)malloc(sizeof(APQ));

  // open tempfiles

  char *tempapq_filename = NULL;

  if( dmat_filename == NULL ) {
    
    if( nn_filename == NULL )
      return;

    tempapq_filename = (char*) malloc(sizeof(char)*(strlen(nn_filename)+strlen(".tempapq.txt")+1));
    sprintf(tempapq_filename,"%s.tempapq.txt",nn_filename);
  }
  else {

    tempapq_filename = (char*) malloc(sizeof(char)*(strlen(dmat_filename)+strlen(".tempapq.txt")+1));
    sprintf(tempapq_filename,"%s.tempapq.txt",dmat_filename);
  }
  FILE* raw_acc_file = fopen(tempapq_filename, "w");

  apq->N = N;
  if( (apq->quadstore = (IFQUAD*)malloc((maxDim+1)*sizeof(IFQUAD))) == NULL ) {
    fprintf(stderr,"Error: malloc failed on quadstore.\n");
  }
  if( (apq->inverted_file = (GArray **)malloc((maxDim+1)*sizeof(GArray*))) == NULL ) {
    fprintf(stderr,"Error: malloc failed on inverted_file.\n");
  }	
  if( (apq->false_lengths = (float*)calloc(N,sizeof(float))) == NULL ) {
    fprintf(stderr,"Error: malloc failed on false_lengths.\n");
  }

  // allocate inverted file arrays and accumulator hashes
  for( i = 0 ; i < maxDim+1; i++ ) {
    apq->inverted_file[i] = g_array_new(FALSE,FALSE,sizeof(IFPAIR));
  }
	
  // fill the inverted file array
  for( i = 0; i < N; i++ ) {
		
    GArray*terms=g_array_index(vecfile->points,GArray*,i);
    for( j = 0; j < (int)(terms->len); j++) {
      IFPAIR ifp = {i,g_array_index(terms,TERM*,j)->val};
      g_array_append_val(apq->inverted_file[g_array_index(terms,TERM*,j)->num],ifp);
    }
  }
	
  // sort into impact order
  for( i = 0 ; i < maxDim+1; i++ ) {
    g_array_sort(apq->inverted_file[i],(GCompareFunc)compare_if);
  }

  int footest = 0;
		
  // perform the APQ runs per point
  for( i = 0; i < N; i++ ) {

    if( i % 1000 == 0 )
	printf("processing point %d\n",i);
		
    apq->accumulator_file_pt = g_hash_table_new(g_int_hash,g_int_equal);
    apq->self_lengths_pt = g_hash_table_new(g_int_hash,g_int_equal);
    apq->impact_file_pt = pqueue_init(100,
				      (pqueue_cmp_pri_f) my_pqueue_cmp_pri_f,
				      (pqueue_get_pri_f) my_pqueue_get_pri_f,
				      (pqueue_set_pri_f) my_pqueue_set_pri_f,
				      (pqueue_get_pos_f) my_pqueue_get_pos_f,
				      (pqueue_set_pos_f) my_pqueue_set_pos_f);
    GArray*terms=g_array_index(vecfile->points,GArray*,i);
    gint quadcount = 0;
     for( j = 0; j < (int)(terms->len); j++) {
       apq->quadstore[quadcount].cv = 0;
       apq->quadstore[quadcount].iv = g_array_index(terms,TERM*,j)->num;
       apq->quadstore[quadcount].vv = g_array_index(terms,TERM*,j)->val;
       apq->quadstore[quadcount].pv = g_array_index(apq->inverted_file[apq->quadstore[quadcount].iv],IFPAIR,0).fv * apq->quadstore[quadcount].vv;
       pqueue_insert(apq->impact_file_pt, &(apq->quadstore[quadcount]));
       quadcount++;
     }

     // compute and output apq data

     pt_accumulate(apq, accumulator_runs, i, raw_acc_file );

     // free up the resources for the point

     GList* keyshead = g_hash_table_get_keys( apq->accumulator_file_pt );
     GList* keys = keyshead;
     while( keys != NULL ) {
       g_slice_free( ACCUM, g_hash_table_lookup(apq->accumulator_file_pt,keys->data) );
       keys = keys->next;
     }
     g_list_free(keyshead);
     g_hash_table_destroy(apq->accumulator_file_pt);
     g_hash_table_destroy(apq->self_lengths_pt);
     pqueue_free(apq->impact_file_pt);
   }
   fclose(raw_acc_file);

   // reopen temp file
   i = 0;
   GError *error = NULL;
   GIOChannel* input_channel = g_io_channel_new_file(tempapq_filename,"r",&error);
   if( input_channel == NULL ) {
     g_error( "Error opening %s\n", tempapq_filename );
     exit(0);
   }
   GIOChannel* output_channel = g_io_channel_new_file( dmat_filename, "w",&error);
   if( output_channel == NULL ) {
     g_error( "Error opening %s\n", dmat_filename );
     exit(0);
   }  

   // read in, correct, and output the accumulators
   gchar* str_return = NULL;
   gsize length = 0;
   gsize terminator_pos = 0;
   GIOStatus g_status = g_io_channel_read_line(	input_channel,
						 &str_return,
						 &length,
						 &terminator_pos,
						 &error);
   ACCUM* inner_accum = NULL;
   gchar file_line[65536];
   while( str_return != NULL ) {

     // parse the line into comma delimited tokens
     GList* vec_row = NULL;    
     gchar** tokens = g_strsplit(	g_strstrip( str_return ),
					 ",",
					 0);
     gint token_idx = 0;
     while( tokens[token_idx] != NULL ) {

       // on even token, create a new accumulated item
       if( token_idx % 2 == 0 ) {
	 inner_accum = (ACCUM*)g_slice_new( ACCUM );
	 inner_accum->pt = (int) g_ascii_strtoll( tokens[token_idx],NULL,10);
      }

      // on odd token, add it to the row
      if( token_idx % 2 == 1 ) {
	inner_accum->acc = (float) g_ascii_strtod( tokens[token_idx], NULL );
	inner_accum->acc /= (apq->false_lengths[i] * apq->false_lengths[inner_accum->pt]);
	vec_row = g_list_prepend(vec_row,inner_accum);
      }
      token_idx++;
    } 
    
    vec_row = g_list_sort(vec_row,(GCompareFunc)compare_acc);
    int inner_k = 0;
    GList* row_traversal = vec_row;
    gsize bytes_written;
    if( row_traversal == NULL ) {

      g_sprintf(file_line, "\n" );
      g_io_channel_write_chars( output_channel, file_line, -1, &bytes_written, &error );
    }
    while(inner_k < k && row_traversal != NULL ) {
      
      inner_accum = (ACCUM*) row_traversal->data;
      g_sprintf(file_line, "(%d,%0.5f)", inner_accum->pt+1,1.f-inner_accum->acc );
      g_io_channel_write_chars( output_channel, file_line, -1, &bytes_written, &error );
      row_traversal = row_traversal->next;
      if( row_traversal != NULL && inner_k < k-1 ) {
	g_sprintf(file_line, " " );
	g_io_channel_write_chars( output_channel, file_line, -1, &bytes_written, &error );
      }
      else {
	footest++;
	g_sprintf(file_line, "\n" );
	g_io_channel_write_chars( output_channel, file_line, -1, &bytes_written, &error );
      }
      inner_k++;
    }
    g_list_free_full( vec_row, acc_destroyer );

    g_strfreev( tokens );
    g_free( str_return );
    g_status = g_io_channel_read_line(	input_channel,
					&str_return,
					&length,
					&terminator_pos,
					&error);
    i++; // increment our point index
  }
  
  // clean up
  
  g_status = g_io_channel_shutdown(	input_channel,
					TRUE,
					&error);
  g_status = g_io_channel_shutdown(	output_channel,
					TRUE,
					&error);
  remove(tempapq_filename);
  free(tempapq_filename);
  
  gint64 apq_stop = g_get_monotonic_time();
  gint64 apq_interval = apq_stop - apq_start;
  fprintf(stdout,"%ld\n", apq_interval/1000);
}

/* initialize the apq */
APQ* build_apq(VECFILE* vecfile) {
	
	int i,j;
	int N = vecfile->points->len;
	int maxDim = vecfile->max_dim;
	
	APQ* apq = (APQ*)malloc(sizeof(APQ));
	
	apq->N = N;
	if( (apq->quadstore = (IFQUAD*)malloc((vecfile->term_count)*sizeof(IFQUAD))) == NULL ) {
	  fprintf(stderr,"Error: malloc failed on quadstore.\n");
	}
	/* if( (apq->accumulatorstore = (ACCUM*)malloc(N*MAX_ACCUMULATORS*sizeof(ACCUM))) == NULL ) { */
	/*   fprintf(stderr,"Error: malloc failed on accumulatorstore.\n"); */
	/* }; */
	if( (apq->inverted_file = (GArray **)malloc((maxDim+1)*sizeof(GArray*))) == NULL ) {
	  fprintf(stderr,"Error: malloc failed on inverted_file.\n");
	}	
	if( (apq->impact_file = (pqueue_t**)malloc(N*sizeof(pqueue_t*))) == NULL ) {
	  fprintf(stderr,"Error: malloc failed on impact_file.\n");
	}	
	if( (apq->accumulator_file = (GHashTable**)malloc(N*sizeof(GHashTable*))) == NULL ) {
	  fprintf(stderr,"Error: malloc failed on accumulator_file.\n");
	}	
	if( (apq->self_lengths = (GHashTable**)malloc(N*sizeof(GHashTable*))) == NULL ) {
	  fprintf(stderr,"Error: malloc failed on self_lengths.\n");
	}	
	if( (apq->false_lengths = (float*)calloc(N,sizeof(float))) == NULL ) {
	  fprintf(stderr,"Error: malloc failed on false_lengths.\n");
	}
	

	// allocate inverted file arrays and accumulator hashes
	for( i = 0 ; i < maxDim+1; i++ ) {
		apq->inverted_file[i] = g_array_new(FALSE,FALSE,sizeof(IFPAIR));
	}
	
	// fill the inverted file array
	for( i = 0; i < N; i++ ) {
		
		apq->accumulator_file[i] = g_hash_table_new(g_int_hash,g_int_equal);
		apq->self_lengths[i] = g_hash_table_new(g_int_hash,g_int_equal);
		
		apq->impact_file[i] = pqueue_init(100,
		            (pqueue_cmp_pri_f) my_pqueue_cmp_pri_f,
		            (pqueue_get_pri_f) my_pqueue_get_pri_f,
		            (pqueue_set_pri_f) my_pqueue_set_pri_f,
		            (pqueue_get_pos_f) my_pqueue_get_pos_f,
					(pqueue_set_pos_f) my_pqueue_set_pos_f);
		GArray*terms=g_array_index(vecfile->points,GArray*,i);
		for( j = 0; j < (int)(terms->len); j++) {
			IFPAIR ifp = {i,g_array_index(terms,TERM*,j)->val};
			g_array_append_val(apq->inverted_file[g_array_index(terms,TERM*,j)->num],ifp);
		}
	}
	
	// sort into impact order
	for( i = 0 ; i < maxDim+1; i++ ) {
		g_array_sort(apq->inverted_file[i],(GCompareFunc)compare_if);
	}
		
	// build priority queue
	int quadcount = 0;
	for( i = 0; i < N; i++ ) {
		GArray*terms=g_array_index(vecfile->points,GArray*,i);
		for( j = 0; j < (int)(terms->len); j++) {
			apq->quadstore[quadcount].cv = 0;
			apq->quadstore[quadcount].iv = g_array_index(terms,TERM*,j)->num;
			apq->quadstore[quadcount].vv = g_array_index(terms,TERM*,j)->val;
			apq->quadstore[quadcount].pv = g_array_index(apq->inverted_file[apq->quadstore[quadcount].iv],IFPAIR,0).fv * apq->quadstore[quadcount].vv;
			pqueue_insert(apq->impact_file[i], &(apq->quadstore[quadcount]));
			quadcount++;
		}
	}
	
	return apq;
}

// void results_iterator(gpointer key, gpointer value, gpointer user_data) {
// 	static int kk = 0;
// 	printf("%d %d %d\n", *(int*)user_data,kk, ((ACCUM*)value)->pt);
// 	kk = (kk + 1) % k;
// }
// 

/* convert apq neighbors to distance matrix format */

DMAT* apq_2_dmat(APQ* apq, int k) {
	
	int i,j;

	GList** results_list = (GList**)malloc(sizeof(GList*)*apq->N);
	unsigned int valcount = 0;
	for( i = 0; i < apq->N; i++ ) {
		results_list[i] = g_hash_table_get_values(apq->accumulator_file[i]);
		if( k > 0)
			results_list[i] = g_list_sort(results_list[i],(GCompareFunc)compare_acc);
		valcount += g_list_length(results_list[i]);
	}
	
	DMAT* dmat = (DMAT*)malloc(sizeof(DMAT));
	dmat->N = apq->N;
	dmat->rows = g_array_new(FALSE,FALSE,sizeof(GHashTable*));
	if( k > 0 ) {
		dmat->vals = (float*)malloc(apq->N*k*sizeof(float));
		dmat->idxs = (int*)malloc(apq->N*k*sizeof(int));
	}
	else {
		dmat->vals = (float*)malloc(apq->N*valcount*sizeof(float));
		dmat->idxs = (int*)malloc(apq->N*valcount*sizeof(int));
	}	
	
	for( i = 0; i < apq->N; i++ ) {
		GHashTable* row = g_hash_table_new(g_int_hash,g_int_equal);
		g_array_append_val(dmat->rows,row);
	}
	
	valcount=0;
	for( i = 0; i < apq->N; i++ ) {
		GHashTable* row = g_array_index(dmat->rows,GHashTable*,i);
		for( j = 0; ((k>0 && j < k)||(k<0)) && results_list[i] != NULL; j++ ) {
			if( k < 0 ) {
				dmat->idxs[valcount] = ((ACCUM*)(results_list[i]->data))->pt;
				float cosangle = ((ACCUM*)(results_list[i]->data))->acc / (sqrtf(apq->false_lengths[i])*sqrtf(apq->false_lengths[dmat->idxs[valcount]]));
				dmat->vals[valcount] = 1.f-cosangle;
				g_hash_table_insert(row,&(dmat->idxs[valcount]),&(dmat->vals[valcount]));
				valcount++;
			}
			else {
				dmat->idxs[i*k+j] = ((ACCUM*)(results_list[i]->data))->pt;
				float cosangle = ((ACCUM*)(results_list[i]->data))->acc / (sqrtf(apq->false_lengths[i])*sqrtf(apq->false_lengths[dmat->idxs[i*k+j]]));
				dmat->vals[i*k+j] = 1.f-cosangle;
				//GHashTable* col = g_array_index(dmat->rows,GHashTable*,dest);
				g_hash_table_insert(row,&(dmat->idxs[i*k+j]),&(dmat->vals[i*k+j]));
				//g_hash_table_insert(col,&(dmat->idxs[i*k+j]),&(dmat->vals[i*k+j]));
			}
			results_list[i] = g_list_next(results_list[i]);
		}
	}
	return dmat;
}

/* output the nearest neighbor format (in octave sparse format)*/
void output_nn( APQ* apq, int k, const char* filename ) {

	FILE* fp;
	int i,j;

	fp = fopen(filename,"w");
	
	GList** results_list = (GList**)malloc(sizeof(GList*)*apq->N);
	for( i = 0; i < apq->N; i++ ) {
		results_list[i] = g_hash_table_get_values(apq->accumulator_file[i]);
		results_list[i] = g_list_sort(results_list[i],(GCompareFunc)compare_acc);
		for( j = 0; j < k && results_list[i] != NULL; j++ ) {
			fprintf(fp,"%d %d %d\n",i,(j+1),((ACCUM*)(results_list[i]->data))->pt);
			results_list[i] = g_list_next(results_list[i]);
		}
	}
	fclose(fp);
}
