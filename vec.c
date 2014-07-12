#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "vec.h"

/* output the distance matrix in vec format */
void output_dmat( DMAT* dmat, const char* filename ) {
	int i,j;
	float w;
	FILE* fp;
	
	fp = fopen(filename,"w");

	for( i = 0; i < dmat->N; i++ ) {
		GHashTable* row = g_array_index(dmat->rows,GHashTable*,i);
		GList* indices = g_hash_table_get_keys(row);
		while(indices != NULL) {
			j = *((int*)(indices->data));
			w = *((float*)g_hash_table_lookup(row,indices->data));
			if( i!=j )
				fprintf(fp, "(%d,%0.4f) ",j+1,w);
			indices = indices->next;
		}
		fprintf(fp,"\n");
	}
	
	fclose(fp);
}

/* Load distance matrix data from the .vec sparse matrix format */
DMAT* load_dmat( const char *filename ) {

	char line[65536];
	char item[256];
	int i = 0, j = 0, k = 0;
	char line_done = 0;
	char paren_done = 0;
	FILE *fp = 0;
	DMAT* dmat = NULL;
	float tempfloat = 0.f;

	// open the file 
	fp = fopen(filename, "r");

	// return if there is an error
	if( fp == NULL ) 
		return NULL;

	// count the number of vals (for every line)
	int tempNNZ = 0;
	int max_nnz = 0;
	int line_num = 0;
	while( fgets( line, 65535, fp ) != NULL ) {
		tempNNZ = 0;
		i = 0;
		while( line[i] != '\0' ) {
			if( line[i] == '(')
				tempNNZ++;
			i++;
		}
		line_num++;
		max_nnz = (tempNNZ > max_nnz)?tempNNZ:max_nnz;
	}
	fclose( fp );

	dmat = (DMAT*)malloc(sizeof(DMAT));
	dmat->N = line_num;
	dmat->vals = (float*)malloc(sizeof(float)*line_num*max_nnz);
	dmat->idxs = (int*)malloc(sizeof(int)*line_num*max_nnz);
	dmat->rows = g_array_new(FALSE,FALSE,sizeof(GHashTable*));

	// read values into data structures
	line_num = 0;
	fp = fopen(filename, "r");
	while( fgets( line, 65535, fp ) != NULL ) {
		
		GHashTable* row = g_hash_table_new(g_int_hash,g_int_equal);
		g_array_append_val(dmat->rows,row);
		k = line_num * max_nnz;
		
		line_done = 0;
		i = 0;
		j = 0;
		while( !line_done ) {
			paren_done = 0;
			if( line[i++] == '(' ) {
				while( !paren_done ) {
					if( line[i] == ',' ) {
						item[j] = '\0';
						dmat->idxs[k] = atoi(item)-1;
						j=0;
					}
					else if( line[i] == ')' ) {
						item[j] = '\0';
						tempfloat = atof(item);
						dmat->vals[k] = tempfloat < 0.f ? 0.f : tempfloat;
						g_hash_table_insert(row,&(dmat->idxs[k]),&(dmat->vals[k]));
						j=0;
						paren_done = 1;
						k++;
					}
					else if( line[i] != '(' && line[i] != ' ' ) {
							item[j++]=line[i];
					}
					i++;
				}
			}
			if( line[i] == '\0' ) {
				line_num++;
				line_done = 1;
			}
		}
	}
	
	return dmat;
}

/*
	Load VECTYPE data from the .vec sparse matrix format
*/
VECFILE* load_vec( const char *filename ) {

	VECFILE* vecfile=NULL;
	char* line = NULL;
	char item[256];
	int i = 0, j = 0, k = 0;
	char line_done = 0;
	char paren_done = 0;
	FILE *fp = 0;
	float tempfloat = 0.f;

	int maxlinelen = (2<<24)-1;
	line = (char*)malloc(sizeof(char)*maxlinelen);
	
	vecfile = (VECFILE*)malloc(sizeof(VECFILE));
	vecfile->max_dim=0;
	vecfile->points = g_array_new(FALSE,FALSE,sizeof(GArray*));
	vecfile->term_count = 0;
	
	// open the file 
	fp = fopen(filename, "r");

	// return if there is an error
	if( fp == NULL ) 
		return NULL;

	// count the number of points (for every line)
	int tempNNZ = 0;
	int line_num = 0;
	while( fgets( line, maxlinelen, fp ) != NULL ) {

		fflush(stdout);
		tempNNZ = 0;
		i = 0;
		while( line[i] != '\0' ) {
			if( line[i] == '(')
				tempNNZ++;
			i++;
		}
		vecfile->max_nnz = (tempNNZ>vecfile->max_nnz)?tempNNZ:vecfile->max_nnz;
		GArray* new_points = g_array_new(FALSE,FALSE,sizeof(TERM*));
		g_array_append_val(vecfile->points,new_points);

		line_done = 0;
		i = 0;
		j = 0;
		while( !line_done ) {
			paren_done = 0;
			if( line[i++] == '(' ) {
				while( !paren_done ) {
					if( line[i] == ',' ) {
						item[j] = '\0';
						vecfile->max_dim = (vecfile->max_dim<atoi(item))?atoi(item):vecfile->max_dim;
						j=0;
					}
					else if( line[i] == ')' ) {
						item[j] = '\0';
						j=0;
						vecfile->term_count++;
						paren_done = 1;
					}
					else if( line[i] != '(' && line[i] != ' ' ) {
							item[j++]=line[i];
					}
					i++;
				}
			}
			if( line[i] == '\0' ) {
				line_num++;
				line_done = 1;
			}
		}
	}
	fclose( fp );

	// allocate data 
	vecfile->terms  = (TERM*)malloc(vecfile->term_count*sizeof(TERM));
	
	// read values into data structures
	line_num = 0;
	fp = fopen(filename, "r");
	while( fgets( line, maxlinelen, fp ) != NULL ) {
		
		line_done = 0;
		i = 0;
		j = 0;
		while( !line_done ) {
			paren_done = 0;
			if( line[i++] == '(' ) {
				while( !paren_done ) {
					if( line[i] == ',' ) {
						item[j] = '\0';
						vecfile->terms[k].num=atoi(item);
						j=0;
					}
					else if( line[i] == ')' ) {
						item[j] = '\0';
						tempfloat = atof(item);
						vecfile->terms[k].val=tempfloat<0.f?0.f:tempfloat;
						TERM* termptr = &(vecfile->terms[k]);
						g_array_append_val(g_array_index(vecfile->points,GArray*,line_num),termptr);
						j=0;
						paren_done = 1;
						k++;
					}
					else if( line[i] != '(' && line[i] != ' ' ) {
							item[j++]=line[i];
					}
					i++;
				}
			}
			if( line[i] == '\0' ) {
				line_num++;
				line_done = 1;
			}
		}
	}
	
	free(line);
	return vecfile;
}
