#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csv.h"

#define MAX_LINE_LEN ((2<<24)-1)
/*
	Load a CSV file to an array of floats
*/
DENSEFILE* load_csv( const char *filename ) {

	char *line = NULL;
	char item[512];

	int num_lines = 0;
	int num_commas = 0;

	line = (char*) malloc( MAX_LINE_LEN * sizeof(char) );
	DENSEFILE* dfile = (DENSEFILE*) malloc( sizeof(DENSEFILE) );
	dfile->M=0;
	
	// open the file 
	
	FILE *fp = fopen( filename, "r" );
	if( fp == NULL ) {
		printf( "ERROR cannot open %s\n", filename );
		exit( 0 );
	}

	// get dataset statistics
	while( fgets( line, MAX_LINE_LEN, fp ) != NULL ) {

		// count the number of points (for every line)
		num_lines++;

		// count the number of dimensions (once)
		if( dfile->M == 0 ) {
			int i = 0;
			while( line[i] != '\0' ) {
				if( line[i] == ',' ) {
					num_commas++;
				}
				i++;
			}
			dfile->M = num_commas+1;
		}
	}
	fclose( fp );
	
	// allocate our data buffer	
	dfile->N = num_lines;
	dfile->data = (float*) malloc( sizeof(float) * (dfile->N) * (dfile->M) );

	// read the data into the buffer
	fp = fopen(filename, "r");
	int k = 0;
	while( fgets( line, MAX_LINE_LEN, fp ) != NULL ) {

		int done = 0;
		int i = 0;
		int j = 0;
		while( !done ) {

			// parse character data
			if( line[i] == ',' ) {
		
				item[j] = '\0';
				dfile->data[k++] = (float) atof( item );
				j = 0;
			}
			else if( line[i] == '\n' || line[i] == '\0' ) {

				item[j] = '\0';
				dfile->data[k++] = (float) atof( item );
				done++;
			}
			else if( line[i] != ' ' ) {

				item[j++] = line[i];
			}
			i++;
		}
	}

	free(line);
	
	return dfile;
}
