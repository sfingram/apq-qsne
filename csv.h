#ifndef CSV_H
#define CSV_H

typedef struct _DENSEFILE {
	int N;
	int M;
	float* data;
} DENSEFILE;

DENSEFILE *load_csv( const char *filename );

#endif