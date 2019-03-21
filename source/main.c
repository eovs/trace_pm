#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include "trace_pm.h"

#define MATRIX "HB.txt"
//#define MATRIX "matr_3_4_27_eug_l.txt"
//#define MATRIX "best34_bd3_l.txt"

//#define PERM_BY_HAND
//#define PERM_R2L
//#define PERM_L2R   

#define GMAX 100
#define HDrow 12
#define HDcol 24
#if 0
int HD[HDrow][HDcol] = 
{
/*
	 0,  -1,  -1,  -1,  0,   0, -1,  -1,  0,  -1,  -1,   0,  1,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	22,   0,  -1,  -1, 17,  -1,  0,   0, 12,  -1,  -1,  -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	 6,  -1,   0,  -1, 10,  -1, -1,  -1, 24,  -1,   0,  -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1,
	 2,  -1,  -1,   0, 20,  -1, -1,  -1, 25,   0,  -1,  -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1,
	23,  -1,  -1,  -1,  3,  -1, -1,  -1,  0,  -1,   9,  11, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1,
	24,  -1,  23,   1, 17,  -1,  3,  -1, 10,  -1,  -1,  -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1,
	25,  -1,  -1,  -1,  8,  -1, -1,  -1,  7,  18,  -1,  -1,  0, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1,
	13,  24,  -1,  -1,  0,  -1,  8,  -1,  6,  -1,  -1,  -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1,
	 7,  20,  -1,  16, 22,  10, -1,  -1, 23,  -1,  -1,  -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1,
	11,  -1,  -1,  -1, 19,  -1, -1,  -1, 13,  -1,   3,  17, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1,
	25,  -1,   8,  -1, 23,  18, -1,  14,  9,  -1,  -1,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,
	 3,  -1,  -1,  -1, 16,  -1, -1,   2, 25,   5,  -1,  -1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0
*/

	  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  0,  -1,  -1,  -1,  0,   0, -1,  -1,  0,  -1,  -1,   0,  
	  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 22,   0,  -1,  -1, 17,  -1,  0,   0, 12,  -1,  -1,  -1,  
	 -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1,  6,  -1,   0,  -1, 10,  -1, -1,  -1, 24,  -1,   0,  -1,  
	 -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1,  2,  -1,  -1,   0, 20,  -1, -1,  -1, 25,   0,  -1,  -1,  
	 -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, 23,  -1,  -1,  -1,  3,  -1, -1,  -1,  0,  -1,   9,  11,  
	 -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, 24,  -1,  23,   1, 17,  -1,  3,  -1, 10,  -1,  -1,  -1,  
	 -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1,  0, 25,  -1,  -1,  -1,  8,  -1, -1,  -1,  7,  18,  -1,  -1,  
	 -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, 13,  24,  -1,  -1,  0,  -1,  8,  -1,  6,  -1,  -1,  -1,  
	 -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1,  7,  20,  -1,  16, 22,  10, -1,  -1, 23,  -1,  -1,  -1,  
	 -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, 11,  -1,  -1,  -1, 19,  -1, -1,  -1, 13,  -1,   3,  17,  
	 -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, 25,  -1,   8,  -1, 23,  18, -1,  14,  9,  -1,  -1,  -1,  
	 -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1,  3,  -1,  -1,  -1, 16,  -1, -1,   2, 25,   5,  -1,  -1,  

/*
	 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  0, -1, -1, -1, -1,  0, -1, -1, -1,  0,  0,  0,
	 0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1, -1,  0,  0, -1, 12, 17, 22,
	-1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1,  0, -1, -1, -1, -1, 24, 10,  6,
	-1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1, -1,  0, 25, 20,  2,
	-1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, 11, -1, -1, -1,  9, -1, -1, -1, -1,  0,  3, 23,
	-1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, 23,  1, -1, -1,  3, -1, -1, 10, 17, 24,
	-1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1,  0, -1, -1, -1, -1, -1, -1, -1, -1, 18,  7,  8, 25,
	-1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, 24, -1, -1, -1, -1,  8, -1, -1,  6,  0, 13,
	-1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, 20, -1, 16, -1, 10, -1, -1, -1, 23, 22,  7,
	-1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, 17, -1, -1, -1,  3, -1, -1, -1, -1, 13, 19, 11,
	-1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1,  8, -1, -1, 18, -1, 14, -1,  9, 23, 25,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1, -1, -1, -1, -1, -1, -1, -1,  2,  5, 25, 16,  3,
*/
/*
	// best of BD
	 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1, -1,  0,  0,  0,  0,
	 0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1,  0,  0, -1, -1, -1, 15, 13,  9,
	-1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1, -1, -1,  0, -1, 16,  5,  7,
	-1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1,  2, 23, 10,
	-1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 22, 16,  7, 15,  8,
	-1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, 26, 15, -1, 13, -1, -1, -1, -1, 19, 21,  1,
	-1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, 18, -1, -1, -1, -1, -1, -1,  3, -1, -1, 13, 12, 20,
	-1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1,  7, -1, -1, -1,  9, -1, -1, -1, -1,  0, 13, 11,
	-1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1,  6, -1,  4,  1, -1, -1, -1, -1, -1,  9, 26,  2,
	-1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1,  5, 13, 15,  3,  0,
	-1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1,  1, -1, 12, -1, 25, -1, -1, -1, 23, 16, 20,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1,  4, 22, -1, -1,  3, 15, 16
*/
};
#endif

int column_weight( ARRAY x, int k )
{
	int i;
	int ncol = get_nrow( x );
	int weight = 0;

	for( i = 0; i < ncol; i++ )
		weight += x.addr[i][k] != -1;
	
	return weight;
}

int main( void )
{
	int i, j;
	int M = 32; 
	int gmax = 12; 
	ARRAY matr;
	FILE *fp;
	int nrow, ncol;
	clock_t t0, t1;
	int SN[GMAX];
	int SA[GMAX];
	int *buf;
	int **ptr;

	for( i = 0; i < GMAX; i++ ) SN[i] = 0;	
	for( i = 0; i < GMAX; i++ ) SA[i] = 0;	


	fp = fopen(MATRIX, "rt" );
	if( fp == NULL )
	{
		printf("There is no file\n");
		return 0;
	}

	fscanf_s( fp, "%d", &nrow );
	fscanf_s( fp, "%d", &ncol );

	put_nrow( &matr, nrow );
	put_ncol( &matr, ncol );
	put_addr( &matr, Alloc2d_int( nrow, ncol ) );
	matr.ndim = 2;
	ptr = get_addr( matr );


	for( i = 0; i < nrow; i++ )
		for( j = 0; j < ncol; j++ )
			fscanf_s( fp, "%d", &matr.addr[i][j] );

#ifdef PERM_R2L
	buf = (int*)malloc( ncol * sizeof(buf[0]) );
	
	for( i = 0; i < nrow; i++ )
	{
		int k = 0;
		for( j = ncol-nrow+1; j < ncol; j++ )
			buf[k++] = ptr[i][j];

		buf[k++] = ptr[i][ncol-nrow];

		for( j = 0; j < ncol-nrow; j++ )
			buf[k++] = ptr[i][j];

		for( j = 0; j < ncol; j++ )
			ptr[i][j] = buf[j];
	}

	free( buf );
#endif

#ifdef PERM_L2R
	buf = (int*)malloc( ncol * sizeof(buf[0]) );

	for( i = 0; i < nrow; i++ )
	{
		int k = 0;
		for( j = nrow; j < ncol; j++ )
			buf[k++] = ptr[i][j];

		buf[k++] = ptr[i][nrow-1];

		for( j = 0; j < nrow-1; j++ )
			buf[k++] = ptr[i][j];

		for( j = 0; j < ncol; j++ )
			ptr[i][j] = buf[j];
	}

	free( buf );
#endif



#ifdef PERM_BY_HAND
	while( 1 )
	{
		int dst, src;
		int i, j;
		int ncol = get_ncol( matr );
		int nrow = get_nrow( matr );

		for( i = 0; i < HDcol; i++ )
			printf("%d ", column_weight( matr, i ) );
		printf("\n");

		printf("permutation: ");
		scanf("%d", &src );
		scanf("%d", &dst );
		
		if( src < 0 || src >= ncol )
			break;

		if( dst < 0 || dst >= ncol )
			break;

		if( dst == src )
			break;
		
		for( i = 0; i < nrow; i++ )
		{
			int tmp = matr.addr[i][dst];
			matr.addr[i][dst] = matr.addr[i][src];
			matr.addr[i][src] = tmp;
		}

	}
#endif

	for( i = 0; i < get_nrow( matr ); i++ )
	{
		for( j = 0; j < get_ncol( matr ); j++ )
			matr.addr[i][j] %= M;
	}


	fp = fopen( "matr_res.txt", "wt" );
	fprintf_s( fp, "%d %d\n", nrow, ncol );
	for( i = 0; i < get_nrow( matr ); i++ )
	{
			for( j = 0; j < get_ncol( matr ); j++ )
#if 0
				fprintf( fp, "%3d ", matr.addr[i][j] );
#else
				fprintf( fp, "%3d,", matr.addr[i][j] );
#endif
			fprintf( fp, "\n" );
	}
	fclose( fp );



	t0 = clock();
	trace_bound_pol_mon_pm( matr, M, gmax, SN, SA );
	t1 = clock();

	printf("time: %d msec\n", (t1 - t0) * 1000 / CLOCKS_PER_SEC );
	for( i = 0; i < gmax; i++ )
//		printf("%6.4f  %2d\n", (float)SN[i] / 1000000, SA[i] );
		printf("%8d  %2d\n", SN[i], SA[i] );

	free( matr.addr );
	return 1;
}