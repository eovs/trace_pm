#ifndef _LIB_H_

#define _LIB_H_

#include "array.h"



typedef struct  
{
	ARRAY array1;
	ARRAY array2;
}TANNER_MON_RES;

typedef struct  
{
	ARRAY array1;
	ARRAY array2;
}HP2A_MON_RES;

typedef struct  
{
	ARRAY ret1;
	ARRAY ret2;
	ARRAY ret3;
	ARRAY ret4;
	int	  ret5;
}MUL_MAT_MAT_MON_RES;

typedef struct 
{
	int degree;
	int *pos;
	int *coef;
	int *ace;
} POLYNOM; 

void trace_bound_pol_mon_pm( ARRAY HD, int M, int gmax, int *SN, int *SA );
TANNER_MON_RES tanner_mon( ARRAY HD );
HP2A_MON_RES hp2a_mon_pm( ARRAY HD, ARRAY W, int M );
MUL_MAT_MAT_MON_RES mul_mat_mat_mon( ARRAY WB, ARRAY DB, ARRAY CB, ARRAY BS, ARRAY A, ARRAY AC, int M);

#endif //_LIB_H_