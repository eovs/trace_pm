#include <stdio.h>
#include <malloc.h>

#include "lib.h"
#include "array.h"


#define my_min( x, y ) ((x) < (y) ? (x) : (y))
#define my_max( x, y ) ((x) < (y) ? (y) : (x))


void sum_sp_pol_mon
( 
	int wx, int *dx, int *cx, int *xs, 
	int wy, int *dy, int *cy, int *ys,
	int *wz, int *dz, int *cz, int *zs 
);

void copy_2D_to_3d( ARRAY dst, ARRAY src )
{
	int i, j;
	int ***y = get_addr( dst );
	int **x  = get_addr( src );
	int nrow_src = get_nrow( src );
	int ncol_src = get_ncol( src );
	int nrow_dst = get_nrow( dst );
	int ncol_dst = get_ncol( dst );

	if( nrow_src != nrow_dst || ncol_src != ncol_dst )
		return;

	for( i = 0; i < nrow_src; i++ )
		for( j = 0; j < ncol_src; j++ )
			y[0][i][j] = x[i][j];
}

void copy_2D_to_3d_( ARRAY dst, ARRAY src )
{
	int i, j;
	int ***y = get_addr( dst );
	int **x  = get_addr( src );
	int nrow_src = get_nrow( src );
	int ncol_src = get_ncol( src );
	int nrow_dst = get_nrow( dst );
	int nsheet_dst = get_nsheet( dst );

	if( nrow_src != nsheet_dst || ncol_src != nrow_dst )
		return;

	for( i = 0; i < nrow_src; i++ )
		for( j = 0; j < ncol_src; j++ )
			y[i][j][0] = x[i][j];
}

void trace_bound_pol_mon_pm( ARRAY HD, int M, int gmax, int *S, int *SA )
{
	int i, j, k;
	int N;
	int d;
	TANNER_MON_RES tmr;
	HP2A_MON_RES hmr;
	MUL_MAT_MAT_MON_RES mmmmr;
	ARRAY HP, HA;
	ARRAY A, AC;
	ARRAY WB, DB, CB, BS;

	for( i = 0; i < gmax; i++ ) S[i] = 0;	// Spectrum
	for( i = 0; i < gmax; i++ ) SA[i] = 0;	// Spectrum
	
	tmr = tanner_mon( HD );
	HP = tmr.array1;
	HA = tmr.array2;


	hmr = hp2a_mon_pm( HP, HA, M );
	A  = hmr.array1;
	AC = hmr.array2;

	//save_ARRAY( "A.txt", A );
	//save_ARRAY( "AC.txt", AC );
	N = get_ncol( A );

	//	WB=double(A>=0);	//% B  is accumulator, WB are Hamming weights, 2D, integer 
	WB = zeros( 2, N, N );
	{
		int **wb = get_addr( WB );
		int **a  = get_addr( A );
		for( i = 0; i < N; i++ )
			for( j = 0; j < N; j++ )
				wb[i][j] = a[i][j] >= 0;
	}
	
	//	DB=zeros(N,N,2);    //% degrees,3D
	DB = zeros( 3, 2, N, N );
	//	DB(:,:,1)=A;        //% degrees, 3D
	copy_2D_to_3d_( DB, A );
	
	//	CB=zeros(N,N,2);    //% coefficients, 3D
	CB = zeros( 3, 2, N, N );
	
	//	CB(:,:,1)=WB;       //%  coefs, 3D
	copy_2D_to_3d_( CB, WB );

	//	BS=zeros(N,N,2); 
	BS = zeros( 3, 2, N, N );
	//	BS(:,:,1)=AC;       //% Accumullator for ACE, 3D
	copy_2D_to_3d_( BS, AC );


	for( d = 3; d < gmax; d += 2 )
	{
		//[WB,DB,CB,BS]=mul_mat_mat_mon(WB,DB,CB,BS,A,AC,M);
		mmmmr = mul_mat_mat_mon( WB, DB, CB, BS, A, AC, M );
		free_ARRAY( WB );
		free_ARRAY( DB );
		free_ARRAY( CB );
		free_ARRAY( BS );
		WB = mmmmr.ret1;
		DB = mmmmr.ret2;
		CB = mmmmr.ret3;
		BS = mmmmr.ret4;

		for( i = 0; i < N; i++ )
		{
			int ***pCB = get_addr( CB );
			int ***pBS = get_addr( BS );
			int ***pDB = get_addr( DB );
			int **pWB  = get_addr( WB );
			int wb = pWB[i][i];

			if( wb > 0 && pDB[i][i][0] == 0 )
			{
				//% Increase spectra
				S[d] = S[d] + pCB[i][i][0];
				if( SA[d] > 0 )
					SA[d] = my_min( SA[d], pBS[i][i][0] );
				else
					SA[d] = pBS[i][i][0];
				
				//% Subtract constant
				wb = wb-1;
				pWB[i][i] = wb;
				//DB(I,I,1:wb)=DB(I,I,2:wb+1);
				for( k = 0; k < wb; k++ )
					pDB[i][i][k] = pDB[i][i][k+1];

				//CB(I,I,1:wb)=CB(I,I,2:wb+1);            
				for( k = 0; k < wb; k++ )
					pCB[i][i][k] = pCB[i][i][k+1];
				//%BS(I,I,1:wb)=BS(I,I,2:wb+1);
				//for( k = 0; k < wb; k++ )
				//	pBS[i][i][k] = pBS[i][i][k+1];
			}
		}
		//if( S[d] > 0 )	break;
	}

	for( i = 3; i < gmax; i += 2 )
		S[i] = S[i] / (i+1);
	
	for( i = 3; i < gmax; i += 2 )
		SA[i] /= 2;

	free_ARRAY( WB );
	free_ARRAY( BS );
	free_ARRAY( CB );
	free_ARRAY( DB );
	free_ARRAY( HP );
	free_ARRAY( HA );
	free_ARRAY( A );
	free_ARRAY( AC );
}

TANNER_MON_RES tanner_mon( ARRAY hd )
{
	int cw[1000];
	int i, j;
	int num;
	int IC;
	int **T, **W;
	TANNER_MON_RES res;
	int **HD = get_addr( hd );
	int nrow = get_nrow( hd );
	int ncol = get_ncol( hd );
	
	//cw - column weights
	num = 0;
	for( i = 0; i < ncol; i++ )
	{
		int s = 0;
		for( j = 0; j < nrow; j++ )
			s += HD[j][i] >= 0;
		cw[i] = s;
		num += cw[i]; 
	}

	res.array1 = minus_ones( 2, num, nrow+ncol, num );
	T = get_addr( res.array1 ); 
	
	res.array2 = zeros( 2, num, nrow+ncol );
	W = get_addr( res.array2 );


	IC = 0; // current column
	for( i = 0; i < ncol; i++ )
	{
		for( j = 0; j < nrow; j++ )
		{
			if( HD[j][i] >= 0 )
			{
				T[i][IC] = 1;
				T[ncol + j][IC] = HD[j][i];
				W[ncol + j][IC] = cw[i];

				IC=IC+1;
			}
		}
	}


	return res;
}


HP2A_MON_RES hp2a_mon_pm( ARRAY hd, ARRAY w, int M )
{
	HP2A_MON_RES res;
	int N    = get_ncol( hd );
	int nrow = get_nrow( hd );
	int **HD = get_addr( hd );
	int **W  = get_addr( w );
	int i, j, h;
	int I, J;
	int num;
	ARRAY a, wa, am, wam, ap, wap; 
	ARRAY e, x, s;
	int **A, **WA;
	int **Am, **WAm;
	int **Ap, **WAp;
	int **E, **X, **S;
	int **Apm, **WApm;
	int *f = (int*)malloc( N * sizeof(int) );
	
	a = minus_ones( 2, N*2, N*2 );	A = get_addr( a );
	wa = zeros( 2, N*2, N*2 );  	WA = get_addr( wa );
	
	// %minus		
	am = minus_ones( 2, N, N ); Am = get_addr( am );
	wam = zeros( 2, N, N );	WAm = get_addr( wam );

	// % plus
	ap = minus_ones( 2, N, N );	Ap = get_addr( ap );
	wap = zeros( 2, N, N );	WAp = get_addr( wap );


	//% Edges: loop over columns of Tanner
	e = zeros( 2, 2, 2*N );		E = get_addr( e );
	x = zeros( 2, 1, 2*N );		X = get_addr( x );
	//% column weights
	s = zeros( 2, 1, 2*N );		S = get_addr( s );
	


	for( i = 0; i < N; i++ )
	{
		for( j = 0, num = 0; j < nrow; j++ ){if( HD[j][i] >= 0 ){f[num++] = j;}}

		E[i][0] = f[0];   //% forward  part 1 to part 2
		E[i][1] = f[1]; 
		E[i+N][0] = f[1]; //% back  from 2 to 1;
		E[i+N][1] = f[0]; 
		X[i][0]   = HD[f[1]][i];
		X[i+N][0] = (M - X[i][0]) % M;
		S[i][0]   = W[f[1]][i];
		S[i+N][0] = W[f[1]][i];
	}

	//% Edge adjacency matrix (Aminus)
	for( i = 0; i < N; i++ )          //% i is head of the current edge
	{
		j = E[i][1];        //% j is tail of the current edge
		//% find all edges starting at j ending anywhere except i+N or i-N 
		//f=find( E(N+1:2*N,1)==j );
		for( h = N, num = 0; h < 2*N; h++ ){if( E[h][0] == j ){f[num++] = h;}}

		for( h = 0; h < num; h++ )
		{
			if( f[h] != i+N )
			{
				A[i][f[h]]     = X[f[h]][0];
				Am[i][f[h]-N]  = X[f[h]][0];
				WA[i][f[h]]    = S[f[h]][0];
				WAm[i][f[h]-N] = S[f[h]][0];
			}
		}
	}

	//%    (Aplus)
	for( i = N; i < 2*N; i++ ) //          % i is head of the current edge
	{
		j = E[i][1];        //% j is tail of the current edge
		//% find all edges starting at j ending anywhere except i+N or i-N 
		//% if i<=N, ib=i+N; else ib=i-N; end
		//f=find(E(1:N,1)==j); 
		for( h = 0, num = 0; h < N; h++ ){if( E[h][0] == j ){f[num++] = h;}}

		for( h = 0; h < num; h++ )
		{
			if( f[h] != i-N )
			{
				A[i][f[h]]     = X[f[h]][0];
				Ap[i-N][f[h]]  = X[f[h]][0];
				WA[i][f[h]]    = S[f[h]][0];
				WAp[i-N][f[h]] = S[f[h]][0];
			}
		}
	}


	//% Multiplication Am*Ap
	res.array1 = minus_ones( 2, N, N ); 	Apm  = get_addr( res.array1 );
	res.array2 = zeros( 2, N, N );			WApm = get_addr( res.array2 );
	
	for( I = 0; I < N; I++ )  //% rows Am
	{
		for( J = 0; J < N; J++ ) //% columns Ap
		{
			//%s=-1; ws=0;  //% accumulator
			for( j = 0; j < N; j++ ) // % along row
			{
				if( Am[I][j] >= 0 && Ap[j][J] >= 0 )
				{
					Apm[I][J]  =  (Am[I][j] +  Ap[j][J]) % M; 
					WApm[I][J] = WAm[I][j] + WAp[j][J];
				}
			}
		}
	}

	free_ARRAY( s );
	free_ARRAY( x );
	free_ARRAY( e );

	free_ARRAY( ap );
	free_ARRAY( wap );

	free_ARRAY( am );
	free_ARRAY( wam );

	free_ARRAY( a );
	free_ARRAY( wa );

	free( f );
	return res;
}


int find_first( int *x, int n, int M )
{
	int i;
	int pos = -1;
	for( i = 0; i < n; i++ )
	{
		if( x[i] >= M )
		{
			pos = i;
			break;
		}
	}
	return pos;
}


void sum_sp_pol_mon
( 
	int wx, int *dx, int *cx, int *xs, 
	int wy, int *dy, int *cy, int *ys,
	int *wz, int *dz, int *cz, int *zs 
)
{
	int cntX, cntY, cntZ;
	
	*wz = 0;
	cntX = 0;
	cntY = 0;
	cntZ   = 0;

	while( 1 )
	{
		if( cntX == wx )
		{
			// dx - done
			while( cntY < wy )
			{
				dz[cntZ] = dy[cntY];
				cz[cntZ] = cy[cntY];
				zs[cntZ] = ys[cntY];
				cntY++;
				cntZ++;
			}
			*wz = cntZ;
			return; 
		}
		else
		{
			if( cntY == wy )
			{
				// dy - done;
				while( cntX < wx )
				{
					dz[cntZ] = dx[cntX];
					cz[cntZ] = cx[cntX];
					zs[cntZ] = xs[cntX];
					cntX++;
					cntZ++;
				}
				*wz = cntZ;
				return; 
			}
			else
			{
				if( dx[cntX] < dy[cntY] )
				{
					dz[cntZ] = dx[cntX];
					cz[cntZ] = cx[cntX];
					zs[cntZ] = xs[cntX];
					cntX++;
					cntZ++;
				}
				else
				{
					if( dx[cntX] > dy[cntY] )
					{
						dz[cntZ] = dy[cntY];
						cz[cntZ] = cy[cntY];
						zs[cntZ] = ys[cntY];
						cntY++;
						cntZ++;
					}
					else
					{
						int t = cx[cntX] + cy[cntY];
						if( t != 0 )
						{
							dz[cntZ] = dx[cntX];
							cz[cntZ] = t;
							zs[cntZ] = my_min( xs[cntX], ys[cntY] );
							cntX++;
							cntY++;
							cntZ++;
						}
					}
				}
			}

		}
	}
}

MUL_MAT_MAT_MON_RES mul_mat_mat_mon( ARRAY WB, ARRAY DB, ARRAY CB, ARRAY BS, ARRAY A, ARRAY AS, int M)
{
	int bb, LB, LX;
	ARRAY WX, CX, DX, XS;
	int I, J;
	int w;
	MUL_MAT_MAT_MON_RES res;

	int **pWB = get_addr( WB );
	int **pA = get_addr( A );
	int **pAS = get_addr( AS );
	int ***pDB = get_addr( DB );
	int ***pCB = get_addr( CB );
	int ***pBS = get_addr( BS );
	int **pWX;
	int ***pDX;
	int ***pCX;
	int ***pXS;
	int wx, rwx;
	int* db;
	int* cb;
	int* bs;
	int* dx;// = malloc(M*sizeof(int));
	int* cx;// = malloc(M*sizeof(int));
	int* xs;// = malloc(M*sizeof(int));
	int* rdx = malloc(M*sizeof(int));
	int* rcx = malloc(M*sizeof(int));
	int* rxs = malloc(M*sizeof(int));
	int tpos;
	int h, i, j;
	int reswb;
	int *resdb = malloc(M*sizeof(int));
	int *rescb = malloc(M*sizeof(int));
	int *resbs = malloc(M*sizeof(int));

	//% A,B are degree matrices
	//	% B is 3D
	//	% A is 2D
	//	% M is modulo

	//	[bb,~,LB]=size(DB);
	bb = get_nrow( DB );
//	LB = get_nsheet( DB );
	LB = get_ncol( DB );

	LX=LB+1;
	//db=zeros(1,LB); % Column degrees 
	
	db = calloc( LB, sizeof(db[0]) );	
	//cb=zeros(1,LB); % coefs
	cb = calloc( LB, sizeof(cb[0]) );
	//bs=zeros(1,LB); % ace
	bs = calloc( LB, sizeof(bs[0]) );
 
	//% ARRAY for result
	//WX=zeros(bb,bb);    % weights
	WX = zeros( 2, bb, bb );	
	//DX=zeros(bb,bb,LX); % degrees
	DX = zeros( 3, M, bb, bb );
	//CX=zeros(bb,bb,LX); % coefs
	CX = zeros( 3, M, bb, bb );
	//XS=zeros(bb,bb,LX); % ACEs
	XS = zeros( 3, M, bb, bb );

	//% main loop
	w=0;  //% max weight

	pWX = get_addr( WX );
	pDX = get_addr( DX );
	pCX = get_addr( CX );
	pXS = get_addr( XS );

	for( I=0; I < bb; I++ ) //% rows of B
	{
		for( J=0; J < bb; J++ ) //% columns of A (square marices)
		{
			dx = &pDX[I][J][0]; //% degrees
			cx = &pCX[I][J][0]; //% coefs
			xs = &pXS[I][J][0]; //% ACEs

			wx = 0; 

			for( j = 0; j < bb; j++ ) //% loop along column
			{
				int a = pA[j][J], as = pAS[j][J];
				int wb = pWB[I][j];

				if( wb > 0 && a >= 0 )
				{
	
					for( h = 0; h < wb; h++ )  //% read from 3D arrays & multiply
					{
						db[h] = pDB[I][j][h] + a;
						cb[h] = pCB[I][j][h];
						bs[h] = pBS[I][j][h] + as;
					}

/*			
					//% Multiplication 
					//bs(1:wb)=bs(1:wb)+as;  % new col weights
					for( h = 0; h < wb; h++ ) 	bs[h] += as;
					//db(1:wb)=db(1:wb)+a;   % new degrees
					for( h = 0; h < wb; h++ )	db[h] += a;
*/

					//t=find(db>=M,1);
					tpos = find_first( db, wb, M );

					//if( ~isempty(t) )
					if( tpos >= 0 )
					{
						for( h = 0, i = tpos; i < wb; i++, h++ ) resdb[h] = db[i] - M;
						for( i = 0; i < tpos; i++, h++ ) resdb[h] = db[i];
						for( i = 0; i < wb; i++ ) db[i] = resdb[i];

						for( h = 0, i = tpos; i < wb; i++, h++ ) rescb[h] = cb[i];
						for( i = 0; i < tpos; i++, h++ ) rescb[h] = cb[i];
						for( i = 0; i < wb; i++ ) cb[i] = rescb[i];

						for( h = 0, i = tpos; i < wb; i++, h++ ) resbs[h] = bs[i];
						for( i = 0; i < tpos; i++, h++ ) resbs[h] = bs[i];
						for( i = 0; i < wb; i++ ) bs[i] = resbs[i];
					}
	

					//% Accumulation
					//[wx,dx,cx,xs]=sum_sp_pol_mon( wx,dx,cx,xs, wb,db,cb,bs );
					sum_sp_pol_mon( wx, dx, cx, xs, wb, db, cb, bs, &rwx, rdx, rcx, rxs );
					wx = rwx;
					for( h = 0; h < wx; h++ )	dx[h] = rdx[h];
					for( h = 0; h < wx; h++ )	cx[h] = rcx[h];
					for( h = 0; h < wx; h++ )	xs[h] = rxs[h];
					h = 0;
				}
			}


			pWX[I][J] = wx;    //% weights
			w = my_max( w, wx );  //% max weight
			
		}
	}

	if( w < my_max(LX,2) )
		LX = w;
	
	free( db );
	free( cb );
	free( bs );

	//free( dx );
	//free( cx );
	//free( xs );

	free( rdx );
	free( rcx );
	free( rxs );

	free( resdb );
	free( rescb );
	free( resbs );

	res.ret1 = WX;
	res.ret2 = DX;
	res.ret3 = CX;
	res.ret4 = XS;
	res.ret5 = LX;
	return res;
}

