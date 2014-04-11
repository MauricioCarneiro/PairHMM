#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "input.h"
//#include "read.h"

#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5
#define MAX(a,b) ((a)>(b)?(a):(b))
#define WARP 32 

template<class T>
struct Context{};

template<>
struct Context<double>
{
	Context()
	{
		for (int x = 0; x < 128; x++)
			ph2pr[x] = pow(10.0, -((double)x) / 10.0);

		INITIAL_CONSTANT = ldexp(1.0, 1020.0);
		LOG10_INITIAL_CONSTANT = log10(INITIAL_CONSTANT);
		RESULT_THRESHOLD = 0.0;
	}

	double LOG10(double v){ return log10(v); }

	static double _(double n){ return n; }
	static double _(float n){ return ((double) n); }
	double ph2pr[128];
	double INITIAL_CONSTANT;
	double LOG10_INITIAL_CONSTANT;
	double RESULT_THRESHOLD;
};

template<>
struct Context<float>
{
	Context()
	{
		for (int x = 0; x < 128; x++)
			ph2pr[x] = powf(10.f, -((float)x) / 10.f);

		INITIAL_CONSTANT = ldexpf(1.f, 120.f);
		LOG10_INITIAL_CONSTANT = log10f(INITIAL_CONSTANT);
		RESULT_THRESHOLD = ldexpf(1.f, -110.f);
	}

	float LOG10(float v){ return log10f(v); }

	static float _(double n){ return ((float) n); }
	static float _(float n){ return n; }
	float ph2pr[128];
	float INITIAL_CONSTANT;
	float LOG10_INITIAL_CONSTANT;
	float RESULT_THRESHOLD;
};

__device__ double __shfl_up(double d,unsigned int i){
   double ret;
   ((float*)&ret)[0] = __shfl_up(((float*)&d)[0], i);
   ((float*)&ret)[1] = __shfl_up(((float*)&d)[1], i);
   return ret;
}

__device__ double __shfl_down(double d,unsigned int i){
   double ret;
   ((float*)&ret)[0] = __shfl_down(((float*)&d)[0], i);
   ((float*)&ret)[1] = __shfl_down(((float*)&d)[1], i);
   return ret;
}

template<class NUMBER>
__global__ void 
__launch_bounds__(128,8)
pairhmm_kernel( NUMBER Yr0, NUMBER* M, NUMBER *X, NUMBER *Y, 
                               NUMBER* p, char* rs, char* hap, 
                               int* q, double* ph2pr, int* offset, int n_mats,
                               NUMBER* output, NUMBER log10_init) {
   NUMBER M_p, M_pp, X_p, X_pp, Y_p, Y_pp, distm, pMM, pGapM, pXX, 
          pMX, pYY, pMY, M_loc, X_loc, Y_loc;
   char _rs;
   int _q;
   int tid = threadIdx.x;
   int wid = blockIdx.x;
   //Fake wid,tid to allow for larger blocks
   tid = (threadIdx.x + blockDim.x * blockIdx.x)%WARP;
   wid = (threadIdx.x + blockDim.x * blockIdx.x)/WARP;
   if (wid > n_mats) return;
   int ROWS = offset[3*wid+4]-offset[3*wid+1];
   int COLS = offset[3*wid+5]-offset[3*wid+2];
   NUMBER M_top, X_top, Y_top;
   NUMBER result=0.0;
   M+=offset[3*wid];
   X+=offset[3*wid];
   Y+=offset[3*wid];
   p+=6*offset[3*wid+1];
   rs+=offset[3*wid+1];
   hap+=offset[3*wid+2];
   q+=offset[3*wid+1];
   ph2pr+=128*wid;
   for (int stripe = 0; stripe < ROWS; stripe+=WARP-1) {
      if ( stripe==0 && tid < 2) {
         M_pp=0.0;
         X_pp = 0.0;//X_pp=Xc0[0];
         Y_pp=Yr0;
         M_p=0.0;
         if (tid==1) X_p = 0.0;//X_p=Xc0[1];
         else X_p=0.0;
         if (tid==0) Y_p=Yr0;
         else Y_p=0.0;
      } else if (tid<2) {
         M_pp = 0.0;
         X_pp = 0.0; //X_pp = Xc0[stripe]; //X[stripe][0];
         Y_pp = 0.0;
         M_p = 0.0;
         if (tid==1) X_p=0.0; //X_p=Xc0[1];
         else X_p = X[stripe/(WARP-1)*COLS+1]; //M[stripe+tid][1-tid];
         if (tid==0) Y_p = Y[stripe/(WARP-1)*COLS+1]; //M[stripe+tid][1-tid];
         else Y_p = 0.0;
      }
      //TODO pad instead
      if (tid>0) {
      _rs = rs[tid-1+stripe];
	   _q = q[tid-1+stripe] & 127;
      } else {
      _rs = rs[tid+stripe];
	   _q = q[tid+stripe] & 127;
      }
      //TODO transpose p for coalesced reads
      pMM = p[6*(tid+stripe)+MM];
      pGapM = p[6*(tid+stripe)+GapM];
      pXX = p[6*(tid+stripe)+XX];
      pMX = p[6*(tid+stripe)+MX];
      pYY = p[6*(tid+stripe)+YY];
      pMY = p[6*(tid+stripe)+MY];
      for (int z = 1; z < COLS+WARP+1;z++) 
      {
         int r = tid+stripe;
         int c = z-tid+1;
         //TODO align at word boundaries
         if (1==z%WARP) {
            M_top = M[stripe/(WARP-1)*COLS+z+tid];
            X_top = X[stripe/(WARP-1)*COLS+z+tid];
            Y_top = Y[stripe/(WARP-1)*COLS+z+tid];
         } else {
            M_top = __shfl_down(M_top,1);
            X_top = __shfl_down(X_top,1);
            Y_top = __shfl_down(Y_top,1);
         }

         if (tid<= z+1 && tid+stripe < ROWS && z-tid < COLS)
         {
            //TODO pad instead
            if (c>0) {
            char _hap = hap[c-1];
			   if (_rs == _hap || _rs == 'N' || _hap == 'N')
			   	distm = double(1.0) - ph2pr[_q];
            else distm = ph2pr[_q];
            }
            else distm = ph2pr[_q];
            if (tid == 0 && stripe == 0) {
               X_p = 0.0; 
               Y_p = Yr0;
               M_p = 0.0;
            } else if (tid == 0 && z > 1) {
               M_p = M_top;
               X_p = X_top;
               Y_p = Y_top;
            } 
            M_loc = distm * (M_pp * pMM + X_pp * pGapM + Y_pp * pGapM);
			   Y_loc = M_p * pMY + Y_p * pYY;
            M_p = __shfl_up(M_p,1);
            Y_p = __shfl_up(Y_p,1);
            X_p = __shfl_up(X_p,1);
			   X_loc = M_p * pMX + X_p * pXX;
            M_pp = M_p;
            X_pp = X_p;
            Y_pp = Y_p;
            if (tid == z+1 && stripe==0) {
               M_p = 0.0;
               Y_p = 0.0;
               X_p = 0.0; //X_p = Xc0[tid];
            } else if (tid == z+1) {
               M_p = 0.0;
               Y_p = 0.0;
               X_p = 0.0; //X_p = Xc0[tid+stripe]; //X[tid+stripe][0]
            } else {
               M_p = M_loc;
               X_p = X_loc;
               Y_p = Y_loc;
            }
            r = tid+stripe;
            c =  z-tid+1;
            //TODO shuffle M_out to write one time
            if (tid>0 && c < COLS && (r==stripe+WARP-1 || r==ROWS-1)) {
               M[((r+WARP-2)/(WARP-1))*COLS+c] = M_loc;
               X[((r+WARP-2)/(WARP-1))*COLS+c] = X_loc;
               Y[((r+WARP-2)/(WARP-1))*COLS+c] = Y_loc;
            }
            if (r==ROWS-1) result += M_loc + X_loc;
         }
#if 0
         NUMBER M_bottom, X_bottom, Y_bottom;
         int write_tid=min(WARP-1, ROWS-stripe-1);
         //shuffle from write_tid to the last thread in the warp
         M_bottom = __shfl_down(M_bottom, 1);
         X_bottom = __shfl_down(X_bottom, 1);
         Y_bottom = __shfl_down(Y_bottom, 1);
         M_loc = __shfl_up(M_loc, WARP-1-write_tid); 
         X_loc = __shfl_up(X_loc, WARP-1-write_tid); 
         Y_loc = __shfl_up(Y_loc, WARP-1-write_tid); 
         if (tid == write_tid) {
            M_bottom = M_loc;
            X_bottom = X_loc;
            Y_bottom = Y_loc;
         }
         if (WARP-1 == (z-write_tid)%WARP || COLS-1 == z-write_tid) {
            r = stripe + write_tid;
            c = z - write_tid + 1;
            M[((r+WARP-2)/(WARP-1))*COLS+c-write_tid+tid] = M_bottom;
            X[((r+WARP-2)/(WARP-1))*COLS+c-write_tid+tid] = X_bottom;
            Y[((r+WARP-2)/(WARP-1))*COLS+c-write_tid+tid] = Y_bottom;
         }
#endif
         
#if 0
         if () {
            M[((r+WARP-2)/(WARP-1))*COLS+c] = M_bottom;
            X[((r+WARP-2)/(WARP-1))*COLS+c] = M_bottom;
            Y[((r+WARP-2)/(WARP-1))*COLS+c] = M_bottom;
         } else {
            M_bottom = __shfl_down(M_bottom,1);
            X_bottom = __shfl_down(X_bottom,1);
            Y_bottom = __shfl_down(Y_bottom,1);
         }
#endif
      }
   }
   if (tid == (ROWS-1)%(WARP-1)) output[wid] = log10(result) - log10_init; 
   //if (tid == (ROWS-1)%(WARP-1)) output[wid] = result;
}
template<class NUMBER>
struct GPUmem {
   int offset[1000][3];
   int* d_offset;
   int index;
   NUMBER* d_M;
   NUMBER* d_X;
   NUMBER* d_Y;
   NUMBER* d_p;
   NUMBER* d_Yr0;
   NUMBER* d_Xc0;
   NUMBER* M;
   NUMBER* X;
   NUMBER* Y;
   NUMBER* p;
   NUMBER* Xc0;
   NUMBER* Yr0;
   char* rs;
   char* hap;
   int* q;
   double* ph2pr;
   char* d_rs;
   char* d_hap;
   int* d_q; 
   double *d_ph2pr;
   GPUmem() {M=0;};
}; 
template<class NUMBER>
int tc2gmem(GPUmem<NUMBER>& gmem, testcase* tc)
{
	int ROWS = tc->rslen + 1;
	int COLS = tc->haplen + 1;

	Context<NUMBER> ctx;

	NUMBER M[ROWS][COLS];
	//NUMBER M2[ROWS][COLS];
	NUMBER X[ROWS][COLS];
	//NUMBER X2[ROWS][COLS];
   NUMBER Xc0[ROWS];
	NUMBER Y[ROWS][COLS];
   NUMBER Yr0[COLS];
	NUMBER p[ROWS][6];
   if (ROWS==1 || COLS==1) return 0;
   extract_tc(&M[0][0],&X[0][0],&Y[0][0],&p[0][0],Xc0,Yr0,tc);
   //TODO check data sizes first
   //      return error if we're out of space
   memcpy(gmem.Xc0+gmem.offset[gmem.index][1],Xc0, sizeof(NUMBER)*ROWS);
   memcpy(gmem.Yr0+gmem.offset[gmem.index][2], Yr0, sizeof(NUMBER)*COLS);
   memcpy(gmem.p+gmem.offset[gmem.index][1]*6, &p[0][0], sizeof(NUMBER)*ROWS*6);
   memcpy(gmem.rs+gmem.offset[gmem.index][1], tc->rs, sizeof(char)*ROWS);
   memcpy(gmem.hap+gmem.offset[gmem.index][2], tc->hap, sizeof(char)*COLS);
   memcpy(gmem.q+gmem.offset[gmem.index][1], tc->q, sizeof(int)*ROWS);
   memcpy(gmem.ph2pr+128*gmem.index, ctx.ph2pr, sizeof(double)*128);
   //gmem.offset[gmem.index+1][0] = gmem.offset[gmem.index][0] + ROWS*COLS;
   gmem.offset[gmem.index+1][0] = gmem.offset[gmem.index][0] + ((ROWS+WARP-2)/(WARP-1))*COLS;
   //printf("ROWS/COLS = %d/%d, offset[%d] = %d\n", ROWS, COLS, gmem.index+1, gmem.offset[gmem.index+1][0]);
   gmem.offset[gmem.index+1][1] = gmem.offset[gmem.index][1] + ROWS;
   gmem.offset[gmem.index+1][2] = gmem.offset[gmem.index][2] + COLS;
   gmem.index++;
   return 0;
}
template<class NUMBER>
int GPUmemAlloc(GPUmem<NUMBER>& gmem) 
{
   cudaDeviceProp deviceProp;
   cudaError_t err = cudaGetDeviceProperties(&deviceProp, 0);
   unsigned long long totalMem = deviceProp.totalGlobalMem/5;
   //TODO 
   //allocations based (very loosely) on ROWS=COLS=10
   //TODO readjust sizes
   cudaMalloc(&gmem.d_M, totalMem/10);
   gmem.M = (NUMBER*)malloc(totalMem/10);
   cudaMalloc(&gmem.d_X, totalMem/10);
   gmem.X = (NUMBER*)malloc(totalMem/10);
   cudaMalloc(&gmem.d_Y, totalMem/10);
   gmem.Y = (NUMBER*)malloc(totalMem/10);
   cudaMalloc(&gmem.d_p, totalMem/11);
   gmem.p = (NUMBER*)malloc(totalMem/11);
   cudaMalloc(&gmem.d_Yr0, totalMem/66);
   gmem.Yr0 = (NUMBER*)malloc(totalMem/66);
   cudaMalloc(&gmem.d_Xc0, totalMem/66);
   gmem.Xc0 = (NUMBER*)malloc(totalMem/66);
   cudaMalloc(&gmem.d_q, totalMem/66);
   gmem.q = (int*)malloc(totalMem/66);
   cudaMalloc(&gmem.d_rs, 700000);
   gmem.rs = (char*)malloc(700000);
   cudaMalloc(&gmem.d_hap, 700000);
   gmem.hap = (char*)malloc(700000);
   cudaMalloc(&gmem.d_ph2pr, totalMem/4);
   gmem.ph2pr = (double*)malloc(totalMem/4);
   err = cudaGetLastError();
   if (err) printf("cudaMalloc error %d: %s\n", err, cudaGetErrorString(err));
   if (err) return 9000+err;
   if (!gmem.M ||
       !gmem.X ||
       !gmem.Y ||
       !gmem.p ||
       !gmem.q ||
       !gmem.rs ||
       !gmem.Yr0 ||
       !gmem.Xc0 ||
       !gmem.hap ||
       !gmem.ph2pr) {
      printf("CPU mem allocation fail\n");
      return 1;
   }
   return 0;
}
template<class NUMBER>
int GPUmemFree(GPUmem<NUMBER>& gmem) 
{
   if (0==gmem.M) {
      return 0;
   }
   gmem.index=0;
   cudaFree(gmem.d_M);
   free(gmem.M);
   gmem.M=0;
   cudaFree(gmem.d_X);
   free(gmem.X);
   gmem.X=0;
   cudaFree(gmem.d_Y);
   free(gmem.Y);
   gmem.Y=0;
   cudaFree(gmem.d_p);
   free(gmem.p);
   gmem.p=0;
   cudaFree(gmem.d_Yr0);
   free(gmem.Yr0);
   gmem.Yr0=0;
   cudaFree(gmem.d_Xc0);
   free(gmem.Xc0);
   gmem.Xc0=0;
   cudaFree(gmem.d_q);
   free(gmem.q);
   gmem.q=0;
   cudaFree(gmem.d_rs);
   free(gmem.rs);
   gmem.rs=0;
   cudaFree(gmem.d_hap);
   free(gmem.hap);
   gmem.hap=0;
   cudaFree(gmem.d_ph2pr);
   free(gmem.ph2pr);
   gmem.ph2pr=0;
   return 0;
}
template<class NUMBER>
void extract_tc(NUMBER* M_in, NUMBER* X_in, NUMBER* Y_in, 
                NUMBER* p_in, NUMBER* Xc0, NUMBER* Yr0, testcase* tc)
{
	int ROWS = tc->rslen + 1;
	int COLS = tc->haplen + 1;
   int r,c;

   NUMBER (*M)[COLS] = (NUMBER (*)[COLS]) &M_in[0];
   NUMBER (*Y)[COLS] = (NUMBER (*)[COLS]) &Y_in[0];
   NUMBER (*X)[COLS] = (NUMBER (*)[COLS]) &X_in[0];
   NUMBER (*p)[6] = (NUMBER (*)[6]) &p_in[0];

	Context<NUMBER> ctx;
	p[0][MM] = ctx._(0.0);
	p[0][GapM] = ctx._(0.0);
	p[0][MX] = ctx._(0.0);
	p[0][XX] = ctx._(0.0);
	p[0][MY] = ctx._(0.0);
	p[0][YY] = ctx._(0.0);
	for (r = 1; r < ROWS; r++)
	{
		int _i = tc->i[r-1] & 127;
		int _d = tc->d[r-1] & 127;
		int _c = tc->c[r-1] & 127;
		p[r][MM] = ctx._(1.0) - ctx.ph2pr[(_i + _d) & 127];
		p[r][GapM] = ctx._(1.0) - ctx.ph2pr[_c];
		p[r][MX] = ctx.ph2pr[_i];
		p[r][XX] = ctx.ph2pr[_c];
		p[r][MY] = (r == ROWS - 1) ? ctx._(1.0) : ctx.ph2pr[_d];
		p[r][YY] = (r == ROWS - 1) ? ctx._(1.0) : ctx.ph2pr[_c];
	}

	for (c = 0; c < COLS; c++)
	{
		M[0][c] = ctx._(0.0);
		X[0][c] = ctx._(0.0);
		Y[0][c] = ctx.INITIAL_CONSTANT / (tc->haplen);
		Yr0[c] = ctx.INITIAL_CONSTANT / (tc->haplen);
	}

   Xc0[0]=X[0][0];
	for (r = 1; r < ROWS; r++)
	{
		M[r][0] = ctx._(0.0);
		X[r][0] = X[r-1][0] * p[r][XX];
		Y[r][0] = ctx._(0.0);
		Xc0[r] = Xc0[r-1] * p[r][XX];
	}
}

template<class NUMBER>
void compute_full_prob_multiple(NUMBER* probs, testcase *tc, int n_tc, 
                                 GPUmem<NUMBER> &gmem, NUMBER *before_last_log = NULL) {
   Context<NUMBER> ctx;
   int err;
   cudaError_t cuerr;
   NUMBER *d_out;

   if (0==n_tc) {
      err = GPUmemFree<NUMBER>(gmem);
      return;
   }
   if (gmem.M==0) {
      err = GPUmemAlloc<NUMBER>(gmem);
   }
   gmem.index=0;
   for (int z=0;z<n_tc;z++)
   {
      err = tc2gmem<NUMBER>(gmem, &tc[z]);
   }
   cudaMalloc(&gmem.d_offset, sizeof(int)*3*(n_tc+1));
   cudaMalloc(&d_out, sizeof(NUMBER)*n_tc);
   cudaMemcpy(gmem.d_offset, &gmem.offset[0][0], sizeof(int)*3*(n_tc+1), cudaMemcpyHostToDevice);
   //cudaMemcpy(gmem.d_Xc0, gmem.Xc0, sizeof(NUMBER)*gmem.offset[n_tc][1], cudaMemcpyHostToDevice);
   cudaMemcpy(gmem.d_Yr0, gmem.Yr0, sizeof(NUMBER)*gmem.offset[n_tc][2], cudaMemcpyHostToDevice);
   cudaMemcpy(gmem.d_p, gmem.p, sizeof(NUMBER)*gmem.offset[n_tc][1]*6, cudaMemcpyHostToDevice);
   cudaMemcpy(gmem.d_rs, gmem.rs, sizeof(char)*gmem.offset[n_tc][1], cudaMemcpyHostToDevice);
   cudaMemcpy(gmem.d_hap, gmem.hap, sizeof(char)*gmem.offset[n_tc][2], cudaMemcpyHostToDevice);
   cudaMemcpy(gmem.d_q, gmem.q, sizeof(int)*gmem.offset[n_tc][1], cudaMemcpyHostToDevice);
   fflush(0);
   cudaMemcpy(gmem.d_ph2pr, gmem.ph2pr, sizeof(NUMBER)*n_tc*128, cudaMemcpyHostToDevice);
   cuerr= cudaGetLastError();
   if (cuerr) printf("Error in memcpy. %d : %s\n", cuerr, cudaGetErrorString(cuerr));
   //One warp handles one matrix
   
   pairhmm_kernel<<<(n_tc+3)/4,WARP*4>>>( gmem.Yr0[0], gmem.d_M, gmem.d_X, 
                                  gmem.d_Y, gmem.d_p, 
                                  gmem.d_rs, gmem.d_hap, gmem.d_q, gmem.d_ph2pr,
                                  gmem.d_offset, n_tc-1, d_out, ctx.LOG10_INITIAL_CONSTANT); 
   if (err) printf("Error!\n");
   cuerr = cudaGetLastError();
   if (cuerr) {
      printf ("Cuda error %d : %s\n", cuerr, cudaGetErrorString(cuerr));
   }
#if 0
   NUMBER *M2;
   NUMBER *X2;
   //TODO Fix this!
   cudaMemcpy(gmem.M, gmem.d_M,
                  sizeof(NUMBER)*gmem.offset[n_tc][0]+1338,
                  cudaMemcpyDeviceToHost);
   cudaMemcpy(gmem.X, gmem.d_X,
                  sizeof(NUMBER)*gmem.offset[n_tc][0]+1338,
                  cudaMemcpyDeviceToHost);
   cuerr = cudaGetLastError();
   if (cuerr) printf ("Memcpy(D2H) error %d : %s\n", cuerr, cudaGetErrorString(cuerr));
   for (int z=0;z<n_tc;z++)
   {
      int ROWS = gmem.offset[z+1][1]-gmem.offset[z][1];
      int COLS = gmem.offset[z+1][2]-gmem.offset[z][2];
      M2 = gmem.M+gmem.offset[z][0]+COLS*(((ROWS-1)+WARP-2)/(WARP-1));
      X2 = gmem.X+gmem.offset[z][0]+COLS*(((ROWS-1)+WARP-2)/(WARP-1));
	   NUMBER result = ctx._(0.0);
	   for (int c = 0; c < COLS; c++)
	   	result += M2[c] + X2[c];

   	if (before_last_log != NULL)
   		*before_last_log = result;	

      probs[z] = ctx.LOG10(out[z]) - ctx.LOG10_INITIAL_CONSTANT;
   }
#else
   cudaMemcpy(probs, d_out, sizeof(NUMBER)*n_tc, cudaMemcpyDeviceToHost);
#endif
   cudaFree(d_out);
} 
template<class NUMBER>
NUMBER compute_full_prob(testcase *tc, NUMBER *before_last_log = NULL)
{
	int r, c;
	int ROWS = tc->rslen + 1;
	int COLS = tc->haplen + 1;

	Context<NUMBER> ctx;

	NUMBER M[ROWS][COLS];
	NUMBER M2[ROWS][COLS];
	NUMBER X[ROWS][COLS];
	NUMBER X2[ROWS][COLS];
   NUMBER Xc0[ROWS];
	NUMBER Y[ROWS][COLS];
   NUMBER Yr0[COLS];
	NUMBER p[ROWS][6];
   extract_tc((NUMBER*)&M[0][0],(NUMBER*)&X[0][0],(NUMBER*)&Y[0][0],
              (NUMBER*)&p[0][0],&Xc0[0],&Yr0[0],tc);
   //     ********     Baseline compute ***********
	for (r = 1; r < ROWS; r++)
   {
		for (c = 1; c < COLS; c++)
		{
			char _rs = tc->rs[r-1];
			char _hap = tc->hap[c-1];
			int _q = tc->q[r-1] & 127;
			NUMBER distm = ctx.ph2pr[_q];
			if (_rs == _hap || _rs == 'N' || _hap == 'N')
				distm = ctx._(1.0) - distm;
			M[r][c] = distm * (M[r-1][c-1] * p[r][MM] + X[r-1][c-1] * p[r][GapM] + Y[r-1][c-1] * p[r][GapM]);
			X[r][c] = M[r-1][c] * p[r][MX] + X[r-1][c] * p[r][XX];
			Y[r][c] = M[r][c-1] * p[r][MY] + Y[r][c-1] * p[r][YY];
		}
   }
#if 0
   //   ********   GPU emulation   ************
   printf("Compute on GPU\n");
   // GPU-esque implementation
   // regs
   NUMBER M_p[WARP+1];
   NUMBER M_pp[WARP+1];
   NUMBER X_p[WARP+1];
   NUMBER X_pp[WARP+1];
   NUMBER Y_p[WARP+1];
   NUMBER Y_pp[WARP+1];
   char _rs[WARP+1];
   int _q[WARP+1];
   NUMBER distm;
   NUMBER pMM[WARP+1];
   NUMBER pGapM[WARP+1];
   NUMBER pXX[WARP+1];
   NUMBER pMX[WARP+1];
   NUMBER pYY[WARP+1];
   NUMBER pMY[WARP+1];
   NUMBER M_loc[WARP+1], X_loc[WARP+1], Y_loc[WARP+1];
   
   for (int stripe = 0; stripe < ROWS; stripe+=WARP-1) 
   {
   //printf("stripe = %d\n", stripe);
   M_pp[0] = M[stripe][0];
   M_pp[1] = M[stripe][0];
   X_pp[0] = X[stripe][0];
   X_pp[1] = X[stripe][0];
   Y_pp[0] = Y[stripe][0];
   Y_pp[1] = Y[stripe][0];
   M_p[0] = M[stripe][1]; 
   M_p[1] = M[stripe+1][0]; 
   X_p[0] = X[stripe][1]; 
   X_p[1] = X[stripe+1][0]; 
   Y_p[0] = Y[stripe][1]; 
   Y_p[1] = Y[stripe+1][0]; 
   for (int tid = 1; tid < WARP+1; tid++)
   {
	   _rs[tid] = tc->rs[tid-1+stripe];
	   _q[tid] = tc->q[tid-1+stripe] & 127;
      pMM[tid] = p[tid+stripe][MM];
      pGapM[tid] = p[tid+stripe][GapM];
      pXX[tid] = p[tid+stripe][XX];
      pMX[tid] = p[tid+stripe][MX];
      pYY[tid] = p[tid+stripe][YY];
      pMY[tid] = p[tid+stripe][MY];
   }
   for (int levi = 1; levi < COLS+WARP+1; levi++) 
   {
      for (int tid = 1; tid < WARP+1; tid++)
      {
         if (tid <= levi && tid+stripe < ROWS && levi-tid+1 < COLS) 
         {  
            r = tid+stripe;
            c = levi-tid+1;
			   char _hap = tc->hap[c-1];
			   if (_rs[tid] == _hap || _rs[tid] == 'N' || _hap == 'N')
			   	distm = double(1.0) - ctx.ph2pr[_q[tid]];
            else distm = ctx.ph2pr[_q[tid]];
            if (tid == 1) {
               //Read from top row
               //To do, read every 32 steps in levi, pass to thread 0 each time
               X_p[0] = X[stripe][c];
               Y_p[0] = Y[stripe][c];
               M_p[0] = M[stripe][c];
            }
            M_loc[tid] = distm * (M_pp[tid] * pMM[tid] + X_pp[tid] * pGapM[tid] + Y_pp[tid] * pGapM[tid]);
			   Y_loc[tid] = M_p[tid] * pMY[tid] + Y_p[tid] * pYY[tid];
			   X_loc[tid] = M_p[tid-1] * pMX[tid] + X_p[tid-1] * pXX[tid];
         }
      }
      for (int tid = WARP; tid > 0; tid--)
      {
         if (tid <= levi && tid+stripe < ROWS && levi-tid+1 < COLS) 
         {  
            M_pp[tid]=M_p[tid-1];
            X_pp[tid]=X_p[tid-1];
            Y_pp[tid]=Y_p[tid-1];
            X_p[tid] = X_loc[tid];
            M_p[tid] = M_loc[tid];
            Y_p[tid] = Y_loc[tid];
            r = tid+stripe;
            c = levi-tid+1;
            if (abs((M_loc[tid]-M[r][c])/M_loc[tid]) > 0.00000000001) printf("M[%d][%d] : %1.16e != %1.16e\n", r,c,M_loc[tid], M[r][c]);
            //TODO write [c][r] for coalescing
            //TODO only write at edges
            M[r][c] = M_loc[tid];
            X[r][c] = X_loc[tid];
            Y[r][c] = Y_loc[tid];
         } else if (tid==levi+1 && tid+stripe < ROWS) {
            M_pp[tid]=M_p[tid-1];
            X_pp[tid]=X_p[tid-1];
            Y_pp[tid]=Y_p[tid-1];
            //Read from col 0
            //TODO make this a row, not a column
            X_p[tid] = X[tid+stripe][0];
            M_p[tid] = M[tid+stripe][0];
            Y_p[tid] = Y[tid+stripe][0];
         }
          
      }
   }
   }
#endif
   NUMBER* d_out;
   NUMBER* d_M;
   NUMBER* d_X;
   NUMBER* d_Y;
   NUMBER* d_p;
   NUMBER* d_Yr0;
   NUMBER* d_Xc0;
   char* d_rs;
   char* d_hap;
   int* d_q; 
   double *d_ph2pr;
   int offsets[100];
   int* d_offsets;
   cudaMalloc(&d_M, sizeof(NUMBER)*ROWS*COLS);
   cudaMalloc(&d_X, sizeof(NUMBER)*ROWS*COLS);
   cudaMalloc(&d_Y, sizeof(NUMBER)*ROWS*COLS);
   cudaMalloc(&d_p, sizeof(NUMBER)*ROWS*6);
   cudaMalloc(&d_rs, sizeof(char)*ROWS);
   cudaMalloc(&d_hap, sizeof(char)*COLS);
   cudaMalloc(&d_q, sizeof(int)*ROWS);
   cudaMalloc(&d_ph2pr, sizeof(double)*128);
   cudaMalloc(&d_Xc0, sizeof(NUMBER)*ROWS);
   cudaMalloc(&d_Yr0, sizeof(NUMBER)*COLS);
   cudaMalloc(&d_offsets, sizeof(int)*100);
   //cudaMemcpy(d_M, &M[0][0], sizeof(NUMBER)*ROWS*COLS,cudaMemcpyHostToDevice);
   //cudaMemcpy(d_X, &X[0][0], sizeof(NUMBER)*ROWS*COLS,cudaMemcpyHostToDevice);
   //cudaMemcpy(d_Y, &Y[0][0], sizeof(NUMBER)*ROWS*COLS,cudaMemcpyHostToDevice);
   cudaMemcpy(d_Xc0, &Xc0, sizeof(NUMBER)*ROWS, cudaMemcpyHostToDevice);
   cudaMemcpy(d_Yr0, &Yr0, sizeof(NUMBER)*COLS, cudaMemcpyHostToDevice);
   cudaMemcpy(d_p, &p[0][0], sizeof(NUMBER)*ROWS*6,cudaMemcpyHostToDevice);
   cudaMemcpy(d_rs, tc->rs, sizeof(char)*ROWS,cudaMemcpyHostToDevice);
   cudaMemcpy(d_hap, tc->hap, sizeof(char)*COLS,cudaMemcpyHostToDevice);
   cudaMemcpy(d_q, tc->q, sizeof(int)*ROWS,cudaMemcpyHostToDevice);
   cudaMemcpy(d_ph2pr, ctx.ph2pr, sizeof(double)*128,cudaMemcpyHostToDevice);
   cudaMalloc(d_out, sizeof(NUMBER));
   offsets[0]=0; offsets[1]=0; offsets[2]=0;
   offsets[3]=ROWS*COLS; offsets[4]=ROWS; offsets[5]=COLS;
   cudaMemcpy(d_offsets, offsets, sizeof(int)*100,cudaMemcpyHostToDevice);
   printf("%d ROWS, %d COLS\n", ROWS, COLS);
  
   pairhmm_kernel<<<1,WARP>>>( Yr0[0], d_M, d_X, d_Y, d_p, 
                   d_rs, d_hap, d_q, d_ph2pr, d_offsets, 1, d_out); 

   cudaMemcpy(M2[ROWS-1], d_M+COLS*(ROWS-1), sizeof(NUMBER)*COLS,cudaMemcpyDeviceToHost);
   cudaMemcpy(X2[ROWS-1], d_X+COLS*(ROWS-1), sizeof(NUMBER)*COLS,cudaMemcpyDeviceToHost);

#ifdef DEBUG_MODE
   printf("Error check...\n");
   cudaMemcpy(&M2[0][0], d_M, sizeof(NUMBER)*ROWS*COLS,cudaMemcpyDeviceToHost);
   cudaMemcpy(&X2[0][0], d_X, sizeof(NUMBER)*ROWS*COLS,cudaMemcpyDeviceToHost);
   for (int a=1;a<ROWS;a++) {
   for (int b=1;b<COLS;b++) {
      if (abs((M[a][b]-M2[a][b])/MAX(M[a][b],M2[a][b]))>0.00000000) 
          printf("M: %d,%d = %e != %e\n", a,b,M[a][b],M2[a][b]);
      if (abs((X[a][b]-X2[a][b])/MAX(X[a][b],X2[a][b]))>0.00000000)  
          printf("X: %d,%d = %e != %e\n", a,b,X[a][b],X2[a][b]);
   }
   }
#endif
   
   cudaFree(d_M);
   cudaFree(d_X);
   cudaFree(d_Y);
   cudaFree(d_p);
   cudaFree(d_rs);
   cudaFree(d_hap);
   cudaFree(d_q);
   cudaFree(d_ph2pr);

	NUMBER result = ctx._(0.0);
	for (c = 0; c < COLS; c++)
		result += M[ROWS-1][c] + X[ROWS-1][c];

	if (before_last_log != NULL)
		*before_last_log = result;	

   cudaFree(d_out);
	return ctx.LOG10(result) - ctx.LOG10_INITIAL_CONSTANT;
}

int main(int argc, char* argv[])
{
	testcase tc[10000];
   int cnt=0;
   int basecnt=0;
   double prob[10000];
   GPUmem<double> gmem;
  
   std::ifstream infile;

   infile.open((char*) argv[1]);
  
#if 1
	while (read_testcase(tc+cnt) == 0)
   {
      if (cnt==799) {
         printf("Computing %d testcases\n", cnt+1);
         compute_full_prob_multiple<double>(prob, tc, cnt+1, gmem);
         for (int q=0;q<cnt+1;q++)
		      printf("%E\n", q+basecnt, prob[q]);
         cnt = -1;
         basecnt+=800;
      }
      cnt++;

   }
   
    //TODO remove
    //cnt=1;
   printf("Computing %d testcases\n", cnt);
   if (cnt>0) compute_full_prob_multiple<double>(prob, tc, cnt, gmem);

   //This call frees memory in gmem
   compute_full_prob_multiple<double>(prob, tc, 0, gmem);

   for (int q=0;q<cnt;q++)
     printf("ans %d: %E\n", q+basecnt, prob[q]);
#else
	while (read_testcase(tc,infile) == 0) {
      printf("%E\n", compute_full_prob<double>(tc));
      //if (1<cnt++) break;
   }
#endif
	return 0;
}

