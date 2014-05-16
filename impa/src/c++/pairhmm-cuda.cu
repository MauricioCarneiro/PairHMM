#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "input.h"
#include "compute_gpu.h"

#include "timing.h"
//#include "read.h"

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

template<class NUMBER>
int tc2gmem(GPUmem<NUMBER>& gmem, testcase* tc, int index)
{
	int ROWS = tc->rslen + 1;
	int COLS = tc->haplen + 1;

	Context<NUMBER> ctx;


	//NUMBER M[ROWS][COLS];
	//NUMBER M2[ROWS][COLS];
	//NUMBER X[ROWS][COLS];
	//NUMBER X2[ROWS][COLS];
   //NUMBER Xc0[ROWS];
	//NUMBER Y[ROWS][COLS];
   //NUMBER Yr0[COLS];
	//NUMBER p[ROWS][6];
   //NUMBER q_new[ROWS];
   if (ROWS==1 || COLS==1) return 0;
   extract_tc(gmem.M,gmem.X,gmem.Y,gmem.p+gmem.offset[index].x*6,gmem.n+gmem.offset[index].x,
              gmem.q+gmem.offset[index].x,tc);
   //TODO check data sizes first
   //      return error if we're out of space
   //memcpy(gmem.Xc0+gmem.offset[gmem.index][1], Xc0, sizeof(NUMBER)*ROWS);
   //memcpy(gmem.Yr0+gmem.offset[gmem.index][2], Yr0, sizeof(NUMBER)*COLS);
   //memcpy(gmem.p+gmem.offset[gmem.index][1]*6, &p[0][0], sizeof(NUMBER)*ROWS*6);
#ifndef __NO_COPY
   memcpy(gmem.rs+gmem.offset[index].x, tc->rs, sizeof(char)*ROWS);
   memcpy(gmem.hap+gmem.offset[index].y, tc->hap, sizeof(char)*COLS);
#endif
   //memcpy(gmem.q+gmem.offset[gmem.index][1], q_new, sizeof(NUMBER)*ROWS);
   //gmem.offset[gmem.index+1][0] = gmem.offset[gmem.index][0] + ROWS*COLS;
   //gmem.offset[gmem.index+1][0] = gmem.offset[gmem.index][0] + ((ROWS+WARP-2)/(WARP-1))*COLS;
   //printf("ROWS/COLS = %d/%d, offset[%d] = %d\n", ROWS, COLS, gmem.index+1, gmem.offset[gmem.index+1][0]);
   //gmem.offset[gmem.index+1][1] = gmem.offset[gmem.index][1] + ROWS;
   //gmem.offset[gmem.index+1][2] = gmem.offset[gmem.index][2] + COLS;
   gmem.index++;
   return 0;
}
template<class NUMBER>
void extract_tc(NUMBER* M_in, NUMBER* X_in, NUMBER* Y_in, 
                NUMBER* p_in, int* n_new, NUMBER* q_new, testcase* tc)
{
	int ROWS = tc->rslen + 1;
   int r;

	Context<NUMBER> ctx;
#ifndef __NO_COPY
	for (r = 1; r < ROWS; r++)
	{
		int _i = tc->i[r-1] & 127;
		int _d = tc->d[r-1] & 127;
		int _c = tc->c[r-1] & 127;
      int _q = tc->q[r-1] & 127;
      //n_new[3*r+II] = _i;
      //n_new[3*r+CC] = _c;
      //n_new[3*r+DD] = _d;
      n_new[r]=_i+128*_d+128*128*_c+128*128*128*_q;
	}
#endif

}

int tc_comp(const void* tc_A, const void* tc_B) {
   int rA = 32*((((testcase*)tc_A)->rslen+31)/32); 
   int rB = 32*((((testcase*)tc_B)->rslen+31)/32); 
   if (rA*((testcase*)tc_A)->haplen + 32*32*(rA/32) <  
       rB*((testcase*)tc_B)->haplen + 32*32*(rB/32)) return 1;
   else return -1;
}
int tc_comp_unsort(const void* tc_A, const void* tc_B) {
   return ((testcase*)tc_A)->index - ((testcase*)tc_B)->index;
}

#define round_up(A,B) (B)*((A+B-1)/(B))
//#define round_up(A,B) A

template<class NUMBER>
void compute_full_prob_multiple(NUMBER* probs, testcase *tc, int n_tc, 
                                 GPUmem<NUMBER> &gmem, NUMBER *before_last_log = NULL) {
   Context<NUMBER> ctx;
   int err;

   Timing GPUAlloc(string("GPU Alloc/Free :  "));
   GPUAlloc.start();
   //qsort(tc, n_tc, sizeof(testcase), tc_comp);
   printf("largest mat: %d x %d\n", tc[0].rslen, tc[0].haplen);
   if (0==n_tc) {
      err = GPUmemFree<NUMBER>(gmem);
      return;
   }
   if (gmem.M==0) {
      err = GPUmemAlloc<NUMBER>(gmem);
   }
   if (err != 0) printf("Error in GPU allocation/deallocation\n");
   GPUAlloc.acc();

   Timing Staging(string("Staging :  "));
   Timing ComputeGPU(string("Compute GPU Time :  "));
   gmem.index=0;
   int total_rows=0;
   int total_cols=0;
   int total_scratch=0;
   int total_cells = 0;
   for (int z=0;z<n_tc;z++)
   {
      //gmem.offset[z][0] = total_scratch;
      gmem.offset[z].x = total_rows;
      gmem.offset[z].y = total_cols;
      total_rows += tc[z].rslen+1;
      total_cols += tc[z].haplen+1;
      total_scratch += ((tc[z].rslen+WARP-1)/(WARP-1))*(tc[z].haplen+1);
      total_cells += tc[z].rslen*tc[z].haplen;
   }
   //TODO clean this up!
   //gmem.offset[n_tc][0] = total_scratch;
   gmem.offset[n_tc].x = total_rows;
   gmem.offset[n_tc].y = total_cols;
   gmem.X = gmem.M + total_cols;
   gmem.d_X = gmem.d_M + total_cols;
   gmem.Y = gmem.X + total_cols;
   gmem.d_Y = gmem.d_X + total_cols;
   //q and n must be aligned to 512 bytes for transfer as textures
   gmem.q = gmem.M + round_up(3*total_cols, 512/sizeof(NUMBER));
   gmem.d_q = gmem.d_M + round_up(3*total_cols, 512/sizeof(NUMBER));
   if (((char*)gmem.d_q-(char*)gmem.d_M)%512 != 0) printf("d_q not aligned\n");
   gmem.n = (int*)(gmem.q + round_up(total_rows,512/sizeof(NUMBER)));
   gmem.d_n = (int*)(gmem.d_q + round_up(total_rows,512/sizeof(NUMBER)));
   if (((char*)gmem.d_n-(char*)gmem.d_M)%512 != 0) printf("d_n not aligned\n");
   gmem.rs = (char*)(gmem.n + total_rows*3);
   gmem.d_rs = (char*)(gmem.d_n + total_rows*3);
   gmem.hap = gmem.rs + total_rows;
   gmem.d_hap = gmem.d_rs + total_rows;
   //Make sure results and d_results are properly aligned
   //TODO put results before the input data and scratch to avoid thise complexity
   gmem.results = (NUMBER*)(gmem.hap + total_cols + sizeof(NUMBER)-(13*total_rows+total_cols)%sizeof(NUMBER));
   gmem.d_results = (NUMBER*)(gmem.d_hap + total_cols + sizeof(NUMBER)-(13*total_rows+total_cols)%sizeof(NUMBER));
   if ((char*)gmem.results+n_tc*sizeof(NUMBER)-(char*)gmem.M > gmem.totalMem) {
      printf("data exceeds GPU memory. Quitting.");
      return;
   }
   int s=0;
   cudaStream_t marker_s;
   cudaStreamCreate(&marker_s);
   for (int start=0;start<n_tc;start+=n_tc/gmem.N_STREAMS) {
      int finish=min(start+n_tc/gmem.N_STREAMS, n_tc);
      Staging.start();
      //CPU_start<<<1,1,0,marker_s>>>();
#pragma omp parallel for shared(gmem, tc) private (z)
      for (int z=start;z<finish;z++)
      {
         err = tc2gmem<NUMBER>(gmem, &tc[z], z);
      }
      //CPU_end<<<1,1,0,marker_s>>>();
      Staging.acc();
      ComputeGPU.start();
      cudaStreamSynchronize(gmem.strm[s]);
      compute_gpu_stream(gmem.offset+start, gmem.rs, gmem.hap, gmem.q, gmem.n, 
                         ctx.INITIAL_CONSTANT, finish-start, gmem, gmem.strm[s], start);
      ComputeGPU.acc();
      s++;
      s %= gmem.N_STREAMS;
   }
   ComputeGPU.start();
   for (s=0;s<gmem.N_STREAMS;s++) cudaStreamSynchronize(gmem.strm[s]);
   //memcpy(probs, gmem.results, sizeof(NUMBER)*n_tc);
   for (int z=0;z<n_tc;z++) probs[tc[z].index] = gmem.results[z];
   ComputeGPU.acc();
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
#if 0
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
                   d_rs, d_hap, d_q, d_offsets, 1, d_out); 

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
   cudaFree(d_out);
#endif

	NUMBER result = ctx._(0.0);
	for (c = 0; c < COLS; c++)
		result += M[ROWS-1][c] + X[ROWS-1][c];

	if (before_last_log != NULL)
		*before_last_log = result;	

	return ctx.LOG10(result) - ctx.LOG10_INITIAL_CONSTANT;
}

int main(int argc, char* argv[])
{
   Timing TotalTime(string("TOTAL: "));
   Timing ComputationTime(string("COMPUTATION: "));
   TotalTime.start();
	testcase *tc = new testcase[MAX_PROBS];
   int cnt=0;
   int basecnt=0;
   double *prob;
   prob = (double*)malloc(MAX_PROBS*sizeof(double));
   GPUmem<double> gmem;
  
   std::ifstream infile;

   if (argc>1) {
      infile.open((char*) argv[1]);
   } 
   
	while (read_testcase(tc+cnt, argc>1 ? infile: std::cin) == 0)
   {
      //printf("In pairhmm-cuda: &tc[%d] = %p\n", cnt, tc+cnt);
      //(tc+cnt)->display();
      tc[cnt].index=cnt;
      if (cnt==MAX_PROBS-1) {
         printf("Computing %d testcases\n", cnt+1);
         fflush(0);
         ComputationTime.start();
         compute_full_prob_multiple<double>(prob, tc, cnt+1, gmem);
         ComputationTime.acc();
         for (int q=0;q<cnt+1;q++) {
		      printf("%E\n", q+basecnt, prob[q]);
         }
         cnt = -1;
         basecnt+=MAX_PROBS;
      }
      cnt++;

   }
   
    //TODO remove
    //cnt=1;
   printf("Computing %d testcases\n", cnt);
   ComputationTime.start();
   if (cnt>0) compute_full_prob_multiple<double>(prob, tc, cnt, gmem);
   ComputationTime.acc();

   //This call frees memory in gmem
   compute_full_prob_multiple<double>(prob, tc, 0, gmem);

   for (int q=0;q<cnt;q++) {
     printf("%E\n", q+basecnt, prob[q]);
   }

   TotalTime.acc();

   delete []tc;
   free(prob);
	return 0;
}

