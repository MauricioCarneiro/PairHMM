//#include "input.h"
#include "compute_gpu.h"
#include "compute_gpu.cuh"
#include "stdio.h"
#include <cfloat>
#define MAT_PER_WARP 8
//#define FOCUS 7046
#define FOCUS 24022
//#define FOCUS 12
#define MIN_ACCEPTED 1e-28f
#define ROUND_UP(A,B) (B)*((A+B-1)/(B))

GPUmem<float> gf_gmem;
GPUmem<double> gd_gmem;

template <class T>
void GPUAlloc();
template <class T>
void GPUFree();
template <class T>
void compute_gpu_head(int* offset_in, char* rs, char* hap, int* n,
                      T init_const, int n_tc, double* results, double FAILED_RUN_RESULT);
template <class T>
void gpu_init(int* offset_in, int n_tc);

template <>
void GPUAlloc<double>() {
   GPUmemAlloc(gd_gmem);
}
template <>
void GPUFree<double>() {
   printf("Free\n");fflush(0);
   GPUmemFree(gd_gmem);
}
template <>
void GPUAlloc<float>() {
   GPUmemAlloc(gf_gmem);
}
template <>
void GPUFree<float>() {
   GPUmemFree(gf_gmem);
}
bool isGPUAllocated() {
   return (gd_gmem.amem || gf_gmem.amem);
}
template <> 
void compute_gpu_head<float>(int* offset_in, char* rs, char* hap, int* n,
                      float init_const, int n_tc, double* results, double FAILED_RUN_RESULT) {
   int2* offset= (int2*) offset_in;
   if (!isGPUAllocated()) GPUAlloc<float>();
   gf_gmem.results = results;
   char* current = (char*)gf_gmem.d_amem;

   gf_gmem.d_results = (double*) current;
   current += sizeof(double) * n_tc;

   gf_gmem.d_offset = (int2*) current;
   current += sizeof(int2) * n_tc;

   gf_gmem.d_n = (int*)current;
   current += sizeof(int) * offset[n_tc].x;

   gf_gmem.d_rs = (char*)current;
   current += offset[n_tc].x;
   current += ROUND_UP(current-gf_gmem.d_amem, 512);

   gf_gmem.d_hap = (char*)current;
   current += offset[n_tc].y;

   if (current-gf_gmem.d_amem > gf_gmem.totalMem) {
     fprintf(stderr, "Error. Problem exceeds GPU memory size. Aborting GPU computation\n"); 
     return;
   }
   compute_gpu_stream((int2*)offset, rs, hap, (float*)NULL, n, init_const, n_tc, gf_gmem, 0, 0, FAILED_RUN_RESULT);
   gf_gmem.results= NULL;
   //if (isGPUAllocated()) GPUFree<float>();
   cudaCheckError(__LINE__, __FILE__);
}
template <> 
void gpu_init<float> (int* offset_in, int n_tc) {
   int2* offset= (int2*) offset_in;
   if (!isGPUAllocated()) {
      GPUAlloc<float>();
   }
   char* current = (char*)gf_gmem.d_amem;

   gf_gmem.d_results = (double*) current;
   current += sizeof(double) * n_tc;

   gf_gmem.d_offset = (int2*) current;
   current += sizeof(int2) * n_tc;

   gf_gmem.d_n = (int*)current;
   current += sizeof(int) * offset[n_tc].x;

   gf_gmem.d_rs = (char*)current;
   current += offset[n_tc].x;
   current += ROUND_UP(current-gf_gmem.d_amem, 512);

   gf_gmem.d_hap = (char*)current;
   current += offset[n_tc].y;

   if (current-gf_gmem.d_amem > gf_gmem.totalMem) {
     fprintf(stderr, "Error. Problem exceeds GPU memory size. Aborting GPU computation\n"); 
     return;
   }
}
template <> 
void gpu_init<double> (int* offset_in, int n_tc) {
   int2* offset= (int2*) offset_in;
   if (!isGPUAllocated()) {
      GPUAlloc<double>();
   }
   char* current = (char*)gd_gmem.d_amem;

   gd_gmem.d_results = (double*) current;
   current += sizeof(double) * n_tc;

   gd_gmem.d_offset = (int2*) current;
   current += sizeof(int2) * n_tc;

   gd_gmem.d_n = (int*)current;
   current += sizeof(int) * offset[n_tc].x;

   gd_gmem.d_rs = (char*)current;
   current += offset[n_tc].x;
   current += ROUND_UP(current-gd_gmem.d_amem, 512);

   gd_gmem.d_hap = (char*)current;
   current += offset[n_tc].y;

   if (current-gd_gmem.d_amem > gd_gmem.totalMem) {
     fprintf(stderr, "Error. Problem exceeds GPU memory size. Aborting GPU computation\n"); 
     return;
   }


}
template <> 
void compute_gpu_head<double>(int* offset_in, char* rs, char* hap, int* n,
                      double init_const, int n_tc, double* results, double FAILED_RUN_RESULT) {
   gpu_init<double>(offset_in, n_tc);
   gd_gmem.results = results;
   compute_gpu_stream((int2*)offset_in, rs, hap, (double*)NULL, n, init_const, n_tc, gd_gmem, 0, 0, FAILED_RUN_RESULT );
   gd_gmem.results= NULL;
   //if (isGPUAllocated()) GPUFree<double>();
   cudaCheckError(__LINE__, __FILE__);
}


void cudaCheckError(int line, const char* file) {
#if 0
   cudaError_t err = cudaGetLastError();
   if (err) {
       fprintf(stderr, "error %d on line %d : %s\n", err, line, cudaGetErrorString(err));
       fflush(0);
   }
#endif
}

void createNewTextureFloat(cudaTextureObject_t& tex, cudaResourceDesc& resDesc, cudaTextureDesc& texDesc, void* devPtr, int size) {

   tex=0;
   memset(&resDesc, 0, sizeof(cudaResourceDesc));
   resDesc.res.linear.devPtr = devPtr;
   resDesc.resType = cudaResourceTypeLinear;
   resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
   resDesc.res.linear.desc.x = 32;
   resDesc.res.linear.sizeInBytes = size*sizeof(float);
   memset(&texDesc, 0, sizeof(cudaTextureDesc));
   texDesc.readMode = cudaReadModeElementType;
   cudaCreateTextureObject(&tex, &resDesc, &texDesc, NULL);
   cudaCheckError(__LINE__, __FILE__);
}

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

template <class PRECISION>
__device__ PRECISION ph2pr(int in) {
   return pow((PRECISION)10.0,-((PRECISION)(in & 127))/(PRECISION)10.0);
}

#define MAX_ROWS 310
template<class NUMBER, int BATCH>
__global__ void
__launch_bounds__(4*WARP,6)
pairhmm_jacopo(
    const char* rs,
    const char* hap,
    const NUMBER*  q, 
    const int* idc,
    cudaTextureObject_t t_n,
    const int2* offsets,    
    const int   n_mats,
    double*     output,
    NUMBER      INITIAL_CONSTANT,
    NUMBER      LOG10_INITIAL_CONSTANT,
    double      FAILED_RUN_RESULT) 
{
    
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    if (pid >= n_mats) return;
    int r, c;
    __shared__ NUMBER ph2pr_save[128];

    //Precompute and store a mapping from ints to floating point qualities
    for (int z=threadIdx.x; z< 128; z+=WARP) 
    {
       ph2pr_save[z] = ph2pr<NUMBER>(z);
    }
    __syncthreads();

    // fetch the row and column offsets
    const int row_off = offsets[pid].x;
    const int col_off = offsets[pid].y;
    const int ROWS    = offsets[pid+1].x - row_off;
    const int COLS    = offsets[pid+1].y - col_off;

    rs+=row_off;
    hap+=col_off;
    idc+=row_off;

    // precompute p[][]
	NUMBER l_p[4][MAX_ROWS];
   

   l_p[MM][0]    = NUMBER(0.0);
	l_p[MX][0]    = NUMBER(0.0);
	l_p[XX][0]    = NUMBER(1.0);
	l_p[MY][0]    = NUMBER(1.0);
   int  l_q[MAX_ROWS];
#if 1
	for (r = 1; r < ROWS; r++)
	{
#ifdef __USE_TEX
        const int _idc = tex1Dfetch<int>(t_n, r + row_off);
#else
        const int _idc = idc[ r ];              //Four quality metrics (0-127) are combined into 1 int
#endif
                                                             // NOTE: not sure if it's worth doing, but this could also be transposed
                                                             // for strided access...
		const int _i = _idc & 127;
		const int _d = (_idc >> 7) & 127;
		const int _c = (_idc >> 14) & 127;
		l_p[MM][r]    = NUMBER(1.0) - ph2pr_save[(_i + _d) & 127];
		l_p[MX][r]    = ph2pr_save[_i];
		l_p[XX][r]    = ph2pr_save[_c];
		l_p[MY][r]    = (r == ROWS - 1) ? NUMBER(1.0) : ph2pr_save[_d];

      l_q[r] = (_idc >> 21) & 127;
	}
#endif

    // local memory copies of data used in the critical path
    char l_rs[MAX_ROWS];
    char l_hap[MAX_ROWS];
#if 1
    for (r = 1; r < ROWS; ++r)
    {
        l_rs[r]  = rs[r-1];       // gmem inputs (NOTE: as above, these 2 vectors could be transposed)
    }

    for (c = 1; c < COLS; ++c)
	    l_hap[c] = hap[c-1];      // gmem inputs (NOTE: as above)
#endif

    NUMBER result = NUMBER(0.0);

    // local memory columns of M, X, Y
    NUMBER M_c[MAX_ROWS];
    NUMBER X_c[MAX_ROWS];
    NUMBER Y_c[MAX_ROWS];

    // initialize the first stripe, i.e. the 0-th column
	M_c[0] = NUMBER(0.0);
	X_c[0] = NUMBER(0.0);
   Y_c[0] = INITIAL_CONSTANT / (COLS-1);
	for (r = 1; r < ROWS; r++)
	{
		M_c[r] = NUMBER(0.0);
		X_c[r] = NUMBER(0.0);
		Y_c[r] = NUMBER(0.0);
	}

#if 1
    c=0;
    // stripe the matrix in BATCH-wide blocks
    for (int stripe = 0; stripe < COLS; stripe += BATCH)
    {
        // register blocks for M, X and Y
	    NUMBER M[BATCH+1];
	    NUMBER X[BATCH+1];
	    NUMBER Y, Y_p;

        // initialize the zero-th row
        for (int bid = 1; bid <= BATCH; ++bid)
        {
		    M[bid] = NUMBER(0.0);
		    X[bid] = NUMBER(0.0);
        }
		  Y = INITIAL_CONSTANT / (COLS-1); 

        // register block for l_hap
        char H[BATCH+1];
        #pragma unroll
        for (int bid = 1; bid <= BATCH; ++bid)
            H[bid] = l_hap[stripe + bid];

        // loop across all rows
        for (r = 1; r < ROWS; r++)
        {
			const char _rs = l_rs[r];
			const int  _q  = l_q[r];

            const NUMBER distm_q = ph2pr_save[_q];  // should use __ldg();

            // load the M,X,Y entries from the previous row of the previous stripe
            M[0] = M_c[r-1];
            X[0] = X_c[r-1];
            Y = Y_c[r-1];
            M_c[r-1] = M[BATCH];
            X_c[r-1] = X[BATCH];

            // prefetch p[r][*] in registers (though the compiler might do it)
            NUMBER p_MM   = l_p[MM][r];
            NUMBER p_MX   = l_p[MX][r];
            NUMBER p_XX   = l_p[XX][r];
            NUMBER p_MY   = l_p[MY][r-1];
            NUMBER p_GapM = NUMBER(1.0) - l_p[XX][r];
            NUMBER p_YY   = r==ROWS ? NUMBER(1.0) : l_p[XX][r-1];

	        NUMBER M_p[BATCH+1];            // save M[r-1]
            #pragma unroll
            for (int bid = 0; bid <= BATCH; ++bid)
                M_p[bid] = M[bid];

            NUMBER X_p = X[0];              // save X[r-1][0]

            // loop across columns in this stripe
            #pragma unroll
            for (int bid = 1; bid <= BATCH; ++bid)
            {
                // fetch the corresponding hap entry
                const char _hap = H[bid];

                const NUMBER distm = (_rs == _hap || _rs == 'N' || _hap == 'N') ?
                      NUMBER(1.0) - distm_q :
                                    distm_q;

                Y_p = Y;
    			    Y = M_p[bid-1] * p_MY + Y_p * p_YY;
                M[bid] = distm * (M_p[bid-1] * p_MM + X_p * p_GapM + Y_p * p_GapM);

                X_p = X[bid]; // save X[r-1][c]
                X[bid] = M_p[bid] * p_MX + X[bid] * p_XX;
            }
            Y_c[r-1] = Y;

        }

        // add up the result
        for (int bid = 1; bid <= BATCH; ++bid)
        {
            // compute the column index
            const int c = stripe + bid;
            if (c < COLS) {
                result += M[bid] + X[bid];
            }
        }
    }
#endif

   if (result < MIN_ACCEPTED) {
      output[pid] = FAILED_RUN_RESULT;
   } else {
	   output[pid] = log10(result) - LOG10_INITIAL_CONSTANT;
   }
}
template<class NUMBER, int BATCH>
__global__ void
__launch_bounds__(4*WARP,6)
pairhmm_jacopo_new(
    const char* rs,
    const char* hap,
    const NUMBER*  q, 
    const int* idc,
    cudaTextureObject_t t_n,
    const int2* offsets,    
    const int   n_mats,
    NUMBER*     output,
    NUMBER      INITIAL_CONSTANT,
    NUMBER      LOG10_INITIAL_CONSTANT) 
{
    
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    if (pid >= n_mats) return;
    int r, c;
    __shared__ NUMBER ph2pr_save[128];

    //Precompute and store a mapping from ints to floating point qualities
    for (int z=threadIdx.x; z< 128; z+=WARP) 
    {
       ph2pr_save[z] = ph2pr<NUMBER>(z);
    }
    __syncthreads();

    // fetch the row and column offsets
    const int row_off = offsets[pid].x;
    const int col_off = offsets[pid].y;
    const int ROWS    = offsets[pid+1].x - row_off;
    const int COLS    = offsets[pid+1].y - col_off;

    rs+=row_off;
    hap+=col_off;
    idc+=row_off;

    // precompute p[][]
	NUMBER l_p[4][MAX_ROWS];
   

   l_p[MM][0]    = NUMBER(0.0);
	l_p[MX][0]    = NUMBER(0.0);
	l_p[XX][0]    = NUMBER(1.0);
	l_p[MY][0]    = NUMBER(1.0);
   int  l_q[MAX_ROWS];
#if 1
	for (r = 1; r < ROWS; r++)
	{
#ifdef __USE_TEX
        const int _idc = tex1Dfetch<int>(t_n, r + row_off);
#else
        const int _idc = idc[ r ];           //Four quality metrics (0-127) are combined into 1 int
#endif
		const int _i = _idc & 127;
		const int _d = (_idc >> 7) & 127;
		const int _c = (_idc >> 14) & 127;
     //TODO (_i+_d)&127 This is here to match results from pairhmm-1-base
		l_p[MM][r]    = NUMBER(1.0) - ph2pr_save[_i + _d & 127];
		l_p[MX][r]    = ph2pr_save[_i];
		l_p[XX][r]    = ph2pr_save[_c];
		l_p[MY][r]    = (r == ROWS - 1) ? NUMBER(1.0) : ph2pr_save[_d];

      l_q[r] = (_idc >> 21) & 127;
	}
#endif

    // local memory copies of data used in the critical path
    char l_rs[MAX_ROWS];
    char l_hap[MAX_ROWS];
#if 1
    for (r = 1; r < ROWS; ++r)
    {
        l_rs[r]  = rs[r-1];    // gmem inputs (NOTE: as above, these 2 vectors could be transposed)
    }

    for (c = 1; c < COLS; ++c)
	    l_hap[c] = hap[c-1];      // gmem inputs (NOTE: as above)
#endif

    NUMBER result = NUMBER(0.0);

    // local memory columns of M, X, Y
    NUMBER M_c[MAX_ROWS];
    NUMBER X_c[MAX_ROWS];
    NUMBER Y_c[MAX_ROWS];

    // initialize the first stripe, i.e. the 0-th column
	M_c[0] = NUMBER(0.0);
	X_c[0] = NUMBER(0.0);
   Y_c[0] = INITIAL_CONSTANT / (COLS-1);
	for (r = 1; r < ROWS; r++)
	{
		M_c[r] = NUMBER(0.0);
		X_c[r] = NUMBER(0.0);
		Y_c[r] = NUMBER(0.0);
	}

#if 1
    c=0;
    // stripe the matrix in BATCH-wide blocks
    for (int stripe = 0; stripe < COLS; stripe += BATCH)
    {
        // register blocks for M, X and Y
	    NUMBER M[BATCH+1];
	    NUMBER X[BATCH+1];
	    NUMBER Y, Y_p;

        // initialize the zero-th row
        for (int bid = 1; bid <= BATCH; ++bid)
        {
		    M[bid] = NUMBER(0.0);
		    X[bid] = NUMBER(0.0);
        }
		  Y = INITIAL_CONSTANT / (COLS-1); 

        // register block for l_hap
        char H[BATCH+1];
        #pragma unroll
        for (int bid = 1; bid <= BATCH; ++bid)
            H[bid] = l_hap[stripe + bid];

        // loop across all rows
        for (r = 1; r < ROWS; r++)
        {
			const char _rs = l_rs[r];
			const int  _q  = l_q[r];

            const NUMBER distm_q = ph2pr_save[_q];  // should use __ldg();

            // load the M,X,Y entries from the previous row of the previous stripe
            M[0] = M_c[r-1];
            X[0] = X_c[r-1];
            Y = Y_c[r-1];
            M_c[r-1] = M[BATCH];
            X_c[r-1] = X[BATCH];

            // prefetch p[r][*] in registers (though the compiler might do it)
            NUMBER p_MM   = l_p[MM][r];
            NUMBER p_MX   = l_p[MX][r];
            NUMBER p_XX   = l_p[XX][r];
            NUMBER p_MY   = l_p[MY][r-1];
            NUMBER p_GapM = NUMBER(1.0) - l_p[XX][r];
            NUMBER p_YY   = r==ROWS ? NUMBER(1.0) : l_p[XX][r-1];

	        NUMBER M_p[BATCH+1];            // save M[r-1]
            #pragma unroll
            for (int bid = 0; bid <= BATCH; ++bid)
                M_p[bid] = M[bid];

            NUMBER X_p = X[0];              // save X[r-1][0]

            // loop across columns in this stripe
            #pragma unroll
            for (int bid = 1; bid <= BATCH; ++bid)
            {
                // fetch the corresponding hap entry
                const char _hap = H[bid];

                const NUMBER distm = (_rs == _hap || _rs == 'N' || _hap == 'N') ?
                      NUMBER(1.0) - distm_q :
                                    distm_q;

                Y_p = Y;
    			    Y = M_p[bid-1] * p_MY + Y_p * p_YY;
                M[bid] = distm * (M_p[bid-1] * p_MM + X_p * p_GapM + Y_p * p_GapM);

                X_p = X[bid]; // save X[r-1][c]
                X[bid] = M_p[bid] * p_MX + X[bid] * p_XX;
            }
            Y_c[r-1] = Y;

        }

        // add up the result
        for (int bid = 1; bid <= BATCH; ++bid)
        {
            // compute the column index
            const int c = stripe + bid;
            if (c < COLS) {
                result += M[bid] + X[bid];
            }
        }
    }

#endif
   if (result < MIN_ACCEPTED) {
      output[pid] = 1.0;
   } else {
	   output[pid] = log10(result) - LOG10_INITIAL_CONSTANT;
   }
}
#define MAX_COLS 310
template<class NUMBER, int BATCH>
__global__
//__launch_bounds__(WARP,1)
void pairhmm_1thread_sweepright(
    const char* rs,
    const char* hap,
    const NUMBER*  q, 
    const int* idc,
    cudaTextureObject_t t_n,
    const int2* offsets,    
    const int   n_mats,
    NUMBER*     output,
    NUMBER      INITIAL_CONSTANT,
    NUMBER      LOG10_INITIAL_CONSTANT) 
{
    
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    if (pid >= n_mats) return;
    int r, c;

    // fetch the row and column offsets
    const int row_off = offsets[pid].x;
    const int col_off = offsets[pid].y;
    const int ROWS    = offsets[pid+1].x - row_off;
    const int COLS    = offsets[pid+1].y - col_off;

    q+=row_off;
    rs+=row_off;
    hap+=col_off;
    idc+=row_off;


	 NUMBER pMM[BATCH+1];
	 NUMBER pMX[BATCH+1];
	 NUMBER pXX[BATCH+1];
	 NUMBER pMY[BATCH+1];
    NUMBER distm_q[BATCH+1];
    NUMBER _rs[BATCH+1];
    char _q[BATCH+1];

    NUMBER result = NUMBER(0.0);

    // local memory columns of M, X, Y
    NUMBER M_r[MAX_COLS];
    NUMBER X_r[MAX_COLS];
    NUMBER Y_r[MAX_COLS];

    // initialize the first stripe, i.e. the 0-th column
	for (c = 0; c < COLS; c++)
	{
		M_r[c] = NUMBER(0.0);
		X_r[c] = NUMBER(0.0);
		Y_r[c] = INITIAL_CONSTANT / (COLS-1);
	}

    // stripe the matrix in BATCH-wide blocks
    for (int stripe = 0; stripe < ROWS; stripe += BATCH)
    {
        // register blocks for M, X and Y
	    NUMBER M[BATCH+1];
	    NUMBER X[BATCH+1];
	    NUMBER Y[BATCH+1];

       //TODO pad by BATCH
       for(int bid = 1; bid <= BATCH && bid + stripe < ROWS; bid++)
       {
          r = bid + stripe;

#ifdef __USE_TEX
           const int _idc = tex1Dfetch<int>(t_n, r + row_off);
#else
           const int _idc = idc[ r ];              //Four quality metrics (0-127) are combined into 1 int
#endif
                                                             // NOTE: not sure if it's worth doing, but this could also be transposed
                                                             // for strided access...
   		const int _i = _idc & 127;
   		const int _d = (_idc >> 7) & 127;
   		const int _c = (_idc >> 14) & 127;
        //TODO (_i+_d)&127 This is here to match results from pairhmm-1-base
		   pMM[bid]    = NUMBER(1.0) - ph2pr<NUMBER>(_i + _d & 127);
   		pMX[bid]    = ph2pr<NUMBER>(_i);
   		pXX[bid]    = ph2pr<NUMBER>(_c);
   		pMY[bid]    = (r == ROWS - 1) ? NUMBER(1.0) : ph2pr<NUMBER>(_d);

         const int _q = (_idc >> 21) & 127;
         distm_q[bid] = ph2pr<NUMBER>(_q);  // should use __ldg();
         _rs[bid] = rs[r-1];
   	}
      // initialize the zero-th col
      for (int bid = 1; bid <= BATCH; ++bid)
      {
		  M[bid] = NUMBER(0.0);
		  X[bid] = NUMBER(0.0);
		  Y[bid] = NUMBER(0.0);
      }

      // loop across all rows
      for (c = 1; c < COLS; c++)
      {

          //TODO read _hap from lmem, not gmem
          const char _hap = hap[c-1];

          // load the M,X,Y entries from the previous row of the previous stripe
          M[0] = M_r[c-1];
          X[0] = X_r[c-1];
          Y[0] = Y_r[c-1];
          M_r[c-1] = M[BATCH];
          X_r[c-1] = X[BATCH];
          Y_r[c-1] = Y[BATCH];

          // prefetch p[r][*] in registers (though the compiler might do it)

          // loop across columns in this stripe
          for (int bid = BATCH; bid >0; --bid)
          {

              const NUMBER distm = (_rs[bid] == _hap || _rs[bid] == 'N' || _hap == 'N') ?
                    NUMBER(1.0) - distm_q[bid] :
                                  distm_q[bid];

              Y[bid] = M[bid] * pMY[bid] + Y[bid] * (r == ROWS-1 ? 1.0 : pXX[bid]); 
              M[bid] = distm * (M[bid-1] * pMM[bid] + (X[bid-1] + Y[bid-1]) * (1.0 - pMX[bid]));
          }

          //TODO remove?
          M[0] = M_r[c];
          Y[0] = Y_r[c];

          for (int bid = 1; bid <= BATCH; bid++)
          {
  			    X[bid] = M[bid-1] * pMX[bid] + X[bid-1] * pXX[bid];
          }

          // add up the result
          if (ROWS-1-stripe <=BATCH) {
             result += X[ROWS-1-stripe] + M[ROWS-1-stripe];
          }
      }

    }

	output[pid] = log10(result) - LOG10_INITIAL_CONSTANT;
//   if (pid==12 && offsets[0].x==0) printf("output[%d] = log10(%e) - %e = %e\n", pid, result, LOG10_INITIAL_CONSTANT, output[pid]);
}
template <class NUMBER>
__global__ void
__launch_bounds__(WARP*4,8)
pairhmm_kernel( NUMBER init_const, NUMBER* M_in, NUMBER *X_in, NUMBER *Y_in,
                          char* rs, char* hap, NUMBER* q, 
                          int *n, cudaTextureObject_t t_n, int2* offset, int n_mats,
                          NUMBER* output, NUMBER log10_init) {
   NUMBER M_p;
   NUMBER M_pp;
   NUMBER X_p;
   NUMBER X_pp;
   NUMBER Y_p;
   NUMBER Y_pp;
   char _rs;
   NUMBER _q;
   NUMBER distm;
   NUMBER pMM;
   NUMBER pXX;
   NUMBER pMX;
   NUMBER pMY;
   NUMBER pGapM;
   NUMBER pYY;
   NUMBER M_loc, X_loc, Y_loc;
#ifdef __SHARED_READ
   __shared__ NUMBER M_top[128], X_top[128], Y_top[128];
#endif
   __shared__ NUMBER ph2pr_save[128];
   for (int z=threadIdx.x;z<128;z+=WARP) {
      ph2pr_save[z] = ph2pr<NUMBER>(z);
   }
   __syncthreads();

   int tid = threadIdx.x + blockDim.x * blockIdx.x;
   int wid = tid/WARP;
   tid %= WARP;
   if (wid > n_mats) return;
   int ROWS = offset[wid+1].x-offset[wid].x;
   int COLS = offset[wid+1].y-offset[wid].y;
   init_const /= COLS-1;
   NUMBER result=0.0;
#ifdef __LMEM
   NUMBER M[500], X[500], Y[500];
#else
   NUMBER *M=M_in+offset[wid].y+ROWS;
   NUMBER *X=X_in+offset[wid].y+ROWS;
   NUMBER *Y=Y_in+offset[wid].y+ROWS;
#endif
   rs+=offset[wid].x; 
   hap+=offset[wid].y; 
   q+=offset[wid].x;
   int n_off=offset[wid].x;
   n+=offset[wid].x;
   int z,r,c,bid;
//   if (wid==FOCUS && tid==0) printf("%s vs %s\n", rs, hap);
   for (int stripe = 1; stripe < ROWS; stripe+=WARP) 
   {
      r = stripe + tid;
	   _rs = rs[r-1];
	   _q = q[r-1];
#ifdef __USE_TEX
      int _n = tex1Dfetch<int>(t_n, r+n_off);
#else
      int _n = n[r];
#endif
      int _i = _n & 127;
      int _d = (_n >> 7) & 127;
      int _c = (_n >> 14) & 127;
      _q = ph2pr_save[(_n >> 21) & 127];
      pMM = 1.0 - ph2pr_save[_i+_d & 127];
      pXX = ph2pr_save[_c];
      pMX = ph2pr_save[_i];
      pMY = r == ROWS-1 ? 1.0 : ph2pr_save[_d];
      pGapM = 1.0 - pXX;
      pYY = r == ROWS-1 ? 1.0 : pXX;
      M_pp = M_p = X_p = X_pp = Y_p = Y_pp = 0.0;
      for (z=1; z<COLS+WARP; z++) {
         c = z - r + stripe;
#ifdef __SHARED_READ
         if (1==((c+r-stripe)%WARP)) {
            M_top[threadIdx.x] = M[c-1+2*r-2*stripe];
            X_top[threadIdx.x] = X[c-1+2*r-2*stripe];
            Y_top[threadIdx.x] = Y[c-1+2*r-2*stripe];
         }
#endif
         if (c>0 && c < COLS) {
			   char _hap = hap[c-1];
			   if (_rs == _hap || _rs == 'N' || _hap == 'N')
	   			distm = NUMBER(1.0) - _q;
            else distm = _q;
//            if (wid==FOCUS) printf("(%d,%d)  _hap: %c  _rs %c\n", r, c, _hap, _rs);
            if (c==1) {
               M_p=Y_p=X_p=0.0;
            }
            if (r-stripe==0) {
               if (stripe==1) {
                  M_pp = 0.0;
                  X_pp = 0.0;
                  Y_pp = init_const;
               } else {
#ifdef __SHARED_READ
                  M_pp = M_top[WARP*(wid%4)+(c-1)%WARP];
                  X_pp = X_top[WARP*(wid%4)+(c-1)%WARP];
                  Y_pp = Y_top[WARP*(wid%4)+(c-1)%WARP];
#else
                  M_pp = M[c-1];
                  X_pp = X[c-1];
                  Y_pp = Y[c-1];
//                  if ((wid==FOCUS || wid==FOCUS) && tid < 31 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): M: reading %e from M[%d]\n", r, c, M_pp, c-1);
//                  if ((wid==FOCUS || wid==FOCUS) && tid < 31 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): X: reading %e from X[%d]\n", r, c, X_pp, c-1);
//                  if ((wid==FOCUS || wid==FOCUS) && tid < 31 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): Y: reading %e from Y[%d]\n", r, c, Y_pp, c-1);
#endif
               }
            }
            M_loc = distm * (M_pp * pMM + (X_pp + Y_pp) * pGapM);
			   Y_loc = M_p * pMY + Y_p * pYY;
//            if ((wid==FOCUS || wid==FOCUS) && tid < 31 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): M: %e * (%e * %e + (%e + %e) * (1-%e)) = %e\n", r, c, distm, M_pp, pMM, X_pp, Y_pp, pXX, M_loc);
//            if ((wid==FOCUS || wid==FOCUS) && tid < 31 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): Y: %e * %e + %e * %e = %e\n", r, c, M_p, pMY, Y_p, (r == ROWS-1 ? 1.0 : pXX), Y_loc);
         } else {
            M_loc = Y_loc = 0.0;
         }
         M_p = __shfl_up(M_p,1);
         X_p = __shfl_up(X_p,1);
         Y_p = __shfl_up(Y_p,1);
         if (r-stripe==0) 
         {
            if (stripe==1) {
               M_p = 0.0;
               X_p = 0.0;
            } else {
#ifdef __SHARED_READ
               M_p = M_top[WARP*(wid%4) + c%WARP];
               X_p = X_top[WARP*(wid%4) + c%WARP];
               Y_p = Y_top[WARP*(wid%4) + c%WARP];
#else
               M_p = M[c]; 
               X_p = X[c];
               Y_p = Y[c];
#endif
//               if ((wid==FOCUS || wid==FOCUS) && tid < 31 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): M: reading %e from M[%d]\n", r, c, M_p, c);
//               if ((wid==FOCUS || wid==FOCUS) && tid < 31 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): X: reading %e from X[%d]\n", r, c, X_p, c);
            }
         }
         if (c > 0 && c < COLS)
         {
			   X_loc = M_p * pMX + X_p * pXX;
//            if ((wid==FOCUS || wid==FOCUS) && tid < 31 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): X: %e * %e + %e * %e = %e\n", r, c, M_p, pMX, X_p, pXX, X_loc);
         }
         M_pp = M_p; M_p = M_loc;
         X_pp = X_p; X_p = X_loc;
         Y_pp = Y_p; Y_p = Y_loc;
         if (r-stripe==WARP-1 && c >=0 && c < COLS) {
            M[c] = M_loc;
            X[c] = X_loc;
            Y[c] = Y_loc;
//            if ((wid==FOCUS || wid==FOCUS) && tid < 31 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): X: writing %e to X[%d]\n", r, c, X_loc, c);
         }
         if (r==ROWS-1 && c > 0 && c < COLS) {
            result += M_loc + X_loc;
//            if ((wid==FOCUS || wid==FOCUS) && offset[0].x==0 && c > 0 && c < COLS) printf("(%d,%d): result += %e + %e = %e\n", r, c, M_loc, X_loc, result);
         }
      }
   }
   if (r==ROWS-1)
   {
      if (result < MIN_ACCEPTED) {
         output[wid] = 1.0;
      } else {
	      output[wid] = log10(result) - log10_init;
      }
   }
}
template <class NUMBER>
__global__ void
__launch_bounds__(WARP*4,10)
pairhmm_kernel_wrap( NUMBER init_const, NUMBER* M_in, NUMBER *X_in, NUMBER *Y_in,
                          char* rs, char* hap, NUMBER* q, 
                          int *n, cudaTextureObject_t t_n, int2* offset, int n_mats,
                          NUMBER* output, NUMBER log10_init) {
   NUMBER M_p;
   NUMBER M_pp;
   NUMBER X_p;
   NUMBER X_pp;
   NUMBER Y_p;
   NUMBER Y_pp;
   char _rs, _rs_next;
   NUMBER _q, _q_next;
   NUMBER distm;
   NUMBER pMM, pMM_next;
   NUMBER pXX, pXX_next;
   NUMBER pMX, pMX_next;
   NUMBER pMY, pMY_next;
   NUMBER pYY, pYY_next;
   NUMBER pGapM, pGapM_next;
   NUMBER M_loc, X_loc, Y_loc;
   NUMBER Y0;
#ifdef __SHARED_READ
   __shared__ NUMBER M_top[128], X_top[128], Y_top[128];
#endif
   __shared__ NUMBER ph2pr_save[128];
   for (int z=threadIdx.x;z<128;z+=WARP) {
      ph2pr_save[z] = ph2pr<NUMBER>(z);
   }
   __syncthreads();

   int tid = threadIdx.x + blockDim.x * blockIdx.x;
   int wid = MAT_PER_WARP*(tid/WARP);
   if (wid > n_mats) return;
#ifdef __LMEM
   NUMBER M[500], X[500], Y[500];
#else
   NUMBER *M=M_in+offset[wid].y;
   NUMBER *X=X_in+offset[wid].y;
   NUMBER *Y=Y_in+offset[wid].y;
#endif
   rs+=offset[wid].x; 
   hap+=offset[wid].y; 
   q+=offset[wid].x;
   int n_off=offset[wid].x;
   n+=offset[wid].x;
   int mid = -1;
   int r=tid%WARP;
   int c=1-r;
   int stripe=1;
   int ROWS=0, COLS=2, nCOLS=0, nROWS=0;
   NUMBER result=0.0;
   nROWS = offset[wid+1].x-offset[wid].x;
   nCOLS = offset[wid+1].y-offset[wid].y;
   {
      int nextr = r+1;
	   _rs_next = rs[nextr-1];
#ifdef __USE_TEX
      int _n = tex1Dfetch<int>(t_n, nextr+n_off);
#else
      int _n = n[nextr];
#endif
      int _i = _n & 127;
      int _d = (_n >> 7) & 127;
      int _c = (_n >> 14) & 127;
      _q_next = ph2pr_save[(_n >> 21) & 127];
//      if (MAT_PER_WARP*((threadIdx.x +blockIdx.x * blockDim.x)/WARP) == FOCUS) printf("threadIdx.x: %d. tid:%d. _q_next = %f\n", tid, tid%WARP, _q_next);
      pMM_next = 1.0 - ph2pr_save[_i+_d & 127];
      pXX_next = ph2pr_save[_c];
      pMX_next = ph2pr_save[_i];
      pMY_next = r == ROWS-1 ? 1.0 : ph2pr_save[_d];
   }
   while(mid < MAT_PER_WARP) 
   {
      while (r>=ROWS) {
         mid++;
         tid = threadIdx.x + blockDim.x * blockIdx.x;
         wid = tid/WARP;
         wid = wid*MAT_PER_WARP + mid;
         tid %= WARP;
         if (wid >=n_mats) return;
         stripe=1;
         r-=ROWS-1;
         //r = tid+1;
         rs+=ROWS;
         hap+=COLS;
         q+=ROWS;
         n_off+=ROWS;
         n+=ROWS;
         ROWS = nROWS;
         if (nCOLS != COLS) init_const *= (1.0*(COLS-1))/(nCOLS-1);
         COLS = nCOLS;
         nCOLS = offset[wid+2].y-offset[wid+1].y;
         nROWS = offset[wid+2].x-offset[wid+1].x;
         //TODO resolve cases where more than 2 matrices are found in a single set of 32 rows
         //nROWS = offset[wid+2].x-offset[wid+1].x;
         //if (tid + ROWS + nROWS - r  < WARP) {
         //   nCOLS = max(nCOLS,offset[wid+3].y-offset[wid+2].y);
         //   nROWS = offset[wid+3].x-offset[wid+2].x;
         //}
//         if (wid==FOCUS || wid == FOCUS-1) printf("thread: %d, block %d, wid = %d. r= %d. n_off = %d, q = %p, ROWS= %d, COLS = %d, rs = %s, hap(%p) = %s\n", threadIdx.x, blockIdx.x, wid, r, n_off, q, ROWS, COLS, rs, (void*)hap, hap);
         result=0.0;
      }
      if (c==1+WARP-tid) { //All threads do this computation at the same time.
         //Compute the next set of q,p and rs
         int nextr = r+WARP;
         //TODO What if we fit more than 3 mats in a single batch of 32 rows?
         if (nextr >= ROWS) nextr++; //skip one for the string termination char (\0)
         if (nextr >= ROWS+nROWS) nextr++;
	      _rs_next = rs[nextr-1];
#ifdef __USE_TEX
         int _n = tex1Dfetch<int>(t_n, nextr+n_off);
#else
         int _n = n[nextr];
#endif
         int _i = _n & 127;
         int _d = (_n >> 7) & 127;
         int _c = (_n >> 14) & 127;
         _q_next = ph2pr_save[(_n >> 21) & 127];
//         if (wid==FOCUS || wid==FOCUS-1) printf("tid:%d. new _q_next = %f\n", tid, _q_next);
         pMM_next = 1.0 - ph2pr_save[_i+_d & 127];
         pXX_next = ph2pr_save[_c];
         pMX_next = ph2pr_save[_i];
         pMY_next = r == ROWS-1 ? 1.0 : ph2pr_save[_d];
         pYY_next = r == ROWS-1 ? 1.0 : pXX;
         pGapM_next = 1.0-pXX;
      }
      if (c==1) {
         //First time on this row, recompute pMM, pXX, pMX, pMY and _q
         //TODO don't do this every time
	      _rs = _rs_next;
	      _q = _q_next;
         pMM = pMM_next;
         pXX = pXX_next;
         pMX = pMX_next;
         pMY = pMY_next;
         M_pp = M_p = X_p = X_pp = Y_p = Y_pp = 0.0;
         pYY = pYY_next;
         pGapM = pGapM_next;
//         if (wid==FOCUS) printf("row: %d tid: %d mid: %d. Switching to _q_next (%e)\n", r, tid, mid, _q);
      }
#ifdef __SHARED_READ
      if (1==((c+r-stripe)%WARP)) {
         M_top[threadIdx.x] = M[c-1+2*r-2*stripe];
         X_top[threadIdx.x] = X[c-1+2*r-2*stripe];
         Y_top[threadIdx.x] = Y[c-1+2*r-2*stripe];
      }
#endif
      if (c>0 && c < COLS) {
		  char _hap = hap[c-1];
		  if (_rs == _hap || _rs == 'N' || _hap == 'N')
	        distm = NUMBER(1.0) - _q;
        else distm = _q;
        if (r==1 && c==1) {
           Y_pp = init_const;
        } 
		  Y_loc = M_p * pMY + Y_p * pYY;
        M_loc = distm * (M_pp * pMM + (X_pp + Y_pp) * pGapM );
//        if ((wid==FOCUS || wid==FOCUS) && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): M: %e * (%e * %e + (%e + %e) * (1-%e)) = %e\n", r,c, distm, M_pp, pMM, X_pp, Y_pp, pXX, M_loc);
//        if ((wid==FOCUS || wid==FOCUS) && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): Y: %e * %e + %e * %e = %e\n", r,c, M_p, pMY, Y_p, (r==ROWS-1 ? 1.0 : pXX), Y_loc);
      } else {
         //TODO do we need this?
         //M_loc = Y_loc = X_loc = 0.0;
      }
      M_p = __shfl_up(M_p,1);
      X_p = __shfl_up(X_p,1);
      Y_p = __shfl_up(Y_p,1);
      if (r == 1)
      { 
         M_p = X_p = 0.0;
         Y_p = init_const;
      } else if (tid==0) {
#ifdef __SHARED_READ
         M_p = M_top[WARP*(wid%4)+c%WARP];
         X_p = X_top[WARP*(wid%4)+c%WARP];
         Y_p = Y_top[WARP*(wid%4)+c%WARP];
#else
         M_p = M[c];
         X_p = X[c];
         Y_p = Y[c];
#endif
      }
		X_loc = M_p * pMX + X_p * pXX;
//      if ((wid==FOCUS || wid==FOCUS) && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): X: %e * %e + %e * %e = %e\n", r,c, M_p, pMX, X_p, pXX, X_loc);
      M_pp = M_p; M_p = M_loc;
      X_pp = X_p; X_p = X_loc;
      Y_pp = Y_p; Y_p = Y_loc;
      if (tid==WARP-1 && c >=0 && c < COLS) {
         M[c] = M_loc;
         X[c] = X_loc;
         Y[c] = Y_loc;
//         if ((wid==FOCUS || wid==FOCUS) && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): writing %e to Y[%d]\n", r, c, Y_loc, c);
      }
      if (r==ROWS-1 && c > 0 && c < COLS) {
         result += M_loc + X_loc;
//         if ((wid==FOCUS || wid==FOCUS) && offset[0].x==0 && c > 0 && c < COLS) printf("(%d,%d): result += %e + %e = %e\n", r, c, M_loc, X_loc, result);
         if (c == COLS-1) {
	         output[wid] = log10(result) - log10_init;
         }
//         if (c == COLS-1 && (wid==FOCUS || wid==FOCUS)) printf("thread: %d, block %d, wid = %d. writing log10(%e) - %e = %e\n", threadIdx.x, blockIdx.x, wid, result, log10_init, output[wid]);
      }
      c++;
      if (c > COLS-1 && c > nCOLS-1 && c >= WARP) { //TODO You don't always have to wait for nCOLS
         //if we're at the end of a column, move down WARP rows
//         if (wid==FOCUS || wid==FOCUS) printf("Row %d (tid %d) reached the end. Changing to row %d (mid=%d)\n", r, tid, r+WARP, mid);
         stripe+=WARP;
         r+=WARP;
         c=1;
      }
   }  //mid
}
template <class NUMBER, unsigned int BATCH>
__global__ void
//__launch_bounds__(WARP*4,3)
pairhmm_kernel_onethread( NUMBER init_const, NUMBER* M_in, NUMBER *X_in, NUMBER *Y_in,
                          char* rs, char* hap, NUMBER* q, 
                          int *n, cudaTextureObject_t t_n, int2* offset, int n_mats,
                          NUMBER* output, NUMBER log10_init) {
   NUMBER M_p[BATCH];
   NUMBER M_pp[BATCH];
   NUMBER X_p[BATCH];
   NUMBER X_pp[BATCH];
   NUMBER Y_p[BATCH];
   NUMBER Y_pp[BATCH];
   char _rs[BATCH];
   NUMBER _q[BATCH];
   NUMBER distm;
   NUMBER pMM[BATCH];
   NUMBER pXX[BATCH];
   NUMBER pMX[BATCH];
   NUMBER pMY[BATCH];
   NUMBER M_loc, X_loc, Y_loc;

   int tid = threadIdx.x + blockDim.x * blockIdx.x;
   if (tid > n_mats) return;
   int ROWS = offset[tid+1].x-offset[tid].x;
   int COLS = offset[tid+1].y-offset[tid].y;
   //if (tid==11) printf("Dimensions : %d x %d\n", ROWS, COLS);
   init_const /= COLS-1;
   NUMBER result=0.0;
   //Make this lmem, not gmem
   NUMBER M[350], X[350], Y[350];
   //M+=offset[3*tid];
   //X+=offset[3*tid];
   //Y+=offset[3*tid];
   rs+=offset[tid].x; 
   hap+=offset[tid].y; 
   q+=offset[tid].x;
   int n_off=3*offset[tid].x;
   n+=3*offset[tid].x;
   int z,r,c,bid;
   for (int stripe = 1; stripe < ROWS; stripe+=BATCH) 
   {
      #pragma unroll
      for (bid=0;bid<BATCH;bid++) {
         r = stripe + bid;
	      _rs[bid] = rs[r-1];
	      _q[bid] = q[r-1];
         int _i = n[r] & 127;
         int _d = (n[r] >> 7) & 127;
         int _c = (n[r] >> 14) & 127;
         //_q = ph2pr<NUMBER>((n[r] >> 21) & 127);
         
 //        if(tid==12 && c < 12) printf("%d: q=%e\n", bid, _q[bid]);
         pMM[bid] = 1.0 - ph2pr<NUMBER>(_i+_d & 127);
         pXX[bid] = ph2pr<NUMBER>(_c);
         pMX[bid] = ph2pr<NUMBER>(_i);
         pMY[bid] = r == ROWS-1 ? 1.0 : ph2pr<NUMBER>(_d);
         M_pp[bid] = M_p[bid] = X_p[bid] = X_pp[bid] = Y_p[bid] = Y_pp[bid] = 0.0;
      }
      for (z=1; z<COLS+BATCH; z++) {
         #pragma unroll
         for (bid=BATCH-1;bid>=0;bid--) {
            r = stripe + bid;
            c = z - bid;
            if (c<1) continue;
			   char _hap = hap[c-1];
			   if (_rs[bid] == _hap || _rs[bid] == 'N' || _hap == 'N')
			   	distm = NUMBER(1.0) - _q[bid];
            else distm = _q[bid];
            if (bid==0) {
               if (stripe==1) {
                  M_loc = distm * init_const * (1.0 - pXX[bid]) ;
			         Y_loc = M_p[bid] * pMY[bid] + Y_p[bid] * (r == ROWS-1 ? 1.0 : pXX[bid]);
			         X_loc = 0.0;
 //                 if (tid==12 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): %e * (%e * %e + (%e + %e) * (1-%e) = %e\n", r,c, distm, 0.0, pMM[bid], 0.0, init_const, pXX[bid], M_loc);
               } else {
                  M_loc = distm * (M[c-1] * pMM[bid] + (X[c-1] + Y[c-1]) * (1.0 - pXX[bid]) );
			         Y_loc = M_p[bid] * pMY[bid] + Y_p[bid] * (r == ROWS-1 ? 1.0 : pXX[bid]);
			         X_loc = M[c] * pMX[bid] + X[c] * pXX[bid];
 //                 if (tid==12 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): %e * (%e * %e + (%e + %e) * (1-%e) = %e\n", r,c, distm, M[c-1], pMM[bid], X[c-1], Y[c-1], pXX[bid], M_loc);
               }
            } else { 
               M_loc = distm * (M_pp[bid-1] * pMM[bid] + (X_pp[bid-1] + Y_pp[bid-1]) * (1.0 - pXX[bid]) );
			      Y_loc = M_p[bid] * pMY[bid] + Y_p[bid] * (r == ROWS-1 ? 1.0 : pXX[bid]);
			      X_loc = M_p[bid-1] * pMX[bid] + X_p[bid-1] * pXX[bid];
 //              if (tid==12 && offset[0].x==0 && c > 0 && c < 12) printf("(%d,%d): %e * (%e * %e + (%e + %e) * (1-%e) = %e\n", r,c, distm, M_pp[bid], pMM[bid], X_pp[bid], Y_pp[bid], pXX[bid], M_loc);
            }
            if (bid==BATCH-1) {
               M[c] = M_loc;
               X[c] = X_loc;
               Y[c] = Y_loc;
            }
            if (r==ROWS-1) {
               result += M_loc + X_loc;
            }
            M_pp[bid] = M_p[bid]; M_p[bid] = M_loc;
            X_pp[bid] = X_p[bid]; X_p[bid] = X_loc;
            Y_pp[bid] = Y_p[bid]; Y_p[bid] = Y_loc;
         } 
      }
   }
   output[tid]=log10(result)-log10_init;
}
#if 1
//TODO Alloc more memory in M,X,Y due to smaller "warps"
//TODO Xc0 and Yr0 are constants
//TODO compute and save result
template <class NUMBER, unsigned int BATCH>
__global__ void
//__launch_bounds__(WARP*4,3)
pairhmm_kernel_onethread_old( NUMBER init_const, NUMBER* M_in, NUMBER *X_in, NUMBER *Y_in,
                          char* rs, char* hap, NUMBER* q, 
                          int *n, cudaTextureObject_t t_n, int2* offset, int n_mats,
                          NUMBER* output, NUMBER log10_init) {
   NUMBER M_p[BATCH+1];
   NUMBER M_pp[BATCH+1];
   NUMBER X_p[BATCH+1];
   NUMBER X_pp[BATCH+1];
   NUMBER Y_p[BATCH+1];
   NUMBER Y_pp[BATCH+1];
   char _rs[BATCH+1];
   NUMBER _q[BATCH+1];
   NUMBER distm;
   NUMBER pMM[BATCH+1];
   NUMBER pXX[BATCH+1];
   NUMBER pMX[BATCH+1];
   NUMBER pMY[BATCH+1];
   NUMBER M_loc[BATCH+1], X_loc[BATCH+1], Y_loc[BATCH+1];
   int r, c;
   
   int tid = threadIdx.x + blockDim.x * blockIdx.x;
   if (tid > n_mats) return;
   int ROWS = offset[tid+1].x-offset[tid].x;
   int COLS = offset[tid+1].y-offset[tid].y;
   if (tid==12) printf("Dimensions : %d x %d\n", ROWS, COLS);
   NUMBER result=0.0;
   //Make this lmem, not gmem
   NUMBER M[250], X[250], Y[250];
   //M+=offset[3*tid];
   //X+=offset[3*tid];
   //Y+=offset[3*tid];
   rs+=offset[tid].x; 
   hap+=offset[tid].y; 
   q+=offset[tid].x;
   int n_off=3*offset[tid].x;
   n+=3*offset[tid].x;
   for (int stripe = 0; stripe < ROWS; stripe+=BATCH) 
   {
   M_pp[0] = 0.0;
   M_pp[1] = 0.0;
   X_pp[0] = 0.0;
   X_pp[1] = 0.0;
   Y_pp[0] = init_const/(COLS-1);
   Y_pp[1] = init_const/(COLS-1);
   M_p[0] = 0.0;
   M_p[1] = 0.0;
   X_p[1] = 0.0;
   Y_p[1] = 0.0;
   if (stripe==0) {
      X_p[0] = 0.0;
      Y_p[0] = init_const/(COLS-1);
   } else {
      X_p[0] = X[1]; 
      Y_p[0] = Y[1]; 
   }
   #pragma unroll 
   for (int bid = 1; bid < BATCH+1; bid++)
   {
      r = bid+stripe;
	   _rs[bid] = rs[bid-1+stripe];
	   _q[bid] = q[bid-1+stripe];
      int _i = n[r] & 127;
      int _d = (n[r] >> 7) & 127;
      int _c = (n[r] >> 14) & 127;
      //_q[bid] = ph2pr<NUMBER>((n[r] >> 21) & 127);
      
      pMM[bid] = 1.0 - ph2pr<NUMBER>(_i+_d & 127);
      pXX[bid] = ph2pr<NUMBER>(_c);
      pMX[bid] = ph2pr<NUMBER>(_i);
      pMY[bid] = r == ROWS-1 ? 1.0 : ph2pr<NUMBER>(_d);
   }
   for (int z = 0; z < COLS+BATCH+1; z++) 
   {
      #pragma unroll 
      for (int bid = 0; bid < BATCH+1; bid++)
      {
         r = bid+stripe;
         c = z-bid+1;
         if (c > 0 && r < ROWS && c < COLS) 
         {  
            if(tid==12 && offset[0].x==0 && c < 12) printf("(%d,%d): ", r, c);
			   char _hap;
            if (c>0) _hap = hap[c-1];
            else _hap = hap[c];
			   if (_rs[bid] == _hap || _rs[bid] == 'N' || _hap == 'N')
			   	distm = NUMBER(1.0) - _q[bid];
            else distm = _q[bid];
            if (r==0) {
               X_p[0] = 0.0;
               Y_p[0] = init_const/(COLS-1);
               M_p[0] = 0.0;
            } else if (c==0) {
               X_p[0] = M_p[0] = X_pp[0] = M_pp[0] = 0.0;
               Y_p[0] = Y_pp[0] = 0.0;
            } else if (bid == 0 && c < 250) {
               if (tid==0 && offset[0].x==0 && c < 12) printf("bid:%d, z:%d  reading %e from column %d\n", bid, z, M[c], c);
               X_p[0] = X[c];
               Y_p[0] = Y[c];
               M_p[0] = M[c];
            }
            M_loc[bid] = distm * (M_pp[bid] * pMM[bid] + (X_pp[bid] + Y_pp[bid]) * (1.0 - pXX[bid]) );
			   Y_loc[bid] = M_p[bid] * pMY[bid] + Y_p[bid] * (r == ROWS-1 ? 1.0 : pXX[bid]);
			   X_loc[bid] = M_p[bid-1] * pMX[bid] + X_p[bid-1] * pXX[bid];
            if (tid==12 && offset[0].x==0 && c > 0 && c < 12) printf("%e * (%e * %e + (%e + %e) * (1-%e) = %e\n", distm, M_pp[bid], pMM[bid], X_pp[bid], Y_pp[bid], pXX[bid], M_loc[bid]);
         }
         if (c < COLS && r-stripe==BATCH-1 && c < 250) {
            if (tid==0 && offset[0].x==0 && c < 12) printf("(%d,%d): writing %e to column %d\n", r, c, M_loc[bid], c);
            M[c] = M_loc[bid];
            X[c] = X_loc[bid];
            Y[c] = Y_loc[bid];
         }
         if (r==ROWS-1 && c > 0 && c < COLS ) {
            result += M_loc[bid] + X_loc[bid]; 
            if (tid==12 && offset[0].x==0 && c < 12) printf("(%d,%d): result += %e + %e = %e\n", r, c, M_loc[bid], X_loc[bid], result);
         }
      }
      #pragma unroll 
      for (int bid = BATCH; bid > 0; bid--)
      {
         if (bid <= z+1 && bid+stripe < ROWS && z-bid+1 < COLS) 
         {  
            r = bid+stripe;
            c = z-bid+1;
            M_pp[bid]=M_p[bid-1];
    //        if (tid==0) printf("putting M_p[%d] (%e) into M_pp[%d]\n", bid-1, M_p[bid-1], bid);
            X_pp[bid]=X_p[bid-1];
            Y_pp[bid]=Y_p[bid-1];
            if (c==0) {
               X_p[bid]=0.0;
               M_p[bid]=0.0;
               Y_p[bid]=0.0;
            } else {
               X_p[bid] = X_loc[bid];
               M_p[bid] = M_loc[bid];
    //           if (tid==0) printf("putting M_loc[%d] (%e) into M_p[%d]\n", bid, M_loc[bid], bid);
               Y_p[bid] = Y_loc[bid];
            }
         } else if (c==0 && r < ROWS) {
            M_pp[bid]=M_p[bid-1];
            X_pp[bid]=X_p[bid-1];
            Y_pp[bid]=Y_p[bid-1];
            //Read from col 0
            X_p[bid] = 0.0;
            M_p[bid] = 0.0;
            Y_p[bid] = 0.0;
         }
          
      }
   }
   
   }
   if (tid==12 && offset[0].x==0)printf("output[%d] = log10(%e) - %e = %e\n", tid, result, log10_init, log10(result)-log10_init);
   output[tid]=log10(result)-log10_init;
}
#endif
__global__ void CPU_start() {}
__global__ void CPU_end() {}
template <unsigned int u>
__global__ void debugMark() {}

template<class NUMBER>
__global__ void
__launch_bounds__(128,9)
repeatedSine(NUMBER* in_out, int size) {
   int bid = blockIdx.x;
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   __shared__ NUMBER shm[256];

   if (tid > size) return;
 
   shm[threadIdx.x] = in_out[tid];
   __syncthreads();

   NUMBER foo = in_out[tid];
   if ((threadIdx.x > 0 && shm[threadIdx.x] > shm[threadIdx.x-1]) && 
       ( threadIdx.x < 127 && shm[threadIdx.x] > shm[threadIdx.x+1]))
      foo = -in_out[tid];
   for (int q=0;q<20;q++) foo = sin(foo);
   if (bid%2==0) for (int q=0;q<20;q++) foo = log(foo);
   else for (int q=0;q<20;q++) foo = sin(foo);
   
   in_out[tid] = foo;
}

template<class NUMBER>
__global__ void 
__launch_bounds__(128,9)
pairhmm_kernel_old( NUMBER init_const, NUMBER* M, NUMBER *X, NUMBER *Y, 
                               char* rs, char* hap, NUMBER* q, 
                               int* n, cudaTextureObject_t t_n, int* offset, int n_mats,
                               NUMBER* output, NUMBER log10_init) {
   NUMBER M_p, M_pp, X_p, X_pp, Y_p, Y_pp, distm, pMM, pGapM, pXX, 
          pMX, pYY, pMY, M_loc, X_loc, Y_loc;
#ifdef __COALESCED_WRITE
   NUMBER M_bottom, X_bottom, Y_bottom;
#endif
#ifdef __COALESCED_READ
   NUMBER M_top, X_top, Y_top;
#endif
#ifdef __SHARED_READ
   __shared__ NUMBER M_top[128], X_top[128], Y_top[128];
#endif
#ifdef __SHARED_WRITE
   __shared__ NUMBER M_bottom[128], X_bottom[128], Y_bottom[128];
#endif
   char _rs;
   NUMBER _q;
   int tid = threadIdx.x;
   int wid = blockIdx.x;
   //Fake wid,tid to allow for larger blocks
   tid = (threadIdx.x + blockDim.x * blockIdx.x)%WARP;
   wid = (threadIdx.x + blockDim.x * blockIdx.x)/WARP;
   if (wid > n_mats) return;
   int ROWS = offset[3*wid+4]-offset[3*wid+1];
   int COLS = offset[3*wid+5]-offset[3*wid+2];
   NUMBER result=0.0;
   M+=offset[3*wid];
   X+=offset[3*wid];
   Y+=offset[3*wid];
   rs+=offset[3*wid+1];
   hap+=offset[3*wid+2];
   q+=offset[3*wid+1];
   //int n_off=3*offset[3*wid+1];
   n+=3*offset[3*wid+1];
   int r;
   for (int stripe = 0; stripe < ROWS-1; stripe+=WARP-1) {
      r = stripe + tid;
      if ( r < 2) {
         M_pp=0.0;
         X_pp = 0.0;//X_pp=Xc0[0];
         Y_pp=init_const/(COLS-1);
         M_p=0.0;
         if (r==1) X_p = 0.0;//X_p=Xc0[1];
         else X_p=0.0;
         if (r==0) Y_p=init_const/(COLS-1);
         else Y_p=0.0;
      } else if (r-stripe<2) {
         M_pp = 0.0;
         X_pp = 0.0; //X_pp = Xc0[stripe]; //X[stripe][0];
         Y_pp = 0.0;
         M_p = 0.0;
         if (r-stripe==1) X_p=0.0; //X_p=Xc0[1];
         else X_p = X[stripe/(WARP-1)*COLS+1]; //M[stripe+tid][1-tid];
         if (r-stripe==0) Y_p = Y[stripe/(WARP-1)*COLS+1]; //M[stripe+tid][1-tid];
         else Y_p = 0.0;
      }
      //TODO pad?
      if (r>0) {
      _rs = rs[r-1];
	   _q = q[r-1];
      } else {
      _rs = rs[r];
	   _q = q[r];
      }
      //TODO transpose n for coalesced reads?
      //TODO pass n as a single integer since only the first 7 bits of each matter
#ifdef DEBUG
      if (wid==1 && 
              tex1Dfetch<int>(t_n, 3*r+DD+n_off) != n[3*r+DD]) {
                   //printf("mismatch : stripe = %d, tid = %d: D = %d, I = %d, C = %d\nFrom (int*)n: D = %d, I = %d, C = %d\n", stripe, tid,
                   //             tex1Dfetch<int>(t_n, 3*r+DD+n_off),
                   //             tex1Dfetch<int>(t_n, 3*r+II+n_off),
                   //             tex1Dfetch<int>(t_n, 3*r+CC+n_off),
                   //             n[3*r+DD],
                   //             n[3*r+II],
                   //             n[3*r+CC]);
      } else if(wid==1) printf("Match : stripe = %d, tid = %d\n", stripe, tid);
#endif
      int _n = n[r];
      int _i = _n & 127;
      int _d = (_n >> 7) & 127;
      int _c = (_n >> 14) & 127;
      //_q = ph2pr<NUMBER>((_n >> 21) & 127);
      pMM = 1.0 - ph2pr<NUMBER>(_i+_d & 127);
      pXX = ph2pr<NUMBER>(_d);
      pMX = ph2pr<NUMBER>(_i);
      pMY = r == ROWS-1 ? 1.0 : ph2pr<NUMBER>(_d);

      pGapM = 1.0 - pXX;
      pYY = r == ROWS-1 ? 1.0 : pXX;
      //TODO get rid of this?
      if (r==0) {
         pMM = pXX = pGapM = pMX = pYY = pMY = 0.0;
      }
      for (int z = 1; z < COLS+WARP+1;z++) 
      {
         int c = z-r+stripe+1;
         //TODO align at word boundaries
#ifdef __COALESCED_READ
         //TODO write/read to/from shared memory instead?
         //TODO don't read if stripe==0
         if (1==z%WARP) {
            M_top = M[stripe/(WARP-1)*COLS+z+r-stripe];
            X_top = X[stripe/(WARP-1)*COLS+z+r-stripe];
            Y_top = Y[stripe/(WARP-1)*COLS+z+r-stripe];
            //if (wid==0) printf("%d. tid: %d. M_top = %e\n", z, tid, M_top);
         } else {
            M_top = __shfl_down(M_top,1);
            X_top = __shfl_down(X_top,1);
            Y_top = __shfl_down(Y_top,1);
         }
#endif
#ifdef __SHARED_READ
         if (1==(c+r-stripe-1)%WARP) {
            M_top[threadIdx.x] = M[stripe/(WARP-1)*COLS+c+2*r-2*stripe-1];
            X_top[threadIdx.x] = X[stripe/(WARP-1)*COLS+c+2*r-2*stripe-1];
            Y_top[threadIdx.x] = Y[stripe/(WARP-1)*COLS+c+2*r-2*stripe-1];
            //__syncthreads(); No need for sync. Warps maintain synchronization
         }
#endif
      

         if (0 <= c && r < ROWS && c <= COLS)
         {
            //TODO pad instead
            if (c>0) {
            char _hap = hap[c-1];
			   if (_rs == _hap || _rs == 'N' || _hap == 'N')
			   	distm = NUMBER(1.0) - _q;
            else distm = _q;
            }
            else distm = _q;
            if (r == 0 && stripe == 0) {
               X_p = 0.0; 
               Y_p = init_const/(COLS-1);
               M_p = 0.0;
            } else if (r-stripe == 0 && c+r-stripe-1 > 1) {
#ifdef __COALESCED_READ
               M_p = M_top;
               X_p = X_top;
               Y_p = Y_top;
#else
#ifdef __SHARED_READ
               M_p = M_top[threadIdx.x + (c+r-stripe-2)%WARP];
               X_p = X_top[threadIdx.x + (c+r-stripe-2)%WARP];
               Y_p = Y_top[threadIdx.x + (c+r-stripe-2)%WARP];
#else
               M_p = M[stripe/(WARP-1)*COLS+z];
               X_p = X[stripe/(WARP-1)*COLS+z];
               Y_p = Y[stripe/(WARP-1)*COLS+z];
#endif
#endif
            } 
            M_loc = distm * (M_pp * pMM + X_pp * pGapM + Y_pp * pGapM);
			   Y_loc = M_p * pMY + Y_p * pYY;
            //if (wid==12 && c > 0 && c < 12 && r > 0) printf("(%d,%d): %e * (%e * %e + (%e + %e) * (1-%e) = %e\n", r, c, distm, M_pp, pMM, X_pp, Y_pp, pXX, M_loc);
            M_p = __shfl_up(M_p,1);
            Y_p = __shfl_up(Y_p,1);
            X_p = __shfl_up(X_p,1);
			   X_loc = M_p * pMX + X_p * pXX;
            M_pp = M_p;
            X_pp = X_p;
            Y_pp = Y_p;
            if (c == 0 && stripe==0) {
               M_p = 0.0;
               Y_p = 0.0;
               X_p = 0.0; //X_p = Xc0[tid];
            } else if (c == 0) {
               M_p = 0.0;
               Y_p = 0.0;
               X_p = 0.0; //X_p = Xc0[tid+stripe]; //X[tid+stripe][0]
            } else {
               M_p = M_loc;
               X_p = X_loc;
               Y_p = Y_loc;
            }
            if (c < COLS && r-stripe==WARP-1) {
#ifdef __COALESCED_WRITE
               M_bottom = M_loc;
               X_bottom = X_loc;
               Y_bottom = Y_loc;
#else
#ifdef __SHARED_WRITE
               M_bottom[(threadIdx.x/WARP)*WARP + c%WARP] = M_loc;
               X_bottom[(threadIdx.x/WARP)*WARP + c%WARP] = X_loc;
               Y_bottom[(threadIdx.x/WARP)*WARP + c%WARP] = Y_loc;
#else
               //Vanilla
               M[((r+WARP-2)/(WARP-1))*COLS+c] = M_loc;
               X[((r+WARP-2)/(WARP-1))*COLS+c] = X_loc;
               Y[((r+WARP-2)/(WARP-1))*COLS+c] = Y_loc;
#endif
#endif
            }
            if (r==ROWS-1) { 
               result += M_loc + X_loc;
               //if (wid==12 && offset[0]==0) printf("(%d,%d): result += %e + %e = %e\n", r, c, M_loc, X_loc, result);
            }
         }
#ifdef __COALESCED_WRITE
         if (z%WARP==0 || z==COLS+WARP-6) {
            int rt = stripe+WARP-1;
            int ct = z-WARP+2;
            M[((rt+WARP-2)/(WARP-1))*COLS+ct-WARP+tid+1] = M_bottom;
            X[((rt+WARP-2)/(WARP-1))*COLS+ct-WARP+tid+1] = X_bottom;
            Y[((rt+WARP-2)/(WARP-1))*COLS+ct-WARP+tid+1] = Y_bottom;
         }
         M_bottom = __shfl_down(M_bottom, 1);
         X_bottom = __shfl_down(X_bottom, 1);
         Y_bottom = __shfl_down(Y_bottom, 1);
#endif
#ifdef __SHARED_WRITE
         if (stripe + WARP < ROWS && z>WARP+2 && (z+3)%WARP==0) {
            int rt = stripe+WARP-1;
            int ct = z-WARP+2;
            M[((rt+WARP-2)/(WARP-1))*COLS+ct-WARP+tid+1] = M_bottom[threadIdx.x];
            X[((rt+WARP-2)/(WARP-1))*COLS+ct-WARP+tid+1] = X_bottom[threadIdx.x];
            Y[((rt+WARP-2)/(WARP-1))*COLS+ct-WARP+tid+1] = Y_bottom[threadIdx.x];
         }
         if (stripe + WARP < ROWS && (z-WARP+2)%WARP >= tid && z-WARP+2==COLS-1) {
            int rt = stripe+WARP-1;
            int ct = z-WARP+2;
            M[((rt+WARP-2)/(WARP-1))*COLS+ct+tid-(z-WARP+2)%WARP] = M_bottom[threadIdx.x];
            X[((rt+WARP-2)/(WARP-1))*COLS+ct+tid-(z-WARP+2)%WARP] = X_bottom[threadIdx.x];
            Y[((rt+WARP-2)/(WARP-1))*COLS+ct+tid-(z-WARP+2)%WARP] = Y_bottom[threadIdx.x];
         }
#endif
      }
   }
   if (tid == (ROWS-2)%(WARP-1)+1) {
      //if (wid==11 && offset[0]==0)printf("output[%d] = log10(%e) - %e = %e\n", wid, result, log10_init, log10(result)-log10_init);
      //printf("output[%d] = LOG10(%1.50E) = %1.50E\n", wid, result, log10(result));
      output[wid] = log10(result) - log10_init; 
   }
   //if (tid == (ROWS-1)%(WARP-1)) output[wid] = result;
}
template<class NUMBER>
int GPUmemAlloc(GPUmem<NUMBER>& gmem) 
{
   cudaDeviceProp deviceProp;
   int devinx;
   cudaGetDevice(&devinx);
   cudaError_t err = cudaGetDeviceProperties(&deviceProp, devinx);
   char *current, *d_current;
   int size = sizeof(NUMBER);
   if (err) printf("Error getting device props\n");
   printf("GPUmemAlloc for device %d\n", devinx);
   gmem.offset = (int2*)malloc((MAX_PROBS+1)*sizeof(int2));
   cudaHostRegister(gmem.offset, (MAX_PROBS+1)*sizeof(int2), cudaHostRegisterMapped);
   cudaMalloc(&gmem.d_offset, (MAX_PROBS+1)*sizeof(int2));
   cudaMalloc(&gmem.d_results, (MAX_PROBS+1)*sizeof(double));
   gmem.totalMem = deviceProp.totalGlobalMem/4;
   //TODO no need to assign d_M, etc.
   cudaMalloc(&gmem.d_amem, gmem.totalMem);
   cudaCheckError(__LINE__, __FILE__);
   gmem.amem = (char*)malloc(gmem.totalMem);
   cudaHostRegister(gmem.amem, gmem.totalMem, cudaHostRegisterMapped);
   cudaCheckError(__LINE__, __FILE__);
   //cudaMallocHost(&gmem.amem, gmem.totalMem); 
   d_current = (char*)gmem.d_amem;
   current = (char*)gmem.amem;

   gmem.d_M = (NUMBER*)d_current;
   gmem.M = (NUMBER*)current;
   current += size*(gmem.totalMem/10/size); d_current += size*(gmem.totalMem/10/size);

   gmem.d_X = (NUMBER*)d_current;
   gmem.X = (NUMBER*)current;
   current += size*(gmem.totalMem/10/size); d_current += size*(gmem.totalMem/10/size);

   gmem.d_Y = (NUMBER*)d_current;
   gmem.Y = (NUMBER*)current;
   current += size*(gmem.totalMem/10/size); d_current += size*(gmem.totalMem/10/size);

   //TODO don't assign these, they get reassigned anyway
   gmem.d_p = (NUMBER*)d_current;
   gmem.p = (NUMBER*)current;
   current += size*(gmem.totalMem/10/size); d_current += size*(gmem.totalMem/10/size);

   gmem.d_q = (NUMBER*)d_current;
   gmem.q = (NUMBER*)current;
   current += size*(gmem.totalMem/10/size); d_current += size*(gmem.totalMem/10/size);

   gmem.d_n = (int*)d_current;
   gmem.n = (int*)current;
   current += 4*(gmem.totalMem/10/4); d_current += 4*(gmem.totalMem/10/4);

   gmem.d_rs = d_current;
   gmem.rs = current;
   current += gmem.totalMem/5; d_current += gmem.totalMem/5;

   gmem.d_hap = d_current;
   gmem.hap = current;
   current += gmem.totalMem/5; d_current += gmem.totalMem/5;

   if( current - (char*)gmem.amem > gmem.totalMem) {
      printf("Error: Requested too much memory on GPU\n");
      return 9000;
   }
   cudaCheckError(__LINE__, __FILE__);
   if (!gmem.amem) {
      printf("CPU mem allocation fail\n");
      return 1;
   }
   if (!gmem.d_amem) {
      printf("GPU mem allocation fail\n");
      printf("gmem.d_amem = %p\n", gmem.d_amem);
      return 1;
   }
   gmem.N_STREAMS=STREAMS;
   gmem.strm = (cudaStream_t*)malloc(sizeof(cudaStream_t)*gmem.N_STREAMS);
   cudaStreamCreate(&gmem.marker_s);
   for (int z=0;z<gmem.N_STREAMS;z++) cudaStreamCreate(&gmem.strm[z]);
   return 0;
}
template<class NUMBER>
int GPUmemFree(GPUmem<NUMBER>& gmem) 
{
   if (NULL==gmem.amem) {
      return 0;
   }
   gmem.index=0;
   gmem.M=0;
   gmem.X=0;
   gmem.Y=0;
   gmem.p=0;
   gmem.q=0;
   gmem.rs=0;
   gmem.hap=0;
   gmem.n=0;
   gmem.d_M=gmem.d_X=gmem.d_Y=gmem.d_p=gmem.d_q=0;
   gmem.d_rs=gmem.d_hap=0;
   gmem.d_n=0;
   cudaFree(gmem.d_amem);
   gmem.d_amem = 0;
   cudaHostUnregister(gmem.amem);
   free(gmem.amem);
   gmem.d_amem = 0;
   gmem.amem = 0;
   free(gmem.offset);gmem.offset=NULL;
   cudaFree(gmem.d_offset); gmem.d_offset = NULL;
   cudaFree(gmem.d_results); gmem.d_results = NULL;
   for (int z=0;z<gmem.N_STREAMS;z++) cudaStreamDestroy(gmem.strm[z]);
   cudaStreamDestroy(gmem.marker_s);
   free(gmem.strm);
   return 0;
}
template int GPUmemAlloc<double>(GPUmem<double>&);
template int GPUmemAlloc<float>(GPUmem<float>&);
template int GPUmemFree<double>(GPUmem<double>&);
template int GPUmemFree<float>(GPUmem<float>&);
template __global__ void debugMark<1>();
template __global__ void debugMark<2>();
template __global__ void debugMark<3>();
template __global__ void debugMark<4>();

#define N_STREAM 4
template <class PRECISION>
void compute_gpu(int2 *offset, char* rs, char* hap, PRECISION* q, 
                 int* n, PRECISION init_const, int n_tc, GPUmem<PRECISION>& gmem, double FAILED_RUN_RESULT)
{
    cudaStream_t strm[N_STREAM];
    printf("%d streams\n", N_STREAM);
    int start[N_STREAM+1];
    int l_n_tc=0;
    for (int z=0;z<N_STREAM;z++,l_n_tc+=n_tc/N_STREAM) {
       cudaStreamCreate(&strm[z]);
       start[z] = l_n_tc;
    }
    start[N_STREAM]=n_tc;
    for (int z=0;z<N_STREAM;z++) {
       compute_gpu_stream(offset+start[z], rs, hap, q, n, init_const, start[z+1]-start[z], 
                          gmem, strm[z], start[z], FAILED_RUN_RESULT);
       //memcpy(&gmem.results[start[z]], gmem.results, sizeof(PRECISION)*(start[z+1]-start[z]));
    }
    for (int z=0;z<N_STREAM;z++) {
       cudaStreamSynchronize(strm[z]);
       cudaStreamDestroy(strm[z]);
    }
}
template <class PRECISION>
void compute_gpu_stream(int2 *offset, char* rs, char* hap, PRECISION* q, 
                           int* n, PRECISION init_const, int n_tc, GPUmem<PRECISION>& gmem, 
                           cudaStream_t strm, int results_inx, double FAILED_RUN_RESULT) 
{
   //GPUmem<PRECISION> gmem;
   //cudaEvent_t start, finish;
   //cudaEventCreate(&start);
   //cudaEventCreate(&finish);

#if 0
   printf("offset = \n");
   for (int z=0;z<n_tc+1;z++) printf("%d:%d,%d\n",z,offset[z].x, offset[z].y);
   printf("rs = ");
   for (int z=0;z<offset[n_tc].x;z++) 
      if (rs[z]=='\0') printf("\n");
      else printf("%c", rs[z]);
   printf("\nhap = ");
   for (int z=0;z<offset[n_tc].y;z++) 
      if (hap[z]=='\0') printf("\n");
      else printf("%c", hap[z]);
   //printf("\nq = \n");
   //for (int z=0;z<offset[n_tc].x;z++) printf("%d:%e\n", z, q[z]);
   printf("init_const = %e\n", init_const);
   printf("n = \n");
   for (int z=0;z<offset[n_tc].x;z++) printf("%d:%d\n", z, n[z]);
#endif
   if (0==gmem.amem) {
      GPUmemAlloc<PRECISION>(gmem);
   }
   cudaMemcpyAsync(gmem.d_offset+results_inx, &offset[0], sizeof(int2)*(n_tc+1), 
                   cudaMemcpyHostToDevice, strm);
   cudaCheckError(__LINE__, __FILE__);
   cudaMemcpyAsync(gmem.d_n+offset[0].x, n+offset[0].x, sizeof(int)*(offset[n_tc].x-offset[0].x), 
                   cudaMemcpyHostToDevice, strm);
   cudaMemcpyAsync(gmem.d_rs+offset[0].x, rs+offset[0].x, sizeof(char)*(offset[n_tc].x-offset[0].x),
                   cudaMemcpyHostToDevice, strm);
   cudaMemcpyAsync(gmem.d_hap+offset[0].y, hap+offset[0].y, sizeof(char)*(offset[n_tc].y-offset[0].y),
                   cudaMemcpyHostToDevice, strm);
   cudaCheckError(__LINE__, __FILE__);
#ifdef __USE_TEX
   printf("Create a texture\n");
   createNewTextureFloat(gmem.n_tex.tex, gmem.n_tex.RD, gmem.n_tex.TD, gmem.d_n, 
                         offset[n_tc].x-offset[0].x);
#endif
   //createNewTextureFloat(gmem.q_tex.tex, gmem.q_tex.RD, gmem.q_tex.TD, gmem.d_q);
   cudaCheckError(__LINE__, __FILE__);
   //One warp handles one matrix
	PRECISION LOG10_INITIAL_CONSTANT = log10(init_const);
   //printf("sizeof(PRECISION) = %d. init_const = %e. log(init_const) = %e\n", 
   //                         sizeof(PRECISION), init_const, LOG10_INITIAL_CONSTANT);
   
#ifdef __WARP_PER_MAT
   printf("One warp/matrix\n");
   //cudaEventRecord(start,strm);
   //pairhmm_kernel<<<(n_tc+3)/4,WARP*4,0,strm>>>( init_const, gmem.d_M, gmem.d_X, 
   pairhmm_kernel_wrap<<<(n_tc+3)/4,WARP*4,0,strm>>>( init_const, gmem.d_M, gmem.d_X, gmem.d_Y, 
                                  gmem.d_rs, gmem.d_hap, gmem.d_q, gmem.d_n, gmem.n_tex.tex,
                                  gmem.d_offset+results_inx, n_tc-1, gmem.d_results+results_inx, 
                                  LOG10_INITIAL_CONSTANT, FAILED_RUN_RESULT); 
#else
   //printf("One thread/matrix\n");
   //cudaEvent_t start, finish;
   //cudaEventCreate(&start);
   //cudaEventCreate(&finish);
   //cudaEventRecord(start,strm);
   //repeatedSine<<<(n_tc+4*WARP-1)/(4*WARP), 4*WARP, 0, strm>>> (gmem.d_q, n_tc-1);

   pairhmm_jacopo<PRECISION,18><<<(n_tc+4*WARP-1)/(4*WARP),4*WARP,0,strm>>>(  gmem.d_rs, gmem.d_hap, 
                                           gmem.d_q, gmem.d_n, gmem.n_tex.tex, 
                                           (int2*)gmem.d_offset+results_inx, n_tc-1, 
                                           gmem.d_results+results_inx, init_const, 
                                           LOG10_INITIAL_CONSTANT, FAILED_RUN_RESULT); 

   //cudaEventRecord(finish,strm);
   //cudaEventSynchronize(finish);
   //float elapsed; 
   //cudaEventElapsedTime(&elapsed, start,finish);
   //printf("Elapsed Time: %f\n", elapsed);
   //cudaEventDestroy(start);
   //cudaEventDestroy(finish);
#endif
   cudaCheckError(__LINE__, __FILE__);
   cudaMemcpyAsync(gmem.results+results_inx, gmem.d_results+results_inx, sizeof(double)*n_tc,
                   cudaMemcpyDeviceToHost,strm);
   //cudaFree(gmem.d_offset);
   //gmem.d_offset=0;
   //GPUmemFree(gmem);
   cudaCheckError(__LINE__, __FILE__);
}
template void compute_gpu<double>(int2*, char*, char*, double*, int*, double, int, GPUmem<double>&, double);
template void compute_gpu<float>(int2*, char*, char*, float*, int*, float, int, GPUmem<float>&, double);
template void compute_gpu_head<float>(int*, char*, char*, int*, float, int, double*, double);
template void compute_gpu_head<double>(int*, char*, char*, int*, double, int, double*, double);
