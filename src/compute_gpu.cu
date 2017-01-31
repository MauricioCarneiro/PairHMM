//#include "input.h"
#include "compute_gpu.h"
#include "compute_gpu.cuh"
#include "stdio.h"
#include <cfloat>
#define MAT_PER_WARP 8
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
template <class T>
void gpu_stream_sync(int s);

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
int get_nstreams() {
   if (gf_gmem.amem != 0) return gf_gmem.N_STREAMS;
   else if (gd_gmem.amem != 0) return gd_gmem.N_STREAMS;
   else return 0;
}
template <class T> 
void compute_gpu_setup(int* offset_in, int n_tc, double* results);
cudaStream_t dbg_strm;
template <> 
void compute_gpu_setup<float>(int* offset_in, int n_tc, double* results) {
   int2* offset= (int2*) offset_in;
   if (!isGPUAllocated()) {
      GPUAlloc<float>();
   }
   cudaCheckError(__LINE__, __FILE__);
   gf_gmem.results = results;
   char* current = (char*)gf_gmem.d_amem;

   gf_gmem.d_results = (double*) current;
   current += sizeof(double) * n_tc;

   gf_gmem.d_offset = (int2*) current;
   current += sizeof(int2) * (n_tc+1);

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
   //if (isGPUAllocated()) GPUFree<float>();
   cudaCheckError(__LINE__, __FILE__);
   cudaStreamCreate(&dbg_strm);
}
template <> 
void compute_gpu_head<float>(int* offset_in, char* rs, char* hap, int* n,
                      float init_const, int n_tc, double* results, double FAILED_RUN_RESULT) {
   compute_gpu_setup<float>(offset_in, n_tc, results);
   compute_gpu_stream(offset_in, rs, hap, n, init_const, n_tc, gf_gmem, 0, 0, FAILED_RUN_RESULT);
   gf_gmem.results= NULL;
   cudaCheckError(__LINE__, __FILE__);
}
//To enable this to be called from a file not compiled with nvcc
template <>
void compute_gpu_stream<float>(int* offset_in, char* rs, char* hap, int* n, float init_const, int n_tc, 
                        int strm, int results_inx, double FAILED_RUN_RESULT) {
   compute_gpu_stream(offset_in, rs, hap, n, init_const, n_tc, gf_gmem, gf_gmem.strm[strm], 
                      results_inx, FAILED_RUN_RESULT);
}
template <>
void compute_gpu_stream<double>(int* offset_in, char* rs, char* hap, int* n, double init_const, int n_tc, 
                        int strm, int results_inx, double FAILED_RUN_RESULT) {
   compute_gpu_stream(offset_in, rs, hap, n, init_const, n_tc, gd_gmem, gd_gmem.strm[strm], 
                      results_inx, FAILED_RUN_RESULT);
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
void gpu_sync() {
   fprintf(stderr,"Sync\n");
   for(int z=0;z<gf_gmem.N_STREAMS;z++)
      cudaStreamSynchronize(gf_gmem.strm[z]);
   cudaMemcpy(gf_gmem.results, gf_gmem.d_results, sizeof(double)*4, cudaMemcpyDeviceToHost);
   //cudaDeviceSynchronize();
}
template <unsigned int I>
void dbgMark() {
    debugMark<I><<<1,1,0,dbg_strm>>>();
}
template <>
void gpu_stream_sync<float>(int s) {
   cudaStreamSynchronize(gf_gmem.strm[s]);
   cudaCheckError(__LINE__, __FILE__);
}
template <>
void gpu_stream_sync<double>(int s) {
   cudaStreamSynchronize(gd_gmem.strm[s]);
   cudaCheckError(__LINE__, __FILE__);
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
   compute_gpu_stream(offset_in, rs, hap, n, init_const, n_tc, gd_gmem, 0, 0, FAILED_RUN_RESULT );
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
    if (pid >= n_mats) {
       return;
    }
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
__global__ void CPU_start() {}
__global__ void CPU_end() {}
template <unsigned int u>
__global__ void debugMark() {}

template<class NUMBER>
int GPUmemAlloc(GPUmem<NUMBER>& gmem) 
{
   cudaDeviceProp deviceProp;
   int devinx;
   cudaGetDevice(&devinx);
   cudaError_t err = cudaGetDeviceProperties(&deviceProp, devinx);
   if (err) printf("Error getting device props\n");
   printf("GPUmemAlloc for device %d\n", devinx);
   gmem.offset = (int2*)malloc((MAX_PROBS+1)*sizeof(int2));
   cudaHostRegister(gmem.offset, (MAX_PROBS+1)*sizeof(int2), cudaHostRegisterMapped);
   gmem.totalMem = deviceProp.totalGlobalMem/4;
   //TODO no need to assign d_M, etc.
   cudaMalloc(&gmem.d_amem, gmem.totalMem);
   cudaCheckError(__LINE__, __FILE__);
   gmem.amem = (char*)malloc(gmem.totalMem);
   cudaHostRegister(gmem.amem, gmem.totalMem, cudaHostRegisterMapped);
   cudaCheckError(__LINE__, __FILE__);
   //cudaMallocHost(&gmem.amem, gmem.totalMem); 
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
   printf("Creating %d streams\n", STREAMS);
   for (int z=0;z<gmem.N_STREAMS;z++) {
      cudaStreamCreate(&gmem.strm[z]);
   }
   cudaCheckError(__LINE__, __FILE__);
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
   gmem.d_results=0;
   gmem.d_offset=0;
   cudaFree(gmem.d_amem);
   gmem.d_amem = 0;
   cudaHostUnregister(gmem.amem);
   free(gmem.amem);
   gmem.d_amem = 0;
   gmem.amem = 0;
   free(gmem.offset);gmem.offset=NULL;
   //cudaFree(gmem.d_offset); gmem.d_offset = NULL;
   //cudaFree(gmem.d_results); gmem.d_results = NULL;
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
       compute_gpu_stream((int*)(offset+start[z]), rs, hap, n, init_const, start[z+1]-start[z], 
                          gmem, strm[z], start[z], FAILED_RUN_RESULT);
    }
    for (int z=0;z<N_STREAM;z++) {
       cudaStreamSynchronize(strm[z]);
       cudaStreamDestroy(strm[z]);
    }
}
void test_copy() {
   cudaCheckError(__LINE__, __FILE__);
   double* junk  = (double*)malloc(sizeof(double)*200000);
   double* d_junk;
   cudaMalloc(&d_junk, sizeof(double)*200000);
   cudaCheckError(__LINE__, __FILE__);
   cudaMemcpy(junk, d_junk, sizeof(double)*1024,
                   cudaMemcpyDeviceToHost);
   cudaCheckError(__LINE__, __FILE__);
}
template <>
void gpu_copy_results<float>(int strm, int start, int finish) {
   cudaCheckError(__LINE__, __FILE__);
   cudaMemcpyAsync(gf_gmem.results+start, gf_gmem.d_results+start, sizeof(double)*(finish-start),
                   cudaMemcpyDeviceToHost, gf_gmem.strm[strm]);
   cudaCheckError(__LINE__, __FILE__);
}
template <>
void gpu_copy_results<double>(int strm, int start, int finish) {
   cudaCheckError(__LINE__, __FILE__);
   cudaMemcpyAsync(gd_gmem.results+start, gd_gmem.d_results+start, sizeof(double)*(finish-start),
                   cudaMemcpyDeviceToHost, gd_gmem.strm[strm]);
   cudaCheckError(__LINE__, __FILE__);
}
template <class PRECISION>
void compute_gpu_stream(int *offset_in, char* rs, char* hap, 
                           int* n, PRECISION init_const, int n_tc, GPUmem<PRECISION>& gmem, 
                           cudaStream_t strm, int results_inx, double FAILED_RUN_RESULT) 
{
   int2* offset = (int2*) offset_in;

   if (0==gmem.amem) {
      GPUmemAlloc<PRECISION>(gmem);
   }
   cudaMemcpyAsync(gmem.d_offset+results_inx, offset, sizeof(int2)*(n_tc+1), 
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
   cudaCheckError(__LINE__, __FILE__);
   //One warp handles one matrix
	PRECISION LOG10_INITIAL_CONSTANT = log10(init_const);
   
   pairhmm_jacopo<PRECISION,18><<<(n_tc+4*WARP-1)/(4*WARP),4*WARP,0,strm>>>(  gmem.d_rs, gmem.d_hap, 
                                           gmem.d_q, gmem.d_n, gmem.n_tex.tex, 
                                           (int2*)gmem.d_offset+results_inx, n_tc-1, 
                                           gmem.d_results+results_inx, init_const, 
                                           LOG10_INITIAL_CONSTANT, FAILED_RUN_RESULT); 
   cudaCheckError(__LINE__, __FILE__);
}
template void compute_gpu<double>(int2*, char*, char*, double*, int*, double, int, GPUmem<double>&, double);
template void compute_gpu<float>(int2*, char*, char*, float*, int*, float, int, GPUmem<float>&, double);
template void compute_gpu_head<float>(int*, char*, char*, int*, float, int, double*, double);
template void compute_gpu_head<double>(int*, char*, char*, int*, double, int, double*, double);
template void compute_gpu_stream<float>(int*, char*, char*, int*, float, int, int, int, double);
template void compute_gpu_stream<double>(int*, char*, char*, int*, double, int, int, int, double);
template void gpu_stream_sync<float>(int);
template void gpu_stream_sync<double>(int);
template void gpu_copy_results<float>(int, int, int);
template void gpu_copy_results<double>(int, int, int);
template void dbgMark<1>();
template void dbgMark<2>();
template void dbgMark<3>();
template void dbgMark<4>();
template void dbgMark<5>();

