#ifndef GPU_COMPUTE_H
#define GPU_COMPUTE_H
#define MM 0
#define MX 1
#define XX 2
#define MY 3
#define YY 4
#define GapM 5
#define II 0
#define CC 1
#define DD 2
#define MAX(a,b) ((a)>(b)?(a):(b))

#define WARP 32 
#ifndef STREAMS
#define STREAMS 4
#endif

struct cudaTextureData {
   cudaResourceDesc RD;
   cudaTextureDesc TD;
   cudaTextureObject_t tex;

};
template<class PRECISION>
struct GPUmem {
   int2* offset;
   int2* d_offset;
   int index;
   int N_STREAMS;
   cudaStream_t *strm;
   cudaStream_t marker_s;
   cudaTextureData n_tex;
   cudaTextureData q_tex;
   char* d_amem; //these two get alloced, all else
   char* amem;   //  are pointers within this allocation
   unsigned long long totalMem; //allocated GPU memory (in sizeofchar))
   /* input */
   PRECISION* p;
   char* rs;
   char* hap;
   int* n;
   PRECISION* q;
   PRECISION* d_p;
   char* d_rs;
   char* d_hap;
   int * d_n;
   PRECISION* d_q; 

   /*  output */
   double* results;
   double* d_results;

   /*  scratch  */
   PRECISION* M;
   PRECISION* X;
   PRECISION* Y;
   PRECISION* d_M;
   PRECISION* d_X;
   PRECISION* d_Y;

   GPUmem() {amem=0;};
}; 
template<class PRECISION>
int GPUmemAlloc(GPUmem<PRECISION>& gmem); 
template<class PRECISION>
int GPUmemFree(GPUmem<PRECISION>& gmem);
void createNewTextureFloat(cudaTextureObject_t& tex, cudaResourceDesc& resDesc, cudaTextureDesc& texDesc, void* devPtr);
__global__ void CPU_start();
__global__ void CPU_end();
template <unsigned int u>
__global__ void debugMark();
template <class PRECISION>
void compute_gpu(int2 *offset, char *rs, char* hap, PRECISION* q, int* n,
                  PRECISION init_const, int n_tc, GPUmem<PRECISION>&, double FAILED_RUN_RESULT);
template <class PRECISION>
void compute_gpu_stream(int *offset, char *rs, char* hap, 
                  int* n, PRECISION init_const, int n_tc, GPUmem<PRECISION>&, cudaStream_t strm, int start,
                  double FAILED_RUN_RESULT);
#endif
