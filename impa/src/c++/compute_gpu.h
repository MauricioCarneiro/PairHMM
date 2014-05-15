#ifndef GPU_COMPUTE_H
#define GPU_COMPUTE_H
#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5
#define II 0
#define CC 1
#define DD 2
#define MAX(a,b) ((a)>(b)?(a):(b))

#define WARP 32 
struct cudaTextureData {
   cudaResourceDesc RD;
   cudaTextureDesc TD;
   cudaTextureObject_t tex;

};
template<class PRECISION>
struct GPUmem {
   int offset[200001][3];
   int* d_offset;
   int index;
   int N_STREAMS;
   cudaStream_t *strm;
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
   PRECISION* results;
   PRECISION* d_results;

   /*  scratch  */
   PRECISION* M;
   PRECISION* X;
   PRECISION* Y;
   PRECISION* d_M;
   PRECISION* d_X;
   PRECISION* d_Y;

   GPUmem() {M=0;};
}; 
template<class PRECISION>
int GPUmemAlloc(GPUmem<PRECISION>& gmem); 
template<class PRECISION>
int GPUmemFree(GPUmem<PRECISION>& gmem);
template <class PRECISION>
void compute_gpu(int offset[][3], char *rs, char* hap, PRECISION* q, int* n,
                  PRECISION init_const, int n_tc, GPUmem<PRECISION>&);
template <class PRECISION>
void compute_gpu_stream(int offset[][3], char *rs, char* hap, PRECISION* q,
                  int* n, PRECISION init_const, int n_tc, GPUmem<PRECISION>&, cudaStream_t strm, int start);
void cudaCheckError(int line, const char* file);
void createNewTextureFloat(cudaTextureObject_t& tex, cudaResourceDesc& resDesc, cudaTextureDesc& texDesc, void* devPtr);
__global__ void CPU_start();
__global__ void CPU_end();
#endif
