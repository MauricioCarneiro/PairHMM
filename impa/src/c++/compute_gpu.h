#ifndef GPU_COMPUTE_H
#define GPU_COMPUTE_H
#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5
#define MAX(a,b) ((a)>(b)?(a):(b))

#define WARP 32 
template<class PRECISION>
struct GPUmem {
   int offset[2000][3];
   int* d_offset;
   int index;
   PRECISION* d_M;
   PRECISION* d_X;
   PRECISION* d_Y;
   PRECISION* d_p;
   PRECISION* d_Yr0;
   PRECISION* d_Xc0;
   PRECISION* M;
   PRECISION* X;
   PRECISION* Y;
   PRECISION* p;
   PRECISION* Xc0;
   PRECISION* Yr0;
   char* rs;
   char* hap;
   PRECISION* q;
   char* d_rs;
   char* d_hap;
   PRECISION* d_q; 
   GPUmem() {M=0;};
}; 
template<class PRECISION>
int GPUmemAlloc(GPUmem<PRECISION>& gmem); 
template<class PRECISION>
int GPUmemFree(GPUmem<PRECISION>& gmem);
template <class PRECISION>
void compute_gpu(int offset[][3], PRECISION *p, char *rs, char* hap, PRECISION* q,
                  PRECISION Yr0, int n_tc, PRECISION* h_out, GPUmem<PRECISION>&);
#endif
