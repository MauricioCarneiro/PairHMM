//#include "input.h"
#include "compute_gpu.h"
#include "stdio.h"

void cudaCheckError(int line, const char* file) {
   cudaError_t err = cudaGetLastError();
   if (err) {
       fprintf(stderr, "error %d on line %d : %s\n", err, line, cudaGetErrorString(err));
       fflush(0);
   }
}

void createNewTextureFloat(cudaTextureObject_t& tex, cudaResourceDesc& resDesc, cudaTextureDesc& texDesc, void* devPtr) {

   tex=0;
   memset(&resDesc, 0, sizeof(cudaResourceDesc));
   resDesc.res.linear.devPtr = devPtr;
   resDesc.resType = cudaResourceTypeLinear;
   resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
   resDesc.res.linear.desc.x = 32;
   resDesc.res.linear.sizeInBytes = 10000*sizeof(float);
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

#if 0
//TODO Alloc more memory in M,X,Y due to smaller "warps"
//TODO Xc0 and Yr0 are constants
//TODO compute and save result
template <class NUMBER, unsigned int BATCH>
__global__ void
__lanuch_bounds__(128,1)
pairhmm_kernel_onethread( NUMBER init_const, NUMBER* M, NUMBER *X, NUMBER *Y
                          NUMBER* p, char* rs, char* hap, 
                          NUMBER* q, int* offset, int n_mats,
                          NUMBER* output, NUMBER log10_init) {
   NUMBER M_p[BATCH+1];
   NUMBER M_pp[BATCH+1];
   NUMBER X_p[BATCH+1];
   NUMBER X_pp[BATCH+1];
   NUMBER Y_p[BATCH+1];
   NUMBER Y_pp[BATCH+1];
   char _rs[BATCH+1];
   int _q[BATCH+1];
   NUMBER distm;
   NUMBER pMM[BATCH+1];
   NUMBER pGapM[BATCH+1];
   NUMBER pXX[BATCH+1];
   NUMBER pMX[BATCH+1];
   NUMBER pYY[BATCH+1];
   NUMBER pMY[BATCH+1];
   NUMBER M_loc[BATCH+1], X_loc[BATCH+1], Y_loc[BATCH+1];
   
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   M+=offset[3*tid]; 
   X+=offset[3*tid]; 
   Y+=offset[3*tid]; 
   p+=offset[3*tid+1]; 
   rs+=offset[3*tid+1]; 
   hap+=offset[3*tid+2]; 
   q+=offset[3*tid+1];
   for (int stripe = 0; stripe < ROWS; stripe+=BATCH-1) 
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
   for (int bid = 1; bid < BATCH+1; bid++)
   {
	   _rs[bid] = tc->rs[bid-1+stripe];
	   _q[bid] = tc->q[bid-1+stripe];
      pMM[bid] = p[6*(bid+stripe)+MM];
      pGapM[bid] = p[6*(bid+stripe)+GapM];
      pXX[bid] = p[6*(bid+stripe)+XX];
      pMX[bid] = p[6*(bid+stripe)+MX];
      pYY[bid] = p[6*(bid+stripe)+YY];
      pMY[bid] = p[6*(bid+stripe)+MY];
   }
   for (int z = 1; z < COLS+BATCH+1; z++) 
   {
      for (int bid = 1; bid < BATCH+1; bid++)
      {
         if (bid <= z && bid+stripe < ROWS && z-bid+1 < COLS) 
         {  
            r = bid+stripe;
            c = z-bid+1;
			   char _hap = tc->hap[c-1];
			   if (_rs[bid] == _hap || _rs[bid] == 'N' || _hap == 'N')
			   	distm = double(1.0) - _q[bid];
            else distm = _q[bid];
            if (bid == 1) {
               X_p[0] = X[stripe/(BATCH-1)*COLS+z+bid][c];
               Y_p[0] = Y[stripe][c];
               M_p[0] = M[stripe][c];
            }
            M_loc[bid] = distm * (M_pp[bid] * pMM[bid] + X_pp[bid] * pGapM[bid] + Y_pp[bid] * pGapM[bid]);
			   Y_loc[bid] = M_p[bid] * pMY[bid] + Y_p[bid] * pYY[bid];
			   X_loc[bid] = M_p[bid-1] * pMX[bid] + X_p[bid-1] * pXX[bid];
         }
      }
      for (int bid = BATCH; bid > 0; bid--)
      {
         if (bid <= z && bid+stripe < ROWS && z-bid+1 < COLS) 
         {  
            M_pp[bid]=M_p[bid-1];
            X_pp[bid]=X_p[bid-1];
            Y_pp[bid]=Y_p[bid-1];
            X_p[bid] = X_loc[bid];
            M_p[bid] = M_loc[bid];
            Y_p[bid] = Y_loc[bid];
            r = bid+stripe;
            c = z-bid+1;
            if (abs((M_loc[bid]-M[r][c])/M_loc[bid]) > 0.00000000001) printf("M[%d][%d] : %1.16e != %1.16e\n", r,c,M_loc[bid], M[r][c]);
            //TODO only write at edges
            M[r][c] = M_loc[bid];
            X[r][c] = X_loc[bid];
            Y[r][c] = Y_loc[bid];
         } else if (bid==z+1 && bid+stripe < ROWS) {
            M_pp[bid]=M_p[bid-1];
            X_pp[bid]=X_p[bid-1];
            Y_pp[bid]=Y_p[bid-1];
            //Read from col 0
            //TODO make this a row, not a column
            X_p[bid] = X[bid+stripe][0];
            M_p[bid] = M[bid+stripe][0];
            Y_p[bid] = Y[bid+stripe][0];
         }
          
      }
   }
   }
}
#endif
__global__ void CPU_start() {}
__global__ void CPU_end() {}

template<class NUMBER>
__global__ void 
__launch_bounds__(128,9)
pairhmm_kernel( NUMBER init_const, NUMBER* M, NUMBER *X, NUMBER *Y, 
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
   int n_off=3*offset[3*wid+1];
   n+=3*offset[3*wid+1];
   for (int stripe = 0; stripe < ROWS-1; stripe+=WARP-1) {
      if ( stripe==0 && tid < 2) {
         M_pp=0.0;
         X_pp = 0.0;//X_pp=Xc0[0];
         Y_pp=init_const/(COLS-1);
         M_p=0.0;
         if (tid==1) X_p = 0.0;//X_p=Xc0[1];
         else X_p=0.0;
         if (tid==0) Y_p=init_const/(COLS-1);
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
      //TODO pad?
      if (tid>0) {
      _rs = rs[tid-1+stripe];
	   _q = q[tid-1+stripe];
      } else {
      _rs = rs[tid+stripe];
	   _q = q[tid+stripe];
      }
      //TODO transpose n for coalesced reads?
      //TODO pass n as a single integer since only the first 7 bits of each matter
#ifdef DEBUG
      if (wid==1 && 
              tex1Dfetch<int>(t_n, 3*(tid+stripe)+DD+n_off) != n[3*(tid+stripe)+DD]) {
                   printf("mismatch : stripe = %d, tid = %d: D = %d, I = %d, C = %d\nFrom (int*)n: D = %d, I = %d, C = %d\n", stripe, tid,
                                tex1Dfetch<int>(t_n, 3*(tid+stripe)+DD+n_off),
                                tex1Dfetch<int>(t_n, 3*(tid+stripe)+II+n_off),
                                tex1Dfetch<int>(t_n, 3*(tid+stripe)+CC+n_off),
                                n[3*(tid+stripe)+DD],
                                n[3*(tid+stripe)+II],
                                n[3*(tid+stripe)+CC]);
      } else if(wid==1) printf("Match : stripe = %d, tid = %d\n", stripe, tid);
#endif
#ifdef __USE_TEX
      pMM = 1.0 - ph2pr<NUMBER>((tex1Dfetch<int>(t_n,3*(tid+stripe)+II+n_off)+tex1Dfetch<int>(t_n,3*(tid+stripe)+DD+n_off)) & 127);
      pXX = ph2pr<NUMBER>(tex1Dfetch<int>(t_n,3*(tid+stripe)+CC+n_off));
      pMX = ph2pr<NUMBER>(tex1Dfetch<int>(t_n,3*(tid+stripe)+II+n_off));
      pMY = tid + stripe == ROWS-1 ? 1.0 : ph2pr<NUMBER>(tex1Dfetch<int>(t_n,3*(tid+stripe)+DD+n_off));
#else
      pMM = 1.0 - ph2pr<NUMBER>(n[3*(tid+stripe)+II]+n[3*(tid+stripe)+DD] & 127);
      pXX = ph2pr<NUMBER>(n[3*(tid+stripe)+CC]);
      pMX = ph2pr<NUMBER>(n[3*(tid+stripe)+II]);
      pMY = tid + stripe == ROWS-1 ? 1.0 : ph2pr<NUMBER>(n[3*(tid+stripe)+DD]);
#endif

      pGapM = 1.0 - pXX;
      pYY = tid + stripe == ROWS-1 ? 1.0 : pXX;
      //TODO get rid of this?
      if (tid+stripe==0) {
         pMM = pXX = pGapM = pMX = pYY = pMY = 0.0;
      }
      for (int z = 1; z < COLS+WARP+1;z++) 
      {
         int r = tid+stripe;
         int c = z-tid+1;
         //TODO align at word boundaries
#ifdef __COALESCED_READ
         //TODO write/read to/from shared memory instead?
         //TODO don't read if stripe==0
         if (1==z%WARP) {
            M_top = M[stripe/(WARP-1)*COLS+z+tid];
            X_top = X[stripe/(WARP-1)*COLS+z+tid];
            Y_top = Y[stripe/(WARP-1)*COLS+z+tid];
            //if (wid==0) printf("%d. tid: %d. M_top = %e\n", z, tid, M_top);
         } else {
            M_top = __shfl_down(M_top,1);
            X_top = __shfl_down(X_top,1);
            Y_top = __shfl_down(Y_top,1);
         }
#endif
#ifdef __SHARED_READ
         if (1==z%WARP) {
            M_top[threadIdx.x] = M[stripe/(WARP-1)*COLS+z+tid];
            X_top[threadIdx.x] = X[stripe/(WARP-1)*COLS+z+tid];
            Y_top[threadIdx.x] = Y[stripe/(WARP-1)*COLS+z+tid];
            //__syncthreads(); No need for sync. Warps maintain synchronization
         }
#endif
      

         if (tid<= z+1 && tid+stripe < ROWS && z-tid < COLS)
         {
            //TODO pad instead
            if (c>0) {
            char _hap = hap[c-1];
			   if (_rs == _hap || _rs == 'N' || _hap == 'N')
			   	distm = double(1.0) - _q;
            else distm = _q;
            }
            else distm = _q;
            if (tid == 0 && stripe == 0) {
               X_p = 0.0; 
               Y_p = init_const/(COLS-1);
               M_p = 0.0;
            } else if (tid == 0 && z > 1) {
#ifdef __COALESCED_READ
               M_p = M_top;
               X_p = X_top;
               Y_p = Y_top;
#else
#ifdef __SHARED_READ
               M_p = M_top[threadIdx.x + (z-1)%WARP];
               X_p = X_top[threadIdx.x + (z-1)%WARP];
               Y_p = Y_top[threadIdx.x + (z-1)%WARP];
#else
               M_p = M[stripe/(WARP-1)*COLS+z];
               X_p = X[stripe/(WARP-1)*COLS+z];
               Y_p = Y[stripe/(WARP-1)*COLS+z];
#endif
#endif
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
            if (c < COLS && tid==WARP-1) {
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
      //printf("output[%d] = LOG10(%1.50E) = %1.50E\n", wid, result, log10(result));
      output[wid] = log10(result) - log10_init; 
   }
   //if (tid == (ROWS-1)%(WARP-1)) output[wid] = result;
}
template<class NUMBER>
int GPUmemAlloc(GPUmem<NUMBER>& gmem) 
{
   cudaDeviceProp deviceProp;
   cudaError_t err = cudaGetDeviceProperties(&deviceProp, 0);
   char *current, *d_current;
   int size = sizeof(NUMBER);
   gmem.totalMem = 3*deviceProp.totalGlobalMem/4;
   //TODO no need to assign d_M, etc.
   cudaMalloc(&gmem.d_amem, gmem.totalMem);
   //gmem.amem = (char*)malloc(gmem.totalMem);
   cudaMallocHost(&gmem.amem, gmem.totalMem); 
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
   err = cudaGetLastError();
   if (err) printf("cudaMalloc error %d: %s\n", err, cudaGetErrorString(err));
   if (err) return 9000+err;
   if (!gmem.M ||
       !gmem.X ||
       !gmem.Y ||
       !gmem.p ||
       !gmem.q ||
       !gmem.rs ||
       !gmem.hap) {
      printf("CPU mem allocation fail\n");
      return 1;
   }
   gmem.N_STREAMS=4;
   gmem.strm = (cudaStream_t*)malloc(sizeof(cudaStream_t)*gmem.N_STREAMS);
   for (int z=0;z<gmem.N_STREAMS;z++) cudaStreamCreate(&gmem.strm[z]);
   return 0;
}
template<class NUMBER>
int GPUmemFree(GPUmem<NUMBER>& gmem) 
{
   if (0==gmem.amem) {
      return 0;
   }
   cudaFree(gmem.d_amem);
   cudaFree(gmem.amem);
   gmem.d_amem = 0;
   gmem.amem = 0;
   gmem.index=0;
   gmem.M=0;
   gmem.X=0;
   gmem.Y=0;
   gmem.p=0;
   gmem.q=0;
   gmem.rs=0;
   gmem.hap=0;
   for (int z=0;z<gmem.N_STREAMS;z++) cudaStreamDestroy(gmem.strm[z]);
   free(gmem.strm);
   return 0;
}
template int GPUmemAlloc<double>(GPUmem<double>&);
template int GPUmemAlloc<float>(GPUmem<float>&);
template int GPUmemFree<double>(GPUmem<double>&);
template int GPUmemFree<float>(GPUmem<float>&);

#define N_STREAM 4
template <class PRECISION>
void compute_gpu(int offset[][3], char* rs, char* hap, PRECISION* q, 
                 int* n, PRECISION init_const, int n_tc, GPUmem<PRECISION>& gmem)
{
    cudaStream_t strm[N_STREAM];
    int start[N_STREAM+1];
    int l_n_tc=0;
    for (int z=0;z<N_STREAM;z++,l_n_tc+=n_tc/N_STREAM) {
       cudaStreamCreate(&strm[z]);
       start[z] = l_n_tc;
    }
    start[N_STREAM]=n_tc;
    for (int z=0;z<N_STREAM;z++) {
       compute_gpu_stream(&offset[start[z]], rs, hap, q, n, init_const, start[z+1]-start[z], gmem, strm[z], start[z]);
       //memcpy(&gmem.results[start[z]], gmem.results, sizeof(PRECISION)*(start[z+1]-start[z]));
    }
    for (int z=0;z<N_STREAM;z++) {
       cudaStreamSynchronize(strm[z]);
       cudaStreamDestroy(strm[z]);
    }
}
template <class PRECISION>
void compute_gpu_stream(int offset[][3], char* rs, char* hap, PRECISION* q, 
                           int* n, PRECISION init_const, int n_tc, GPUmem<PRECISION>& gmem, 
                           cudaStream_t strm, int results_inx) 
{
   //GPUmem<PRECISION> gmem;
   cudaError_t cuerr;

#if 0
   printf("offset = \n");
   for (int z=0;z<n_tc+1;z++) printf("%d:%d,%d,%d\n",z,offset[z][0], offset[z][1],offset[z][2]); 
   printf("p = \n");
   for (int z=0;z<offset[n_tc][1]*6;z++) printf("%d:%1.13f\n", z, p[z]);
   printf("rs = ");
   for (int z=0;z<offset[n_tc][1];z++) printf("%c", rs[z]);
   printf("\nhap = ");
   for (int z=0;z<offset[n_tc][2];z++) printf("%c", hap[z]);
   printf("\nq = \n");
   for (int z=0;z<offset[n_tc][1];z++) printf("%d:%e\n", z, q[z]);
   printf("init_const = %e\n", init_const);
   printf("n = \n");
   for (int z=0;z<offset[n_tc][1]*3;z++) printf("%d:%d\n", z/3, n[z]);
#endif
   if (0==gmem.M) {
      GPUmemAlloc<PRECISION>(gmem);
   }
   cudaMalloc(&gmem.d_offset, sizeof(int)*3*(n_tc+1));
   cudaMemcpyAsync(gmem.d_offset, &offset[0][0], sizeof(int)*3*(n_tc+1), cudaMemcpyHostToDevice, strm);
#ifdef __CONDENSE_MEM
   cudaMemcpyAsync(gmem.d_p, p, sizeof(PRECISION)*offset[n_tc][1]*6 +
                           sizeof(PRECISION)*offset[n_tc][1] +
                           sizeof(int)*offset[n_tc][1]*3 + 
                           sizeof(char)*offset[n_tc][1] +
                           sizeof(char)*offset[n_tc][2], cudaMemcpyHostToDevice, strm);
#else
   cudaMemcpyAsync(gmem.d_q+offset[0][1], q+offset[0][1], sizeof(PRECISION)*(offset[n_tc][1]-offset[0][1]), cudaMemcpyHostToDevice, strm);
   cudaMemcpyAsync(gmem.d_n+offset[0][1]*3, n+offset[0][1]*3, sizeof(int)*(offset[n_tc][1]-offset[0][1])*3, cudaMemcpyHostToDevice, strm);
   cudaMemcpyAsync(gmem.d_rs+offset[0][1], rs+offset[0][1], sizeof(char)*(offset[n_tc][1]-offset[0][1]), cudaMemcpyHostToDevice, strm);
   cudaMemcpyAsync(gmem.d_hap+offset[0][2], hap+offset[0][2], sizeof(char)*(offset[n_tc][2]-offset[0][2]), cudaMemcpyHostToDevice, strm);
#endif
   createNewTextureFloat(gmem.n_tex.tex, gmem.n_tex.RD, gmem.n_tex.TD, gmem.d_n);
   createNewTextureFloat(gmem.q_tex.tex, gmem.q_tex.RD, gmem.q_tex.TD, gmem.d_q);
   cuerr= cudaGetLastError();
   if (cuerr) printf("Error in memcpy. %d : %s\n", cuerr, cudaGetErrorString(cuerr));
   //One warp handles one matrix
	PRECISION INITIAL_CONSTANT = ldexp(1.0, 1020.0);
	PRECISION LOG10_INITIAL_CONSTANT = log10(INITIAL_CONSTANT);
   
   pairhmm_kernel<<<(n_tc+3)/4,WARP*4,0,strm>>>( init_const, gmem.d_M, gmem.d_X, 
                                  gmem.d_Y, gmem.d_rs, gmem.d_hap, gmem.d_q, gmem.d_n, gmem.n_tex.tex,
                                  gmem.d_offset, n_tc-1, gmem.d_results+results_inx, LOG10_INITIAL_CONSTANT); 
   cuerr = cudaGetLastError();
   if (cuerr) {
      printf ("Cuda error %d : %s\n", cuerr, cudaGetErrorString(cuerr));
   }
   cudaMemcpyAsync(gmem.results+results_inx, gmem.d_results+results_inx, sizeof(PRECISION)*n_tc,
                   cudaMemcpyDeviceToHost,strm);
   //GPUmemFree(gmem);
}
template void compute_gpu<double>(int [][3], char*, char*, double*, int*, double, int, GPUmem<double>&);
template void compute_gpu<float>(int [][3], char*, char*, float*, int*, float, int, GPUmem<float>&);
