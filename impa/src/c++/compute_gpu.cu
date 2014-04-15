//#include "input.h"
#include "compute_gpu.h"
#include "stdio.h"

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
                               NUMBER* q, int* offset, int n_mats,
                               NUMBER* output, NUMBER log10_init) {
   NUMBER M_p, M_pp, X_p, X_pp, Y_p, Y_pp, distm, pMM, pGapM, pXX, 
          pMX, pYY, pMY, M_loc, X_loc, Y_loc;
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
   NUMBER M_top, X_top, Y_top;
   NUMBER result=0.0;
   M+=offset[3*wid];
   X+=offset[3*wid];
   Y+=offset[3*wid];
   p+=6*offset[3*wid+1];
   rs+=offset[3*wid+1];
   hap+=offset[3*wid+2];
   q+=offset[3*wid+1];
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
	   _q = q[tid-1+stripe];
      } else {
      _rs = rs[tid+stripe];
	   _q = q[tid+stripe];
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
			   	distm = double(1.0) - _q;
            else distm = _q;
            }
            else distm = _q;
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
            if (r==ROWS-1) { 
               result += M_loc + X_loc;
            }
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
         
      }
   }
   if (tid == (ROWS-1)%(WARP-1)) {
      output[wid] = log10(result) - log10_init; 
   }
   //if (tid == (ROWS-1)%(WARP-1)) output[wid] = result;
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
   gmem.q = (NUMBER*)malloc(totalMem/66);
   cudaMalloc(&gmem.d_rs, 700000);
   gmem.rs = (char*)malloc(700000);
   cudaMalloc(&gmem.d_hap, 700000);
   gmem.hap = (char*)malloc(700000);
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
       !gmem.hap) {
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
   return 0;
}
template int GPUmemAlloc<double>(GPUmem<double>&);
template int GPUmemAlloc<float>(GPUmem<float>&);
template int GPUmemFree<double>(GPUmem<double>&);
template int GPUmemFree<float>(GPUmem<float>&);

template <class PRECISION>
void compute_gpu(int offset[][3], PRECISION* p, char* rs, char* hap, PRECISION* q, 
                           PRECISION Yr0, int n_tc, PRECISION* h_out, GPUmem<PRECISION>& gmem) 
{
   //GPUmem<PRECISION> gmem;
   PRECISION *d_out;
   cudaError_t cuerr;

   if (0==gmem.M) {
      GPUmemAlloc<PRECISION>(gmem);
   }
   cudaMalloc(&d_out,  sizeof(PRECISION)*n_tc);
   cudaMalloc(&gmem.d_offset, sizeof(int)*3*(n_tc+1));
   cudaMemcpy(gmem.d_offset, &offset[0][0], sizeof(int)*3*(n_tc+1), cudaMemcpyHostToDevice);
   cudaMemcpy(gmem.d_p, p, sizeof(PRECISION)*offset[n_tc][1]*6, cudaMemcpyHostToDevice);
   cudaMemcpy(gmem.d_rs, rs, sizeof(char)*offset[n_tc][1], cudaMemcpyHostToDevice);
   cudaMemcpy(gmem.d_hap, hap, sizeof(char)*offset[n_tc][2], cudaMemcpyHostToDevice);
   cudaMemcpy(gmem.d_q, q, sizeof(PRECISION)*offset[n_tc][1], cudaMemcpyHostToDevice);
   fflush(0);
   cuerr= cudaGetLastError();
   if (cuerr) printf("Error in memcpy. %d : %s\n", cuerr, cudaGetErrorString(cuerr));
   //One warp handles one matrix
	PRECISION INITIAL_CONSTANT = ldexp(1.0, 1020.0);
	PRECISION LOG10_INITIAL_CONSTANT = log10(INITIAL_CONSTANT);
   
   pairhmm_kernel<<<(n_tc+3)/4,WARP*4>>>( Yr0, gmem.d_M, gmem.d_X, 
                                  gmem.d_Y, gmem.d_p, 
                                  gmem.d_rs, gmem.d_hap, gmem.d_q,
                                  gmem.d_offset, n_tc-1, d_out, LOG10_INITIAL_CONSTANT); 
   cuerr = cudaGetLastError();
   if (cuerr) {
      printf ("Cuda error %d : %s\n", cuerr, cudaGetErrorString(cuerr));
   }
#if 0
   PRECISION *M2;
   PRECISION *X2;
   //TODO Fix this!
   cudaMemcpy(gmem.M, gmem.d_M,
                  sizeof(PRECISION)*gmem.offset[n_tc][0]+1338,
                  cudaMemcpyDeviceToHost);
   cudaMemcpy(gmem.X, gmem.d_X,
                  sizeof(PRECISION)*gmem.offset[n_tc][0]+1338,
                  cudaMemcpyDeviceToHost);
   cuerr = cudaGetLastError();
   if (cuerr) printf ("Memcpy(D2H) error %d : %s\n", cuerr, cudaGetErrorString(cuerr));
   for (int z=0;z<n_tc;z++)
   {
      int ROWS = gmem.offset[z+1][1]-gmem.offset[z][1];
      int COLS = gmem.offset[z+1][2]-gmem.offset[z][2];
      M2 = gmem.M+gmem.offset[z][0]+COLS*(((ROWS-1)+WARP-2)/(WARP-1));
      X2 = gmem.X+gmem.offset[z][0]+COLS*(((ROWS-1)+WARP-2)/(WARP-1));
	   PRECISION result = 0.0;
	   for (int c = 0; c < COLS; c++)
	   	result += M2[c] + X2[c];

   	if (before_last_log != NULL)
   		*before_last_log = result;	

      probs[z] = log10(out[z]) - LOG10_INITIAL_CONSTANT;
   }
#else
   cudaMemcpy(h_out, d_out, sizeof(PRECISION)*n_tc, cudaMemcpyDeviceToHost);
#endif
   cudaFree(d_out);
   //GPUmemFree(gmem);
}
template void compute_gpu<double>(int [][3], double*, char*, char*, double*, double, int, double*, GPUmem<double>&);
template void compute_gpu<float>(int [][3], float*, char*, char*, float*, float, int, float*, GPUmem<float>&);
