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
__launch_bounds__(128,6)
pairhmm_kernel( NUMBER init_const, NUMBER* M, NUMBER *X, NUMBER *Y, 
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
               Y_p = init_const/(COLS-1);
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
            if (tid>0 && c < COLS && (r==stripe+WARP-1 /*|| r==ROWS-1*/)) {
               M[((r+WARP-2)/(WARP-1))*COLS+c] = M_loc;
               X[((r+WARP-2)/(WARP-1))*COLS+c] = X_loc;
               Y[((r+WARP-2)/(WARP-1))*COLS+c] = Y_loc;
            }
            if (r==ROWS-1) { 
               result += M_loc + X_loc;
               //printf("stripe = %d t:%d b:%d result = %1.50E\n", stripe, threadIdx.x, blockIdx.x, result);
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
   if (tid == (ROWS-2)%(WARP-1)+1) {
      //printf("t : %d  LOG10(%1.50E) = %1.50E\n", tid, result, log10(result));
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
   gmem.totalMem = 3*deviceProp.totalGlobalMem/4;
   //TODO no need to assign d_M, etc.
   //TODO remove Xc0
   cudaMalloc(&gmem.d_amem, gmem.totalMem);
   //gmem.amem = (char*)malloc(gmem.totalMem);
   cudaMallocHost(&gmem.amem, gmem.totalMem); 
   d_current = (char*)gmem.d_amem;
   current = (char*)gmem.amem;

   gmem.d_M = (NUMBER*)d_current;
   gmem.M = (NUMBER*)current;
   current += 8*(gmem.totalMem/10/8); d_current += 8*(gmem.totalMem/10/8);

   gmem.d_X = (NUMBER*)d_current;
   gmem.X = (NUMBER*)current;
   current += 8*(gmem.totalMem/10/8); d_current += 8*(gmem.totalMem/10/8);

   gmem.d_Y = (NUMBER*)d_current;
   gmem.Y = (NUMBER*)current;
   current += 8*(gmem.totalMem/10/8); d_current += 8*(gmem.totalMem/10/8);

   //TODO don't assign these, they get reassigned anyway
   gmem.d_p = (NUMBER*)d_current;
   gmem.p = (NUMBER*)current;
   current += 8*(gmem.totalMem/10/8); d_current += 8*(gmem.totalMem/10/8);

   gmem.d_q = (NUMBER*)d_current;
   gmem.q = (NUMBER*)current;
   current += 8*(gmem.totalMem/10/8); d_current += 8*(gmem.totalMem/10/8);

   gmem.d_rs = d_current;
   gmem.rs = current;
   current += 700000; d_current += 700000;

   gmem.d_hap = d_current;
   gmem.hap = current;
   current += 700000; d_current += 700000;

   gmem.d_Yr0 = (NUMBER*)d_current;
   gmem.Yr0 = (NUMBER*)current;
   current += gmem.totalMem/10; d_current += gmem.totalMem/10;

   gmem.d_Xc0 = (NUMBER*)d_current;
   gmem.Xc0 = (NUMBER*)current;
   current += gmem.totalMem/10; d_current += gmem.totalMem/10;

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
       !gmem.Yr0 ||
       !gmem.Xc0 ||
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
   gmem.Yr0=0;
   gmem.Xc0=0;
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
void compute_gpu(int offset[][3], PRECISION* p, char* rs, char* hap, PRECISION* q, 
                           PRECISION init_const, int n_tc, GPUmem<PRECISION>& gmem)
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
       compute_gpu_stream(&offset[start[z]], p, rs, hap, q, init_const, start[z+1]-start[z], gmem, strm[z], start[z]);
       //memcpy(&gmem.results[start[z]], gmem.results, sizeof(PRECISION)*(start[z+1]-start[z]));
    }
    for (int z=0;z<N_STREAM;z++) {
       cudaStreamSynchronize(strm[z]);
       cudaStreamDestroy(strm[z]);
    }
}
template <class PRECISION>
void compute_gpu_stream(int offset[][3], PRECISION* p, char* rs, char* hap, PRECISION* q, 
                           PRECISION init_const, int n_tc, GPUmem<PRECISION>& gmem, cudaStream_t strm,
                           int results_inx) 
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
#endif
   if (0==gmem.M) {
      GPUmemAlloc<PRECISION>(gmem);
   }
   cudaMalloc(&gmem.d_offset, sizeof(int)*3*(n_tc+1));
   cudaMemcpyAsync(gmem.d_offset, &offset[0][0], sizeof(int)*3*(n_tc+1), cudaMemcpyHostToDevice, strm);
#ifdef __CONDENSE_MEM
   cudaMemcpyAsync(gmem.d_p, p, sizeof(PRECISION)*offset[n_tc][1]*6 +
                           sizeof(PRECISION)*offset[n_tc][1] +
                           sizeof(char)*offset[n_tc][1] +
                           sizeof(char)*offset[n_tc][2], cudaMemcpyHostToDevice, strm);
#else
   cudaMemcpyAsync(gmem.d_p+offset[0][1]*6, p+offset[0][1]*6, sizeof(PRECISION)*(offset[n_tc][1]-offset[0][1])*6, 
                          cudaMemcpyHostToDevice, strm);
   cudaMemcpyAsync(gmem.d_q+offset[0][1], q+offset[0][1], sizeof(PRECISION)*(offset[n_tc][1]-offset[0][1]), cudaMemcpyHostToDevice, strm);
   cudaMemcpyAsync(gmem.d_rs+offset[0][1], rs+offset[0][1], sizeof(char)*(offset[n_tc][1]-offset[0][1]), cudaMemcpyHostToDevice, strm);
   cudaMemcpyAsync(gmem.d_hap+offset[0][2], hap+offset[0][2], sizeof(char)*(offset[n_tc][2]-offset[0][2]), cudaMemcpyHostToDevice, strm);
#endif
   cuerr= cudaGetLastError();
   if (cuerr) printf("Error in memcpy. %d : %s\n", cuerr, cudaGetErrorString(cuerr));
   //One warp handles one matrix
	PRECISION INITIAL_CONSTANT = ldexp(1.0, 1020.0);
	PRECISION LOG10_INITIAL_CONSTANT = log10(INITIAL_CONSTANT);
   
   pairhmm_kernel<<<(n_tc+3)/4,WARP*4,0,strm>>>( init_const, gmem.d_M, gmem.d_X, 
                                  gmem.d_Y, gmem.d_p, 
                                  gmem.d_rs, gmem.d_hap, gmem.d_q,
                                  gmem.d_offset, n_tc-1, gmem.d_results+results_inx, LOG10_INITIAL_CONSTANT); 
   cuerr = cudaGetLastError();
   if (cuerr) {
      printf ("Cuda error %d : %s\n", cuerr, cudaGetErrorString(cuerr));
   }
   cudaMemcpyAsync(gmem.results+results_inx, gmem.d_results+results_inx, sizeof(PRECISION)*n_tc, cudaMemcpyDeviceToHost,strm);
   //GPUmemFree(gmem);
}
template void compute_gpu<double>(int [][3], double*, char*, char*, double*, double, int, GPUmem<double>&);
template void compute_gpu<float>(int [][3], float*, char*, char*, float*, float, int, GPUmem<float>&);
