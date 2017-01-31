#define MAX_PROBS 512*1024
template <class PRECISION>
void compute_gpu_head(int* offset, char* rs, char* hap, int* n,
                   PRECISION init_const, int n_tc, double* results, double FAILED_RUN_RESULT);

template <class PRECISION>
void compute_gpu_setup(int *offset, int n_tc, double* results); 

template <class PRECISION>
void compute_gpu_stream(int *offset, char *rs, char* hap, 
                  int* n, PRECISION init_const, int n_tc, int strm, int start,
                  double FAILED_RUN_RESULT);
template <class PRECISION>
void gpu_copy_results(int strm, int start, int finish);
template <class PRECISION>
void gpu_stream_sync(int s);
template <class T>
void GPUAlloc();
template <class T>
void GPUFree();
bool isGPUAllocated();
int get_nstreams();
void test_copy();
template <unsigned int I>
void dbgMark();

void cudaCheckError(int line, const char* file);
void gpu_sync();
