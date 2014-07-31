#define MAX_PROBS 512*1024
template <class PRECISION>
void compute_gpu_head(int* offset, char* rs, char* hap, int* n,
                   PRECISION init_const, int n_tc, double* results, double FAILED_RUN_RESULT);

template <class T>
void GPUAlloc();
template <class T>
void GPUFree();
bool isGPUAllocated();

void cudaCheckError(int line, const char* file);
