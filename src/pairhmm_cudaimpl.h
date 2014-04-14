#ifndef PAIRHMM_CUDAIMPL_H
#define PAIRHMM_CUDAIMPL_H

#include "pairhmm_impl.h" 
#include "stdio.h"
#include "compute_gpu.h"
#include <string.h>

template <class PRECISION>
class PairhmmCudaImpl: public PairhmmImpl<PRECISION, Diagonals3<PRECISION>, Constants<PRECISION>, 1>  {
  using Base =  PairhmmImpl<PRECISION,Diagonals3<PRECISION>, Constants<PRECISION>, 1>;
public:
  PairhmmCudaImpl(const size_t initial_size = Base::INITIAL_SIZE): Base {initial_size} { }
  virtual ~PairhmmCudaImpl() { }
  std::vector<double> calculate (const Testcase& testcase) {
    return Base::calculate(testcase);
  }
  std::vector<PRECISION> calculate (const std::vector<Read<PRECISION,char>>& padded_reads,
                                    const std::vector<Haplotype<PRECISION>>& padded_haplotypes) {
     int r=0, c=0, rc=0, n_mat=0;
     std::vector<int>offsets;
     auto results = std::vector<PRECISION>{};
     results.reserve(padded_reads.size() * padded_haplotypes.size());
     GPUmem<PRECISION> gmem;

     //first count the rows, cols and r*c/WARP_SIZE
     //offsets is are total number of rows*cols/WARP, rows and cols
     //in all matrices up to this one
     offsets.push_back(0); offsets.push_back(0); offsets.push_back(0);
     for (const auto& read : padded_reads) {
       Base::m_constants.update(read);
       for (const auto& hap : Base::m_padded_haplotypes) {
         int this_r = read.original_length;
         int this_c = hap.original_length;
         int this_rc;
         r += this_r;
         c += this_c;
         this_rc = ((this_r + this->GPU_WARP_SIZE - 2)/(this->GPU_WARP_SIZE - 1)) * this_c;
         rc += this_rc;
         offsets.push_back(rc);
         offsets.push_back(r);
         offsets.push_back(c);
         printf("rs = %s\n", read.bases);
         printf("hap = %s\n", hap.bases);
         ++n_mat;
       }
     }

     printf("r=%d. c=%d, rc=%d", r,c,rc);
     //GPUAlloc(r,c,rc);
     
     // single allocation for these for an easier copy to the GPU
     //      PRECISION p[6*offsets[3*n_mat+1]];
     //      char rs[offsets[3*n_mat+1]];
     //      char haps[offsets[3*n_mat+2]];
     //      PRECISION q[offsets[3*n_mat+1]];
     PRECISION* p = (PRECISION*) malloc(sizeof(PRECISION)*7*r+ sizeof(char)*(r+c));
     char *rs = (char*)(p+6*r);
     char *haps = (char*)(rs+r);
     PRECISION *q = (PRECISION*)(haps+c);

     int inx=0;
     auto &diags = this->m_diagonals;
     auto &consts = this->m_constants;
     for (const auto& read : padded_reads) {
       Base::m_constants.update(read);
       for (const auto& hap : Base::m_padded_haplotypes) {
          Base::m_diagonals.update(Base::INITIAL_CONSTANT / hap.original_length);
          size_t r,c;
          r = offsets[3*(inx+1)+1]-offsets[3*inx+1];
          c = offsets[3*(inx+1)+2]-offsets[3*inx+2];
          memcpy(p+offsets[3*inx+1]*6+MM, &consts.mm[0], sizeof(PRECISION)*r);
          memcpy(p+offsets[3*inx+1]*6+GapM, &consts.gm[0], sizeof(PRECISION)*r);
          memcpy(p+offsets[3*inx+1]*6+MX, &consts.mx[0], sizeof(PRECISION)*r);
          memcpy(p+offsets[3*inx+1]*6+XX, &consts.xx[0], sizeof(PRECISION)*r);
          memcpy(p+offsets[3*inx+1]*6+YY, &consts.yy[0], sizeof(PRECISION)*r);
          memcpy(p+offsets[3*inx+1]*6+MY, &consts.my[0], sizeof(PRECISION)*r);
          memcpy(rs+offsets[3*inx+1], &read.bases[0], sizeof(char)*r);
          memcpy(haps+offsets[3*inx+2], &hap.bases[0], sizeof(char)*c);
          memcpy(q+offsets[3*inx+1], &read.base_quals[0], sizeof(PRECISION)*r);
          results.push_back(0.0);
          printf("base_qual.size = %lu\n", read.base_quals.size());
         printf("\n");
     
       }
     }
     //int* foo[3] = (int (*)[3])&offsets[0];

     
     compute_gpu<PRECISION>((int (*)[3])&offsets[0], p, rs, haps, q,
                            diags.y[0], n_mat, &results[0], gmem);
     return results;
  }
  double do_compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) override {
     return 0.0;
  }
};





#endif
