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
  std::vector<PRECISION> calculate (const std::vector<Read<PRECISION,PRECISION>>& padded_reads,
                                    const std::vector<Haplotype<PRECISION>>& padded_haplotypes) {
     int r=0, c=0, rc=0, n_mat=0;
     std::vector<int>offsets;
     auto results = std::vector<PRECISION>{};
     results.reserve(padded_reads.size() * padded_haplotypes.size());
     GPUmem<PRECISION> gmem;
     auto hap_reverse = std::vector<char>{};
     auto rs_copy = std::vector<char>{};

     //first count the rows, cols and r*c/WARP_SIZE
     //offsets is are total number of rows*cols/WARP, rows and cols
     //in all matrices up to this one
     offsets.push_back(0); offsets.push_back(0); offsets.push_back(0);
     for (const auto& read : padded_reads) {
       Base::m_constants.update(read);
       for (const auto& hap : padded_haplotypes) {
         int this_r = read.original_length+1;
         int this_c = hap.original_length+1;
         int this_rc;
         r += this_r;
         c += this_c;
         this_rc = ((this_r + this->GPU_WARP_SIZE - 2)/(this->GPU_WARP_SIZE - 1)) * this_c;
         rc += this_rc;
         offsets.push_back(rc);
         offsets.push_back(r);
         offsets.push_back(c);
         //printf("rs (size %d)= ", (int)read.bases.size());
         for (int q=read.left_padding; q<read.bases.size();q++) {
         //   printf("%c", (char)read.bases[q]);
            rs_copy.push_back(read.bases[q]);
         }
         rs_copy.push_back('\0');
         //printf("\n");
         //printf("hap (size %d) = ", (int)hap.original_length);
         for (int q=hap.original_length+hap.original_length; q>hap.original_length;q--) {
         //   printf("%c", (char)hap.bases[q]);
            hap_reverse.push_back(hap.bases[q]);
         }
         hap_reverse.push_back('\0');
         //printf("\n");
         ++n_mat;
       }
     }

     printf("r=%d. c=%d, rc=%d\n", r,c,rc);
     //GPUAlloc(r,c,rc);
     
     // single allocation for these for an easier copy to the GPU
     //      PRECISION p[6*offsets[3*n_mat+1]];
     //      PRECISION q[offsets[3*n_mat+1]];
     //      char rs[offsets[3*n_mat+1]];
     //      char haps[offsets[3*n_mat+2]];
 
     /*** Change this order at your peril  ***/
     PRECISION* p = (PRECISION*) malloc(sizeof(PRECISION)*7*r+ sizeof(PRECISION)*(r+c));
     PRECISION *q = p+6*r;
     char *rs = (char*)(q+r);
     char *haps = (char*)(rs+r);

     int inx=0;
     //auto &diags = this->m_diagonals;
     auto &consts = this->m_constants;
     for (const auto& read : padded_reads) {
       Base::m_constants.update(read);
       for (const auto& hap : padded_haplotypes) {
          //Base::m_diagonals.update(Base::INITIAL_CONSTANT / hap.original_length);
          size_t r,c;
          r = offsets[3*(inx+1)+1]-offsets[3*inx+1];
          c = offsets[3*(inx+1)+2]-offsets[3*inx+2];
          //printf("r= %d\n",(int)r);
          //printf("c= %d\n",(int)c);
          //printf("init = %e\n", Base::INITIAL_CONSTANT);
        //memcpy(p+6*offsets[3*inx+1]+MM, &consts.mm[0], sizeof(PRECISION)*r);
        //memcpy(p+6*offsets[3*inx+1]+GapM, &consts.gm[0], sizeof(PRECISION)*r);
        //memcpy(p+6*offsets[3*inx+1]+MX, &consts.mx[0], sizeof(PRECISION)*r);
        //memcpy(p+6*offsets[3*inx+1]+XX, &consts.xx[0], sizeof(PRECISION)*r);
        //memcpy(p+6*offsets[3*inx+1]+YY, &consts.yy[0], sizeof(PRECISION)*r);
        //memcpy(p+6*offsets[3*inx+1]+MY, &consts.my[0], sizeof(PRECISION)*r);
        int cnt=0;
        for (int z=offsets[3*inx+1];z<offsets[3*inx+4];z++,cnt++) {
           p[z*6+MM] = consts.mm[cnt];
           p[z*6+GapM] = consts.gm[cnt];
           p[z*6+MX] = consts.mx[cnt];
           p[z*6+XX] = consts.xx[cnt];
           p[z*6+MY] = consts.my[cnt];
           p[z*6+YY] = consts.yy[cnt];
        }
          //printf("rs_copy = ");
          //for (int z=offsets[3*inx+1];z<offsets[3*inx+4];z++) printf("%c", rs_copy[z]);
          //printf("\n offset = %d\n", offsets[3*inx+1]);
          memcpy(rs+offsets[3*inx+1], &rs_copy[offsets[3*inx+1]], sizeof(char)*r);
          //printf("rs[%d] = ", offsets[3*inx+1]);
          //for (int z=offsets[3*inx+1];z<offsets[3*inx+4];z++) printf("%c", rs[z]);
          //printf("\n");
          
          memcpy(haps+offsets[3*inx+2], &hap_reverse[offsets[3*inx+2]], sizeof(char)*c);
          memcpy(q+offsets[3*inx+1], &read.base_quals[1], sizeof(PRECISION)*r);
          q[offsets[3*inx+1]+r-1] = 1.00000; // possibly unnecessary
          results.push_back(0.0);
          //`for (int q=offsets[3*inx+1]; q<offsets[3*inx+4];q++) printf("%c", rs_copy[q]);
          //printf(" vs. \n");
          //for (int q=offsets[3*inx+2]; q<offsets[3*inx+5];q++) printf("%c", hap_reverse[q]);
          //printf("\n");
          ++inx;
       }
     }
     
     compute_gpu<PRECISION>((int (*)[3])&offsets[0], p, rs, haps, q,
                            Base::INITIAL_CONSTANT, n_mat, &results[0], gmem);
     if (gmem.M != 0) GPUmemFree(gmem);
     return results;
  }
  void calculate_failed (const Testcase& testcase, std::vector<double>& results) {
  }
  double do_compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) override {
     return 0.0;
  }
};





#endif
