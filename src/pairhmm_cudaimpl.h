#ifndef PAIRHMM_CUDAIMPL_H
#define PAIRHMM_CUDAIMPL_H

#include "pairhmm_impl.h" 
#include "stdio.h"
#include "compute_gpu.h"
#include <string.h>
#include "chronos.h"

struct int2 {
   int x; int y;
};

struct pair {
   Read<uint8_t,uint8_t> read;
   Haplotype<uint8_t> hap;
   int inx;
   double result;

   public:
      pair(Read<uint8_t,uint8_t> r, Haplotype<uint8_t> h, int i) {
         this->read=r;this->hap=h;this->inx=i; result = 0.0;
      }
};

bool pair_comp (pair* pA, pair* pB) {
     int rA = pA->read.bases.size() * pA->hap.bases.size();
     int rB = pB->read.bases.size() * pB->hap.bases.size();
     return rA<rB;
}

template <class PRECISION>
class PairhmmCudaImpl: public PairhmmImpl<PRECISION, Diagonals3<PRECISION>, Constants<PRECISION>, 1>  {
  using Base =  PairhmmImpl<PRECISION,Diagonals3<PRECISION>, Constants<PRECISION>, 1>;
  std::vector<int2> offsets;
  std::vector<char> hap_reverse;
  std::vector<char> rs_copy;
  std::vector<int> n;
  size_t next = 0;
  int2 this_offset;
  std::vector<pair*> pair_list;
  double* output;
  double sow_time = 0.f, setup_time = 0.f, compute_time = 0.f, sort_time = 0.f, extract_time = 0.f;
  void clear() {
     offsets.clear();
     hap_reverse.clear();
     rs_copy.clear();
     n.clear();
  }
  void one_pair_extract(pair p) {
     int this_r = p.read.bases.size()+1;
     int this_c = p.hap.bases.size()+1;
     this_offset.x += this_r;
     this_offset.y += this_c;
     offsets.push_back(this_offset);
     for (int z=0; z<p.hap.bases.size(); z++) {
        //printf("%c", (char)hap.bases[q]);
        hap_reverse.push_back((char)p.hap.bases[z]);
     }
     hap_reverse.push_back('\0');
     //printf("\n");
     for (int z=0;z<p.read.bases.size();z++) {
        //TODO try this at 256 (w/ unsigned)
        int i = (int)p.read.ins_quals[z] & 127;
        int d = (int)p.read.del_quals[z] & 127;
        int c = (int)p.read.gcp_quals[z] & 127;
        int q = (int)p.read.base_quals[z] & 127;
        n.push_back(i + d * 128 + c * 128 * 128 + q * 128 * 128 * 128);
        //printf("push_back %d,%d,%d,%d => %d\n", i, d, c, q, (((( q * 128 ) + c) * 128 + d) * 128) + i);
        rs_copy.push_back((char)p.read.bases[z]);
     }
     rs_copy.push_back('\0');
     n.push_back(0);
     results.push_back(0.0); //this will get replaced with the actual value
  }
  void one_pair_extract_parallel(pair p, int2 this_off) {
  //   this version assumes rs_copy and hap_reverse are "reserved" and 
  //   that the offset list is already built
     for (int z=0; z<p.hap.bases.size(); z++) {
        //printf("%c", (char)hap.bases[q]);
        hap_reverse[z+this_off.y]=(char)p.hap.bases[z];
     }
     hap_reverse[p.hap.bases.size()+this_off.y]='\0';
     //printf("\n");
     for (int z=0;z<p.read.bases.size();z++) {
        //TODO try this at 256 (w/ unsigned)
        int i = (int)p.read.ins_quals[z] & 127;
        int d = (int)p.read.del_quals[z] & 127;
        int c = (int)p.read.gcp_quals[z] & 127;
        int q = (int)p.read.base_quals[z] & 127;
        n[z+this_off.x]= i + d * 128 + c * 128 * 128 + q * 128 * 128 * 128;
        rs_copy[z+this_off.x]=(char)p.read.bases[z];
     }
     rs_copy[p.read.bases.size()+this_off.x]='\0';
     n[p.read.bases.size()+this_off.x]=0;
  }
  void one_testcase_extract( const std::vector<Read<uint8_t,uint8_t>>& padded_reads,
                                    const std::vector<Haplotype<uint8_t>>& padded_haplotypes) {
     for (const auto& read : padded_reads) {
       //Base::m_constants.update(read);
       for (const auto& hap : padded_haplotypes) {
          pair p;
          p.read = read;
          p.hap = hap;
          one_pair_extract(p);
       }
     }
  }
#ifdef __TIMING_DETAIL
#define CHECKPT(ABC) ABC+=time.elapsed();\
                         fprintf(stderr,\
                         "%s %f ms\n",\
                         #ABC,time.elapsed());\
                         time.reset(); 
#else 
#define CHECKPT(ABC) ;;
#endif
public:
  std::vector<double> results;
  PairhmmCudaImpl(const size_t initial_size = Base::INITIAL_SIZE): Base {initial_size} { }
  virtual ~PairhmmCudaImpl() { }
  void sow(const Testcase& testcase) {
    Chronos time;
    time.reset();
    sow(testcase.reads, testcase.haplotypes);
    sow_time += time.elapsed();
  }
  std::vector<double> calculate (const Testcase& testcase) {
    //const auto padded_haplotypes = this->pad_haplotypes<PRECISION>(testcase); // updates m_max_original_read_length (important!)
    //const auto padded_reads = this->pad_reads_noconvert(testcase.reads);
    sow(testcase.reads, testcase.haplotypes);
    return reap();
  }
  void gpu_compute( ) {
#ifdef __TIMING_DETAIL
     fprintf(stderr, "Total sow time: %f ms\n\n", sow_time);
#endif
     fprintf(stderr,"Solving %lu matrices\n", pair_list.size());
     Chronos time; time.reset();
     output = (double*)realloc(output, sizeof(double) * pair_list.size());
     std::sort(pair_list.begin(), pair_list.end(), pair_comp);
     //std::sort(pair_list.begin(), pair_list.end());

     CHECKPT(sort_time);

     int total_rows = 0, total_cols = 0, total_cells = 0;
     offsets.reserve(pair_list.size()+1);
     for (int z = 0; z < pair_list.size(); z++) {
        offsets[z].x = total_rows;
        offsets[z].y = total_cols;
        total_rows += pair_list[z]->read.bases.size()+1;
        total_cols += pair_list[z]->hap.bases.size()+1;
        total_cells += pair_list[z]->read.bases.size()*pair_list[z]->hap.bases.size();
        results.push_back(0.0);
     }
     offsets[pair_list.size()].x = total_rows;
     offsets[pair_list.size()].y = total_cols;
     n.reserve(total_rows);
     rs_copy.reserve(total_rows);
     hap_reverse.reserve(total_cols);
     #pragma omp parallel for 
     for (int z = 0; z < pair_list.size(); z++) {
        one_pair_extract_parallel(*pair_list[z], offsets[z]);
     }
     CHECKPT(setup_time);

     compute_gpu_head<PRECISION>((int*)&offsets[0], (char*)&rs_copy[0], (char*)&hap_reverse[0], (int*)&n[0],
                         (PRECISION)Base::INITIAL_CONSTANT, results.size()-next, output, 
                         (double)this->FAILED_RUN_RESULT );
     CHECKPT(compute_time);

     for (int z=0; z < results.size()-next; z++) {
        results[pair_list[z]->inx] = output[z];
     }

     clear();
     next=results.size();
     for (auto const& pt : pair_list) delete(pt);
     pair_list.clear();
     CHECKPT(extract_time);
  }
  void sow (const std::vector<Read<uint8_t,uint8_t>>& padded_reads,
                                    const std::vector<Haplotype<uint8_t>>& padded_haplotypes) {
     if (pair_list.size() + padded_reads.size() * padded_haplotypes.size() > MAX_PROBS) {
        Chronos time; time.reset();
        gpu_compute();
        fprintf(stderr, " Total gpu_compute time: %f ms\n\n", time.elapsed());
     } 
     //If results have just been harvested, clear results
     if (0==next && 0 == rs_copy.size()) results.clear();
     if (0==offsets.size()) {
        this_offset.x=0;this_offset.y=0;offsets.push_back(this_offset);
     }
     for (const auto& read : padded_reads) {
        for (const auto& hap : padded_haplotypes) {
        pair_list.push_back(new pair(read, hap, pair_list.size()+next));
        }
     }
  }
  std::vector<double> reap() {
     printf("Before reap. Solving %lu matrices\n", pair_list.size());
     if (0 != pair_list.size() ) { 
        gpu_compute();
     }
     next=0;
#ifdef __TIMING_DETAIL
     fprintf(stderr, "Totals:\n");
     fprintf(stderr, "    sow_time : %f\n", sow_time);
     fprintf(stderr, "    sort_time : %f ( %f \%)\n", sort_time, 100.f*sort_time/sow_time);
     fprintf(stderr, "    setup_time : %f ( %f \%)\n", setup_time, 100.f*setup_time/sow_time);
     fprintf(stderr, "    compute_time : %f ( %f \%)\n", compute_time, 100.f*compute_time/sow_time);
     fprintf(stderr, "    extract_time : %f ( %f \%)\n", extract_time, 100.f*extract_time/sow_time);
     sow_time = sort_time = setup_time = compute_time = extract_time = 0.f;
#endif
     return results;
  }
  void calculate_failed (const Testcase& testcase, std::vector<double>& results) {
  }
  double do_compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) override {
     return 0.0;
  }
};





#endif
