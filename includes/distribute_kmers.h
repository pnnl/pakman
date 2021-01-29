// ***********************************************************************
//
// PaKman: Algorithm for generating genomic contigs on distributed-memory machines
// 
// Priyanka Ghosh (Pacific Northwest National Laboratory)
// Sriram Krishnamoorthy (Pacific Northwest National Laboratory)
// Ananth Kalyanaraman (Washington State University)
//               
//
// ***********************************************************************
//
//       Copyright (2020) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************

#ifndef DISTRIBUTE_KMERS_H
#define DISTRIBUTE_KMERS_H

#include <stdio.h>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <queue>          // std::priority_queue
#include <string>
#include <assert.h>
#include <cstring>
#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <atomic>
#include <string.h>
#include <algorithm>
#include <numeric>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include<cstdint>
//#include "serial.h"

#ifdef USE_CEREAL
#define CEREAL_THREAD_SAFE 1
#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/archives/portable_binary.hpp>
#endif


#define KMER_LENGTH                    (WINDW_SIZE+1)
#define LMER_SIZE                      (pow(4, LMER_LENGTH))
#define MN_LENGTH                      (KMER_LENGTH-1)

#define MOD 2147483647
#define HL 31


#ifdef EXTEND_KMER
//static_assert(std::is_same_v<__uint128_t, unsigned __int128>);
//typedef unsigned __int128 uint128_t;
//typedef uint128_t kmer_t;
typedef __uint128_t kmer_t;
#else
typedef uint64_t kmer_t;
#endif

using BasePair = uint8_t;
typedef BasePair ElType;
typedef uint64_t lmer_t;

#ifdef EXTEND_KMER
#define KMER_MASK            ((__uint128_t)(~0L)) >> ((sizeof(kmer_t)*8) - (2*KMER_LENGTH))
#define BASE_MASK            ((__uint128_t)(~0L)) >> ((sizeof(kmer_t)*8) - (2*64))
#define LMER_MASK            ((__uint128_t)(~0L)) >> ((sizeof(lmer_t)*8) - (2*LMER_LENGTH))
#define MN_MASK              ((__uint128_t)(~0L)) >> ((sizeof(kmer_t)*8) - (2*MN_LENGTH))
#define SUCC_MASK            ((__uint128_t)(~0L)) >> ((sizeof(kmer_t)*8) - (2*(KMER_LENGTH-1)))
#else
#define KMER_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*KMER_LENGTH))
#define BASE_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*32))
#define LMER_MASK           ((~0UL)) >> ((sizeof(lmer_t)*8) - (2*LMER_LENGTH))
#define MN_MASK             ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*MN_LENGTH))
#define SUCC_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*(KMER_LENGTH-1)))
#endif

//enum ElType { A, C, T, G };
enum KTerminal {N, Y};
enum KdirType {P, S};
extern int rank, size;

class BasePairVector {
 public:
 #ifdef EXTEND_KMER
  using T = __uint128_t;
#else
  using T = uint64_t;
#endif
  static const size_t bsize = 2;
  //static const T bsize_mask = (T{1}<<bsize) - 1;
  static const size_t nels_per_value = (sizeof(T)*8)/bsize;

  //assumes sizeof(T) % bsize == 0
  
  //BasePairVector() = default;
  BasePairVector() { size_ = 0 ;}
  ~BasePairVector() = default;

  BasePairVector(const BasePairVector& bv)
      : vec_{bv.vec_},
        size_ {bv.size_} {}

  BasePair operator [] (size_t pos) const {
    return get(pos);
  };

  BasePair get(size_t pos) const {
    assert(pos < size_);
    size_t idx = pos / nels_per_value;
    size_t idx_p = size_ / nels_per_value;
    size_t off_size = idx < idx_p ? nels_per_value : size_ % nels_per_value;
    size_t off_pos = pos % nels_per_value;

    return ((vec_[idx]) >> ((off_size-(off_pos)-1)*bsize)) & (0x3);
    
  };

  void set(BasePair val, size_t pos) {
    assert(pos < size_);
    assert(val < (1<<bsize));
    size_t idx = pos / nels_per_value;

    //vec_[idx] = ((vec_[idx]<<bsize) | val) & static_cast<T>(KMER_MASK);
    vec_[idx] = ((vec_[idx]<<bsize) | val) & static_cast<T>(BASE_MASK);
  };


#ifdef USE_CEREAL
  template<class Archive>
  void serialize(Archive & archive)
  {
       archive(vec_, size_); // serialize things by passing them to the archive
  }
#endif


  void push_back(BasePair val) {
    assert(val < (1<<bsize));
    size_ += 1;
    if((size_%nels_per_value) == 1) {
      vec_.push_back(T{0});
    }
    set(val, size_-1);
  };

  size_t size() const {
    return size_;
  }
  
  T extract_pred(T pos, size_t dist, size_t length)
  {
    assert (dist < nels_per_value);

    size_t rem = length - dist;
    //if (!rem)
    //    return pos;

    kmer_t pmask = mask(dist) << (rem*bsize);

    return ((pmask & pos) >> (rem*bsize));
  }

  void resize(size_t nsize) {
    if (size_) {
       size_t new_size = nsize/nels_per_value + (nsize%nels_per_value ? 1 : 0);
       size_t idx_vsize = nsize/nels_per_value;
       size_t idx_size = size_/nels_per_value;
       size_t size_length = idx_vsize < idx_size ? nels_per_value : size_ % nels_per_value;
       size_t new_vsize_off = nsize%nels_per_value;
   
       if(idx_vsize >= vec_.size())
          printf("In rank: %d, size_: %lu, vec_ size: %lu, idx_vsize: %lu, nsize: %lu\n", 
                 rank, size_, vec_.size(), idx_vsize, nsize);
       assert(idx_vsize < vec_.size());   
       vec_[idx_vsize] = extract_pred(vec_[idx_vsize], new_vsize_off, size_length);
       vec_.resize(new_size);
       size_ = nsize;
    }
    else {
      vec_.resize(nsize/nels_per_value + (nsize%nels_per_value ? 1 : 0));
      size_=nsize;
    }
  }

  void mnode_extract_pred(const BasePairVector& bpv, size_t len)
  {
      //assert(bpv.size() >= nels_per_value);
      assert(bpv.size() > MN_LENGTH);
      assert (size_ == 0);
      size_t new_size = len/nels_per_value + (len%nels_per_value ? 1 : 0);
      vec_.resize(new_size);
      size_ = len;
      size_t ori_size = bpv.retrieve_size();
      //size_t ori_off = bpv.size()%nels_per_value;
      //size_t num_itr = ori_off > 0 ? ori_size-2 : ori_size-1;
      size_t idx=0;
     
      size_t move_len = bpv.size() >= nels_per_value ? (nels_per_value-MN_LENGTH) : len;
      size_t num_itr = (len-move_len)/nels_per_value;
      size_t num_off = (len-move_len)%nels_per_value;

      //vec_[idx] = mask(1) & bpv.retrieve_val(idx);
      vec_[idx] = mask(move_len) & bpv.retrieve_val(idx);
      for (size_t i=1; i<=num_itr; i++)
      {
          vec_[idx]=((vec_[idx] << (MN_LENGTH*bsize)) | 
                  extract_pred(bpv.retrieve_val(i), MN_LENGTH, nels_per_value));
          idx++;
          //vec_[idx] = mask(1) & bpv.retrieve_val(i);
          vec_[idx] = mask(nels_per_value-MN_LENGTH) & bpv.retrieve_val(i);
      }
      if (num_off) {
          if (num_off > MN_LENGTH) {
              vec_[idx] = ((vec_[idx] << (MN_LENGTH*bsize)) |
                            extract_pred(bpv.retrieve_val(ori_size-1), MN_LENGTH, num_off));
              idx++;
              vec_[idx] = mask(num_off%MN_LENGTH) & bpv.retrieve_val(ori_size-1);
          }
          else 
              vec_[idx] = ((vec_[idx] << (num_off*bsize)) |
                     extract_pred(bpv.retrieve_val(ori_size-1), num_off, num_off));
      }
      //if (num_off)
      //    vec_[idx] = ((vec_[idx] << (num_off*bsize)) |
      //            extract_pred(bpv.retrieve_val(ori_size-1), num_off, num_off));
      /*
      if (ori_off)
          vec_[idx] = ((vec_[idx] << (ori_off*bsize)) |
                  extract_pred(bpv.retrieve_val(ori_size-1), ori_off, ori_off));
      */
  }

  void mnode_extract_succ(const BasePairVector& bpv, size_t len)
  {
      //assert(bpv.size() >= nels_per_value);
      assert(bpv.size() > MN_LENGTH);
      assert (size_ == MN_LENGTH);
      size_t new_len = len + size_;
      size_t new_size = new_len/nels_per_value + (new_len%nels_per_value ? 1 : 0);
      vec_.resize(new_size);
      size_ = new_len;
      size_t idx=0, i=1;

      size_t move_len = new_len >= nels_per_value ? (nels_per_value-MN_LENGTH) : len;
      size_t size_length = bpv.size() > nels_per_value ? nels_per_value : bpv.size();
      vec_[idx] = (vec_[idx] << (move_len*bsize)) | extract_pred(bpv.retrieve_val(idx), move_len, size_length);
      //vec_[idx] = (vec_[idx] << (move_len*bsize)) | extract_pred(bpv.retrieve_val(idx), move_len, nels_per_value);
      //vec_[idx] = (vec_[idx] << bsize) | extract_pred(bpv.retrieve_val(idx), 1, nels_per_value);
      size_t num_itr = (len-move_len)/nels_per_value;
      size_t num_off = (len-move_len)%nels_per_value;
      //size_t num_itr = (len-1)/nels_per_value;
      //size_t num_off = (len-1)%nels_per_value;

      for (i=1; i<=num_itr; i++)
      {
          vec_[i] = mask(MN_LENGTH) & bpv.retrieve_val(idx);
          idx++;
          vec_[i] = (vec_[i] << ((nels_per_value-MN_LENGTH)*bsize)) | extract_pred(bpv.retrieve_val(idx), (nels_per_value-MN_LENGTH), nels_per_value);
          //vec_[i] = (vec_[i] << bsize) | extract_pred(bpv.retrieve_val(idx), 1, nels_per_value);
      }
      //if (num_off)
      //    vec_[i] = extract_pred(bpv.retrieve_val(idx), num_off, MN_LENGTH);
      if (num_off) {
         if (num_off > MN_LENGTH) {
             vec_[i] = mask(MN_LENGTH) & bpv.retrieve_val(idx);
             idx++;
             vec_[i] = (vec_[i] << ((num_off%MN_LENGTH)*bsize)) | extract_pred(bpv.retrieve_val(idx), (num_off%MN_LENGTH), num_off);
         }
         else
            vec_[i] = extract_pred(bpv.retrieve_val(idx), num_off, MN_LENGTH);
      }

  }

  void append(const BasePairVector& bpv) {
    for(size_t i=0; i<bpv.size(); i++) {
      push_back(bpv[i]);
    }
  }

  T extract(size_t dist) {
      assert (dist < nels_per_value);

      if (size_ < nels_per_value)
          return extract_pred(vec_[0], dist, size_);
      else
      {
      size_t idx_size = size_/nels_per_value;
      size_t idx_pos = dist/nels_per_value;
      size_t size_length = idx_pos < idx_size ? nels_per_value : size_ % nels_per_value;
      return extract_pred(vec_[idx_pos], dist, size_length);
      }
  } 

  T extract_succ(size_t dist) {
      assert (dist < nels_per_value);

      if (size_ < nels_per_value)
          return mask(dist) & vec_[0];
      else
      {
      size_t idx_size = size_/nels_per_value;
      size_t idx_off = size_%nels_per_value;
      //size_t rem = dist - idx_off;
      int rem = dist - idx_off;
      if (rem < 0)
          return mask(dist) & vec_[idx_size];
      
      if (rem > 0) {

         if (rem > nels_per_value) {
             printf("rem: %d, size_: %lu, dist: %lu\n", rem, size_, dist);
         } 
          T partial_kmer = mask(rem) & vec_[idx_size-1];
          if (idx_off)
              return ((partial_kmer << (idx_off*bsize)) | vec_[idx_size]);
          else
              return partial_kmer;
      }
      else
          return vec_[idx_size];
      }
  }

  void populate(T val, size_t len) {
      assert (len < nels_per_value);
      size_ = len;
      vec_.push_back(val);
  }

 friend bool operator==(const BasePairVector& k1, const BasePairVector& k2);

 const std::vector<T>& vec() const {
     return vec_;
     }

 std::vector<T>& vec() {
     return vec_;
     }

  T retrieve_val(size_t idx) const {
      assert (idx < vec_.size());

      return vec_[idx];
  }

  size_t retrieve_size() const {
      return vec_.size();
  }

 protected:
  T mask(size_t off) const {
    assert(off  < nels_per_value);
    return ((T)1 << (off*bsize))-1;
  }

  unsigned int shift(size_t off) const {
    return off*bsize;
  }
  
  std::vector<T> vec_;
  size_t size_;  
};  // BasePairVector

inline bool operator==(const BasePairVector& k1, const BasePairVector& k2)
{
    if (k1.size() != k2.size())
        return false;
    else {
          for (size_t i=0; i<k1.size(); i++) {
               if (k1[i] != k2[i])
                   return false;
          }
    }

    return true;
}

template<typename T>
inline bool
operator == (const std::vector<T>& v1,
             const std::vector<T>& v2) {
  if(v1.size() != v2.size()) {
    return false;
  }
  return std::equal(v1.begin(), v1.end(), v2.begin());
}

/*inline bool
operator == (const BasePairVector& lhs,
             const BasePairVector& rhs) {
  return lhs.size() == rhs.size() &&
      lhs.vec() == rhs.vec();
}
*/

inline char el_to_char(unsigned x) {
      static char symbols[] = {'A', 'C', 'T', 'G', '\0', '\0', 'N'};

        return symbols[(ElType)x];
}

inline kmer_t mnmer_shift(kmer_t kmer_in, 
                        ElType el) {

  //assert(el>=A && el<=G);
  
  return (kmer_t)((kmer_in<<2) | (kmer_t)el) & (kmer_t)MN_MASK;
  //return ((kmer_in<<2) | (kmer_t)el) & KMER_MASK;
}

inline kmer_t mn_extract_pred(kmer_t kmer_in, size_t dist) {

    kmer_t mask = ((kmer_t)1 << (dist*2))-1;
    size_t rem = MN_LENGTH - dist;
    kmer_t pmask = mask << (rem*2);

    return ((pmask & kmer_in) >> (rem*2));
}

inline kmer_t mn_extract_succ(kmer_t kmer_in, size_t dist) {

    kmer_t mask = ((kmer_t)1 << (dist*2))-1;

    return (mask & kmer_in);
}

struct KeyHasher
{
      std::size_t operator()(const BasePairVector& k) const
      {
          using std::size_t;
          using std::hash;
          size_t kmer = 17;
          for (size_t i=0; i<k.size(); i++) {
               kmer += kmer * 31 + hash<uint8_t>()( k[i] );
          }
          return kmer;
      }
};

/*struct KeyEqual
{
      std::size_t operator()(const BasePairVector& k1, const BasePairVector& k2) const
      {
          kmer_t kmer1 = 0, kmer2 = 0;
          for (kmer_t i=0; i<k1.size(); i++) {
               kmer1 = mnmer_shift(kmer1, k1[i]);
          }

          for (kmer_t i=0; i<k2.size(); i++) {
               kmer2 = mnmer_shift(kmer2, k2[i]);
          }
          
          if (kmer1 == kmer2)
              return true;
          else
              return false;
      }
};
*/

class Comp_rev{
    const std::vector<std::pair<int, int>> & _v;
  public:
    Comp_rev(const std::vector<std::pair<int, int>> & v) : _v(v) {}
    bool operator()(size_t i, size_t j){
         return ((_v[i].second > _v[j].second) ||
                 ((_v[i].second == _v[j].second) &&
                  ( _v[i].first > _v[j].first))
                );

   }
};

struct ContigInfo {
    BasePairVector contig_data;
    std::string contig_name;

    ContigInfo()=default;

};

struct WireInfo {
  int suffix_id;
  int offset_in_suffix;
  int count;

  WireInfo()=default;

#ifdef USE_CEREAL
  template<class Archive>
  void serialize(Archive & archive)
  {
      archive(suffix_id, offset_in_suffix, count);
  }
#endif

};

struct PrefixInfo {
  int prefix_pos;
  int num_wires;

  PrefixInfo()=default;
  
  #ifdef USE_CEREAL
  template<class Archive>
  void serialize(Archive & archive)
  {
      archive(prefix_pos, num_wires);
  }
  #endif

};

struct MacroNode {
  std::vector<BasePairVector> suffixes;
  std::vector<bool> suffixes_terminal;
  std::vector<std::pair<int, int>> suffix_count;
  std::vector<BasePairVector> prefixes;
  std::vector<bool> prefixes_terminal;
  std::vector<std::pair<int, int>> prefix_count;
  BasePairVector k_1_mer;
  std::vector<WireInfo> wiring_info;
  std::vector<PrefixInfo> prefix_begin_info;

  //MacroNode() : suffixes(0), suffixes_terminal(0), suffix_count(0), 
  //              prefixes(0), prefixes_terminal(0), prefix_count(0),
  //              k_1_mer(BasePairVector()), wiring_info(0) {}
  MacroNode() = default;

#ifdef USE_CEREAL
  template<class Archive>
  void serialize(Archive & archive)
  {
      archive(suffixes, suffixes_terminal, suffix_count,
              prefixes, prefix_count, prefixes_terminal,
              k_1_mer, 
              wiring_info,
              prefix_begin_info);
  }
#endif


  void setup_wiring() {
    int pc = 0, sc = 0;
    int null_suffix_id = -1, null_prefix_id = -1;

    for(size_t i=0; i<suffixes.size(); i++) {
      if(suffixes[i].size() != 0) {
        sc += suffix_count[i].second;
      }
    }
    for(size_t i=0; i<prefixes.size(); i++) {
      if(prefixes[i].size() != 0) {
        pc += prefix_count[i].second;
      }
    }
    for(size_t i=0; i<suffixes.size(); i++) {
      if(suffixes[i].size() == 0) {
        assert(suffixes_terminal[i]);
        assert(null_suffix_id == -1);
        null_suffix_id = i;
        if (suffix_count[i].first == -1) {
            suffix_count[i] = std::make_pair(1, std::max(pc-sc,0));
            //printf("num_suffix: %lu, null_suffix_id: %d\n", suffixes.size(), null_suffix_id); 
            assert(null_suffix_id == suffixes.size()-1);
        }
      }
    }
    for(size_t i=0; i<prefixes.size(); i++) {
      if(prefixes[i].size() == 0) {
        assert(prefixes_terminal[i]);
        assert(null_prefix_id == -1);
        null_prefix_id = i;
        if (prefix_count[i].first == -1) {
            prefix_count[i] = std::make_pair(1, std::max(sc-pc,0));
            //printf("num_prefix: %lu, null_prefix_id: %d\n", prefixes.size(), null_prefix_id); 
            assert(null_prefix_id == prefixes.size()-1);
        }
      }
    }
 
#ifdef DEBUG_SERMOD
    if (sc + suffix_count[null_suffix_id].second != pc + prefix_count[null_prefix_id].second)
    {
         kmer_t test_mn=0;
         for (int l=0; l<MN_LENGTH; l++)
              test_mn = mnmer_shift(test_mn, k_1_mer[l]);

         std::string out_mn;
         for (size_t l=0; l<k_1_mer.size(); l++)
              out_mn += el_to_char(k_1_mer[l]);
 
        printf("Wiring error for node: %s, node_id: %lu\n", out_mn.c_str(), test_mn);
    }
#endif    
    assert(sc + suffix_count[null_suffix_id].second == pc + prefix_count[null_prefix_id].second);
    assert(suffix_count[null_suffix_id].second==0 || prefix_count[null_prefix_id].second==0);

    std::vector<int> indices_p(prefix_count.size());
    std::vector<int> indices_s(suffix_count.size());
    std::iota(indices_p.begin(), indices_p.end(), 0);
    std::iota(indices_s.begin(), indices_s.end(), 0);

    std::sort(indices_p.begin(), indices_p.end(), Comp_rev(prefix_count));
    std::sort(indices_s.begin(), indices_s.end(), Comp_rev(suffix_count));

    int leftover = sc + suffix_count[null_suffix_id].second;
    std::vector<int> offset_in_suffix(suffixes.size(), 0);
    int top_p=0, top_s=0;
    int var_p=0, var_s=0;
    int wire_idx=0, p_size=0, last_largest_pid=-1, prefix_begin_pos=-1;

    while(leftover > 0) {
      int largest_sid = indices_s[top_s];
      int largest_pid = indices_p[top_p];

      int count = std::min((prefix_count[largest_pid].second-var_p),
                           (suffix_count[largest_sid].second-var_s));
      
      wiring_info[wire_idx] = WireInfo{largest_sid,
              offset_in_suffix[largest_sid],
              count};

      if (last_largest_pid != largest_pid) {
          prefix_begin_pos = wire_idx;
          last_largest_pid = largest_pid;
      }

      p_size++;
      wire_idx++;
      leftover -= count;
      var_p += count;
      var_s += count;
      offset_in_suffix[largest_sid] += count;

      if (var_p == prefix_count[largest_pid].second) {
          var_p = 0;
	  top_p++;
          prefix_begin_info[largest_pid] = PrefixInfo{prefix_begin_pos, p_size};
          p_size=0;
      }

      if (var_s == suffix_count[largest_sid].second) {
          var_s = 0;
	  top_s++;
      }

    }
  }

};

/*
inline bool
operator == (const MacroNode& mn1, const MacroNode& mn2) {
  return mn1.k_1_mer == mn2.k_1_mer  && 
      mn1.suffixes == mn2.suffixes && 
      mn1.suffixes_terminal == mn2.suffixes_terminal && 
      mn1.suffix_count == mn2.suffix_count &&
      mn1.prefixes == mn2.prefixes && 
      mn1.prefix_count == mn2.prefix_count && 
      mn1.prefixes_terminal == mn2.prefixes_terminal
      ;
}
*/

typedef struct begin_kmer_id
{
    kmer_t node;
    int terminal_prefix_id;

} BeginMN;

struct MnodeInfo {
    kmer_t search_mn;
    BasePairVector search_ext;

    MnodeInfo()=default;

#ifdef USE_CEREAL
    template<class Archive>
    void serialize(Archive & archive)
    {
       archive(search_mn, search_ext);
    }
#endif

};

struct TransferNode {
  kmer_t search_mn;
  BasePairVector search_ext;
  BasePairVector mn_ext;
  std::pair<int, int> mn_count;
  //KTerminal mn_terminal;
  bool mn_terminal;
  KdirType direction;

  TransferNode()=default;

};

typedef struct __attribute__ ((__packed__)) kmer_pairs
{ 
  kmer_t seq;
  int k_count;

} KmerPairs;
static_assert(sizeof(KmerPairs) == (sizeof(kmer_t)+sizeof(int)), "struct KmerPairs shouldn't be padded");

//data structure for storing Reads
typedef struct Rd_Sequence
{
	// Read data
	char *read_data;
	size_t read_data_size;

} input_read_data;

typedef struct each_contig_entry
{
       long int my_contig_id;
       std::string contig_data;

} contig_entry;

typedef struct contig_series
{
       int contig_count;
       std::vector < contig_entry> c_series;

} contig_thrd_list;

inline bool Comp_pair(const std::pair<kmer_t,MacroNode> &a, 
               const std::pair<kmer_t,MacroNode> &b) {
        return a.first < b.first;

}


class Comp{
    const std::vector<kmer_t> & _v;
  public:
    Comp(const std::vector<kmer_t> & v) : _v(v) {}
    bool operator()(size_t i, size_t j){
         return _v[i] < _v[j];
   }
};

inline ElType kmerel(kmer_t kmer, unsigned pos) {
  assert(pos < KMER_LENGTH);

  return ElType(((kmer) >> ((KMER_LENGTH-(pos)-1)*2)) & (0x3));
}

inline ElType lmerel(lmer_t kmer, unsigned pos) {
  assert(pos < LMER_LENGTH);

  return ElType(((kmer) >> ((LMER_LENGTH-(pos)-1)*2)) & (0x3));
}

inline ElType kmerel_mn(kmer_t kmer, unsigned pos) {
      assert(pos < MN_LENGTH);

        return ElType(((kmer) >> ((MN_LENGTH-(pos)-1)*2)) & (0x3));
}


inline lmer_t kmer_to_lmer(kmer_t kmer_in, unsigned pos, lmer_t kmer_out) {
  assert(pos < KMER_LENGTH);

  //ElType int_el = ElType(((kmer_in) >> ((KMER_LENGTH-(pos)-1)*2)) & (0x3));
  //assert(int_el>=A && int_el<=G);
  //return (kmer_t)((kmer_out<<2) | (kmer_t)int_el) & (kmer_t)LMER_MASK;

  return (lmer_t)((kmer_out<<2) | (lmer_t)(ElType(((kmer_in) >> ((KMER_LENGTH-(pos)-1)*2)) & (0x3)))) & (lmer_t)LMER_MASK;
}

inline ElType char_to_el(char ch) {
              return (ElType)((((ElType)ch)>>1) & 0x7);
}


inline kmer_t kmer_cons(kmer_t kmer_in, 
                        unsigned pos,
                        ElType el) {
  //assert(el>=A && el<=G);
  assert(pos < KMER_LENGTH);
 
  return kmer_in | ((kmer_t)el << ((KMER_LENGTH-pos-1)*2));
}

inline kmer_t kmer_shift(kmer_t kmer_in, 
                        ElType el) {

  //assert(el>=A && el<=G);
  
  return (kmer_t)((kmer_in<<2) | (kmer_t)el) & (kmer_t)KMER_MASK;
  //return ((kmer_in<<2) | (kmer_t)el) & KMER_MASK;
}

inline lmer_t lmer_shift(lmer_t kmer_in, 
                        ElType el) {

  //assert(el>=A && el<=G);
  
  return (lmer_t)((kmer_in<<2) | (lmer_t)el) & (lmer_t)LMER_MASK;
  //return ((kmer_in<<2) | (kmer_t)el) & KMER_MASK;
}

inline kmer_t mn_shift(kmer_t kmer_in)
{
    return (kmer_t)((kmer_in>>2)) & (kmer_t)MN_MASK;
}

inline kmer_t tokmer(const char *kmer_str,
                     int kmer_len) {
  int i;
  assert(kmer_len == KMER_LENGTH);
  assert(kmer_str != NULL);
  assert(kmer_str[kmer_len] == '\0');

  kmer_t km = 0;
  for(i=0; i<kmer_len; i++) {
    ElType el = char_to_el(kmer_str[i]);
    km = kmer_cons(km, i, el);
  }  
  return km;
}

template <typename T> 
inline long uhash31( uint64_t a, uint64_t b, T x)
{

  T result;
  long lresult;  

  // return a hash of x using a and b mod (2^31 - 1)
// may need to do another mod afterwards, or drop high bits
// depending on d, number of bad guys
// 2^31 - 1 = 2147483647

  //  result = ((long long) a)*((long long) x)+((long long) b);
  result=(a * x) + b;
  result = ((result >> HL) + result) & MOD;
  lresult=(long) result; 
  
  return(lresult);
}

inline long hash31(long long a, long long b, long long x)
{

  long long result;
  long lresult;  

  // return a hash of x using a and b mod (2^31 - 1)
// may need to do another mod afterwards, or drop high bits
// depending on d, number of bad guys
// 2^31 - 1 = 2147483647

  //  result = ((long long) a)*((long long) x)+((long long) b);
  result=(a * x) + b;
  result = ((result >> HL) + result) & MOD;
  lresult=(long) result; 
  
  return(lresult);
}

/*kmer_t cvt_inv(char* lmer) {
    kmer_t l_num=0;

    for (int i = 0; i<LMER_LENGTH; i++) {
         l_num |= ( (kmer_t)char_to_el(lmer[i]) << ((LMER_LENGTH-i-1)*2 ) );
    }
    return l_num;
}*/
 
int hash_fcn (const char* word, unsigned M);
void parse_alphabet_file (FILE * fp);
void free_reads (int num_reads);
bool is_sigma (char c);
void change_to_num (char *pred_kmer, int *num);
void change_to_char (char *pred_kmer, int num);
int find_max (int array[]);
int find_min (int array[], int *min);
int compare (const void * a, const void * b);
int FindIndex( const int a[], int size, int value);
char* DivideReads(MPI_File *in, const int rank, const int size, 
        const int overlap, uint64_t *nlines, size_t *data_size);
//input_read_data perform_input_reading (const int rank, 
//        const int size, char *argv[]);
input_read_data perform_input_reading (const int rank, const int size,
                                       std::string &fileName, int read_length);
void perform_kmer_counting(const char *ptr, size_t length);
//void SortAndAggregate(std::vector<kmer_t>& arr, std::vector<int>& count);
void SortAndAggregate(std::vector<KmerPairs>& arr);
void sort_recv_buffer(std::vector<KmerPairs>& kmer_recv_buf, 
                      std::vector<int>& rcounts_kmer, 
                      std::vector<int>& rdisp_kmer);
void Sliding_window_l (const char *ptr, size_t length);
void Sliding_window (const char *ptr, size_t length, int *n_kmers, 
                     std::vector<std::vector<kmer_t>> &partial_kmer_counts);
//void Sliding_window (const char *ptr, size_t length, int *n_kmers, 
//                     std::vector<std::vector<kmer_t>> &kmers_per_proc, 
//                     std::vector<std::vector<int>> &kmer_cnt_tmp_buf);
//void transfer_kmers (std::vector<int>& scounts_kmer, 
//                     std::vector<std::vector<KmerPairs>> &kmers_per_proc); 
void transfer_kmers (std::vector<int>& scounts_kmer, 
                     std::vector<KmerPairs> &kmers_per_proc); 
//void transfer_kmers (int *scounts_kmer, 
//                     std::vector<std::vector<kmer_t>> &kmers_per_proc, 
//                     std::vector<std::vector<int>> &kmer_cnt_tmp_buf);
//
void process_remaining_kmers(
                     std::vector<std::vector<kmer_t>> &partial_kmer_counts); 
//void process_remaining_kmers(
//                     std::vector<std::vector<kmer_t>> &kmers_per_proc, 
//                     std::vector<std::vector<int>> &kmer_cnt_tmp_buf);
//
void free_kmer_count_buffers();
int convert_hash_fcn (const char* word);
unsigned int hash_str(const char *str);
std::vector <KmerPairs> construct_kmer_histogram (std::vector< std::vector< std::vector<KmerPairs> > > &global_mn_list);
std::vector<KmerPairs> transfer_macro_nodes(std::vector<KmerPairs> &mn_send_buf, 
                                            std::vector<int> &scounts);
                                            //size_t ssize);
void begin_mnode_construction(std::vector<std::pair<kmer_t,MacroNode>> &MN_map);
void construct_macro_nodes (std::vector<std::pair<kmer_t,MacroNode>> &MN_map, 
        std::vector<KmerPairs> &all_kmers_from_procs);
void debug_wired_mnodes(std::vector<std::pair<kmer_t,MacroNode>>& MN_map);
int retrieve_proc_id (lmer_t min_lmer);
int retrieve_proc_id (kmer_t min_lmer);
//template <typename T> 
//int retrieve_proc_id (T min_lmer);
void initiate_mnode_wiring(std::vector<std::pair<kmer_t,MacroNode>> &MN_map);
void check(MacroNode& mn);
size_t begin_iterative_compaction (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                                   std::vector<BasePairVector> &partial_contig_list);
void identify_begin_kmers (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                               std::vector<BeginMN> &list_of_begin_kmers);
void generate_compacted_pakgraph(std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                                 std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map);
void gather_all_macro_nodes (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
       std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map);
int  find_mnode_exists (kmer_t key, std::vector<std::pair<kmer_t,MacroNode>> &MN_map);
void walk(BasePairVector &partial_contig, int freq,
        int offset_in_prefix, int prefix_id, const MacroNode &mn,
        std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map,
        std::vector<BasePairVector> &local_contig_list);
        //std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map, FILE * pFile);

void traverse_pakgraph (std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map,
                        std::vector<BeginMN> &list_of_begin_kmers,
                        std::vector<BasePairVector> &partial_contig_list);
void check_for_self_loops (kmer_t search_pkmer, kmer_t search_skmer,
        kmer_t node, std::array<bool, 2>& self_loop_info);
//void generate_id_set (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
//        std::vector<BeginMN> &id_list_nodes);
void generate_id_set (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
        std::vector<size_t> &id_list_nodes);
MnodeInfo retrieve_mn_sinfo(BasePairVector &succ_ext, BasePairVector &key);
MnodeInfo retrieve_mn_pinfo(BasePairVector &succ_ext, BasePairVector &key);
/*
void push_to_pred (ModNodeInfo &new_mnode, int pos, 
        std::vector<std::pair<kmer_t,MacroNode>> &MN_map);
void push_to_succ (ModNodeInfo &new_mnode, int pos, 
        std::vector<std::pair<kmer_t,MacroNode>> &MN_map);
void iterate_and_pack_mn (std::vector<size_t>& id_list_nodes,
                          std::vector< std::vector<ModNodeInfo> >& mn_nodes_per_proc,
                          std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                          std::vector<BasePairVector> &local_contig_list);
*/
void push_to_pred (TransferNode &new_mnode, int pos, 
        std::vector<std::pair<kmer_t,MacroNode>> &MN_map);
void push_to_succ (TransferNode &new_mnode, int pos, 
        std::vector<std::pair<kmer_t,MacroNode>> &MN_map);
void output(BasePairVector &partial_contig);

void print_id_list (std::vector<size_t>& id_list_nodes , const std::string& func, std::vector<std::pair<kmer_t,MacroNode>>& MN_map);
void iterate_and_pack_mn (std::vector<size_t>& id_list_nodes,
                          std::vector< std::vector<TransferNode> >& mn_nodes_per_proc,
                          std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                          std::vector<BasePairVector> &local_contig_list);
//std::vector<ModNodeInfo> serialize_and_transfer 
//                         (std::vector< std::vector<ModNodeInfo> >& mn_nodes_per_proc);
/*void serialize_and_transfer 
                         (std::vector< std::vector<ModNodeInfo> >& mn_nodes_per_proc,
                          std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                          std::vector<size_t> &rewire_pos_list,
                          int num_itr);
*/
void serialize_and_transfer 
                         (std::vector< std::vector<TransferNode> >& mn_nodes_per_proc,
                          std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                          std::vector<size_t> &rewire_pos_list,
                          int num_itr);

void print_output_contigs(std::vector<BasePairVector>& partial_contig_list,
                  const std::string& func);
void print_map_content_itr (std::vector<std::pair<kmer_t,MacroNode>>& MN_map, const std::string& func, int num_itr);
void print_rewire_list (std::vector<size_t>& rewire_pos_list, 
                        std::vector<std::pair<kmer_t,
                        MacroNode>>& MN_map, 
                        const std::string& func, int num_itr);
//////////////////////////////////////////////////////////////
#endif
