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

#ifndef SERIALIZE_H
#define SERIALIZE_H

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
#include <iostream>
#include <sstream>
#include <random>
#include <utility>
#include "distribute_kmers.h"
#include "serialize.h"

extern int rank, size;

template<typename POD>
std::ostream& serialize(std::ostream& os, const POD& v)
{
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
        "Can only serialize POD types with this function");

  //printf("in POD, sizeof(POD): %lu\n", sizeof(v));
  os.write(reinterpret_cast<const char*>(&v), sizeof(POD));
  return os;
}

template<typename POD1, typename POD2>
std::ostream& serialize(std::ostream& os, const std::pair<POD1, POD2>& p)
{
    static_assert(std::is_trivial<POD1>::value && std::is_standard_layout<POD1>::value,
        "Can only serialize POD types with this function");

    static_assert(std::is_trivial<POD2>::value && std::is_standard_layout<POD2>::value,
        "Can only serialize POD types with this function");

   serialize(os, std::get<0>(p));
   serialize(os, std::get<1>(p));
   return os;

}

template<typename POD>
std::ostream& serialize(std::ostream& os, std::vector<POD> const& vec)
{
    // this only works on built in data types (PODs)
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
        "Can only serialize POD types with this function");

    auto size = vec.size();
    os.write(reinterpret_cast<char const*>(&size), sizeof(size));
    for (auto& v:  vec) {
         serialize(os, v);
    }
    //os.write(reinterpret_cast<char const*>(v.data()), v.size() * sizeof(POD));
    return os;
}

std::ostream& serialize(std::ostream& os, const BasePairVector &bp)
{
  auto size=bp.size();
  os.write(reinterpret_cast<char const*>(&size), sizeof(size));
  return serialize(os, bp.vec());
  //return os;
}

std::ostream& serialize(std::ostream& os, std::vector<BasePairVector> const& bvec)
{
  auto size = bvec.size();
  os.write(reinterpret_cast<char const*>(&size), sizeof(size));
  for (auto& v:  bvec) {
       serialize(os, v);
  }
  return os;

}

std::ostream& serialize(std::ostream& os, std::vector<std::pair<int, int>> const& pr)
{
  auto size = pr.size();
  os.write(reinterpret_cast<char const*>(&size), sizeof(size));
  for (auto& v:  pr) {
       serialize(os, v);
  }
  return os;

}

std::ostream& serialize(std::ostream& os, std::vector<bool> const& tvec)
{
  auto size = tvec.size();
  os.write(reinterpret_cast<char const*>(&size), sizeof(size));
  //os.write(reinterpret_cast<char const*>(tvec.data()), tvec.size() * sizeof(bool));
  /*for (const auto& v: tvec) {
       serialize(os, v);
  }
  */
  for(std::vector<bool>::size_type i = 0; i < size;)
  {
        unsigned char aggr = 0;
        for(unsigned char mask = 1; mask > 0 && i < size; ++i, mask <<= 1)
            if(tvec.at(i))
                aggr |= mask;
        os.write(reinterpret_cast<char const*>(&aggr), sizeof(unsigned char));
  }
  return os;

}

std::ostream& serialize(std::ostream& os, const WireInfo &winfo)
{
  serialize(os, winfo.suffix_id);
  serialize(os, winfo.offset_in_suffix);
  serialize(os, winfo.count);
  return os;

}


std::ostream& serialize(std::ostream& os, std::vector<WireInfo> const& wi)
{
  auto size = wi.size();
  os.write(reinterpret_cast<char const*>(&size), sizeof(size));
  for (auto& v: wi) {
       serialize(os, v);
  }
  return os;

}

std::ostream& serialize(std::ostream& os, const PrefixInfo &pi)
{
  serialize(os, pi.prefix_pos);
  serialize(os, pi.num_wires);
  return os;

}

std::ostream& serialize(std::ostream& os, std::vector<PrefixInfo> const& vp)
{
  auto size = vp.size();
  os.write(reinterpret_cast<char const*>(&size), sizeof(size));
  for (auto& v:  vp) {
       serialize(os, v);
  }
  return os;

}

std::ostream& serialize(std::ostream& os, const TransferNode& tn)
{
  serialize(os, tn.search_mn);
  serialize(os, tn.search_ext);
  serialize(os, tn.mn_ext);
  serialize(os, tn.mn_count);
  serialize(os, tn.mn_terminal);
  serialize(os, tn.direction);
  return os;

}

std::ostream& serialize(std::ostream& os, const MacroNode& mn)
{
  serialize(os, mn.suffixes);
  serialize(os, mn.suffixes_terminal);
  serialize(os, mn.suffix_count);
  serialize(os, mn.prefixes);
  serialize(os, mn.prefixes_terminal);
  serialize(os, mn.prefix_count);
  serialize(os, mn.k_1_mer);
  serialize(os, mn.wiring_info);
  serialize(os, mn.prefix_begin_info);
  return os;

}

template<typename POD>
std::istream& deserialize(std::istream& is, POD& v)
{
  static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
        "Can only deserialize POD types with this function");

  is.read(reinterpret_cast<char*>(&v), sizeof(POD));
  return is;
}

template<typename POD1, typename POD2>
std::istream& deserialize(std::istream& is, std::pair<POD1, POD2>& p)
{
  static_assert(std::is_trivial<POD1>::value && std::is_standard_layout<POD1>::value,
        "Can only deserialize POD types with this function");

  static_assert(std::is_trivial<POD2>::value && std::is_standard_layout<POD2>::value,
        "Can only deserialize POD types with this function");

  deserialize(is, std::get<0>(p));
  deserialize(is, std::get<1>(p));
  return is;

}

template<typename POD>
std::istream& deserialize(std::istream& is, std::vector<POD>& vec)
{
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
        "Can only deserialize POD types with this function");

    decltype(vec.size()) size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    vec.resize(size);
    for (auto& v: vec) {
         deserialize(is, v);
    }
    //is.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(POD));
    return is;
}

std::istream& deserialize(std::istream& is, BasePairVector &bp)
{
  decltype(bp.size()) size;
  is.read(reinterpret_cast<char*>(&size), sizeof(size));
  bp.resize(size);
  deserialize(is, bp.vec());
  return is;
}

std::istream& deserialize(std::istream& is, std::vector<BasePairVector>& bvec)
{
  decltype(bvec.size()) size;
  is.read(reinterpret_cast<char*>(&size), sizeof(size));
  bvec.resize(size);
  for (auto& v: bvec) {
       deserialize(is, v);
  }
  return is;

}

std::istream& deserialize(std::istream& is, std::vector<std::pair<int, int>>& pr)
{
  decltype(pr.size()) size;
  is.read(reinterpret_cast<char*>(&size), sizeof(size));
  pr.resize(size);
  for (auto& v: pr) {
       deserialize(is, v);
  }
  return is;

}

std::istream& deserialize(std::istream& is, std::vector<bool>& tvec)
{
  decltype(tvec.size()) size;
  is.read(reinterpret_cast<char*>(&size), sizeof(size));
  tvec.resize(size);
  //is.read(reinterpret_cast<char*>(tvec.data()), tvec.size() * sizeof(bool));
  /*
  for (const auto& v: tvec) {
       deserialize(is, v);
  }
  */
  for(std::vector<bool>::size_type i = 0; i < size;)
  {
        unsigned char aggr;
        is.read(reinterpret_cast<char*>(&aggr), sizeof(unsigned char));
        for(unsigned char mask = 1; mask > 0 && i < size; ++i, mask <<= 1)
            tvec.at(i) = aggr & mask;
  }
  return is;

}

std::istream& deserialize(std::istream& is, WireInfo &winfo)
{
 deserialize(is, winfo.suffix_id);
 deserialize(is, winfo.offset_in_suffix);
 deserialize(is, winfo.count);
 return is;

}

std::istream& deserialize(std::istream& is, std::vector<WireInfo>& wi)
{
  decltype(wi.size()) size;
  is.read(reinterpret_cast<char*>(&size), sizeof(size));
  wi.resize(size);
  for (auto& v: wi) {
       deserialize(is, v);
  }
  return is;

}

std::istream& deserialize(std::istream& is, PrefixInfo &pi)
{
  deserialize(is, pi.prefix_pos);
  deserialize(is, pi.num_wires);
  return is;

}

std::istream& deserialize(std::istream& is, std::vector<PrefixInfo>& vp)
{
  decltype(vp.size()) size;
  is.read(reinterpret_cast<char*>(&size), sizeof(size));
  vp.resize(size);
  for (auto& v: vp) {
       deserialize(is, v);
  }
  return is;

}

std::istream& deserialize(std::istream& is, TransferNode& tn)
{
  deserialize(is, tn.search_mn);
  deserialize(is, tn.search_ext);
  deserialize(is, tn.mn_ext);
  deserialize(is, tn.mn_count);
  deserialize(is, tn.mn_terminal);
  deserialize(is, tn.direction);
  return is;

}

std::istream& deserialize(std::istream& is, MacroNode& mn)
{
  deserialize(is, mn.suffixes);
  deserialize(is, mn.suffixes_terminal);
  deserialize(is, mn.suffix_count);
  deserialize(is, mn.prefixes);
  deserialize(is, mn.prefixes_terminal);
  deserialize(is, mn.prefix_count);
  deserialize(is, mn.k_1_mer);
  deserialize(is, mn.wiring_info);
  deserialize(is, mn.prefix_begin_info);
  return is;

}

#ifdef CHAR_SERIALIZE

template<typename T>
size_t serialize_size_pod(const T &v) {
   printf("In size_pod rank: %d\n", rank);
       return sizeof(T);
}

template<typename T>
size_t serialize_size_vec(const std::vector<T>& vec) {
   printf("In size_vec rank: %d\n", rank);
       size_t ret=sizeof(vec.size());
       for (auto& v:  vec) {
            ret += serialize_size_pod(v);
       }
       return ret;
}

template<typename T1,
         typename T2>
size_t serialize_size_pair(const std::pair<T1, T2>& p) {
   printf("In size_pair rank: %d\n", rank);
  return serialize_size_pod(std::get<0>(p)) + serialize_size_pod(std::get<1>(p));
}

size_t serialize_size_bp(const BasePairVector &bp) {
   printf("In size_bp rank: %d\n", rank);
       return sizeof(bp.size()) + serialize_size_vec(bp.vec());
}

size_t serialize_size_TransferNode(const TransferNode& tn) {
   printf("In size_TransferNode rank: %d\n", rank);
   return serialize_size_pod(tn.search_mn) +
          serialize_size_bp(tn.search_ext)+
          serialize_size_bp(tn.mn_ext) +
          serialize_size_pair(tn.mn_count) +
          serialize_size_pod(tn.mn_terminal) +
          serialize_size_pod(tn.direction);
}

/*
inline size_t
ser_size(const MacroNode& mn) {
  return ser_size(mn.k_1_mer) +
      ser_size(mn.prefixes.size()) + // num prefixes
      ser_size(mn.suffixes.size()) + //num suffixes
      ser_size(mn.prefixes_terminal) +
      ser_size(mn.prefix_count) +
      ser_size(mn.suffixes_terminal) +
      ser_size(mn.suffix_count) +
      ser_size(mn.prefixes) +
      ser_size(mn.suffixes);
}    
*/

/*
size_t serialize_size_Info(const MnodeInfo& mn) {
    return  serialize_size(mn.search_mn) +
        serialize_size(mn.search_ext);
}

size_t serialize_size_Tuple(const MNnodeTuple& mn) {
    return serialize_size(mn.mn_ext) +
           serialize_size(mn.mn_count) +
           serialize_size(mn.mn_terminal);
}

size_t serialize_size_Node(const ModNodeInfo& mn) {
    return serialize_size_Info(mn.modifyMN) +
           serialize_size_Tuple(mn.modifyMN_info) +
           serialize_size(mn.direction);
}
*/

/*
size_t serialize_size_TransferNode(const TransferNode& tn) {
   return serialize_size(tn.search_mn) +
          serialize_size(tn.search_ext) +
          serialize_size(tn.mn_ext) +
          serialize_size(tn.mn_count) +
          serialize_size(tn.mn_terminal) +
          serialize_size(tn.direction);
}
*/
//////////////////////////////////////////////////////////////


template<typename T>
char* serialize_pod(char* ptr1, const T& val) {
   printf("In POD rank: %d, bVal ptr: %lu\n", rank, ptr1);
   if (ptr1!=nullptr)
       *(T*) ptr1 = val;
   printf("rank: %d, aVal ptr: %lu\n", rank, ptr1);
   return ptr1 + sizeof(val);
}

template<typename T1,
          typename T2>
char* serialize_pair(char *ptr2, const std::pair<T1, T2>& p) {
   printf("In Pair rank: %d, before ptr: %lu\n", rank, ptr2);
   ptr2 = serialize_pod(ptr2, std::get<0>(p));
   ptr2 = serialize_pod(ptr2, std::get<1>(p));
   return ptr2;
}

template<typename T>
char* serialize_vec(char* ptr3, const std::vector<T>& vec) {
   printf("In vector rank: %d, before ptr: %lu, vec size: %lu\n", rank, ptr3, vec.size());
      ptr3 = serialize_pod(ptr3, vec.size());  
      for (auto& v:  vec) {
           ptr3 = serialize_pod(ptr3, v);
      }
      return ptr3;
}

char* serialize_bp(char* ptr4, const BasePairVector &bp) {
      printf("In BasePairVec rank: %d, size: %lu, before ptr: %lu\n", rank, bp.size(), ptr4);
      ptr4 = serialize_pod(ptr4, bp.size());
      return serialize_vec(ptr4, bp.vec());
}

char* serialize_TransferNode(char* ptr, const TransferNode& tn) {
   printf("In TransferNode rank: %d, before ptr: %lu\n", rank, ptr);
   ptr = serialize_pod(ptr, tn.search_mn);
   ptr = serialize_bp(ptr, tn.search_ext);
   ptr = serialize_bp(ptr, tn.mn_ext);
   ptr = serialize_pair(ptr, tn.mn_count);
   ptr = serialize_pod(ptr, tn.mn_terminal);
   ptr = serialize_pod(ptr, tn.direction);
   return ptr;
}

/*
char* serialize_Info(char* ptr,
    const MnodeInfo& mnode) {
   printf("In ser Info rank: %d, before ptr: %lu\n", rank, ptr);
   ptr = serialize(ptr, mnode.search_mn);
   //printf("rank: %d, after ptr: %lu\n", rank, ptr);
   ptr = serialize(ptr, mnode.search_ext);
   printf("rank: %d, Ext ptr: %lu\n", rank, ptr);
   return ptr;
}

char* serialize_Tuple(char* ptr,
    const MNnodeTuple& mtup) {
       printf("In ser Tuple rank: %d, before ptr: %lu\n", rank, ptr);
   ptr = serialize(ptr, mtup.mn_ext);
   ptr = serialize(ptr, mtup.mn_count);
   ptr = serialize(ptr, mtup.mn_terminal);
   return ptr;
}

char* serialize_Node(char* ptr,
    const ModNodeInfo& minfo) {
       printf("In ser Node rank: %d, before ptr: %lu\n", rank, ptr);
   ptr = serialize_Info(ptr, minfo.modifyMN);
   ptr = serialize_Tuple(ptr, minfo.modifyMN_info);
   ptr = serialize(ptr, minfo.direction);
   return ptr;
}

char* serialize_TransferNode(char* ptr, const TransferNode& tn) {
   printf("In TransferNode rank: %d, before ptr: %lu\n", rank, ptr);
   ptr = serialize(ptr, tn.search_mn);
   ptr = serialize(ptr, tn.search_ext);
   ptr = serialize(ptr, tn.mn_ext);
   ptr = serialize(ptr, tn.mn_count);
   ptr = serialize(ptr, tn.mn_terminal);
   ptr = serialize(ptr, tn.direction);
   return ptr;
}
*/

//////////////////////////////////////////////////////////////

template<typename T>
char* deserialize_pod(char* ptr, T& val) {
   printf("In DES POD rank: %d, bVal ptr: %lu\n", rank, ptr);
    val = *(T*) ptr;
    return ptr + sizeof(val);
}

template<typename T1,
          typename T2>
char* deserialize_pair(char *ptr, std::pair<T1, T2>& p) {
   printf("In DES Pair rank: %d, before ptr: %lu\n", rank, ptr);
   ptr = deserialize_pod(ptr, std::get<0>(p));
   ptr = deserialize_pod(ptr, std::get<1>(p));
   return ptr;
}

template<typename T>
char* deserialize_vec(char* ptr, std::vector<T>& vec) {
     decltype(vec.size()) sz;
     printf("In DES vector rank: %d, before ptr: %lu, vec size: %lu\n", rank, ptr, sz);
     ptr = deserialize_pod(ptr, sz);
     //printf("de-sz: %d\n", sz);     
     vec.resize(sz);
     for (auto& v: vec) {
          ptr = deserialize_pod(ptr, v);
     }
     return ptr;
}

char* deserialize_bp(char* ptr, BasePairVector &bp) {
     decltype(bp.size()) sz;
     printf("In DES BasePairVec rank: %d, before ptr: %lu\n", rank, ptr);
     ptr = deserialize_pod(ptr, sz);
     bp.resize(sz);
     printf("bp_size: %lu\n", sz);
     ptr = deserialize_vec(ptr, bp.vec());
     return ptr;
}

char* deserialize_TransferNode(char* ptr,
      TransferNode& tn) {
   printf("In DES TransferNode rank: %d, before ptr: %lu\n", rank, ptr);
   ptr = deserialize_pod(ptr, tn.search_mn);
   ptr = deserialize_bp(ptr, tn.search_ext);
   ptr = deserialize_bp(ptr, tn.mn_ext);
   ptr = deserialize_pair(ptr, tn.mn_count);
   ptr = deserialize_pod(ptr, tn.mn_terminal);
   ptr = deserialize_pod(ptr, tn.direction);
   return ptr;
}

/*
char*
deserialize_Info(char* ptr,
      MnodeInfo& mnode) {
   ptr = deserialize(ptr, mnode.search_mn);
   ptr = deserialize(ptr, mnode.search_ext);
   return ptr;
}

char*
deserialize_Tuple(char* ptr,
      MNnodeTuple& mtup) {
   ptr = deserialize(ptr, mtup.mn_ext);
   ptr = deserialize(ptr, mtup.mn_count);
   ptr = deserialize(ptr, mtup.mn_terminal);
   return ptr;
}

char*
deserialize_Node(char* ptr,
      ModNodeInfo& minfo) {
   ptr = deserialize_Info(ptr, minfo.modifyMN);
   ptr = deserialize_Tuple(ptr, minfo.modifyMN_info);
   ptr = deserialize(ptr, minfo.direction);
   return ptr;
}

char* deserialize_TransferNode(char* ptr,
      TransferNode& tn) {
   ptr = deserialize(ptr, tn.search_mn);
   ptr = deserialize(ptr, tn.search_ext);
   ptr = deserialize(ptr, tn.mn_ext);
   ptr = deserialize(ptr, tn.mn_count);
   ptr = deserialize(ptr, tn.mn_terminal);
   ptr = deserialize(ptr, tn.direction);
   return ptr;
}
*/

//////////////////////////////////////////////////////////////

#endif //CHAR_SERIALIZE

#endif
