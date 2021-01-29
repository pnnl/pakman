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
#include <iostream>
#include <sstream>
#include <random>
#include <utility>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include "distribute_kmers.h"

//extern int rank, size;

template<typename POD>
std::ostream& serialize(std::ostream& os, const POD& v);


template<typename POD1, typename POD2>
std::ostream& serialize(std::ostream& os, const std::pair<POD1, POD2>& p);

template<typename POD>
std::ostream& serialize(std::ostream& os, std::vector<POD> const& vec);

std::ostream& serialize(std::ostream& os, std::vector<std::pair<int, int>> const& pr);

std::ostream& serialize(std::ostream& os, std::vector<bool> const& tvec);

std::ostream& serialize(std::ostream& os, const BasePairVector &bp);

std::ostream& serialize(std::ostream& os, std::vector<BasePairVector> const& bvec);

std::ostream& serialize(std::ostream& os, const WireInfo &winfo);
std::ostream& serialize(std::ostream& os, std::vector<WireInfo> const& wi);

std::ostream& serialize(std::ostream& os, const PrefixInfo &pi);
std::ostream& serialize(std::ostream& os, std::vector<PrefixInfo> const& vp);

std::ostream& serialize(std::ostream& os, const TransferNode& tn);

std::ostream& serialize(std::ostream& os, const MacroNode& mn);




template<typename POD>
std::istream& deserialize(std::istream& is, POD& v);

template<typename POD1, typename POD2>
std::istream& deserialize(std::istream& is, std::pair<POD1, POD2>& p);

template<typename POD>
std::istream& deserialize(std::istream& is, std::vector<POD>& vec);

std::istream& deserialize(std::istream& is, BasePairVector &bp);
std::istream& deserialize(std::istream& is, std::vector<BasePairVector>& bvec);

std::istream& deserialize(std::istream& is, std::vector<std::pair<int, int>>& pr);

std::istream& deserialize(std::istream& is, std::vector<bool>& tvec);

std::istream& deserialize(std::istream& is, WireInfo &winfo);
std::istream& deserialize(std::istream& is, std::vector<WireInfo>& wi);

std::istream& deserialize(std::istream& is, PrefixInfo &pi);
std::istream& deserialize(std::istream& is, std::vector<PrefixInfo>& vp);

std::istream& deserialize(std::istream& is, TransferNode& tn);

std::istream& deserialize(std::istream& is, MacroNode& mn);

#ifdef CHAR_SERIALIZE
template<typename T>
size_t serialize_size(const T &v);

template<typename T>
size_t serialize_size(const std::vector<T>& vec);

template<typename T1,
         typename T2>
size_t serialize_size(const std::pair<T1, T2>& p);

size_t serialize_size(const BasePairVector &bp);

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

size_t serialize_size_Info(const MnodeInfo& mn);

size_t serialize_size_Tuple(const MNnodeTuple& mn);

size_t serialize_size_Node(const ModNodeInfo& mn);

size_t serialize_size_TransferNode(const TransferNode& tn);

//////////////////////////////////////////////////////////////


template<typename T>
char* serialize(char* ptr, const T& val);

template<typename T1,
          typename T2>
char* serialize(char *ptr, const std::pair<T1, T2>& p);

template<typename T>
char* serialize(char* ptr, const std::vector<T>& vec);

char* serialize(char* ptr, const BasePairVector &bp);

char* serialize_Info(char* ptr,
    const MnodeInfo& mnode);

char* serialize_Tuple(char* ptr,
    const MNnodeTuple& mtup);

char* serialize_Node(char* ptr,
    const ModNodeInfo& minfo);

char* serialize_TransferNode(char* ptr, const TransferNode& tn);

//////////////////////////////////////////////////////////////

template<typename T>
char* deserialize(char* ptr, T& val);

template<typename T1,
          typename T2>
char*
deser(char *ptr, std::pair<T1, T2>& p);

template<typename T>
char* deserialize(char* ptr, std::vector<T>& vec);

char* deserialize(char* ptr, BasePairVector &bp);

char*
deserialize_Info(char* ptr,
      MnodeInfo& mnode);

char*
deserialize_Tuple(char* ptr,
      MNnodeTuple& mtup);

char*
deserialize_Node(char* ptr,
      ModNodeInfo& minfo);

char* deserialize_TransferNode(char* ptr,
      TransferNode& tn);

//////////////////////////////////////////////////////////////
#endif //CHAR_SERIALIZE

#endif
