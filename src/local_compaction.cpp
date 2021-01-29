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

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <atomic>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <parallel/algorithm>
#include <numeric>
#include <unordered_map>
#include "distribute_kmers.h"
//#include "bigmpi.h"
#include <iostream>     // std::cout
#include <sstream>      // std::istringstream
#include <string>
#include <utility> 
#include <fstream>
#include <omp.h>
#include "serialize.h"

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

#define FREQ_CUTOFF 1
#define BEGIN_CUTOFF 1

extern long int this_contig_id;
extern int rank, size;
extern int coverage;
extern int num_threads;

/*inline kmer_t pkmer_shift(kmer_t kmer_in, 
                        ElType el) {

  assert(el>=A && el<=G);
  
  return (kmer_t)((kmer_in>>2) | (kmer_t)((kmer_t)el<<(2*(KMER_LENGTH-1)))) & (kmer_t)KMER_MASK;
  //return ((kmer_in<<2) | (kmer_t)el) & KMER_MASK;
}*/

char num_to_char (int i) {

     char c;
     assert (i>=0 && i<=3);

     switch (i) {
   
     case 0:
           c = 'A';
           break;

     case 1:
           c = 'C';
           break;
          
     case 2:
           c = 'T';
           break;
     
     case 3:
           c = 'G';
           break;
      
    default: 
           c = 'X';
           break;
  }
  return c;
}


/*
void check(MacroNode& mn) {
  //wiring info properties:
  assert(mn.wiring_info.size() == mn.prefixes.size());
  std::vector<size_t> suffix_sz(mn.suffixes.size(), 0);
  for(size_t i=0; i<mn.wiring_info.size(); i++) {
    size_t sz = 0;
    for(size_t j=0; j<mn.wiring_info[i].size(); j++) {
      assert(mn.wiring_info[i][j].count > 0);
      sz += mn.wiring_info[i][j].count;
      assert(mn.wiring_info[i][j].suffix_id >= 0);
      assert(mn.wiring_info[i][j].suffix_id < mn.suffixes.size());
      suffix_sz[mn.wiring_info[i][j].suffix_id] += mn.wiring_info[i][j].count;
    }
    if (sz != mn.prefix_count[i].second) {
        BasePairVector key = mn.k_1_mer; 
        std::string test_data;
        for (size_t l=0; l<key.size(); l++)
             test_data += el_to_char(key[l]);
        printf("i:%lu, key: %s, sz: %lu, prefix_count: %d\n", i, test_data.c_str(), sz, mn.prefix_count[i].second);
        printf("num prefixes: %lu, num_suffixes: %lu\n", mn.prefixes.size(), mn.suffixes.size());
        assert(mn.prefixes.size() == mn.prefix_count.size());
        for (size_t l=0; l<mn.prefixes.size(); l++) {
             BasePairVector x = mn.prefixes[l];
             printf("x size: %lu, count: (%d,%d)\n", x.size(), mn.prefix_count[l].first, mn.prefix_count[l].second);
             for (size_t t=0; t<x.size(); t++)
                  printf("%c", el_to_char(x[t]));
             printf("\n");
        }
        for(size_t p=0; p<mn.wiring_info[i].size(); p++) {
            printf("wire_count: %d, wire_sid: %d, wire_off: %d\n", mn.wiring_info[i][p].count, 
                    mn.wiring_info[i][p].suffix_id, mn.wiring_info[i][p].offset_in_suffix);
        }
    }
    assert(sz == mn.prefix_count[i].second);
  }
  for(size_t i=0; i<mn.suffixes.size(); i++) {
    assert(mn.suffix_count[i].second == suffix_sz[i]);
  }
}
*/


/*
Macronode *find_kmer_exists(kmer_t key, std::vector<MacroNode> &MN_map, int *indx)
{
     std::vector<MacroNode>::iterator it;
    
     it = std::lower_bound(MN_map.begin(), MN_map.end(), key,
                      [] (const MacroNode& lhs, kmer_t rhs) {
                             return (lhs.macronode_seq < rhs);
                      });

     int pos = it - MN_map.begin();
     if (MN_map[pos].macronode_seq == key) {

         *indx = pos;
         return &(MN_map[pos]);
 
     } else
       return NULL;

}
*/

std::vector<int> flatten_k(const std::vector<std::vector<int>> &orig)
{
        std::vector<int> ret;
            for(const auto &v: orig)
                        ret.insert(ret.end(), v.begin(), v.end());
                return ret;
}

std::vector<std::pair<kmer_t,MacroNode>> flatten_mn(const std::vector<std::vector<std::pair<kmer_t,MacroNode>>> &orig)
{
        std::vector<std::pair<kmer_t,MacroNode>> ret;
            for(const auto &v: orig)
                        ret.insert(ret.end(), v.begin(), v.end());
                return ret;
}

void gather_all_macro_nodes (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
       std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map)
{

    double c1 = MPI_Wtime();
    // serialize 

    int num_threads=0;
#pragma omp parallel
{
    num_threads=omp_get_num_threads();
}
    if (rank==0)
        printf ("Number of threads: %d \n", num_threads);
 
    std::string send_data;
    int num_nodes=0;
    uint64_t pad_width=1000;

    // global buffers
    std::string os_thd;
    std::ostringstream os(std::stringstream::binary);
    //cereal::BinaryOutputArchive oarchive(os);
    uint64_t send_data_size=0;
    int scounts_dd=0;

    {
      for(size_t i=0; i<MN_map.size(); i++)
      {
          MacroNode &mn = MN_map[i].second;
          //oarchive(mn);
          serialize(os, mn);
          /*
          std::ostringstream tmp_os(std::stringstream::binary);
          cereal::BinaryOutputArchive tmp_archive(tmp_os);
          tmp_archive(mn);
          std::string data = tmp_os.str();
          per_node_size.push_back(data.size());
          */
      }
      os_thd = os.str();
      if (os_thd.size() > 0)
      {
          size_t nonPaddedSize = os_thd.size();
          size_t padded_size = ((nonPaddedSize/pad_width) + (nonPaddedSize%pad_width!=0))*pad_width;

          os_thd.resize(padded_size, 0);
          /*
          size_t pad_size = ceil((nonPaddedSize/pad_width)+0.5)*pad_width;
          const size_t numPaddingElements = (pad_size - nonPaddedSize % pad_size) % pad_size;
          if(numPaddingElements > 0)
             os_thd.resize(nonPaddedSize + numPaddingElements, 0);

          assert(pad_size == os_thd.size());
          */
      }

      send_data_size = os_thd.length();
      scounts_dd = send_data_size/pad_width;

      send_data.append(os_thd);

      num_nodes = MN_map.size();

    }
    MN_map.clear(); //deleting the local copy of macro_nodes
    MN_map.shrink_to_fit(); //deleting the local copy of macro_nodes

    double c2 = MPI_Wtime();
    double ser_time=0.0, global_ser_time=0.0;
    ser_time = c2 - c1;
    
    MPI_Reduce(&ser_time, &global_ser_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for serializing across all procs (secs): %f \n", 
                            (double)global_ser_time/(double)size);


    double c3 = MPI_Wtime();

    //buffers for Allgather 1 operation
    std::vector<uint64_t> rcounts (size,0);
    std::vector<uint64_t> rdisp (size,0);
    std::vector<int> rcounts_dd(size,0);
    std::vector<int> rdisp_dd(size,0);
    std::string recv_buffer;
    char* recv_ptr=nullptr;

    //buffers for Allgather 2 operation
    std::vector<int> rcounts_node (size,0);
    //std::vector<int> rdisp_node (size,0);
    //std::vector<int> recv_data_buffer;

    //Allgather 1 operation
    MPI_Allgather(&send_data_size, 1, MPI_UINT64_T, 
            rcounts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);

    uint64_t rsize=0;
    for (int t=0; t<size; t++) rsize += rcounts[t];

    for (int t=0; t<size; t++) {
         rdisp[t] = (t>0) ? (rcounts[t-1] + rdisp[t-1]) : 0;
         rcounts_dd[t] = rcounts[t]/pad_width;
         rdisp_dd[t] = (t>0) ? (rcounts_dd[t-1] + rdisp_dd[t-1]) : 0;
    }

    recv_buffer.resize(rsize, 'F');
    recv_ptr=&recv_buffer[0];
    MPI_Datatype rowtype;

    //create contiguous derived data type
    MPI_Type_contiguous(pad_width, MPI_BYTE, &rowtype);
    MPI_Type_commit(&rowtype);

    //Allgather 2 operation
    MPI_Allgather(&num_nodes, 1, MPI_INT, 
            rcounts_node.data(), 1, MPI_INT, MPI_COMM_WORLD);

    uint64_t rsize2=0;
    for (int t=0; t<size; t++) rsize2 += rcounts_node[t];

    //for (int i=0; i<size; i++)
    //     rdisp_node[i] = (i>0) ? (rcounts_node[i-1] + rdisp_node[i-1]) : 0;

    //recv_data_buffer.resize(rsize2);

#ifdef CHECK_ALLGATHERV
     uint64_t global_sfound_cnt=0, sfound_cnt=0;
    
     std::size_t sfound = send_data.find("F");
     if (sfound!=std::string::npos)
         sfound_cnt = std::count(send_data.begin(), send_data.end(), 'F');

     MPI_Allreduce(&sfound_cnt, &global_sfound_cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
         //fprintf (stderr,"rank: %d, Error in alltoallv recv_buffer, F was found at: %lu\n", rank, found);
#endif 
    
    //Allgatherv 1 operation
    int result1 = MPI_Allgatherv(send_data.c_str(), scounts_dd, rowtype, 
            recv_ptr, rcounts_dd.data(), rdisp_dd.data(), rowtype, MPI_COMM_WORLD);

    if (result1 != MPI_SUCCESS) {
        printf("rank: %d, MPI_Allgatherv for string failed with return value: %d\n", rank, result1);
        MPI_Finalize();
        exit(2);
    }

#ifdef CHECK_ALLGATHERV
     uint64_t global_rfound_cnt=0, rfound_cnt=0;

     std::size_t rfound = recv_buffer.find("F");
     if (rfound!=std::string::npos)
         rfound_cnt = std::count(recv_buffer.begin(), recv_buffer.end(), 'F');
         //fprintf (stderr,"rank: %d, Error in alltoallv recv_buffer, F was found at: %lu\n", rank, found);
     
     /*
     MPI_Allreduce(&rfound_cnt, &global_rfound_cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

     if (global_rfound_cnt != global_sfound_cnt)
         fprintf (stderr,"rank: %d, Error in alltoallv, F counts dont match! global_sfound_cnt: %lu, global_rfound_cnt: %lu\n",
                  rank, global_sfound_cnt, global_rfound_cnt);
     
     assert(global_rfound_cnt == global_sfound_cnt);
     */
     assert(rfound_cnt == global_sfound_cnt);

#endif

    os_thd.clear();
    //rcounts.clear();
    //rdisp.clear();
    send_data.clear();
    rcounts_dd.clear();
    rdisp_dd.clear();

     // free datatype
     MPI_Type_free(&rowtype);

    //Allgatherv 2 operation
    /*
    int result2 = MPI_Allgatherv(send_data_buffer.data(), num_nodes, MPI_INT, 
            recv_data_buffer.data(), rcounts_node.data(), rdisp_node.data(), MPI_INT, MPI_COMM_WORLD);

    if (result2 != MPI_SUCCESS) {
        printf("rank: %d, MPI_Allgatherv for node size failed with return value: %d\n", rank, result2);
        MPI_Finalize();
        exit(2);
    }

    rcounts_node.clear();
    rdisp_node.clear();
    send_data_buffer.clear();
    */

    double c4 = MPI_Wtime();
    double allgather_time=0.0, global_allgather_time=0.0;
    allgather_time = c4 - c3;
    
    MPI_Reduce(&allgather_time, &global_allgather_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for Allgather across all procs (secs): %f \n", 
                            (double)global_allgather_time/(double)size);

    MPI_Barrier(MPI_COMM_WORLD);

    //deserialize

    double c5 = MPI_Wtime();
    {

      std::istringstream is(std::ios::binary | std::ios::out | std::ios::in);
      int num=0;

      for (int i=0; i<size; i++)
      {
          is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()+rdisp[i]), rcounts[i]);
          //cereal::BinaryInputArchive iarchive(is); // Create an input archive
          int nnodes = rcounts_node[i];

          //try {
              while(nnodes)
              {
               deserialize(is, global_MN_map[num].second);
               //iarchive(global_MN_map[num].second);
               BasePairVector key = global_MN_map[num].second.k_1_mer;
               kmer_t kmer_1=0;
               for (size_t i=0; i<key.size(); i++)
                   kmer_1 = mnmer_shift(kmer_1, key[i]);

               global_MN_map[num].first = kmer_1;
               nnodes--;
               num++;
              }
          //}
          //catch (cereal::Exception& e) {
          //      std::cout << e.what() << std::endl;
          //}

      }

      assert (rsize2 == num);
      assert (rsize2 == global_MN_map.size());

    }
    
    recv_buffer.clear();
    rcounts.clear();
    rdisp.clear();
    rcounts_node.clear();
    //recv_data_buffer.clear();

    MPI_Barrier(MPI_COMM_WORLD);

    double c6 = MPI_Wtime();
    double deser_time=0.0, global_deser_time=0.0;
    deser_time = c6 - c5;
    
    MPI_Reduce(&deser_time, &global_deser_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for de-serializing across all procs (secs): %f \n", 
                            (double)global_deser_time/(double)size);

    double c7 = MPI_Wtime();
    //sort the vector with all the macro nodes
    __gnu_parallel::sort(global_MN_map.begin(), global_MN_map.end(), Comp_pair);

    double c8 = MPI_Wtime();
    double msort_time=0.0, global_msort_time=0.0;
    msort_time = c8 - c7;
    
    MPI_Reduce(&msort_time, &global_msort_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for sorting entire set of macro nodes across all procs (secs): %f \n", 
                            (double)global_msort_time/(double)size);

}



// void gather_all_macro_nodes (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
//        std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map)
// {

//     double c1 = MPI_Wtime();
//     // serialize 

//     int num_threads=0;
// #pragma omp parallel
// {
//     num_threads=omp_get_num_threads();
// }
//     if (rank==0)
//         printf ("Number of threads: %d \n", num_threads);
 
//     std::string send_data;
//     int num_nodes=0;
//     uint64_t pad_width=1000;

//     // global buffers
//     std::string os_thd;
//     std::ostringstream os(std::stringstream::binary);
//     cereal::BinaryOutputArchive oarchive(os);
//     uint64_t send_data_size=0;
//     int scounts_dd=0;

//     {
//       for(size_t i=0; i<MN_map.size(); i++)
//       {
//           MacroNode &mn = MN_map[i].second;
//           oarchive(mn);
//           /*
//           std::ostringstream tmp_os(std::stringstream::binary);
//           cereal::BinaryOutputArchive tmp_archive(tmp_os);
//           tmp_archive(mn);
//           std::string data = tmp_os.str();
//           per_node_size.push_back(data.size());
//           */
//       }
//       os_thd = os.str();
//       if (os_thd.size() > 0)
//       {
//           size_t nonPaddedSize = os_thd.size();
//           size_t pad_size = ceil((nonPaddedSize/pad_width)+0.5)*pad_width;
//           const size_t numPaddingElements = (pad_size - nonPaddedSize % pad_size) % pad_size;
//           if(numPaddingElements > 0)
//              os_thd.resize(nonPaddedSize + numPaddingElements, 0);

//           assert(pad_size == os_thd.size());
//       }

//       send_data_size = os_thd.length();
//       scounts_dd = send_data_size/pad_width;

//       send_data.append(os_thd);

//       num_nodes = MN_map.size();

//     }
//     MN_map.clear(); //deleting the local copy of macro_nodes
//     MN_map.shrink_to_fit(); //deleting the local copy of macro_nodes

//     double c2 = MPI_Wtime();
//     double ser_time=0.0, global_ser_time=0.0;
//     ser_time = c2 - c1;
    
//     MPI_Reduce(&ser_time, &global_ser_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//     if (rank == 0) printf ("Average time for serializing across all procs (secs): %f \n", 
//                             (double)global_ser_time/(double)size);


//     double c3 = MPI_Wtime();

//     //buffers for Allgather 1 operation
//     std::vector<uint64_t> rcounts (size,0);
//     std::vector<uint64_t> rdisp (size,0);
//     std::vector<int> rcounts_dd(size,0);
//     std::vector<int> rdisp_dd(size,0);
//     std::string recv_buffer;
//     char* recv_ptr=nullptr;

//     //buffers for Allgather 2 operation
//     std::vector<int> rcounts_node (size,0);
//     //std::vector<int> rdisp_node (size,0);
//     //std::vector<int> recv_data_buffer;

//     //Allgather 1 operation
//     MPI_Allgather(&send_data_size, 1, MPI_UINT64_T, 
//             rcounts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);

//     uint64_t rsize=0;
//     for (int t=0; t<size; t++) rsize += rcounts[t];

//     for (int t=0; t<size; t++) {
//          rdisp[t] = (t>0) ? (rcounts[t-1] + rdisp[t-1]) : 0;
//          rcounts_dd[t] = rcounts[t]/pad_width;
//          rdisp_dd[t] = (t>0) ? (rcounts_dd[t-1] + rdisp_dd[t-1]) : 0;
//     }

//     recv_buffer.resize(rsize, 'F');
//     recv_ptr=&recv_buffer[0];
//     MPI_Datatype rowtype;

//     //create contiguous derived data type
//     MPI_Type_contiguous(pad_width, MPI_BYTE, &rowtype);
//     MPI_Type_commit(&rowtype);

//     //Allgather 2 operation
//     MPI_Allgather(&num_nodes, 1, MPI_INT, 
//             rcounts_node.data(), 1, MPI_INT, MPI_COMM_WORLD);

//     uint64_t rsize2=0;
//     for (int t=0; t<size; t++) rsize2 += rcounts_node[t];

//     //for (int i=0; i<size; i++)
//     //     rdisp_node[i] = (i>0) ? (rcounts_node[i-1] + rdisp_node[i-1]) : 0;

//     //recv_data_buffer.resize(rsize2);

// #ifdef CHECK_ALLGATHERV
//      uint64_t global_sfound_cnt=0, sfound_cnt=0;
    
//      std::size_t sfound = send_data.find("F");
//      if (sfound!=std::string::npos)
//          sfound_cnt = std::count(send_data.begin(), send_data.end(), 'F');

//      MPI_Allreduce(&sfound_cnt, &global_sfound_cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
//          //fprintf (stderr,"rank: %d, Error in alltoallv recv_buffer, F was found at: %lu\n", rank, found);
// #endif 
    
//     //Allgatherv 1 operation
//     int result1 = MPI_Allgatherv(send_data.c_str(), scounts_dd, rowtype, 
//             recv_ptr, rcounts_dd.data(), rdisp_dd.data(), rowtype, MPI_COMM_WORLD);

//     if (result1 != MPI_SUCCESS) {
//         printf("rank: %d, MPI_Allgatherv for string failed with return value: %d\n", rank, result1);
//         MPI_Finalize();
//         exit(2);
//     }

// #ifdef CHECK_ALLGATHERV
//      uint64_t global_rfound_cnt=0, rfound_cnt=0;

//      std::size_t rfound = recv_buffer.find("F");
//      if (rfound!=std::string::npos)
//          rfound_cnt = std::count(recv_buffer.begin(), recv_buffer.end(), 'F');
//          //fprintf (stderr,"rank: %d, Error in alltoallv recv_buffer, F was found at: %lu\n", rank, found);
     
//      /*
//      MPI_Allreduce(&rfound_cnt, &global_rfound_cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

//      if (global_rfound_cnt != global_sfound_cnt)
//          fprintf (stderr,"rank: %d, Error in alltoallv, F counts dont match! global_sfound_cnt: %lu, global_rfound_cnt: %lu\n",
//                   rank, global_sfound_cnt, global_rfound_cnt);
     
//      assert(global_rfound_cnt == global_sfound_cnt);
//      */
//      assert(rfound_cnt == global_sfound_cnt);

// #endif

//     os_thd.clear();
//     //rcounts.clear();
//     //rdisp.clear();
//     send_data.clear();
//     rcounts_dd.clear();
//     rdisp_dd.clear();

//      // free datatype
//      MPI_Type_free(&rowtype);

//     //Allgatherv 2 operation
//     /*
//     int result2 = MPI_Allgatherv(send_data_buffer.data(), num_nodes, MPI_INT, 
//             recv_data_buffer.data(), rcounts_node.data(), rdisp_node.data(), MPI_INT, MPI_COMM_WORLD);

//     if (result2 != MPI_SUCCESS) {
//         printf("rank: %d, MPI_Allgatherv for node size failed with return value: %d\n", rank, result2);
//         MPI_Finalize();
//         exit(2);
//     }

//     rcounts_node.clear();
//     rdisp_node.clear();
//     send_data_buffer.clear();
//     */

//     double c4 = MPI_Wtime();
//     double allgather_time=0.0, global_allgather_time=0.0;
//     allgather_time = c4 - c3;
    
//     MPI_Reduce(&allgather_time, &global_allgather_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//     if (rank == 0) printf ("Average time for Allgather across all procs (secs): %f \n", 
//                             (double)global_allgather_time/(double)size);

//     MPI_Barrier(MPI_COMM_WORLD);

//     //deserialize

//     double c5 = MPI_Wtime();
//     {

//       std::istringstream is(std::ios::binary | std::ios::out | std::ios::in);
//       int num=0;

//       for (int i=0; i<size; i++)
//       {
//           is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()+rdisp[i]), rcounts[i]);
//           cereal::BinaryInputArchive iarchive(is); // Create an input archive
//           int nnodes = rcounts_node[i];

//           try {
//               while(nnodes)
//               {
//                iarchive(global_MN_map[num].second);
//                BasePairVector key = global_MN_map[num].second.k_1_mer;
//                kmer_t kmer_1=0;
//                for (size_t i=0; i<key.size(); i++)
//                    kmer_1 = mnmer_shift(kmer_1, key[i]);

//                global_MN_map[num].first = kmer_1;
//                nnodes--;
//                num++;
//               }
//           }
//           catch (cereal::Exception& e) {
//                 std::cout << e.what() << std::endl;
//           }

//       }

//       assert (rsize2 == num);
//       assert (rsize2 == global_MN_map.size());

//     }
    
//     recv_buffer.clear();
//     rcounts.clear();
//     rdisp.clear();
//     rcounts_node.clear();
//     //recv_data_buffer.clear();

//     MPI_Barrier(MPI_COMM_WORLD);

//     double c6 = MPI_Wtime();
//     double deser_time=0.0, global_deser_time=0.0;
//     deser_time = c6 - c5;
    
//     MPI_Reduce(&deser_time, &global_deser_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//     if (rank == 0) printf ("Average time for de-serializing across all procs (secs): %f \n", 
//                             (double)global_deser_time/(double)size);

//     double c7 = MPI_Wtime();
//     //sort the vector with all the macro nodes
//     __gnu_parallel::sort(global_MN_map.begin(), global_MN_map.end(), Comp_pair);

//     double c8 = MPI_Wtime();
//     double msort_time=0.0, global_msort_time=0.0;
//     msort_time = c8 - c7;
    
//     MPI_Reduce(&msort_time, &global_msort_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//     if (rank == 0) printf ("Average time for sorting entire set of macro nodes across all procs (secs): %f \n", 
//                             (double)global_msort_time/(double)size);

// }


#ifdef MULTITHREADED
void gather_all_macro_nodes_nodd (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
       std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map)
{

    double c1 = MPI_Wtime();
    // serialize 

    int num_threads=0;
#pragma omp parallel
{
    num_threads=omp_get_num_threads();
}
    if (rank==0)
        printf ("Number of threads: %d \n", num_threads);
 
    std::vector<int> send_data_buffer(MN_map.size());
    std::string send_data;
    int num_nodes=0;
    int send_node_size_all=0;

    {
      // global buffers
      std::vector<std::string> os_thd (num_threads);
      std::vector< std::vector<int> > node_size;
      for (int i=0; i<num_threads; i++)
           node_size.push_back(std::vector<int> ());
      
#pragma omp parallel
{
      // local buffers 
      std::ostringstream os(std::stringstream::binary);
      cereal::BinaryOutputArchive oarchive(os);
      std::vector<int> per_node_size;

#pragma omp for schedule(guided)
      for(size_t i=0; i<MN_map.size(); i++)
      {
          //std::pair<kmer_t, MacroNode> ref = MN_map[i];
          MacroNode &mn = MN_map[i].second;
          oarchive(mn);

          std::ostringstream tmp_os(std::stringstream::binary);
          cereal::BinaryOutputArchive tmp_archive(tmp_os);
          tmp_archive(mn);
          std::string data = tmp_os.str();
          per_node_size.push_back(data.size());

      }
      os_thd[omp_get_thread_num()] = os.str();
      node_size[omp_get_thread_num()] = per_node_size;
      per_node_size.clear();

}
#pragma omp parallel for reduction(+:num_nodes)
      for (int i=0; i<num_threads; i++)
           num_nodes += node_size[i].size();

      for (int i=0; i<num_threads; i++)
           send_data.append(os_thd[i]);

      for (int i=0; i<num_threads; i++)
           for (int j=0; j<node_size[i].size(); j++)
                send_node_size_all += node_size[i][j];

      //printf("rank: %d, size of send stream buffer: %lu, num_nodes: %d\n", rank, send_data.size(), num_nodes);
      //send_data_buffer = flatten_k(node_size);
      //send_data_buffer.resize(num_nodes);
      assert(num_nodes == MN_map.size());

      size_t off=0;
      for (int i=0; i<num_threads; i++) {
           std::copy(node_size[i].begin(), node_size[i].end(), send_data_buffer.begin()+off);
           off += node_size[i].size();
      }
      node_size.clear();

    }
    MN_map.clear(); //deleting the local copy of macro_nodes
    MN_map.shrink_to_fit(); //deleting the local copy of macro_nodes

    double c2 = MPI_Wtime();
    double ser_time=0.0, global_ser_time=0.0;
    ser_time = c2 - c1;
    
    MPI_Reduce(&ser_time, &global_ser_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for serializing across all procs (secs): %f \n", 
                            (double)global_ser_time/(double)size);


    double c3 = MPI_Wtime();

    int send_data_size = send_data.length();
    assert(send_node_size_all == send_data_size);

    //buffers for Allgather 1 operation
    std::vector<int> rcounts (size,0);
    std::vector<int> rdisp (size,0);
    std::string recv_buffer;
    char* recv_ptr=nullptr;

    //buffers for Allgather 2 operation
    std::vector<int> rcounts_node (size,0);
    std::vector<int> rdisp_node (size,0);
    std::vector<int> recv_data_buffer;

    //Allgather 1 operation
    MPI_Allgather(&send_data_size, 1, MPI_INT, 
            rcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    uint64_t rsize=0;
    for (int t=0; t<size; t++) rsize += rcounts[t];

    for (int i=0; i<size; i++)
         rdisp[i] = (i>0) ? (rcounts[i-1] + rdisp[i-1]) : 0;

    recv_buffer.resize(rsize, 'F');
    recv_ptr=&recv_buffer[0];

    //Allgather 2 operation
    MPI_Allgather(&num_nodes, 1, MPI_INT, 
            rcounts_node.data(), 1, MPI_INT, MPI_COMM_WORLD);

    uint64_t rsize2=0;
    for (int t=0; t<size; t++) rsize2 += rcounts_node[t];

    for (int i=0; i<size; i++)
         rdisp_node[i] = (i>0) ? (rcounts_node[i-1] + rdisp_node[i-1]) : 0;

    recv_data_buffer.resize(rsize2);

#ifdef CHECK_ALLGATHERV
     uint64_t global_sfound_cnt=0, sfound_cnt=0;
    
     std::size_t sfound = send_data.find("F");
     if (sfound!=std::string::npos)
         sfound_cnt = std::count(send_data.begin(), send_data.end(), 'F');

     MPI_Allreduce(&sfound_cnt, &global_sfound_cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
         //fprintf (stderr,"rank: %d, Error in alltoallv recv_buffer, F was found at: %lu\n", rank, found);
#endif 
    
    //Allgatherv 1 operation
    int result1 = MPI_Allgatherv(send_data.c_str(), send_data_size, MPI_BYTE, 
            recv_ptr, rcounts.data(), rdisp.data(), MPI_BYTE, MPI_COMM_WORLD);

    if (result1 != MPI_SUCCESS) {
        printf("rank: %d, MPI_Allgatherv for string failed with return value: %d\n", rank, result1);
        MPI_Finalize();
        exit(2);
    }

#ifdef CHECK_ALLGATHERV
     uint64_t global_rfound_cnt=0, rfound_cnt=0;

     std::size_t rfound = recv_buffer.find("F");
     if (rfound!=std::string::npos)
         rfound_cnt = std::count(recv_buffer.begin(), recv_buffer.end(), 'F');
         //fprintf (stderr,"rank: %d, Error in alltoallv recv_buffer, F was found at: %lu\n", rank, found);
     
     /*
     MPI_Allreduce(&rfound_cnt, &global_rfound_cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

     if (global_rfound_cnt != global_sfound_cnt)
         fprintf (stderr,"rank: %d, Error in alltoallv, F counts dont match! global_sfound_cnt: %lu, global_rfound_cnt: %lu\n",
                  rank, global_sfound_cnt, global_rfound_cnt);
     
     assert(global_rfound_cnt == global_sfound_cnt);
     */
     assert(rfound_cnt == global_sfound_cnt);

#endif

    rcounts.clear();
    rdisp.clear();
    send_data.clear();

    //Allgatherv 2 operation
    int result2 = MPI_Allgatherv(send_data_buffer.data(), num_nodes, MPI_INT, 
            recv_data_buffer.data(), rcounts_node.data(), rdisp_node.data(), MPI_INT, MPI_COMM_WORLD);

    if (result2 != MPI_SUCCESS) {
        printf("rank: %d, MPI_Allgatherv for node size failed with return value: %d\n", rank, result2);
        MPI_Finalize();
        exit(2);
    }

    rcounts_node.clear();
    rdisp_node.clear();
    send_data_buffer.clear();

    double c4 = MPI_Wtime();
    double allgather_time=0.0, global_allgather_time=0.0;
    allgather_time = c4 - c3;
    
    MPI_Reduce(&allgather_time, &global_allgather_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for Allgather across all procs (secs): %f \n", 
                            (double)global_allgather_time/(double)size);

    MPI_Barrier(MPI_COMM_WORLD);

    //deserialize

    double c5 = MPI_Wtime();
    {
      std::vector< std::vector< std::pair<kmer_t,MacroNode> > > mn_nodes;
      for (int i=0; i<num_threads; i++)
           mn_nodes.push_back(std::vector<std::pair<kmer_t,MacroNode>> ());

      int node_chunk[num_threads];
      size_t node_offset[num_threads];

      //node_chunk[i] corresponds to chunk of nodes read by each thread
      for (int i=0; i<num_threads; i++){
           int chunk_size = rsize2/num_threads;
           node_chunk[i] = (i==(num_threads-1)) ? (chunk_size + (rsize2%num_threads)) : chunk_size;
      }

      //node_offset[i] corresponds to the offset (position in recv_buffer) from which to extract the chunk
      for (int i=0; i<num_threads; i++) {
           node_offset[i] = (i>0) ? (node_chunk[i-1] + node_offset[i-1]) : 0;
      }

#pragma omp parallel shared(node_offset, node_chunk, mn_nodes, recv_data_buffer, recv_buffer)
{
      int chunk_size = node_chunk[omp_get_thread_num()];
      std::vector<std::pair<kmer_t,MacroNode>> local_buf(chunk_size);

      size_t start_offset=0, start_range=0;
      for (int i=0; i<node_offset[omp_get_thread_num()]; i++)
           start_offset += recv_data_buffer[i];

      for (int i=node_offset[omp_get_thread_num()]; i<(node_offset[omp_get_thread_num()]+node_chunk[omp_get_thread_num()]); i++)
           start_range += recv_data_buffer[i];

      std::istringstream is;
      //is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()), 22112);
      //is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()+22112), 14048);
      is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()+start_offset), start_range);
      cereal::BinaryInputArchive iarchive(is); // Create an input archive
      int num=0;

    //try {
      while(chunk_size)
      {
          iarchive(local_buf[num].second);
          BasePairVector key = local_buf[num].second.k_1_mer;
          kmer_t kmer_1=0;
          for (size_t i=0; i<key.size(); i++)
               kmer_1 = mnmer_shift(kmer_1, key[i]);

          local_buf[num].first = kmer_1;
          chunk_size--;
          num++;
      }
   // }
   //   catch (cereal::Exception& e) {
   //       std::cout << e.what() << std::endl;
   //   }

      mn_nodes[omp_get_thread_num()] = local_buf;
      local_buf.clear();
}
      //global_MN_map = flatten_mn(mn_nodes);
      size_t offset=0;
      for (int i=0; i<num_threads; i++) {
           std::copy(mn_nodes[i].begin(), mn_nodes[i].end(), global_MN_map.begin()+offset);
           offset += mn_nodes[i].size();
      }
      mn_nodes.clear();
    }
    
    recv_buffer.clear();
    recv_data_buffer.clear();

    double c6 = MPI_Wtime();
    double deser_time=0.0, global_deser_time=0.0;
    deser_time = c6 - c5;
    
    MPI_Reduce(&deser_time, &global_deser_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for de-serializing across all procs (secs): %f \n", 
                            (double)global_deser_time/(double)size);

    MPI_Barrier(MPI_COMM_WORLD);

    double c7 = MPI_Wtime();
    //sort the vector with all the macro nodes
    __gnu_parallel::sort(global_MN_map.begin(), global_MN_map.end(), Comp_pair);

    double c8 = MPI_Wtime();
    double msort_time=0.0, global_msort_time=0.0;
    msort_time = c8 - c7;
    
    MPI_Reduce(&msort_time, &global_msort_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for sorting entire set of macro nodes across all procs (secs): %f \n", 
                            (double)global_msort_time/(double)size);

}
#endif

#ifdef BIGMPI
void gather_all_macro_nodes_mt (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
       std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map)
{

    double c1 = MPI_Wtime();
    // serialize 

    int num_threads=0;
#pragma omp parallel
{
    num_threads=omp_get_num_threads();
}
    if (rank==0)
        printf ("Number of threads: %d \n", num_threads);
 
    std::vector<int> send_data_buffer(MN_map.size());
    std::string send_data;
    uint64_t num_nodes=0;

    {
      // global buffers
      std::vector<std::string> os_thd (num_threads);
      std::vector< std::vector<int> > node_size;
      for (int i=0; i<num_threads; i++)
           node_size.push_back(std::vector<int> ());
      
#pragma omp parallel
{
      // local buffers 
      std::ostringstream os(std::stringstream::binary);
      cereal::BinaryOutputArchive oarchive(os);
      std::vector<int> per_node_size;

#pragma omp for schedule(guided)
      for(size_t i=0; i<MN_map.size(); i++)
      {
          std::pair<kmer_t, MacroNode> ref = MN_map[i];
          MacroNode mn = ref.second;
          oarchive(mn);

          std::ostringstream tmp_os(std::stringstream::binary);
          cereal::BinaryOutputArchive tmp_archive(tmp_os);
          tmp_archive(mn);
          std::string data = tmp_os.str();
          per_node_size.push_back(data.size());

      }
      os_thd[omp_get_thread_num()] = os.str();
      node_size[omp_get_thread_num()] = per_node_size;
      per_node_size.clear();

}
#pragma omp parallel for reduction(+:num_nodes)
      for (int i=0; i<num_threads; i++)
           num_nodes += node_size[i].size();

      for (int i=0; i<num_threads; i++)
           send_data.append(os_thd[i]);

      //printf("rank: %d, size of send stream buffer: %lu, num_nodes: %d\n", rank, send_data.size(), num_nodes);
      //send_data_buffer = flatten_k(node_size);
      //send_data_buffer.resize(num_nodes);
      assert(num_nodes == MN_map.size());

      size_t off=0;
      for (int i=0; i<num_threads; i++) {
           std::copy(node_size[i].begin(), node_size[i].end(), send_data_buffer.begin()+off);
           off += node_size[i].size();
      }
      node_size.clear();

    }
    MN_map.clear(); //deleting the local copy of macro_nodes
    MN_map.shrink_to_fit(); //deleting the local copy of macro_nodes

    double c2 = MPI_Wtime();
    double ser_time=0.0, global_ser_time=0.0;
    ser_time = c2 - c1;
    
    MPI_Reduce(&ser_time, &global_ser_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for serializing across all procs (secs): %f \n", 
                            (double)global_ser_time/(double)size);


    double c3 = MPI_Wtime();

    uint64_t send_data_size = send_data.length();

    //buffers for Allgather 1 operation
    std::vector<MPI_Count> rcounts (size,0);
    std::vector<MPI_Aint> rdisp (size,0);
    std::string recv_buffer;
    char* recv_ptr=nullptr;

    //buffers for Allgather 2 operation
    std::vector<MPI_Count> rcounts_node (size,0);
    std::vector<MPI_Aint> rdisp_node (size,0);
    std::vector<int> recv_data_buffer;

    //Allgather 1 operation
    MPIX_Allgather_x(&send_data_size, 1, MPI_UINT64_T, 
            rcounts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);

    uint64_t rsize=0;
    for (int t=0; t<size; t++) rsize += rcounts[t];

    for (int i=0; i<size; i++)
         rdisp[i] = (i>0) ? (rcounts[i-1] + rdisp[i-1]) : 0;

    recv_buffer.resize(rsize);
    recv_ptr=&recv_buffer[0];

    //Allgather 2 operation
    MPIX_Allgather_x(&num_nodes, 1, MPI_UINT64_T, 
            rcounts_node.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);

    uint64_t rsize2=0;
    for (int t=0; t<size; t++) rsize2 += rcounts_node[t];

    for (int i=0; i<size; i++)
         rdisp_node[i] = (i>0) ? (rcounts_node[i-1] + rdisp_node[i-1]) : 0;

    recv_data_buffer.resize(rsize2);

    //Allgatherv 1 operation
    int result1 = MPIX_Allgatherv_x(send_data.c_str(), send_data_size, MPI_BYTE, 
            recv_ptr, rcounts.data(), rdisp.data(), MPI_BYTE, MPI_COMM_WORLD);

    if (result1 != MPI_SUCCESS) {
        printf("rank: %d, MPI_Allgatherv for string failed with return value: %d\n", rank, result1);
        MPI_Finalize();
        exit(2);
    }

    rcounts.clear();
    rdisp.clear();
    send_data.clear();

    //Allgatherv 2 operation
    int result2 = MPIX_Allgatherv_x(send_data_buffer.data(), num_nodes, MPI_INT, 
            recv_data_buffer.data(), rcounts_node.data(), rdisp_node.data(), MPI_INT, MPI_COMM_WORLD);

    if (result2 != MPI_SUCCESS) {
        printf("rank: %d, MPI_Allgatherv for node size failed with return value: %d\n", rank, result2);
        MPI_Finalize();
        exit(2);
    }

    rcounts_node.clear();
    rdisp_node.clear();
    send_data_buffer.clear();

    double c4 = MPI_Wtime();
    double allgather_time=0.0, global_allgather_time=0.0;
    allgather_time = c4 - c3;
    
    MPI_Reduce(&allgather_time, &global_allgather_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for Allgather across all procs (secs): %f \n", 
                            (double)global_allgather_time/(double)size);

    MPI_Barrier(MPI_COMM_WORLD);

    //deserialize

    double c5 = MPI_Wtime();
    {
      std::vector< std::vector< std::pair<kmer_t,MacroNode> > > mn_nodes;
      for (int i=0; i<num_threads; i++)
           mn_nodes.push_back(std::vector<std::pair<kmer_t,MacroNode>> ());

      int node_chunk[num_threads];
      size_t node_offset[num_threads];

      //node_chunk[i] corresponds to chunk of nodes read by each thread
      for (int i=0; i<num_threads; i++){
           int chunk_size = rsize2/num_threads;
           node_chunk[i] = (i==(num_threads-1)) ? (chunk_size + (rsize2%num_threads)) : chunk_size;
      }

      //node_offset[i] corresponds to the offset (position in recv_buffer) from which to extract the chunk
      for (int i=0; i<num_threads; i++) {
           node_offset[i] = (i>0) ? (node_chunk[i-1] + node_offset[i-1]) : 0;
      }

#pragma omp parallel shared(node_offset, node_chunk, mn_nodes, recv_data_buffer, recv_buffer)
{
      int chunk_size = node_chunk[omp_get_thread_num()];
      std::vector<std::pair<kmer_t,MacroNode>> local_buf(chunk_size);

      size_t start_offset=0, start_range=0;
      for (int i=0; i<node_offset[omp_get_thread_num()]; i++)
           start_offset += recv_data_buffer[i];

      for (int i=node_offset[omp_get_thread_num()]; i<(node_offset[omp_get_thread_num()]+node_chunk[omp_get_thread_num()]); i++)
           start_range += recv_data_buffer[i];

      std::istringstream is;
      //is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()), 22112);
      //is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()+22112), 14048);
      is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()+start_offset), start_range);
      cereal::BinaryInputArchive iarchive(is); // Create an input archive
      int num=0;

    //try {
      while(chunk_size)
      {
          iarchive(local_buf[num].second);
          BasePairVector key = local_buf[num].second.k_1_mer;
          kmer_t kmer_1=0;
          for (size_t i=0; i<key.size(); i++)
               kmer_1 = mnmer_shift(kmer_1, key[i]);

          local_buf[num].first = kmer_1;
          chunk_size--;
          num++;
      }
   // }
   //   catch (cereal::Exception& e) {
   //       std::cout << e.what() << std::endl;
   //   }

      mn_nodes[omp_get_thread_num()] = local_buf;
      local_buf.clear();
}
      //global_MN_map = flatten_mn(mn_nodes);
      size_t offset=0;
      for (int i=0; i<num_threads; i++) {
           std::copy(mn_nodes[i].begin(), mn_nodes[i].end(), global_MN_map.begin()+offset);
           offset += mn_nodes[i].size();
      }
      mn_nodes.clear();
    }
    
    recv_buffer.clear();
    recv_data_buffer.clear();

    double c6 = MPI_Wtime();
    double deser_time=0.0, global_deser_time=0.0;
    deser_time = c6 - c5;
    
    MPI_Reduce(&deser_time, &global_deser_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for de-serializing across all procs (secs): %f \n", 
                            (double)global_deser_time/(double)size);

    MPI_Barrier(MPI_COMM_WORLD);

    double c7 = MPI_Wtime();
    //sort the vector with all the macro nodes
    __gnu_parallel::sort(global_MN_map.begin(), global_MN_map.end(), Comp_pair);

    double c8 = MPI_Wtime();
    double msort_time=0.0, global_msort_time=0.0;
    msort_time = c8 - c7;
    
    MPI_Reduce(&msort_time, &global_msort_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for sorting entire set of macro nodes across all procs (secs): %f \n", 
                            (double)global_msort_time/(double)size);

}
#endif

