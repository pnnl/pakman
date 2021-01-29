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
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "distribute_kmers.h"
//#include "bigmpi.h"

extern int rank, size;
extern int coverage;
extern int num_threads;

///* Vector containing k-mers allocated to each process */
//extern std::vector<kmer_t> kmer_proc_buf;
//extern std::vector<int> kmer_cnt_proc_buf;

//timers
double mn_transfer_time=0.0, global_mn_transfer_time=0.0;
//double mn_alltoall_time=0.0, mn_global_alltoall_time=0.0;
//double mn_alltoallv_time=0.0, mn_global_alltoallv_time=0.0;

int visit_value (int num)                                               
{
        double ceil_val=(double)num/coverage;
        return (int)ceil(ceil_val);
}


bool pairCompare(const KmerPairs &a, const KmerPairs &b) {
        return a.seq < b.seq;

}

bool Comp_suffix(const KmerPairs &v1, const KmerPairs &v2) {

         kmer_t succ_k1 = (kmer_t)v1.seq & (kmer_t)SUCC_MASK;

         kmer_t succ_k2 = (kmer_t)v2.seq & (kmer_t)SUCC_MASK;

         return succ_k1 < succ_k2;

}


std::vector<KmerPairs> transfer_macro_nodes(std::vector<KmerPairs> &mn_send_buf, 
                                            std::vector<int> &scounts)
                                            //size_t ssize)
{

    std::vector<int> rcounts (size,0);
    std::vector<int> rdisp (size,0);
    std::vector<int> sdisp (size,0);
    uint64_t rsize=0;
    //kmer_t ssize=0, rsize=0;
 
    /*for (int t=0; t<size; t++) {
         scounts[t] = mn_node_list[t].size();
         ssize += scounts[t];
    }*/

     //create contiguous derived data type
     MPI_Datatype rowtype;
     MPI_Type_contiguous(sizeof(KmerPairs), MPI_BYTE, &rowtype);
     MPI_Type_commit(&rowtype);

    double comm1 = MPI_Wtime ();
    MPI_Alltoall(scounts.data(), 1, MPI_INT, rcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    double comm2 = MPI_Wtime ();
    mn_transfer_time += (comm2 - comm1);
    //alltoall_time += (comm2 - comm1);

    /*kmer_t p=0,k=0;
    std::vector<KmerPairs> mn_send_buf(ssize);
    for (int t=0; t<size; t++)
    {
        for(k=0; k<mn_node_list[t].size(); k++) {
            mn_send_buf[p] = mn_node_list[t][k];
            p++;
        }
        assert(k == scounts[t]);
    }
    assert(p == ssize);
    */

    for (int t=0; t<size; t++) rsize += rcounts[t];
    std::vector<KmerPairs> mn_recv_buf(rsize);

    for (int t=0; t<size; t++) {
        sdisp[t] = (t>0) ? (scounts[t-1] + sdisp[t-1]) : 0;
        rdisp[t] = (t>0) ? (rcounts[t-1] + rdisp[t-1]) : 0;
    }

    double comm3 = MPI_Wtime ();
    MPI_Alltoallv(mn_send_buf.data(), scounts.data(), sdisp.data(), rowtype, 
            mn_recv_buf.data(), rcounts.data(), rdisp.data(), rowtype, MPI_COMM_WORLD);
    double comm4 = MPI_Wtime ();
    mn_transfer_time += (comm4 - comm3);
    //alltoallv_time += (comm4 - comm3);

    rcounts.clear();
    scounts.clear();
    rdisp.clear();
    sdisp.clear();
    mn_send_buf.clear();
    mn_send_buf.shrink_to_fit();

     // free datatype
     MPI_Type_free(&rowtype);

    MPI_Reduce(&mn_transfer_time, &global_mn_transfer_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for communication in Transfer of MN nodes (secs): %f \n", 
                            (double)global_mn_transfer_time/(double)size);

    //MPI_Reduce(&alltoallv_time, &global_alltoallv_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //if (rank == 0) printf ("Average time for AlltoallV in Transfer of MN nodes (secs): %f \n", 
    //                        (double)global_alltoallv_time/(double)size);

    return mn_recv_buf;

}

size_t bsearch_pkmer_exists (kmer_t mn_key, std::vector<KmerPairs> &kmer_list)
{
    std::vector<KmerPairs>::iterator it;
    size_t pos = -1;
    
    it = std::lower_bound(kmer_list.begin(), kmer_list.end(), mn_key,
                      [] (const KmerPairs& lhs, kmer_t rhs) {
                          kmer_t pred_k1 = (kmer_t)lhs.seq ^ 1;
                          pred_k1 = mn_shift(pred_k1);
                          return (pred_k1 < rhs);
                      });

    pos = it - kmer_list.begin();
    //kmer_t pred_k1 = (kmer_t)(kmer_list[pos].seq) ^ 1;
    //pred_k1 = mn_shift(pred_k1);
    //if (pred_k1 == key)    
    return pos;

}

size_t bsearch_skmer_exists (kmer_t mn_key, std::vector<KmerPairs> &kmer_list)
{
    std::vector<KmerPairs>::iterator it;
    size_t pos = -1;
    
    it = std::lower_bound(kmer_list.begin(), kmer_list.end(), mn_key,
                      [] (const KmerPairs& lhs, kmer_t rhs) {
                          kmer_t succ_k1 = (kmer_t)lhs.seq & (kmer_t)SUCC_MASK;
                          return (succ_k1 < rhs);
                      });

    pos = it - kmer_list.begin();

    return pos;

}

void construct_macro_nodes (std::vector<std::pair<kmer_t,MacroNode>> &MN_map, 
        std::vector<KmerPairs> &all_kmers_from_procs)
{

    double mn1 = MPI_Wtime ();
   
    std::vector< std::vector<kmer_t> > glist_of_mn_nodes(num_threads);

#pragma omp parallel shared(glist_of_mn_nodes, all_kmers_from_procs)
{
    std::vector<kmer_t> tmp_mn_nodes;
    kmer_t pred_k_1=0, succ_k_1=0;

#pragma omp for
    for (size_t f=0; f<all_kmers_from_procs.size(); f++) 
    {
        /* check if k-mer is candidate for macro_node creation */
        kmer_t full_kmer = all_kmers_from_procs[f].seq;

        pred_k_1 = (kmer_t)full_kmer ^ 1;
        pred_k_1 = mn_shift(pred_k_1);
        if (retrieve_proc_id(pred_k_1) == rank)
            tmp_mn_nodes.push_back(pred_k_1);

        succ_k_1 = (kmer_t)full_kmer & (kmer_t)SUCC_MASK;
        if (retrieve_proc_id(succ_k_1) == rank)
            tmp_mn_nodes.push_back(succ_k_1);
    }

    //sort and unique to obtain the list of distinct macro nodes for each thread
    //sort(tmp_mn_nodes.begin(), tmp_mn_nodes.end());
    //tmp_mn_nodes.erase( unique(tmp_mn_nodes.begin(), tmp_mn_nodes.end() ), tmp_mn_nodes.end() ); 

    glist_of_mn_nodes[omp_get_thread_num()] = tmp_mn_nodes;
    tmp_mn_nodes.clear();

}
    
    std::vector<kmer_t> list_of_mn_nodes;
    for (int i=0; i<num_threads; i++) {
         std::copy(glist_of_mn_nodes[i].begin(), glist_of_mn_nodes[i].end(), std::back_inserter(list_of_mn_nodes));
                 
    }
    
    //sort and unique to obtain the list of distinct macro nodes
    __gnu_parallel::sort(list_of_mn_nodes.begin(), list_of_mn_nodes.end());
    list_of_mn_nodes.erase( unique(list_of_mn_nodes.begin(), list_of_mn_nodes.end() ), list_of_mn_nodes.end() ); 

    MN_map.reserve(list_of_mn_nodes.size());
    {
      struct get_pair {
          std::pair<kmer_t,MacroNode> operator() (const kmer_t &p) {
              return std::make_pair( p, MacroNode() ); }
      };
      std::transform( list_of_mn_nodes.begin(), list_of_mn_nodes.end(), back_inserter(MN_map), get_pair() );
    }

    assert (list_of_mn_nodes.size() == MN_map.size());

    double mn2 = MPI_Wtime ();
    double mn_contruct_time=0.0, global_mn_contruct_time=0.0;

    mn_contruct_time = mn2 - mn1;
    MPI_Reduce(&mn_contruct_time, &global_mn_contruct_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for identifying the macro_nodes (secs): %f \n", 
                            (double)global_mn_contruct_time/(double)size);

    double mn3 = MPI_Wtime ();
    //sort kmer list based on prefixes
    __gnu_parallel::sort(all_kmers_from_procs.begin(), all_kmers_from_procs.end(), pairCompare);

    //sort kmer list based on suffixes
    std::vector<KmerPairs> all_kmers_suff = all_kmers_from_procs;
    __gnu_parallel::sort(all_kmers_suff.begin(), all_kmers_suff.end(), Comp_suffix);

    double mn4 = MPI_Wtime ();
    double sort_time=0.0, global_sort_time=0.0;

    sort_time = mn4 - mn3;
    MPI_Reduce(&sort_time, &global_sort_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for sorting the pref/suff kmer buffers (secs): %f \n", 
                            (double)global_sort_time/(double)size);

#ifdef DEBUG_SORT
  //FILE *f = fopen("sort_by_suffixes.txt", "w");
  std::string outfile1="sort_by_pref_"+std::to_string(rank)+".log";
  FILE *f = fopen(outfile1.c_str(), "w"); 
  //FILE *f = fopen("sort_by_prefixes.txt", "w"); 

  char kmer_str[KMER_LENGTH+1];

  for (size_t i=0; i<all_kmers_from_procs.size(); i++) {

      for (int k=0; k<KMER_LENGTH; k++) {
           kmer_str[k] = el_to_char(kmerel(all_kmers_from_procs[i].seq, k));
      }
      kmer_str[KMER_LENGTH] = '\0';

      fprintf(f, "%s, %lu\n", kmer_str, all_kmers_from_procs[i].seq);
  }

  fclose(f);

  //f = fopen("mn_list.txt", "w");
  std::string outfile2="mn_list_"+std::to_string(rank)+".log";
  f = fopen(outfile2.c_str(), "w"); 

  char mn_str[MN_LENGTH+1];

  for (size_t i=0; i<list_of_mn_nodes.size(); i++) {

      for (int k=0; k<MN_LENGTH; k++) {
           mn_str[k] = el_to_char(kmerel_mn(list_of_mn_nodes[i], k));
      }
      mn_str[MN_LENGTH] = '\0';

      fprintf(f, "%s, %lu\n", mn_str, list_of_mn_nodes[i]);
  }

  fclose(f);

#endif

    double mn5 = MPI_Wtime ();
    //begin populating the prefixes and suffixes of nodes in MN_map
    //we block partition and assign distinct macro nodes to each of the threads
    
#pragma omp parallel
{
    kmer_t last_mn;
    kmer_t pred_k=0, succ_k=0;
    size_t chunk = (size_t)(list_of_mn_nodes.size()/num_threads);
    size_t lo = chunk*omp_get_thread_num();
    size_t hi = std::min(chunk*(omp_get_thread_num()+1), list_of_mn_nodes.size());

    if (omp_get_thread_num() == 0)
        lo = 0;

    if (omp_get_thread_num() == num_threads-1)
        hi = list_of_mn_nodes.size()-1;

    size_t lo_s = lo;
    size_t hi_s = hi;

    size_t k_low_pos = bsearch_pkmer_exists (list_of_mn_nodes[lo], all_kmers_from_procs);
    assert(k_low_pos >= 0);
    size_t k_hi_pos = bsearch_pkmer_exists (list_of_mn_nodes[hi], all_kmers_from_procs);
    assert(k_hi_pos >= 0);

    //kmer_t tmp_pred_k = 0;
    if (omp_get_thread_num() == num_threads-1){
        //printf("rank: %d, num_kmers: %lu, num_mn: %lu\n", rank, all_kmers_from_procs.size(), MN_map.size());

        /*printf("rank: %d, before: hi: %lu, mn_at_hi: %lu, k_hi_pos: %lu, kmer_at_hi: %lu\n", 
                rank, hi, list_of_mn_nodes[hi], k_hi_pos, all_kmers_from_procs[k_hi_pos].seq);
        tmp_pred_k = (kmer_t)all_kmers_from_procs[k_hi_pos].seq ^ 1;
        tmp_pred_k = mn_shift(tmp_pred_k);
        printf("rank: %d, predk at hi: %lu\n", rank, tmp_pred_k);
        */
        k_hi_pos = all_kmers_from_procs.size();

        //printf("rank: %d, after: hi: %lu, mn_at_hi: %lu, k_hi_pos: %lu, kmer_at_high: %lu\n", 
        //        rank, hi, list_of_mn_nodes[hi], k_hi_pos, all_kmers_from_procs[k_hi_pos].seq);
    }

    //iterate over kmer list to populate the suffixes
    last_mn = MN_map[lo].first;
    for (size_t f=k_low_pos; f<k_hi_pos; f++) 
    {
        pred_k=0, succ_k=0;
        kmer_t full_kmer = all_kmers_from_procs[f].seq;
        int kmer_count = all_kmers_from_procs[f].k_count;

        pred_k = (kmer_t)full_kmer ^ 1;
        pred_k = mn_shift(pred_k);

        if (retrieve_proc_id(pred_k) == rank)
        {

            if (pred_k != last_mn) {
                /*while (1) {
                    printf("rank: %d, in_while: thread: %d, pred_k: %lu, last_mn: %lu, LO: %lu, k_lo: %lu\n", 
                            rank, omp_get_thread_num(), pred_k, last_mn, lo, f);
                    if (MN_map[lo].first == pred_k)
                        break;
                    else
                        lo++;
                }*/
                while (MN_map[lo].first != pred_k) {
                       lo++;
                }
                if (lo > hi)
                    break;
                else{
                    last_mn = MN_map[lo].first;
                }
            }

            assert(pred_k == last_mn);

            //if (pred_k == tmp_pred_k)
            //    printf("in rank: %d, found last mn: %lu, lo: %lu, mn_at_lo: %lu\n", rank, tmp_pred_k, lo, MN_map[lo].first);

            MacroNode &mn_val = MN_map[lo].second;
            if (mn_val.k_1_mer.size() == 0) 
            {
                BasePairVector mybvp;
                for (int k=0; k<MN_LENGTH; k++){
                     BasePair x = kmerel_mn(pred_k, k);
                     mybvp.push_back(x);
                }
                mn_val.k_1_mer = mybvp;
            }

            BasePairVector suff;
            suff.push_back(kmerel(full_kmer, MN_LENGTH));
            mn_val.suffixes.push_back(suff);
            mn_val.suffixes_terminal.push_back(false);
            mn_val.suffix_count.push_back(std::make_pair(kmer_count, visit_value(kmer_count)));
        }

    }

//#pragma omp barrier

    lo = lo_s;
    hi = hi_s;
    k_low_pos = bsearch_skmer_exists (list_of_mn_nodes[lo], all_kmers_suff);
    assert(k_low_pos >= 0);
    k_hi_pos = bsearch_skmer_exists (list_of_mn_nodes[hi], all_kmers_suff);
    assert(k_hi_pos >= 0);

    if (omp_get_thread_num() == num_threads-1)
        k_hi_pos = all_kmers_suff.size();

    //iterate over kmer list to populate the prefixes
    last_mn = MN_map[lo].first;
    for (size_t f=k_low_pos; f<k_hi_pos; f++)
    {
        pred_k=0, succ_k=0;
        kmer_t full_kmer = all_kmers_suff[f].seq;
        int kmer_count = all_kmers_suff[f].k_count;

        succ_k = (kmer_t)full_kmer & (kmer_t)SUCC_MASK;

        if (retrieve_proc_id(succ_k) == rank)
        {

            if (succ_k != last_mn) {
                //printf("Suffix::thread: %d, lo: %lu, hi: %lu, succ_k: %lu, last_mn: %lu\n", 
                 //omp_get_thread_num(),lo,hi,succ_k,MN_map[lo].first);
                /*while (1) {
                    if (MN_map[lo].first == succ_k)
                        break;
                    else
                        lo++;
                }*/
                while (MN_map[lo].first != succ_k) {
                    lo++;
                }
                if (lo > hi)
                    break;
                else {
                    last_mn = MN_map[lo].first;
                }
            }

            assert (succ_k == last_mn);
            MacroNode &mn_val = MN_map[lo].second;
            if (mn_val.k_1_mer.size() == 0)
            {
                BasePairVector mybvp;
                for (int k=0; k<MN_LENGTH; k++) {
                         BasePair x = kmerel_mn(succ_k, k);
                         mybvp.push_back(x);
                    }
                    mn_val.k_1_mer = mybvp;
            }
        
            BasePairVector pred;
            pred.push_back(kmerel(full_kmer, 0));
            mn_val.prefixes.push_back(pred);
            mn_val.prefixes_terminal.push_back(false);
            mn_val.prefix_count.push_back(std::make_pair(kmer_count, visit_value(kmer_count)));
        }

    }

}// end of parallel region

    double mn6 = MPI_Wtime ();
    double populate_time=0.0, global_populate_time=0.0;

    populate_time = mn6 - mn5;
    MPI_Reduce(&populate_time, &global_populate_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for populating the suff and pref buffers (secs): %f \n", 
                            (double)global_populate_time/(double)size);

    list_of_mn_nodes.clear();
    // delete glist_of_mn_nodes
    for (int i=0; i<num_threads; i++)
         glist_of_mn_nodes[i].clear();

    all_kmers_suff.clear();
    all_kmers_suff.shrink_to_fit();

}

void begin_mnode_construction(std::vector<std::pair<kmer_t,MacroNode>> &MN_map)
{
    double start_time = MPI_Wtime ();

    std::vector< std::vector< std::vector<KmerPairs> > > global_mn_list (num_threads);
    std::vector <KmerPairs> all_kmers_from_procs = construct_kmer_histogram(global_mn_list);

    if (rank==0)
        fprintf(stderr, "Begin macro node generation\n");
        
//#ifdef CONSTRUCT_MNODES

    double mn3 = MPI_Wtime ();
    // Iterate through k-mers and contruct the macro_nodes
    if (rank==0)
        fprintf(stderr, "Starting construct_macro_nodes\n");
        
    construct_macro_nodes(MN_map, all_kmers_from_procs);
    all_kmers_from_procs.clear();
    all_kmers_from_procs.shrink_to_fit();
    double mn4 = MPI_Wtime ();

    double map_time=0.0, global_map_time=0.0;
    map_time = mn4 - mn3;
    MPI_Reduce(&map_time, &global_map_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for constructing MN map across all procs (secs): %f \n", 
                            (double)global_map_time/(double)size);


    uint64_t MN_map_size = MN_map.size();
    uint64_t global_MN_map_size=0;
    //printf("Rank: %d, number of macro_nodes created: %lu\n", rank, MN_map_size);

    MPI_Allreduce(&MN_map_size, &global_MN_map_size, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0)
        printf("Total number of macro_nodes across all proc's: %lu\n", global_MN_map_size);


    double end_time = MPI_Wtime ();
    double mn_contruct_time=0.0, global_mn_contruct_time=0.0;

    mn_contruct_time = end_time - start_time;
    MPI_Reduce(&mn_contruct_time, &global_mn_contruct_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for MN node construction across all procs (secs): %f \n", 
                            (double)global_mn_contruct_time/(double)size);

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0)
        fprintf(stderr, "Completed MN node construction. Starting wiring\n");
        
//#endif // end of CONSTRUCT_MNODES

}

