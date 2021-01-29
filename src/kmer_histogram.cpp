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
#include <set>
#include <algorithm>
#include <unordered_map>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "distribute_kmers.h"
#include "timers.h"

extern long int MAX_KMER_COUNT;
extern int rank, size;
extern int coverage;
extern int num_buckets;
extern int num_threads;

extern std::vector<KmerPairs> kmer_proc_buf;


std::vector <KmerPairs> construct_kmer_histogram (std::vector< std::vector< std::vector<KmerPairs> > > &global_mn_list)
{
    size_t local_kmer_count = kmer_proc_buf.size();
   
    if (rank==0) printf("Number of threads: %d\n", num_threads);
    std::vector<std::vector<uint64_t>> global_bucket_list (num_threads);
    int min_bucket=0;

    if (num_buckets)
    {
    std::vector<uint64_t> all_hist_bucket(num_buckets,0);

#pragma omp parallel shared(local_kmer_count)
{
    // Iterate once to plot the histogram to calculate min_bucket
    std::vector<uint64_t> hist_bucket(num_buckets,0);

#pragma omp for
    for (size_t f=0; f<local_kmer_count; f++) {
         int val = kmer_proc_buf[f].k_count;
         if (val < num_buckets) {
             hist_bucket[val] += 1;
         }
    }
    global_bucket_list[omp_get_thread_num()] = hist_bucket;
    hist_bucket.clear();
#pragma omp barrier

#pragma omp for
    for (int i=0; i<num_buckets; i++) {
        for (int j=0; j<num_threads; j++) {
             all_hist_bucket[i] += global_bucket_list[j][i];
        }
    }

}

    std::vector<uint64_t> global_bucket(num_buckets,0);

    MPI_Allreduce(all_hist_bucket.data(), global_bucket.data(), num_buckets, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(hist_bucket.data(), global_bucket.data(), num_buckets, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

    std::vector<uint64_t>::iterator result = std::min_element(std::begin(global_bucket)+1, std::end(global_bucket));
    min_bucket = std::distance(std::begin(global_bucket), result);

    all_hist_bucket.clear();
    global_bucket.clear();

    if (rank==0)
        fprintf(stderr, "Identified min_bucket\n");
        
    } // if condition checking num_buckets>0

    std::vector<uint64_t> valid_entries_per_thrd (num_threads);
    uint64_t valid_entries=0, all_valid_entries=0;

    // Iterate once again to obtain the number of valid_entries, given the min_bucket cutoff
    double mn1 = MPI_Wtime ();

#pragma omp parallel
{
    std::vector< std::vector<KmerPairs> > mn_list(size);
    uint64_t valid_entries_per_t=0;
    kmer_t pred_k=0, succ_k=0;
    int home_p=0, home_s=0;

#pragma omp for
    for (size_t f=0; f<local_kmer_count; f++) {
         if (kmer_proc_buf[f].k_count >= min_bucket) 
         {
             valid_entries_per_t++;
             KmerPairs k_entry;
             kmer_t full_kmer = kmer_proc_buf[f].seq;
             //k_entry.seq = kmer_proc_buf[f];
             //k_entry.k_count = kmer_cnt_proc_buf[f];
             k_entry.seq = kmer_proc_buf[f].seq;
             k_entry.k_count = kmer_proc_buf[f].k_count;

             pred_k = (kmer_t)full_kmer ^ 1;
             pred_k = mn_shift(pred_k);
             home_p = retrieve_proc_id(pred_k);

             succ_k = (kmer_t)full_kmer & (kmer_t)SUCC_MASK;
             home_s = retrieve_proc_id(succ_k);
             if (home_p == home_s)
                 mn_list[home_p].push_back(k_entry);
             else {
                 mn_list[home_p].push_back(k_entry);
                 mn_list[home_s].push_back(k_entry);
             }
         }
    }

    global_mn_list[omp_get_thread_num()] = mn_list;
    valid_entries_per_thrd[omp_get_thread_num()] = valid_entries_per_t;

    for (int i=0; i<size; i++)
         mn_list[i].clear();

} // end of parallel region

    #pragma omp parallel for reduction(+:valid_entries)
        for (int t=0; t<num_threads; t++)
             valid_entries += valid_entries_per_thrd[t];

    MPI_Allreduce(&valid_entries, &all_valid_entries, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Reduce(&valid_entries, &all_valid_entries, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("Total valid (>= min_bucket(%d)) k-mer entries across all proc's: %lu \n", min_bucket, all_valid_entries); 
        printf("**********************************************************************\n\n");
    }

    valid_entries_per_thrd.clear();
    free_kmer_count_buffers();
    MPI_Barrier(MPI_COMM_WORLD);

    if(!all_valid_entries) {
        //cleanup macro-node setup metadata
        for (int i=0; i<num_threads; i++) {
           for (int t=0; t<size; t++) {
             global_mn_list[i][t].clear();
           }
        global_mn_list[i].clear();
       }

       std::cerr << "Number of k-mers post pruning based on frequency cuttoff at: "<< min_bucket << " are null. Contig generation cannot proceed." << std::endl;
       MPI_Abort(MPI_COMM_WORLD, -99);
    }

    std::vector<int> scounts (size,0);
    size_t ssize=0;
    size_t off=0;
    size_t num_kentries=0;

    for (int i=0; i<size; i++) {
        num_kentries=0;
    #pragma omp parallel for reduction(+:num_kentries)
        for (int t=0; t<num_threads; t++)
            num_kentries += global_mn_list[t][i].size();

        scounts[i] = num_kentries;
        ssize += scounts[i];
    }

    std::vector<KmerPairs> mn_send_buf(ssize);

    for (int i=0; i<size; i++) {
        for (int t=0; t<num_threads; t++) {
            std::copy(global_mn_list[t][i].begin(), global_mn_list[t][i].end(), mn_send_buf.begin()+off);
            off += global_mn_list[t][i].size();
        }
    }

    for (int i=0; i<num_threads; i++) {
        for (int t=0; t<size; t++) {
             global_mn_list[i][t].clear();
        }
        global_mn_list[i].clear();
    }

    if (rank==0)
        fprintf(stderr, "Starting transfer_macro_nodes \n");
        
    std::vector <KmerPairs> all_kmers_from_procs = transfer_macro_nodes (mn_send_buf, scounts);

    double mn2 = MPI_Wtime ();
    double regroup_kmer_time=0.0, global_regroup_kmer_time=0.0;

    regroup_kmer_time = mn2 - mn1;
    MPI_Reduce(&regroup_kmer_time, &global_regroup_kmer_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for Re-grouping for k-mers across all procs (secs): %f \n", 
                            (double)global_regroup_kmer_time/(double)size);


    return all_kmers_from_procs;

}
