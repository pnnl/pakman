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

#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "distribute_kmers.h"
#include "timers.h"
//#include "serial.h"

long int MAX_KMER_COUNT=0;
int rank, size;
int coverage=0;
std::string inputFileName;
int read_length=0;
int num_buckets=0;
int num_threads=0;
int num_contigs=0;
int node_threashold=0;

int num_batch_transfers=0;
long int this_contig_id = 0;

std::vector<lmer_t> lmer_frequency(LMER_SIZE,0);
std::vector<lmer_t> global_lmer_frequency(LMER_SIZE,0);

/* Vector containing k-mers allocated to each process */
std::vector<KmerPairs> kmer_proc_buf;

//kmer_t *hasha, *hashb;
uint64_t hasha = 68111;
uint64_t hashb = 105929;

void parseCommandLine(const int argc, char * const argv[]);

int retrieve_proc_id (lmer_t min_lmer)
{
    //determine which proc holds the l-mer bucket based on a block-cyclic
    //distribution of buckets across the processes

    //int proc_num = (int)((int)(floor(min_lmer/4)) % (int)size); // static block-cyclic
    //int proc_num = (int)((int)(hash31(hasha[min_lmer], hashb[min_lmer], min_lmer)) % (int)size); // randomized
    int proc_num = (int)((int)(hash31(hasha, hashb, min_lmer)) % (int)size); // randomized

    return proc_num;
}
#ifdef EXTEND_KMER
int retrieve_proc_id (kmer_t min_lmer)
{
    //determine which proc holds the l-mer bucket based on a block-cyclic
    //distribution of buckets across the processes

    //int proc_num = (int)((int)(floor(min_lmer/4)) % (int)size); // static block-cyclic
    //int proc_num = (int)((int)(hash31(hasha[min_lmer], hashb[min_lmer], min_lmer)) % (int)size); // randomized

    //int proc_num = (int)((int)(hash31(hasha, hashb, min_lmer)) % (int)size); // randomized
    int proc_num = (int)((int)(uhash31(hasha, hashb, min_lmer)) % (int)size);

    return proc_num;
}
#endif

void set_num_threads()
{

#pragma omp parallel
{
    num_threads = omp_get_num_threads();
}

}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    parseCommandLine(argc, argv);
    set_num_threads();

    
    input_read_data rdata = perform_input_reading(rank, size, inputFileName, read_length);


    perform_kmer_counting (rdata.read_data, rdata.read_data_size);
    free(rdata.read_data);

#ifndef PERFORM_CONTIG_GENERATION

    free_kmer_count_buffers();

#endif

#ifdef PERFORM_CONTIG_GENERATION

    /************ CONTIG GENERATION BEGINS *******************/
    /* Phase 1) Macro Node construction and Initial Wiring
     * Phase 2) Construction of Independent Set. Merge and delete all candidate nodes from Independent set.
     * Phase 3) Gather all remaining nodes
     * Phase 4) Enumerate/Generate the contigs with the Walk algorithm
     */


    /* PHASE 1: MACRO NODE CONSTRUCTION */
    std::vector<std::pair<kmer_t,MacroNode>> MN_map;
    begin_mnode_construction(MN_map);


    /* PHASE 1: MACRO NODE WIRING */
    initiate_mnode_wiring(MN_map);

#ifdef DEBUG_WIRE_INIT
    debug_wired_mnodes(MN_map);
#endif

    /* PHASE 2: INDEPENDENT SET CONSTRUCTION */
    
    // temp vector for storing all partial contigs generated during Phase 2
    std::vector<BasePairVector> partial_contig_list;
    size_t global_num_nodes = begin_iterative_compaction(MN_map, partial_contig_list);

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef COMPACT_PGRAGH
    //retain a list of terminal prefixes for each individual process, potential begin k-mers
    std::vector<BeginMN> list_of_begin_kmers;
    identify_begin_kmers (MN_map, list_of_begin_kmers);

#ifdef DEBUG_BKMERS
    debug_begin_kmers_list (MN_map, list_of_begin_kmers);
#endif

    /* Perform an Allgather such that all macro_nodes are accessible to all procs */

    std::vector<std::pair<kmer_t,MacroNode>> global_MN_map(global_num_nodes); // new map for storing all the macro_nodes
    generate_compacted_pakgraph(MN_map, global_MN_map);

    traverse_pakgraph(global_MN_map, list_of_begin_kmers, partial_contig_list);
#endif // end of COMPACT_PGRAGH

#endif // end of PERFORM_CONTIG_GENERATION

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

void parseCommandLine(const int argc, char * const argv[])
{
  int ret;

  while ((ret = getopt(argc, argv, "f:b:r:c:t:n:")) != -1) {
    switch (ret) {
    case 'f':
       inputFileName.assign(optarg);
       //std::cout << inputFileName << std::endl;
       break;
    case 'b':
       MAX_KMER_COUNT = atol(optarg);
       //std::cout << MAX_KMER_COUNT << std::endl;
       break;
    case 'r':
       read_length = atoi(optarg);
       //std::cout << read_length << std::endl;
       break;
    case 'c':
       coverage = atoi(optarg);
       //std::cout << coverage << std::endl;
       break;
    case 't':
       num_buckets = atoi(optarg);
       //std::cout << num_buckets << std::endl;
       break;
    case 'n':
       node_threashold = atoi(optarg);
       //std::cout << node_threashold << std::endl;
       break;
    default:
       assert(0 && "Should not reach here!!");
       break;
    }
  }

  if (rank==0 && inputFileName.empty()) {
      std::cerr << "Must specify an input FASTA file name with -f" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }

//  if (rank==0 && !MAX_KMER_COUNT) {
//      std::cerr << "Must specify a batch size for k-mer counting with -b" << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, -99);
//  }

  if (!MAX_KMER_COUNT) {
      MAX_KMER_COUNT=100000000;
      if (rank==0) std::cout << "batch size not specified with -b, set to default value MAX_KMER_COUNT=100000000" << std::endl;
  }

  if (rank==0 && !read_length) {
      std::cerr << "Must specify read_length with -r" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }

  if (rank==0 && !coverage) {
      std::cerr << "Must specify coverage of the input data with -c" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }

  if (rank==0 && (read_length>250 || read_length<100)) {
      std::cerr << "Must provide short reads of length >100 and <=250 with -r" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }

  if (!num_buckets) {
      num_buckets=21;
      if (rank==0) std::cout << "num_buckets not specified with -t, set to default value num_buckets=21" << std::endl;
  }

  if (!node_threashold) {
      node_threashold=100000;
      if (rank ==0 )std::cout << "node_threashold not specified with -n, set to default value node_threashold=100K" << std::endl;
  }

  if (rank == 0) {
           printf("K-mer size: %d, L-mer size: %d, Number of Processes: %d, MAX_KMER_COUNT: %ld Coverage: %d,",
           WINDW_SIZE+1, LMER_LENGTH, size, MAX_KMER_COUNT, coverage);
           printf("Avg read length: %d, Num buckets: %d, Node threshold: %d\n", 
           read_length, num_buckets, node_threashold);
  }


} // parseCommandLine


