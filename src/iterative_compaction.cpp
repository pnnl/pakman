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
#include <fstream>
#include "distribute_kmers.h"
#include "timers.h"

extern int rank, size;
extern int coverage;
extern int num_threads;
extern int node_threashold;

void identify_begin_kmers (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                               std::vector<BeginMN> &list_of_begin_kmers)
{

//#ifdef CAL_BKMERS
    //for (std::unordered_map<BasePairVector,MacroNode,KeyHasher>::const_iterator it = MN_map.begin();
    //for (std::vector<std::pair<kmer_t,MacroNode>>::iterator it = MN_map.begin();
    //        it != MN_map.end(); it++) {
    for (size_t it=0; it<MN_map.size(); it++)
    {
        MacroNode &mn = MN_map[it].second;
        for (size_t l=0; l<mn.prefixes.size(); l++) {
            if (mn.prefix_count[l].second > 0 && mn.prefixes_terminal[l]) {
                list_of_begin_kmers.push_back(BeginMN{MN_map[it].first, (int)l});
            }
        }
    }

#ifdef SANITY_CHECK
    printf("rank: %d, number of begin kmers: %lu\n", rank, list_of_begin_kmers.size());
#endif

//#endif //end of CAL_BKMERS


}

void generate_compacted_pakgraph(std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                                 std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map)
{

//#ifdef SERIALIZE

    if (rank==0)
        fprintf(stderr, "Begin Serialize and gather\n");
        
    double mn7 = MPI_Wtime();
    //std::vector<std::pair<kmer_t,MacroNode>> global_MN_map(global_MN_map_size); // new map for storing all the macro_nodes
    gather_all_macro_nodes (MN_map, global_MN_map);

    double mn8 = MPI_Wtime ();
    double cereal_time=0.0, global_cereal_time=0.0;
    cereal_time = mn8 - mn7;
    
    MPI_Reduce(&cereal_time, &global_cereal_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for cereal and Allgather across all procs (secs): %f \n\n", 
                            (double)global_cereal_time/(double)size);

    MPI_Barrier(MPI_COMM_WORLD);

    // sanity check
#ifdef SANITY_CHECK
    printf("rank: %d, number of macro_nodes: %lu\n", rank, global_MN_map.size());
#endif

    if (rank==0)
        fprintf(stderr, "Completed Serialize and gather. Starting Walk phase\n");
        
//#endif //end of SERIALIZE

}


size_t begin_iterative_compaction (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                                   std::vector<BasePairVector> &partial_contig_list)
{

//#ifdef PHASE2
    size_t num_nodes=0, global_num_nodes=0;
#ifdef DEBUG_P2TIME
    size_t id_set_size=0, global_id_set_size=0;

#ifdef PRINT_P2STATS
    FILE *fc1=NULL, *fc2=NULL, *fc3=NULL, *fc4=NULL, *fc5=NULL;

    if (rank==0)
    {
    std::string output_filename1("p2_idset.dat");
    std::string output_filename2("p2_itrpack.dat");
    std::string output_filename3("p2_sermod.dat");
    std::string output_filename4("p2_MN.dat");
    std::string output_filename5("p2_ID.dat");

    fc1 = fopen(output_filename1.c_str(), "w");
    fc2 = fopen(output_filename2.c_str(), "w");
    fc3 = fopen(output_filename3.c_str(), "w");
    fc4 = fopen(output_filename4.c_str(), "w");
    fc5 = fopen(output_filename5.c_str(), "w");
    }
#endif // end of PRINT_P2STATS

#endif
    int num_itr=0;
    std::vector<size_t> id_list_nodes;
    std::vector<size_t> rewire_pos_list;

    /* PHASE 2 TIMERS */
    double merge_mn=0.0, global_merge_mn=0.0;
    double idset_create=0.0, global_idset_create=0.0;
    double pack_mn=0.0, global_pack_mn=0.0;
    double update_mn=0.0, global_update_mn=0.0;

    double mn11 = MPI_Wtime ();

    if (rank==0)
        fprintf(stderr, "Starting phase 2 iterative phase\n");
        
    while (1)
    { 
     num_nodes = MN_map.size();
     ++num_itr; 
     id_list_nodes.clear();
     rewire_pos_list.clear();

     MPI_Allreduce(&num_nodes, &global_num_nodes, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#ifdef DEBUG_P2TIME
     if (rank == 0)
         printf("Itr: %d, Total number of macro_nodes across all proc's: %lu\n", num_itr, global_num_nodes);

#ifdef PRINT_P2STATS
     if (rank == 0)
         fprintf(fc4,"%d, %lu\n", num_itr, global_num_nodes);
#endif //end of PRINT_P2STATS

#endif

     if (global_num_nodes <= node_threashold)
     {
        if (rank==0)
            fprintf(stderr, "Completed iterative phase 2\n");
      
        #ifdef PRINT_P2STATS
        if (rank==0) {
            fclose(fc1);  
            fclose(fc2);
            fclose(fc3);
            fclose(fc4);
            fclose(fc5);
        }
        #endif // end PRINT_P2STATS
 
        double mn12 = MPI_Wtime ();
        merge_mn = mn12 - mn11;

        //print_output_contigs(partial_contig_list, std::string("contigs_par_r"));
        //partial_contig_list.clear();
    
        MPI_Reduce(&idset_create, &global_idset_create, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) printf ("Average time for part1: id set creation across all procs (secs): %f \n", 
                            (double)global_idset_create/(double)size);

        MPI_Reduce(&pack_mn, &global_pack_mn, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) printf ("Average time for part2: packing nodes across all procs (secs): %f \n", 
                            (double)global_pack_mn/(double)size);

        MPI_Reduce(&update_mn, &global_update_mn, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) printf ("Average time for part3: update nodes across all procs (secs): %f \n", 
                            (double)global_update_mn/(double)size);

        MPI_Reduce(&merge_mn, &global_merge_mn, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) printf ("Average time for phase to merge and dropping of nodes across all procs (secs): %f \n\n", 
                            (double)global_merge_mn/(double)size);

        MPI_Reduce(&p2alltoall_time, &p2global_alltoall_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) printf ("Average time for Alltoall in Serialize and Transfer (secs): %f \n", 
                            (double)p2global_alltoall_time/(double)size);

        MPI_Reduce(&p2alltoallv_time, &p2global_alltoallv_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) printf ("Average time for AlltoallV in Serialize and Transfer (secs): %f \n", 
                            (double)p2global_alltoallv_time/(double)size);
         
        return global_num_nodes;
         //break;
     }

#ifdef DEBUG_WIRE
    const std::string s0("map_all_mn_out_");
    print_map_content_itr(MN_map, s0, num_itr);
    
#endif

     double p1 = MPI_Wtime ();
     generate_id_set(MN_map, id_list_nodes);
     double p2 = MPI_Wtime ();
     idset_create += (p2 - p1);
     //print_id_list(id_list_nodes, std::string("id_list"), MN_map);
     //print_id_list(id_list_nodes, std::string("id_list"));

#ifdef DEBUG_P2TIME
     double temp_id = (p2 - p1);
     double g_temp_id=0.0;

     MPI_Reduce(&temp_id, &g_temp_id, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     if (rank == 0) printf ("part1: id set creation (itr: %d)  (secs): %f \n", num_itr, 
                            (double)g_temp_id/(double)size);

#ifdef PRINT_P2STATS                              
     if (rank == 0) 
         fprintf (fc1, "%d, %f\n", num_itr, (double)g_temp_id/(double)size);
#endif // end of PRINT_P2STATS
                              
     //if (rank==0)
     //   fprintf(stderr, "part1: id set creation for itr: %d\n", num_itr);
        
#endif

    /* PHASE 2: ITERATIVE PHASE TO MERGE MACRO_NODES */

#ifdef DEBUG_P2TIME
    id_set_size = id_list_nodes.size();
    MPI_Allreduce(&id_set_size, &global_id_set_size, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0)
        printf("Itr: %d, global ID_set size across all proc's: %lu\n", num_itr, global_id_set_size);

#ifdef PRINT_P2STATS
    if (rank == 0)
        fprintf(fc5,"%d, %lu\n", num_itr, global_id_set_size);
#endif // end of PRINT_P2STATS

#endif


//#endif //end of PHASE2

//#ifdef MERGE

    p1 = MPI_Wtime ();

    //std::vector< std::vector<ModNodeInfo> > mn_nodes_per_proc(size);
    std::vector< std::vector<TransferNode> > mn_nodes_per_proc(size);
    iterate_and_pack_mn (id_list_nodes, mn_nodes_per_proc, MN_map, partial_contig_list);

    //printf("Itr: %d, mn_nodes_per_proc size for proc: %d is %lu\n", num_itr, rank, mn_nodes_per_proc[rank].size());

     size_t new_vec_size = MN_map.size() - id_list_nodes.size();
     size_t last_id = MN_map.size()-1;
     for (int64_t i=id_list_nodes.size()-1; i>=0; i--)
     {
          size_t del_key = id_list_nodes[i];
          //if (!(del_key >= 0 && del_key < MN_map.size())) 
          //   printf("rank: %d, del_key: %lu, id: %lu, Map_size: %lu\n", rank, del_key, i, MN_map.size());
          assert(del_key >= 0 && del_key < MN_map.size());

          std::swap(MN_map[del_key], MN_map[last_id]);
          last_id--;
     }

     MN_map.resize(new_vec_size);
     //MN_map.shrink_to_fit();
     __gnu_parallel::sort(MN_map.begin(), MN_map.end(), Comp_pair); 

     MPI_Barrier(MPI_COMM_WORLD);

     p2 = MPI_Wtime ();
     pack_mn += (p2 - p1);

     //if (rank==0)
     //    printf("Itr: %d, iterate and pack (secs): %f \n", num_itr, pack_mn);


#ifdef DEBUG_P2TIME
     temp_id = (p2 - p1);
     g_temp_id=0.0;

     MPI_Reduce(&temp_id, &g_temp_id, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     if (rank == 0) printf ("part2: iterate and pack  (itr: %d)  (secs): %f \n", num_itr,
                            (double)g_temp_id/(double)size);

#ifdef PRINT_P2STATS
     if (rank == 0) 
         fprintf (fc2, "%d, %f\n", num_itr, (double)g_temp_id/(double)size);
#endif // end of PRINT_P2STATS

     //if (rank==0)
     //   fprintf(stderr, "part2: iterate and pack for itr: %d\n", num_itr);
        
#endif

//#endif // end of MERGE

#ifdef SER_PACK_V1

     p1 = MPI_Wtime ();

     std::vector<ModNodeInfo> mn_nodes_to_modify = serialize_and_transfer(mn_nodes_per_proc);

    //Modifying the selected MN nodes
    //rewire_pos_list.reserve(mn_nodes_to_modify.size());
 
    for (size_t i=0; i<mn_nodes_to_modify.size(); i++)
    {
        kmer_t search_key = mn_nodes_to_modify[i].modifyMN.search_mn;

        int pos = find_mnode_exists(search_key, MN_map);
        if (MN_map[pos].first != search_key) {
            printf("rank: %d, MN node key: %lu, was not found in map at the time of pushing to: %d\n", 
                    rank, search_key, mn_nodes_to_modify[i].direction);
            MPI_Finalize();
            exit(2);
        }

        if (mn_nodes_to_modify[i].direction==P)
           push_to_pred(mn_nodes_to_modify[i], pos, MN_map);
        else
           push_to_succ(mn_nodes_to_modify[i], pos, MN_map);

        rewire_pos_list.push_back(pos); 
    }

#endif //end of SER_PACK_V1

//#ifdef MOD_REWIRE
#ifdef DEBUG_WIRE
    const std::string s("map_mn_out_");
    print_map_content_itr(MN_map, s, num_itr);
    
#endif

    p1 = MPI_Wtime ();

    id_list_nodes.clear();
    serialize_and_transfer(mn_nodes_per_proc, MN_map, rewire_pos_list, num_itr);

    std::sort(rewire_pos_list.begin(), rewire_pos_list.end());
    rewire_pos_list.erase( unique( rewire_pos_list.begin(), rewire_pos_list.end() ), rewire_pos_list.end() );

#ifdef DEBUG_WIRE
    const std::string st("map_rewire_");
    print_rewire_list(rewire_pos_list, MN_map, st, num_itr);
    
#endif
    for (size_t i=0; i<rewire_pos_list.size(); i++) {
         MacroNode &mn = MN_map[rewire_pos_list[i]].second;

         mn.wiring_info.clear();
         mn.prefix_begin_info.clear();

         int num_p = mn.prefixes.size();
         int num_s = mn.suffixes.size();

         mn.wiring_info.resize(num_p+num_s+1);                                                                                 
         mn.prefix_begin_info.resize(num_p);

         mn.setup_wiring();
    }

#ifdef SER_PACK_V1
    mn_nodes_to_modify.clear();
#endif
    rewire_pos_list.clear();

     p2 = MPI_Wtime ();
     update_mn += (p2 - p1);

     //if (rank==0)
     //    printf("Itr: %d, serialize and update nodes (secs): %f \n", num_itr, update_mn);

#ifdef DEBUG_P2TIME
     temp_id = (p2 - p1);
     g_temp_id=0.0;

     MPI_Reduce(&temp_id, &g_temp_id, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     if (rank == 0) printf ("part3: serialize and update nodes (itr: %d)  (secs): %f \n\n", num_itr,
                            (double)g_temp_id/(double)size);

#ifdef PRINT_P2STATS
     if (rank == 0) 
         fprintf (fc3, "%d, %f\n", num_itr, (double)g_temp_id/(double)size);
#endif //end of PRINT_P2STATS

     //if (rank==0)
     //   fprintf(stderr, "part3: serialize and update nodes for itr: %d\n", num_itr);
        
#endif

//#endif // end of MOD_REWIRE

 } // end of while loop

    

}
 
