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
#include <iostream>
#include <sstream>
#include <utility>
#include "distribute_kmers.h"
#include "serialize.h"
//#include "bigmpi.h"

extern int rank, size;
extern int coverage;
extern int num_threads;

//Timers
extern double p2alltoall_time;
extern double p2alltoallv_time;

void generate_id_set (std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
        std::vector<size_t> &id_list_nodes)
{

#ifdef DEBUG_IDSET
    char proc_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    strcpy(output_file_name,"debug_idset_");
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
#endif

    for (size_t it=0; it<MN_map.size(); it++)
    {
          kmer_t key = MN_map[it].first;
          kmer_t tmp_kmer=0;
          kmer_t max_kmer=key;
          MacroNode &mn = MN_map[it].second;
          //BasePairVector key_node = mn.k_1_mer;
          bool key_notin_idset=false;
#ifdef DEBUG_IDSET
        char mn_str[MN_LENGTH+1];
        for (int k=0; k<MN_LENGTH; k++) {
             mn_str[k] = el_to_char(kmerel_mn(key, k));
        }
        mn_str[MN_LENGTH] = '\0';
	    //fprintf(f, "pos: %lu, Key_id: %lu, key: %s\n", it, key, mn_str);
        fprintf(f, "pos: %lu, key: %s\n", it, mn_str);
#endif

          for (size_t i=0; i<mn.prefixes.size(); i++)
          {
             tmp_kmer=0;
             if (!mn.prefixes_terminal[i]) 
             { 
             int len_pref = mn.prefixes[i].size();
             if (len_pref >= MN_LENGTH) {

                 /*extract length=MN_LENGTH from the prefix BPV extension */
                 tmp_kmer = mn.prefixes[i].extract(MN_LENGTH);
                 //for (int k=0; k<MN_LENGTH; k++)
                 //     tmp_kmer = mnmer_shift(tmp_kmer, mn.prefixes[i][k]);
             }
             else if (len_pref > 0)
             {
                 size_t remainder = MN_LENGTH - len_pref;
                 tmp_kmer = mn.prefixes[i].extract(len_pref);
                 kmer_t temp_ext = mn_extract_pred(key, remainder);
                 tmp_kmer = ((tmp_kmer << (remainder*2)) | temp_ext);

                 //for (int k=0; k<len_pref; k++)
                 //     tmp_kmer = mnmer_shift(tmp_kmer, mn.prefixes[i][k]);
                 //for (int k=0; k<remainder; k++)
                 //     tmp_kmer = mnmer_shift(tmp_kmer, key_node[k]);
             }
             }
#ifdef DEBUG_IDSET
        if (tmp_kmer) {
        for (int k=0; k<MN_LENGTH; k++) {
             mn_str[k] = el_to_char(kmerel_mn(tmp_kmer, k));
        }
        mn_str[MN_LENGTH] = '\0';
	    //fprintf(f, "pred_id: %lu, p_mn: %s\n", tmp_kmer, mn_str);
        fprintf(f, "pred_id: p_mn: %s\n", mn_str);
        }
#endif

             if (tmp_kmer > key) {
                 max_kmer=tmp_kmer;
                 key_notin_idset=true;
                 break;
             }
          }

	  tmp_kmer=0;
          if (!key_notin_idset) {
              for (size_t i=0; i<mn.suffixes.size(); i++)
              {
                  tmp_kmer=0;
                  if (!mn.suffixes_terminal[i])
                  {
                  int len_suff = mn.suffixes[i].size();
                  if (len_suff >= MN_LENGTH) {
                      tmp_kmer = mn.suffixes[i].extract_succ(MN_LENGTH);
                      //int start = len_suff-MN_LENGTH;
                      //for (int k=start; k<len_suff; k++)
                      //    tmp_kmer = mnmer_shift(tmp_kmer, mn.suffixes[i][k]);
                  }
                  else if (len_suff > 0)
                  {
                      size_t remainder = MN_LENGTH - len_suff;
                      tmp_kmer = mn_extract_succ(key, remainder);
                      kmer_t temp_ext = mn.suffixes[i].extract(len_suff);
                      tmp_kmer = ((tmp_kmer << (len_suff*2)) | temp_ext);

                      //for (int k=len_suff; k<key_node.size(); k++)
                      //     tmp_kmer = mnmer_shift(tmp_kmer, key_node[k]);
                      //for (int k=0; k<len_suff; k++)
                      //     tmp_kmer = mnmer_shift(tmp_kmer, mn.suffixes[i][k]);
                  }
                  }
#ifdef DEBUG_IDSET
        if (tmp_kmer) {
        for (int k=0; k<MN_LENGTH; k++) {
             mn_str[k] = el_to_char(kmerel_mn(tmp_kmer, k));
        }
        mn_str[MN_LENGTH] = '\0';
	    //fprintf(f, "succ_id: %lu, s_mn: %s\n", tmp_kmer, mn_str);
        fprintf(f, "succ_id: s_mn: %s\n", mn_str);
        }
#endif

                  if (tmp_kmer > key) {
                      max_kmer=tmp_kmer;
                      key_notin_idset=true;
                      break;
                  }
              }
          }

          if (!key_notin_idset) {
              assert(max_kmer == key);
              //id_list_nodes.push_back((int)it);
              id_list_nodes.push_back(it);
              assert(it >= 0 && it < MN_map.size());
              //id_list_nodes.push_back(BeginMN{key,(int)it});
              //id_list_nodes.push_back(BeginMN{key,(int)std::distance(MN_map.begin()-it)});
          }
     } // end of for loop
}


void iterate_and_pack_mn (std::vector<size_t>& id_list_nodes,
                          std::vector< std::vector<TransferNode> >& mn_nodes_per_proc,
                          std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                          std::vector<BasePairVector> &local_contig_list)
{
 
     int itr_p=0, itr_s=0;

#ifdef OPT_ITR

     std::vector<int> prealloc_per_proc(size, 0);

     for (size_t i=0; i<id_list_nodes.size(); i++)
     {
         itr_p=0;
         itr_s=0;

         kmer_t &del_key = MN_map[id_list_nodes[i]].first;
         MacroNode &del_mn = MN_map[id_list_nodes[i]].second;

         for (int k=0; k<del_mn.prefix_begin_info.size(); k++)
         {
             if (del_mn.prefix_begin_info[k].num_wires)
             {
                 itr_p=k;

                 MnodeInfo search_pred_param;
                 if ((del_mn.prefixes[itr_p].size() > 0) && (!del_mn.prefixes_terminal[itr_p]))
                      search_pred_param = retrieve_mn_pinfo(del_mn.prefixes[itr_p], del_mn.k_1_mer);

                 for (int t=0; t<del_mn.prefix_begin_info[k].num_wires; t++)
                 {
                     itr_s=del_mn.wiring_info[del_mn.prefix_begin_info[k].prefix_pos+t].suffix_id;

                     if (!(del_mn.prefixes_terminal[itr_p] && del_mn.suffixes_terminal[itr_s])) {
                         MnodeInfo search_succ_param;
                         if ((del_mn.suffixes[itr_s].size() > 0) && (!del_mn.suffixes_terminal[itr_s]))
                              search_succ_param = retrieve_mn_sinfo(del_mn.suffixes[itr_s], del_mn.k_1_mer);

                         std::array<bool, 2> self_loop_info = {false,false};
                         check_for_self_loops(search_pred_param.search_mn,
                                              search_succ_param.search_mn,
                                              del_key,
                                              self_loop_info);

                         if (!del_mn.prefixes_terminal[itr_p]) {
                             if (!self_loop_info[0]) {
                                 prealloc_per_proc[retrieve_proc_id(search_pred_param.search_mn)]++;
                             }
                         }

                         if (!del_mn.suffixes_terminal[itr_s]) {
                             if (!self_loop_info[1]) {
                                 prealloc_per_proc[retrieve_proc_id(search_succ_param.search_mn)]++;
                             }
                         }
                      }
                   } // end of for loop over num_wires

              } // end of if condition
          }
      } // end of for over list of id_set nodes

      for (int i=0; i<size; i++)
           mn_nodes_per_proc[i].resize(prealloc_per_proc[i]);

      prealloc_per_proc.clear(); 
      std::vector <int> proc_idx (size,0);

#endif

#ifdef DEBUG_KSIZE
    char proc_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    strcpy(output_file_name, "merge_itr");
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
#endif

     //for (size_t i=0; i<id_list_nodes.size(); i++)
     //for (size_t i=0; i<5; i++)
     for (size_t i=0; i<id_list_nodes.size(); i++)
     {
         itr_p=0;
         itr_s=0;

         kmer_t &del_key = MN_map[id_list_nodes[i]].first;
         MacroNode &del_mn = MN_map[id_list_nodes[i]].second;

#ifdef DEBUG_MERGE
         BasePairVector mn_key = del_mn.k_1_mer;

         kmer_t test_mn=0;
         for (int l=0; l<MN_LENGTH; l++)
              test_mn = mnmer_shift(test_mn, mn_key[l]);
/*
         std::string out_mn;
         for (size_t l=0; l<mn_key.size(); l++)
              out_mn += el_to_char(mn_key[l]);

        char mn_str[MN_LENGTH+1];
        for (int l=0; l<MN_LENGTH; l++) {
             mn_str[l] = el_to_char(kmerel_mn(test_mn, l));
        }
        mn_str[MN_LENGTH] = '\0';
*/
        //printf("Rank: %d is processing key: %lu, basepair: %s, mn_id: %lu, mn: %s, pos: %d\n", rank, del_key,
        //                                         out_mn.c_str(), test_mn, mn_str, id_list_nodes[i].terminal_prefix_id);
        assert(test_mn == del_key);
#endif

         for (int k=0; k<del_mn.prefix_begin_info.size(); k++)
         {
             if (del_mn.prefix_begin_info[k].num_wires)
             {
                 itr_p=k;

                 //BasePairVector p_ext = del_mn.prefixes[itr_p];
                 MnodeInfo search_pred_param;
                 if ((del_mn.prefixes[itr_p].size() > 0) && (!del_mn.prefixes_terminal[itr_p]))
                      search_pred_param = retrieve_mn_pinfo(del_mn.prefixes[itr_p], del_mn.k_1_mer);

                 for (int t=0; t<del_mn.prefix_begin_info[k].num_wires; t++)
                 {
                     itr_s=del_mn.wiring_info[del_mn.prefix_begin_info[k].prefix_pos+t].suffix_id;
                     int count=del_mn.wiring_info[del_mn.prefix_begin_info[k].prefix_pos+t].count;

                     if (del_mn.prefixes_terminal[itr_p] && del_mn.suffixes_terminal[itr_s]) {
                         BasePairVector partial_contig = del_mn.prefixes[itr_p];
                         partial_contig.append(del_mn.k_1_mer);
                         partial_contig.append(del_mn.suffixes[itr_s]);
                         local_contig_list.push_back(partial_contig);
                         //output(partial_contig);
                     }
                     else
                     {
                         //BasePairVector s_ext = del_mn.suffixes[itr_s];
                         MnodeInfo search_succ_param;
                         if ((del_mn.suffixes[itr_s].size() > 0) && (!del_mn.suffixes_terminal[itr_s]))
                              search_succ_param = retrieve_mn_sinfo(del_mn.suffixes[itr_s], del_mn.k_1_mer);

                         std::array<bool, 2> self_loop_info = {false,false};
                         check_for_self_loops(search_pred_param.search_mn,
                                              search_succ_param.search_mn,
                                              del_key,
                                              self_loop_info);

                         if (!del_mn.prefixes_terminal[itr_p]) {
                             if (!self_loop_info[0]) {

                                 //determine new value terminal
                                 bool new_pnode_type=false;
                                 if (del_mn.suffixes_terminal[itr_s])
                                     new_pnode_type=true;
                                 else{
                                     if(self_loop_info[1])
                                         new_pnode_type=true;
                                     else
                                         new_pnode_type=false;
                                 }
                              //int is_terminal = (new_pnode_type) ? (KTerminal)Y : (KTerminal)N;
                                 /*ModNodeInfo send_pred_data;
                                 send_pred_data.modifyMN = search_pred_param;
                                 send_pred_data.modifyMN_info = MNnodeTuple{del_mn.suffixes[itr_s],
                                     std::make_pair(std::min(del_mn.prefix_count[itr_p].first,
                                                 del_mn.suffix_count[itr_s].first), count),
                                     new_pnode_type};
                                 send_pred_data.direction=(KdirType)P;
                                 */
#ifdef OPT_ITR
                                 int node_id = retrieve_proc_id(search_pred_param.search_mn);
                                 assert(proc_idx[node_id] < mn_nodes_per_proc[node_id].size());
                                 mn_nodes_per_proc[node_id][proc_idx[node_id]] = send_pred_data;
                                 proc_idx[node_id]++;
#else               
                                 mn_nodes_per_proc[retrieve_proc_id(search_pred_param.search_mn)].push_back(TransferNode{
                                     search_pred_param.search_mn, search_pred_param.search_ext, del_mn.suffixes[itr_s],
                                     std::make_pair(std::min(del_mn.prefix_count[itr_p].first, del_mn.suffix_count[itr_s].first), count),
                                     new_pnode_type, (KdirType)P}
                                     );
#ifdef DEBUG_KSIZE
                                 char mn_str1[MN_LENGTH+1], mn_str2[MN_LENGTH+1];
                                 for (int k=0; k<MN_LENGTH; k++) {
                                 mn_str1[k] = el_to_char(kmerel_mn(del_key, k));
                                 mn_str2[k] = el_to_char(kmerel_mn(search_pred_param.search_mn, k));
                                 }
                                 mn_str1[MN_LENGTH] = '\0';
                                 mn_str2[MN_LENGTH] = '\0';
#ifdef EXTEND_KMER
                                fprintf(f, "Pred: prefix_size: %lu, itr_p: %d, key_str: %s, new_key_str: %s dest rank: %d\n",
                                             del_mn.prefixes[itr_p].size(), itr_p, mn_str1, mn_str2, retrieve_proc_id(search_pred_param.search_mn));
#else

                                 fprintf(f, "Pred: prefix_size: %lu, itr_p: %d, key:%lu, key_str: %s, new_key: %lu, new_key_str: %s\n",
                                             del_mn.prefixes[itr_p].size(), itr_p, del_key, mn_str1, search_pred_param.search_mn, mn_str2);
#endif

#endif

#endif
                                 
                             }
                         }

                         if (!del_mn.suffixes_terminal[itr_s]) {
                             if (!self_loop_info[1]) {

                                 //determine new value terminal
                                 bool new_snode_type=false;
                                 if (del_mn.prefixes_terminal[itr_p])
                                     new_snode_type=true;
                                 else{
                                     if(self_loop_info[0])
                                         new_snode_type=true;
                                     else
                                         new_snode_type=false;
                                 }
                                 /*
                                 ModNodeInfo send_succ_data;
                                 send_succ_data.modifyMN = search_succ_param;
                                 send_succ_data.modifyMN_info = MNnodeTuple{del_mn.prefixes[itr_p],
                                     std::make_pair(std::min(del_mn.prefix_count[itr_p].first,
                                                 del_mn.suffix_count[itr_s].first), count),
                                     new_snode_type};
                                 send_succ_data.direction=(KdirType)S;
                                 */
#ifdef OPT_ITR
                                 int node_id = retrieve_proc_id(search_succ_param.search_mn);
                                 assert(proc_idx[node_id] < mn_nodes_per_proc[node_id].size());
                                 mn_nodes_per_proc[node_id][proc_idx[node_id]] = send_succ_data;
                                 proc_idx[node_id]++;
#else
                                 //mn_nodes_per_proc[retrieve_proc_id(search_succ_param.search_mn)].push_back(send_succ_data);
                                 mn_nodes_per_proc[retrieve_proc_id(search_succ_param.search_mn)].push_back(TransferNode{
                                    search_succ_param.search_mn, search_succ_param.search_ext, del_mn.prefixes[itr_p],
                                    std::make_pair(std::min(del_mn.prefix_count[itr_p].first, del_mn.suffix_count[itr_s].first), count),
                                    new_snode_type, (KdirType)S}
                                 );
#ifdef DEBUG_KSIZE

#ifdef EXTEND_KMER
                                char mn_succ1[MN_LENGTH+1], mn_succ2[MN_LENGTH+1];
                                 for (int k=0; k<MN_LENGTH; k++) {
                                 mn_succ1[k] = el_to_char(kmerel_mn(del_key, k));
                                 mn_succ2[k] = el_to_char(kmerel_mn(search_succ_param.search_mn, k));
                                 }
                                 mn_succ1[MN_LENGTH] = '\0';
                                 mn_succ2[MN_LENGTH] = '\0';
                                 fprintf(f, "Succ: key:%s, new_key: %s, dest rank: %d\n",mn_succ1, mn_succ2, retrieve_proc_id(search_succ_param.search_mn));
#else
                                 fprintf(f, "Succ: key:%lu, new_key: %lu\n",del_key, search_succ_param.search_mn);
#endif

#endif
#endif
                             }
                         }

                     } // end of else

                 }// end of for loop across num_wires
             }// end of if condition

         } // end of for loop for a node
         //MN_map.erase(MN_map.begin()+id_list_nodes[i].terminal_prefix_id);
         //MN_map.erase(std::remove(MN_map.begin(), MN_map.end(), del_key), MN_map.end());

                     //[](auto& elem){ return elem.first == del_key;} ), 
     }// end of for for list of nodes in id_set

#ifdef DEBUG_KSIZE
     fclose(f);
#endif

#ifdef OPT_ITR
     proc_idx.clear();
#endif

}


void check_for_self_loops (kmer_t search_pkmer, kmer_t search_skmer,
        kmer_t node, std::array<bool, 2>& self_loop_info)
{
    if (search_pkmer == node)
        self_loop_info[0]=true;

    if (search_skmer == node)
        self_loop_info[1]=true;

}

/*
//std::vector<ModNodeInfo> serialize_and_transfer 
void serialize_and_transfer_by_cereal 
                         (std::vector< std::vector<ModNodeInfo> >& mn_nodes_per_proc,
                          std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                          std::vector<size_t> &rewire_pos_list,
                          int num_itr)
{

#ifdef DEBUG_P2TIME
     double st_time=0.0, global_st_time=0.0;
#endif

     // Perform Alltoallv to update the nodes 
     std::string send_buffer;
     std::string tmp_buffer;
     // send and recv buffers for obtaining the actual number of macro_nodes 
     std::vector<uint64_t> send_count_buf(size,0);
     std::vector<uint64_t> recv_count_buf(size,0);
#ifdef BIGMPI
     // send and recv buffers for obtaining the serialized data 
     std::vector<MPI_Count> scounts(size,0); //sending serialized data in bytes
     std::vector<MPI_Count> rcounts (size,0);
     std::vector<MPI_Aint> sdisp (size,0);
     std::vector<MPI_Aint> rdisp (size,0);
#else 
     // send and recv buffers for obtaining the serialized data 
     std::vector<uint64_t> scounts(size,0); //sending serialized data in bytes
     std::vector<uint64_t> rcounts (size,0);
     std::vector<uint64_t> rdisp (size,0);
     std::vector<int> scounts_dd(size,0);
     std::vector<int> rcounts_dd(size,0);
     std::vector<int> sdisp_dd(size,0);
     std::vector<int> rdisp_dd(size,0);
#endif

#ifdef DEBUG_KSIZE
    char proc_id[3];
    char itr_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    sprintf(itr_id, "%d", num_itr); 
    strcpy(output_file_name, "pre-ser_");
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)], "_itr");
    strcpy(&output_file_name[strlen(output_file_name)],itr_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (int i=0; i<size; i++)
    {
     for (size_t j=0; j<mn_nodes_per_proc[i].size(); j++) {
          kmer_t key = mn_nodes_per_proc[i][j].modifyMN.search_mn;
          BasePairVector ext = mn_nodes_per_proc[i][j].modifyMN.search_ext;
          char mn_str[MN_LENGTH+1];
          for (int k=0; k<MN_LENGTH; k++) {
              mn_str[k] = el_to_char(kmerel_mn(key, k));
          }
         mn_str[MN_LENGTH] = '\0';
 
         std::string test_ext;
         for (size_t l=0; l<ext.size(); l++)
              test_ext += el_to_char(ext[l]);

         BasePairVector tup_ext = mn_nodes_per_proc[i][j].modifyMN_info.mn_ext;
         std::string test_tup_ext;
         for (size_t l=0; l<tup_ext.size(); l++)
              test_tup_ext += el_to_char(tup_ext[l]);

#ifdef EXTEND_KMER
        fprintf(f, "key_str: %s, ext: %s, tup_ext: %s, tup_dir: %d\n",
                    mn_str, test_ext.c_str(), test_tup_ext.c_str(), mn_nodes_per_proc[i][j].direction);
#else
         fprintf(f, "key:%lu, key_str: %s, ext: %s, tup_ext: %s, tup_dir: %d\n",
                     key, mn_str, test_ext.c_str(), test_tup_ext.c_str(), mn_nodes_per_proc[i][j].direction);
#endif
     }
    }
   
    fclose(f);    
#endif

     std::string recv_buffer;
     char* recv_ptr=nullptr;

     std::ostringstream os_cnt(std::ios::binary | std::ios::out | std::ios::in);

     //serialize
     uint64_t ssize=0,rsize=0,r_mnodes=0;
     uint64_t pad_width=1000;

     {
         cereal::BinaryOutputArchive oarchive(os_cnt);

#ifdef ASSERT_CHECK
         size_t cum_node_sizes=0;
#endif

         for (int i=0; i<size; i++)
         {
             os_cnt.str("");
             tmp_buffer.clear();
#ifdef ASSERT_CHECK
             cum_node_sizes=0;
#endif

             send_count_buf[i] = mn_nodes_per_proc[i].size();
             for (size_t j=0; j<mn_nodes_per_proc[i].size(); j++) {
                  oarchive(mn_nodes_per_proc[i][j]);
#ifdef ASSERT_CHECK
                  cum_node_sizes += ser_size(mn_nodes_per_proc[i][j]);
#endif
             }
            
             tmp_buffer = os_cnt.str();
             size_t numPaddingElements=0;
	     if (tmp_buffer.size() > 0)
	     {
             size_t nonPaddedSize = tmp_buffer.size();
             size_t pad_size = ceil((nonPaddedSize/pad_width)+0.5)*pad_width;
             numPaddingElements = (pad_size - nonPaddedSize % pad_size) % pad_size;
             //printf("nonPaddedSize: %lu, pad_size: %lu, numPaddingElements: %lu\n", nonPaddedSize, pad_size, numPaddingElements);

             if(numPaddingElements > 0)
                tmp_buffer.resize(nonPaddedSize + numPaddingElements, 0);

             assert(pad_size == tmp_buffer.size());
	     }

             scounts[i] = tmp_buffer.length();
             scounts_dd[i] = scounts[i]/pad_width;
 
             send_buffer.append(tmp_buffer);

#ifdef ASSERT_CHECK
             assert (scounts[i] == (cum_node_sizes+numPaddingElements));
#endif

         }

         
     }

     // clear the buffers
     tmp_buffer.clear();
     //node_offset.clear();
     for (int i=0; i<size; i++)
          mn_nodes_per_proc[i].clear();
     mn_nodes_per_proc.clear();

     for (int t=0; t<size; t++) ssize += scounts[t];
     assert (ssize == send_buffer.length());

     MPI_Barrier(MPI_COMM_WORLD);

     double comm1 = MPI_Wtime ();
#ifdef BIGMPI
     //Sending the counts
     MPIX_Alltoall_x(send_count_buf.data(), 1, MPI_UINT64_T, recv_count_buf.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);

     //Sending the serialized data buffer
     MPIX_Alltoall_x(scounts.data(), 1, MPI_UINT64_T, rcounts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
#else
     //Sending the counts
     MPI_Alltoall(send_count_buf.data(), 1, MPI_UINT64_T, recv_count_buf.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);

     //Sending the serialized data buffer
     MPI_Alltoall(scounts.data(), 1, MPI_UINT64_T, rcounts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
#endif
     double comm2 = MPI_Wtime ();
     p2alltoall_time += (comm2 - comm1);
#ifdef DEBUG_P2TIME
     st_time += (comm2 - comm1);
#endif

     // r_mnodes denotes the number of MN entries to deserialize and modify per rank
     for (int t=0; t<size; t++) r_mnodes += recv_count_buf[t];

     for (int t=0; t<size; t++) rsize += rcounts[t];

     for (int t=0; t<size; t++) {
          sdisp_dd[t] = (t>0) ? (scounts_dd[t-1] + sdisp_dd[t-1]) : 0;
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

#ifdef DEBUG_CEREAL
    char proc_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    strcpy(output_file_name,"debug_cereal_");
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    fprintf(f,"%s\n", send_buffer.c_str());
 
#endif

#ifdef CHECK_ALLTOALLV
     uint64_t global_sfound_cnt=0, sfound_cnt=0;
    
     std::size_t sfound = send_buffer.find("F");
     if (sfound!=std::string::npos)
         sfound_cnt = std::count(send_buffer.begin(), send_buffer.end(), 'F');

     MPI_Allreduce(&sfound_cnt, &global_sfound_cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
         //fprintf (stderr,"rank: %d, Error in alltoallv recv_buffer, F was found at: %lu\n", rank, found);
#endif 
    
     MPI_Barrier(MPI_COMM_WORLD);

     double comm3 = MPI_Wtime ();
#ifdef BIGMPI
     int result = MPIX_Alltoallv_x(send_buffer.c_str(), scounts.data(), sdisp.data(), MPI_BYTE,
             recv_ptr, rcounts.data(), rdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
#else
     //int result = MPI_Alltoallv(send_buffer.c_str(), scounts.data(), sdisp.data(), MPI_BYTE,
     //        recv_ptr, rcounts.data(), rdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
     int result = MPI_Alltoallv(send_buffer.c_str(), scounts_dd.data(), sdisp_dd.data(), rowtype,
             recv_ptr, rcounts_dd.data(), rdisp_dd.data(), rowtype, MPI_COMM_WORLD);
#endif

     if (result != MPI_SUCCESS) {
         printf("rank: %d, MPI_Alltoallv in Phase 2 failed with return value: %d\n", rank, result);
         MPI_Finalize();
         exit(2);
     }

     double comm4 = MPI_Wtime ();
     p2alltoallv_time += (comm4 - comm3);
#ifdef DEBUG_P2TIME
     st_time += (comm4 - comm3);
#endif

     // free datatype
     MPI_Type_free(&rowtype);

#ifdef DEBUG_CEREAL
    
    fprintf(f,"%s\n", recv_buffer.c_str());

    fclose(f);
 
#endif
     
#ifdef CHECK_ALLTOALLV
     uint64_t global_rfound_cnt=0, rfound_cnt=0;

     std::size_t rfound = recv_buffer.find("F");
     if (rfound!=std::string::npos)
         rfound_cnt = std::count(recv_buffer.begin(), recv_buffer.end(), 'F');
         //fprintf (stderr,"rank: %d, Error in alltoallv recv_buffer, F was found at: %lu\n", rank, found);

     MPI_Allreduce(&rfound_cnt, &global_rfound_cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

     if (global_rfound_cnt != global_sfound_cnt)
         fprintf (stderr,"rank: %d, Error in alltoallv, F counts dont match! global_sfound_cnt: %lu, global_rfound_cnt: %lu\n",
                  rank, global_sfound_cnt, global_rfound_cnt);

     assert(global_rfound_cnt == global_sfound_cnt);

#endif

     os_cnt.str("");
     scounts.clear();
     scounts_dd.clear();
     sdisp_dd.clear();
     send_buffer.clear();
     send_count_buf.clear();
     rcounts_dd.clear();
     rdisp_dd.clear();

#ifdef DESER_V1
     //deserialize
     std::vector<ModNodeInfo> mn_nodes_to_modify (r_mnodes);

     {
         std::istringstream is(recv_buffer, std::ios::binary);
         cereal::BinaryInputArchive iarchive(is);
         size_t k=0; 

         try {
             while(r_mnodes)
             {
                 iarchive(mn_nodes_to_modify[k]);
                 k++;
                 r_mnodes--;
             }
         }
         catch (cereal::Exception& e) {
             std::cout << e.what() << std::endl;
         }

         assert(k==mn_nodes_to_modify.size());
     }

     MPI_Barrier(MPI_COMM_WORLD);

     recv_buffer.clear();

#endif // end of DESER_V1
     
     {
         assert(recv_buffer.size() == rsize);
         std::istringstream is(std::ios::binary | std::ios::out | std::ios::in);
         //is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()), rsize);
         //std::istringstream is(recv_buffer, std::ios::binary);
         //cereal::BinaryInputArchive iarchive(is);
         uint64_t k=0;

         for (int i=0; i<size; i++)
         {
             is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()+rdisp[i]), rcounts[i]);
             cereal::BinaryInputArchive iarchive(is);
             uint64_t nnodes = recv_count_buf[i];

             try {
                 while(nnodes)
                 {
                     ModNodeInfo mi;
                     iarchive(mi);

                     kmer_t search_key = mi.modifyMN.search_mn;
                     int pos = find_mnode_exists(search_key, MN_map);
                     if (MN_map[pos].first != search_key) {
                         printf("rank: %d, MN node key: %lu, was not found in map at the time of pushing to: %d\n", 
                             rank, search_key, mi.direction);
                         MPI_Finalize();
                         exit(2);
                     }
                     assert(MN_map[pos].first == search_key);

                     if (mi.direction==P)
                         push_to_pred(mi, pos, MN_map);
                     else
                         push_to_succ(mi, pos, MN_map);

                     rewire_pos_list.push_back(pos); 

                     k++;
                     nnodes--;
                 }
             }
             catch (cereal::Exception& e) {
                 std::cout << e.what() << rank << std::endl;
                 //std::cout << e.what() << std::endl;
             }
         }

         if (k!=r_mnodes)
             fprintf(stderr, "Error!! for rank: %d, Not matching, k:%lu, r_mnodes:%lu, ssize:%lu, rsize:%lu\n", 
                                               rank, k, r_mnodes, ssize, rsize);

         assert(k==r_mnodes);
     }
             
     MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG_P2TIME
     MPI_Reduce(&st_time, &global_st_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     if (rank == 0) printf ("part3: serialize and update communication time  (secs): %f \n\n",
                            (double)global_st_time/(double)size);
#endif

     recv_buffer.clear();
     rcounts.clear();
     rdisp.clear();
     recv_count_buf.clear();
     
#ifdef DESER_V1
     return mn_nodes_to_modify;
#endif

}
*/


void serialize_and_transfer 
                         (std::vector< std::vector<TransferNode> >& mn_nodes_per_proc,
                          std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                          std::vector<size_t> &rewire_pos_list,
                          int num_itr)
{

#ifdef DEBUG_P2TIME
     double st_time=0.0, global_st_time=0.0;
#endif

     /* Perform Alltoallv to update the nodes */
     std::string send_buffer;
     std::string tmp_buffer;
     /* send and recv buffers for obtaining the actual number of macro_nodes */
     std::vector<uint64_t> send_count_buf(size,0);
     std::vector<uint64_t> recv_count_buf(size,0);
#ifdef BIGMPI
     /* send and recv buffers for obtaining the serialized data */
     std::vector<MPI_Count> scounts(size,0); //sending serialized data in bytes
     std::vector<MPI_Count> rcounts (size,0);
     std::vector<MPI_Aint> sdisp (size,0);
     std::vector<MPI_Aint> rdisp (size,0);
#else 
     /* send and recv buffers for obtaining the serialized data */
     std::vector<uint64_t> scounts(size,0); //sending serialized data in bytes
     std::vector<uint64_t> rcounts (size,0);
     std::vector<uint64_t> rdisp (size,0);
     std::vector<int> scounts_dd(size,0);
     std::vector<int> rcounts_dd(size,0);
     std::vector<int> sdisp_dd(size,0);
     std::vector<int> rdisp_dd(size,0);
#endif

#ifdef DEBUG_KSIZE
    char proc_id[3];
    char itr_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    sprintf(itr_id, "%d", num_itr); 
    strcpy(output_file_name, "pre-ser_");
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)], "_itr");
    strcpy(&output_file_name[strlen(output_file_name)],itr_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (int i=0; i<size; i++)
    {
     for (size_t j=0; j<mn_nodes_per_proc[i].size(); j++) {
          kmer_t key = mn_nodes_per_proc[i][j].search_mn;
          BasePairVector ext = mn_nodes_per_proc[i][j].search_ext;
          char mn_str[MN_LENGTH+1];
          for (int k=0; k<MN_LENGTH; k++) {
              mn_str[k] = el_to_char(kmerel_mn(key, k));
          }
         mn_str[MN_LENGTH] = '\0';
 
         std::string test_ext;
         for (size_t l=0; l<ext.size(); l++)
              test_ext += el_to_char(ext[l]);

         BasePairVector tup_ext = mn_nodes_per_proc[i][j].mn_ext;
         std::string test_tup_ext;
         for (size_t l=0; l<tup_ext.size(); l++)
              test_tup_ext += el_to_char(tup_ext[l]);


         fprintf(f, "i: %d, key_str: %s, ext: %s, tup_ext: %s, tup_dir: %d\n",
                     i, mn_str, test_ext.c_str(), test_tup_ext.c_str(), mn_nodes_per_proc[i][j].direction);
     }
    }
   
    fclose(f);    
#endif

     std::string recv_buffer;
     char* recv_ptr=nullptr;

     std::ostringstream os_cnt(std::ios::binary | std::ios::out | std::ios::in);

     //serialize
     uint64_t ssize=0,rsize=0,r_mnodes=0;
     uint64_t pad_width=1000;

/*
#ifdef ASSERT_CHECK
         size_t cum_node_sizes=0;
#endif
*/
         for (int i=0; i<size; i++)
         {
             os_cnt.str("");
             tmp_buffer.clear();
/*
#ifdef ASSERT_CHECK
             cum_node_sizes=0;
#endif
*/

             send_count_buf[i] = mn_nodes_per_proc[i].size();
             for (size_t j=0; j<mn_nodes_per_proc[i].size(); j++) {
                  serialize(os_cnt, mn_nodes_per_proc[i][j]);
/*
#ifdef ASSERT_CHECK
                  cum_node_sizes += ser_size(mn_nodes_per_proc[i][j]);
#endif
*/
             }
/*
#ifdef ASSERT_CHECK
             cum_node_sizes=0;
#endif
             for (size_t j=node_offset[i]; j<(node_offset[i]+send_count_buf[i]); j++)
             {
#ifdef ASSERT_CHECK
                 cum_node_sizes += ser_size(mn_nodes_per_proc_1d[j]);
#endif
                 oarchive(mn_nodes_per_proc_1d[j]);
             }
*/             
             tmp_buffer = os_cnt.str();
	     if (tmp_buffer.size() > 0)
	     {
             size_t nonPaddedSize = tmp_buffer.size();
             size_t padded_size = ((nonPaddedSize/pad_width) + (nonPaddedSize%pad_width!=0))*pad_width;

             tmp_buffer.resize(padded_size, 0);
             //printf("rank: %d, nonPaddedSize: %lu, buffer_size: %lu\n", rank, nonPaddedSize, tmp_buffer.length());
             /*
             size_t nonPaddedSize = tmp_buffer.size();
             size_t pad_size = ceil((nonPaddedSize/pad_width)+0.5)*pad_width;
             const size_t numPaddingElements = (pad_size - nonPaddedSize % pad_size) % pad_size;
             //printf("nonPaddedSize: %lu, pad_size: %lu, numPaddingElements: %lu\n", nonPaddedSize, pad_size, numPaddingElements);

             if(numPaddingElements > 0)
                tmp_buffer.resize(nonPaddedSize + numPaddingElements, 0);

             assert(pad_size == tmp_buffer.size());
             */
	     }

             scounts[i] = tmp_buffer.length();
             scounts_dd[i] = scounts[i]/pad_width;
 
             send_buffer.append(tmp_buffer);

//#ifdef ASSERT_CHECK
//             assert (scounts[i] == (cum_node_sizes+numPaddingElements));
//#endif


         /*
         cereal::BinaryOutputArchive oarchive2(os_buf);

         for (int i=0; i<size; i++) {
             for (size_t j=0; j<mn_nodes_per_proc[i].size(); j++)
                  oarchive2(mn_nodes_per_proc[i][j]);
         }

         send_buffer = os_buf.str();
         */
     }

     // clear the buffers
     tmp_buffer.clear();
     //node_offset.clear();
     for (int i=0; i<size; i++)
          mn_nodes_per_proc[i].clear();
     mn_nodes_per_proc.clear();

     for (int t=0; t<size; t++) ssize += scounts[t];
     assert (ssize == send_buffer.length());

     MPI_Barrier(MPI_COMM_WORLD);

     double comm1 = MPI_Wtime ();
#ifdef BIGMPI
     //Sending the counts
     MPIX_Alltoall_x(send_count_buf.data(), 1, MPI_UINT64_T, recv_count_buf.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);

     //Sending the serialized data buffer
     MPIX_Alltoall_x(scounts.data(), 1, MPI_UINT64_T, rcounts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
#else
     //Sending the counts
     MPI_Alltoall(send_count_buf.data(), 1, MPI_UINT64_T, recv_count_buf.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);

     //Sending the serialized data buffer
     MPI_Alltoall(scounts.data(), 1, MPI_UINT64_T, rcounts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
#endif
     double comm2 = MPI_Wtime ();
     p2alltoall_time += (comm2 - comm1);
#ifdef DEBUG_P2TIME
     st_time += (comm2 - comm1);
#endif

     /* r_mnodes denotes the number of MN entries to deserialize and modify per rank */
     for (int t=0; t<size; t++) r_mnodes += recv_count_buf[t];

     for (int t=0; t<size; t++) rsize += rcounts[t];

     for (int t=0; t<size; t++) {
          sdisp_dd[t] = (t>0) ? (scounts_dd[t-1] + sdisp_dd[t-1]) : 0;
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


#ifdef CHECK_ALLTOALLV
     uint64_t global_sfound_cnt=0, sfound_cnt=0;
    
     std::size_t sfound = send_buffer.find("F");
     if (sfound!=std::string::npos)
         sfound_cnt = std::count(send_buffer.begin(), send_buffer.end(), 'F');

     MPI_Allreduce(&sfound_cnt, &global_sfound_cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
         //fprintf (stderr,"rank: %d, Error in alltoallv recv_buffer, F was found at: %lu\n", rank, found);
#endif 
    
     MPI_Barrier(MPI_COMM_WORLD);

     double comm3 = MPI_Wtime ();
#ifdef BIGMPI
     int result = MPIX_Alltoallv_x(send_buffer.c_str(), scounts.data(), sdisp.data(), MPI_BYTE,
             recv_ptr, rcounts.data(), rdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
#else
     //int result = MPI_Alltoallv(send_buffer.c_str(), scounts.data(), sdisp.data(), MPI_BYTE,
     //        recv_ptr, rcounts.data(), rdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
     int result = MPI_Alltoallv(send_buffer.c_str(), scounts_dd.data(), sdisp_dd.data(), rowtype,
             recv_ptr, rcounts_dd.data(), rdisp_dd.data(), rowtype, MPI_COMM_WORLD);
#endif

     if (result != MPI_SUCCESS) {
         printf("rank: %d, MPI_Alltoallv in Phase 2 failed with return value: %d\n", rank, result);
         MPI_Finalize();
         exit(2);
     }

     double comm4 = MPI_Wtime ();
     p2alltoallv_time += (comm4 - comm3);
#ifdef DEBUG_P2TIME
     st_time += (comm4 - comm3);
#endif

     // free datatype
     MPI_Type_free(&rowtype);

#ifdef CHECK_ALLTOALLV
     uint64_t global_rfound_cnt=0, rfound_cnt=0;

     std::size_t rfound = recv_buffer.find("F");
     if (rfound!=std::string::npos)
         rfound_cnt = std::count(recv_buffer.begin(), recv_buffer.end(), 'F');
         //fprintf (stderr,"rank: %d, Error in alltoallv recv_buffer, F was found at: %lu\n", rank, found);

     MPI_Allreduce(&rfound_cnt, &global_rfound_cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

     if (global_rfound_cnt != global_sfound_cnt)
         fprintf (stderr,"rank: %d, Error in alltoallv, F counts dont match! global_sfound_cnt: %lu, global_rfound_cnt: %lu\n",
                  rank, global_sfound_cnt, global_rfound_cnt);

     assert(global_rfound_cnt == global_sfound_cnt);

#endif

     os_cnt.str("");
     scounts.clear();
     scounts_dd.clear();
     sdisp_dd.clear();
     send_buffer.clear();
     send_count_buf.clear();
     rcounts_dd.clear();
     rdisp_dd.clear();

     
     {
         assert(recv_buffer.size() == rsize);
         std::istringstream is(std::ios::binary | std::ios::out | std::ios::in);
         //is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()), rsize);
         //std::istringstream is(recv_buffer, std::ios::binary);
         //cereal::BinaryInputArchive iarchive(is);
         uint64_t k=0;

         for (int i=0; i<size; i++)
         {
             is.rdbuf()->pubsetbuf(const_cast<char*>(recv_buffer.c_str()+rdisp[i]), rcounts[i]);
             uint64_t nnodes = recv_count_buf[i];

             //try {
                 while(nnodes)
                 {
                     TransferNode mi;
                     deserialize(is, mi);

                     kmer_t search_key = mi.search_mn;
                     int pos = find_mnode_exists(search_key, MN_map);
                     if (MN_map[pos].first != search_key) {
                         //printf("rank: %d, MN node key: %lu, was not found in map at the time of pushing to: %d\n", 
                         printf("rank: %d, MN node key was not found in map at the time of pushing to: %d\n", 
                             rank, mi.direction);
                         MPI_Finalize();
                         exit(2);
                     }
                     assert(MN_map[pos].first == search_key);

                     if (mi.direction==P)
                         push_to_pred(mi, pos, MN_map);
                     else
                         push_to_succ(mi, pos, MN_map);

                     rewire_pos_list.push_back(pos); 

                     k++;
                     nnodes--;
                 }
             //}
             //catch (cereal::Exception& e) {
             //    std::cout << e.what() << rank << std::endl;
             //}
         }

         if (k!=r_mnodes)
             fprintf(stderr, "Error!! for rank: %d, Not matching, k:%lu, r_mnodes:%lu, ssize:%lu, rsize:%lu\n", 
                                               rank, k, r_mnodes, ssize, rsize);

         assert(k==r_mnodes);
     }
             
     MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG_P2TIME
     MPI_Reduce(&st_time, &global_st_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     if (rank == 0) printf ("part3: serialize and update communication time  (secs): %f \n\n",
                            (double)global_st_time/(double)size);
#endif

     recv_buffer.clear();
     rcounts.clear();
     rdisp.clear();
     recv_count_buf.clear();
     
#ifdef DESER_V1
     return mn_nodes_to_modify;
#endif

}


//void push_to_pred (std::vector<ModNodeInfo> &new_mnode, int pos, 
void push_to_pred (TransferNode &new_mnode, int pos, 
        std::vector<std::pair<kmer_t,MacroNode>> &MN_map)
{

    MacroNode &mn_found = MN_map[pos].second;
    assert(mn_found.suffixes.size() == mn_found.suffix_count.size());
    assert(mn_found.suffixes.size() == mn_found.suffixes_terminal.size());
    
    //BasePairVector search_ext = new_mnode.modifyMN.search_ext;
    BasePairVector search_ext = new_mnode.search_ext;
    BasePairVector new_pext = search_ext;
    //new_pext.append(new_mnode.modifyMN_info.mn_ext);
    new_pext.append(new_mnode.mn_ext);
    std::pair<int,int> new_kcount = new_mnode.mn_count;
    bool new_ntype = new_mnode.mn_terminal;
    //bool new_ntype =  new_mnode.mn_terminal==Y ? true : false;

    std::vector<BasePairVector>::iterator viter =
        std::find(mn_found.suffixes.begin(), mn_found.suffixes.end(), BasePairVector(search_ext));
    if (viter == mn_found.suffixes.end()) // not found
    {
        mn_found.suffixes.push_back(new_pext);
        mn_found.suffix_count.push_back(new_kcount);
        mn_found.suffixes_terminal.push_back(new_ntype);
    }
    else // found 
    {
        int idx = std::distance(mn_found.suffixes.begin(), viter);
        if (mn_found.suffixes_terminal[idx]) {
            mn_found.suffixes.push_back(new_pext);
            mn_found.suffix_count.push_back(new_kcount);
            mn_found.suffixes_terminal.push_back(new_ntype);
        }
        else
        {
            mn_found.suffixes[idx]=new_pext;
            mn_found.suffix_count[idx]=new_kcount;
            mn_found.suffixes_terminal[idx]=new_ntype;
        }
    }

//#ifdef DEBUG_SERMOD
    //mn_found.setup_wiring();

}


void push_to_succ (TransferNode &new_mnode, int pos, 
        std::vector<std::pair<kmer_t,MacroNode>> &MN_map)
{

    MacroNode &mn_found = MN_map[pos].second;
    assert(mn_found.prefixes.size() == mn_found.prefix_count.size());
    assert(mn_found.prefixes.size() == mn_found.prefixes_terminal.size());
    
    //BasePairVector search_ext = new_mnode.modifyMN.search_ext;
    BasePairVector search_ext = new_mnode.search_ext;
    //BasePairVector new_sext = new_mnode.modifyMN_info.mn_ext;
    BasePairVector new_sext = new_mnode.mn_ext;
    new_sext.append(search_ext);
    std::pair<int,int> new_kcount = new_mnode.mn_count;
    bool new_ntype = new_mnode.mn_terminal;
    //bool new_ntype =  new_mnode.mn_terminal==Y ? true : false; 

    std::vector<BasePairVector>::iterator viter =
        std::find(mn_found.prefixes.begin(), mn_found.prefixes.end(), BasePairVector(search_ext));
    if (viter == mn_found.prefixes.end()) // not found
    {
        mn_found.prefixes.push_back(new_sext);
        mn_found.prefix_count.push_back(new_kcount);
        mn_found.prefixes_terminal.push_back(new_ntype);
    }
    else // found 
    {
        int idx = std::distance(mn_found.prefixes.begin(), viter);
        if (mn_found.prefixes_terminal[idx]) {
            mn_found.prefixes.push_back(new_sext);
            mn_found.prefix_count.push_back(new_kcount);
            mn_found.prefixes_terminal.push_back(new_ntype);
        }
        else
        {
            mn_found.prefixes[idx]=new_sext;
            mn_found.prefix_count[idx]=new_kcount;
            mn_found.prefixes_terminal[idx]=new_ntype;
        }
    }

    //mn_found.setup_wiring();

}

