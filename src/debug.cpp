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

extern long int MAX_KMER_COUNT;
extern int rank, size;
extern int coverage;
extern int num_buckets;
extern int num_threads;



void print_id_list (std::vector<size_t>& id_list_nodes , const std::string& func, std::vector<std::pair<kmer_t,MacroNode>>& MN_map)
{
    char proc_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    strcpy(output_file_name, func.c_str());
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (size_t i=0; i<id_list_nodes.size(); i++) {
         kmer_t &del_key = MN_map[id_list_nodes[i]].first;
         char mn_str[MN_LENGTH+1];
         for (int k=0; k<MN_LENGTH; k++) {
              mn_str[k] = el_to_char(kmerel_mn(del_key, k));
         }
         mn_str[MN_LENGTH] = '\0';

         fprintf(f, "key_str: %s\n", mn_str);
         //fprintf(f, "%lu, key:%lu, key_str: %s\n", id_list_nodes[i], del_key, mn_str);
    }

    fclose(f);

}

void print_map_content (std::vector<std::pair<kmer_t,MacroNode>>& MN_map, const std::string& func)
{
    char proc_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    //strcpy(output_file_name,"wired_mn_out_");
    strcpy(output_file_name, func.c_str());
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (std::vector<std::pair<kmer_t,MacroNode>>::const_iterator it = MN_map.begin();
         it != MN_map.end(); it++) {
        
         kmer_t kmer_1 = it->first; 
         MacroNode mn = it->second;
         BasePairVector key = mn.k_1_mer;
         std::string test_data;

         for (size_t l=0; l<key.size(); l++)
              test_data += el_to_char(key[l]);
         fprintf(f,"key_id: %lu, key: %s, num_prefix: %lu, num_suffix: %lu\n", kmer_1, 
                                 test_data.c_str(), mn.prefixes.size(), mn.suffixes.size());
         assert(mn.prefixes.size() == mn.prefix_count.size());
         assert(mn.suffixes.size() == mn.suffix_count.size());
         for (size_t l=0; l<mn.prefixes.size(); l++) {
              BasePairVector x = mn.prefixes[l];
              fprintf(f,"x size: %lu, count: (%d,%d), type: %d\n", x.size(), mn.prefix_count[l].first, mn.prefix_count[l].second, static_cast<int>(mn.prefixes_terminal[l]));
              for (size_t t=0; t<x.size(); t++)
                   fprintf(f,"%c", el_to_char(x[t]));
              fprintf(f,"\n");
         }
         for (size_t l=0; l<mn.suffixes.size(); l++) {
              BasePairVector x = mn.suffixes[l];
              fprintf(f,"x size: %lu, count: (%d,%d), type: %d\n", x.size(), mn.suffix_count[l].first, mn.suffix_count[l].second, static_cast<int>(mn.suffixes_terminal[l]));
              for (size_t t=0; t<x.size(); t++)
                   fprintf(f,"%c", el_to_char(x[t]));
              fprintf(f,"\n");
         }
         for(size_t i=0; i<mn.wiring_info.size(); i++) {
                 fprintf(f,"wire_sid: %d, wire_off: %d, wire_count: %d\n", mn.wiring_info[i].suffix_id, 
                    mn.wiring_info[i].offset_in_suffix, mn.wiring_info[i].count);
             }
         for(size_t i=0; i<mn.prefix_begin_info.size(); i++) {
                 fprintf(f,"pid: %lu, pos_in_wire: %d, num_wires: %d\n", i, mn.prefix_begin_info[i].prefix_pos, 
                    mn.prefix_begin_info[i].num_wires);
             }
    }
    fclose(f);

}
void print_map_content_itr (std::vector<std::pair<kmer_t,MacroNode>>& MN_map, const std::string& func, int num_itr)
{
    char proc_id[3];
    char itr_id[3];
    char output_file_name[30];

    sprintf(proc_id, "%d", rank); 
    sprintf(itr_id, "%d", num_itr); 
    //strcpy(output_file_name,"wired_mn_out_");
    strcpy(output_file_name, func.c_str());
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],"_itr");
    strcpy(&output_file_name[strlen(output_file_name)],itr_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    size_t itr=0;
    for (std::vector<std::pair<kmer_t,MacroNode>>::const_iterator it = MN_map.begin();
         it != MN_map.end(); it++) {
        
         kmer_t kmer_1 = it->first; 
         MacroNode mn = it->second;
         BasePairVector key = mn.k_1_mer;
         std::string test_data;

         for (size_t l=0; l<key.size(); l++)
              test_data += el_to_char(key[l]);

#ifdef EXTEND_KMER
          fprintf(f,"itr: %lu, key: %s, num_prefix: %lu, num_suffix: %lu\n",
                                 itr, test_data.c_str(), mn.prefixes.size(), mn.suffixes.size());
#else
         fprintf(f,"key_id: %lu, key: %s, num_prefix: %lu, num_suffix: %lu\n", kmer_1, 
                                 test_data.c_str(), mn.prefixes.size(), mn.suffixes.size());
#endif

         assert(mn.prefixes.size() == mn.prefix_count.size());
         assert(mn.suffixes.size() == mn.suffix_count.size());
         for (size_t l=0; l<mn.prefixes.size(); l++) {
              BasePairVector x = mn.prefixes[l];
              fprintf(f,"x size: %lu, count: (%d,%d), type: %d\n", x.size(), mn.prefix_count[l].first, mn.prefix_count[l].second, static_cast<int>(mn.prefixes_terminal[l]));
              for (size_t t=0; t<x.size(); t++)
                   fprintf(f,"%c", el_to_char(x[t]));
              fprintf(f,"\n");
         }
         for (size_t l=0; l<mn.suffixes.size(); l++) {
              BasePairVector x = mn.suffixes[l];
              fprintf(f,"x size: %lu, count: (%d,%d), type: %d\n", x.size(), mn.suffix_count[l].first, mn.suffix_count[l].second, static_cast<int>(mn.suffixes_terminal[l]));
              for (size_t t=0; t<x.size(); t++)
                   fprintf(f,"%c", el_to_char(x[t]));
              fprintf(f,"\n");
         }
         for(size_t i=0; i<mn.wiring_info.size(); i++) {
                 fprintf(f,"wire_sid: %d, wire_off: %d, wire_count: %d\n", mn.wiring_info[i].suffix_id, 
                    mn.wiring_info[i].offset_in_suffix, mn.wiring_info[i].count);
             }
         for(size_t i=0; i<mn.prefix_begin_info.size(); i++) {
                 fprintf(f,"pid: %lu, pos_in_wire: %d, num_wires: %d\n", i, mn.prefix_begin_info[i].prefix_pos, 
                    mn.prefix_begin_info[i].num_wires);
             }
         itr++;
    }
    fclose(f);

}
void print_rewire_list (std::vector<size_t>& rewire_pos_list, 
                        std::vector<std::pair<kmer_t,
                        MacroNode>>& MN_map, 
                        const std::string& func, int num_itr)
{
    char proc_id[3];
    char itr_id[3];
    char output_file_name[30];

    sprintf(proc_id, "%d", rank); 
    sprintf(itr_id, "%d", num_itr); 
    //strcpy(output_file_name,"wired_mn_out_");
    strcpy(output_file_name, func.c_str());
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],"_itr");
    strcpy(&output_file_name[strlen(output_file_name)],itr_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (size_t i=0; i<rewire_pos_list.size(); i++) {
         MacroNode &mn = MN_map[rewire_pos_list[i]].second;
         kmer_t kmer_1 = MN_map[rewire_pos_list[i]].first;

         BasePairVector key = mn.k_1_mer;
         std::string test_data;

         for (size_t l=0; l<key.size(); l++)
              test_data += el_to_char(key[l]);
         fprintf(f,"key: %s, num_prefix: %lu, num_suffix: %lu\n",
                                 test_data.c_str(), mn.prefixes.size(), mn.suffixes.size());
         assert(mn.prefixes.size() == mn.prefix_count.size());
         assert(mn.suffixes.size() == mn.suffix_count.size());
         for (size_t l=0; l<mn.prefixes.size(); l++) {
              BasePairVector x = mn.prefixes[l];
              fprintf(f,"x size: %lu, count: (%d,%d), type: %d\n", x.size(), mn.prefix_count[l].first, mn.prefix_count[l].second, static_cast<int>(mn.prefixes_terminal[l]));
              for (size_t t=0; t<x.size(); t++)
                   fprintf(f,"%c", el_to_char(x[t]));
              fprintf(f,"\n");
         }
         fprintf(f,"----\n");
         for (size_t l=0; l<mn.suffixes.size(); l++) {
              BasePairVector x = mn.suffixes[l];
              fprintf(f,"x size: %lu, count: (%d,%d), type: %d\n", x.size(), mn.suffix_count[l].first, mn.suffix_count[l].second, static_cast<int>(mn.suffixes_terminal[l]));
              for (size_t t=0; t<x.size(); t++)
                   fprintf(f,"%c", el_to_char(x[t]));
              fprintf(f,"\n");
         }
         fprintf(f,"*****\n");
         for(size_t i=0; i<mn.wiring_info.size(); i++) {
                 fprintf(f,"wire_sid: %d, wire_off: %d, wire_count: %d\n", mn.wiring_info[i].suffix_id, 
                    mn.wiring_info[i].offset_in_suffix, mn.wiring_info[i].count);
             }
         for(size_t i=0; i<mn.prefix_begin_info.size(); i++) {
                 fprintf(f,"pid: %lu, pos_in_wire: %d, num_wires: %d\n", i, mn.prefix_begin_info[i].prefix_pos, 
                    mn.prefix_begin_info[i].num_wires);
             }
    }
    fclose(f);

}

void debug_wired_mnodes(std::vector<std::pair<kmer_t,MacroNode>>& MN_map)
{
    const std::string s("wired_mn_out_");
    print_map_content(MN_map, s);

}

void debug_begin_kmers_list (std::vector<std::pair<kmer_t,MacroNode>>& MN_map,
                             std::vector<BeginMN> &list_of_begin_kmers)
{

    char proc_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    strcpy(output_file_name,"bkmers_");
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    char mn_str[MN_LENGTH+1];
    for (size_t i=0; i<list_of_begin_kmers.size(); i++) {
        kmer_t key = list_of_begin_kmers[i].node;
        for (int k=0; k<MN_LENGTH; k++) {
             mn_str[k] = el_to_char(kmerel_mn(key, k));
        }
        mn_str[MN_LENGTH] = '\0';

        int terminal_prefix_id = list_of_begin_kmers[i].terminal_prefix_id;
        fprintf(f, "key_id: %lu, key: %s, t_prefix_id: %d\n", key, mn_str, terminal_prefix_id);
    }
  
    fclose(f);

}

void debug_global_pakgraph(std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map)
{

#ifdef DEBUG_MN_MAP
    
    if (rank == 0)
    {
    char proc_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    strcpy(output_file_name,"mn_out_");
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");
    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    //int num_itr=0;
    for (std::vector<std::pair<kmer_t,MacroNode>>::const_iterator it = global_MN_map.begin();
         it != global_MN_map.end(); it++) 
    {
        
         kmer_t kmer_1 = it->first; 
         MacroNode mn = it->second;
         BasePairVector key = mn.k_1_mer;
         std::string test_data;

         //if (mn.prefixes.size() > 2)
         //{
         for (size_t l=0; l<key.size(); l++)
              test_data += el_to_char(key[l]);
         fprintf(f,"key_id: %lu, key: %s, num_prefix: %lu, num_suffix: %lu\n", kmer_1, test_data.c_str(), mn.prefixes.size(), mn.suffixes.size());
         assert(mn.prefixes.size() == mn.prefix_count.size());
         assert(mn.suffixes.size() == mn.suffix_count.size());
         for (size_t l=0; l<mn.prefixes.size(); l++) {
              BasePairVector x = mn.prefixes[l];
              fprintf(f,"x size: %lu, count: (%d,%d)\n", x.size(), mn.prefix_count[l].first, mn.prefix_count[l].second);
              for (size_t t=0; t<x.size(); t++)
                   fprintf(f,"%c", el_to_char(x[t]));
              fprintf(f,"\n");
         }
         for (size_t l=0; l<mn.suffixes.size(); l++) {
              BasePairVector x = mn.suffixes[l];
              fprintf(f,"x size: %lu, count: (%d,%d)\n", x.size(), mn.suffix_count[l].first, mn.suffix_count[l].second);
              for (size_t t=0; t<x.size(); t++)
                   fprintf(f,"%c", el_to_char(x[t]));
              fprintf(f,"\n");
         }
         for(size_t i=0; i<mn.wiring_info.size(); i++) {
                 fprintf(f,"wire_sid: %d, wire_off: %d, wire_count: %d\n", mn.wiring_info[i].suffix_id, 
                    mn.wiring_info[i].offset_in_suffix, mn.wiring_info[i].count);
         }
         for(size_t i=0; i<mn.prefix_begin_info.size(); i++) {
                 fprintf(f,"prefix_id: %lu, pos_in_wire: %d, num_wires: %d\n", i, mn.prefix_begin_info[i].prefix_pos, 
                    mn.prefix_begin_info[i].num_wires);
         }
         //num_itr++;
         //if (num_itr == 20)
         //    break;
         //}
    }
    fclose(f);
    }

    MPI_Barrier(MPI_COMM_WORLD);

#endif
}
   
/*
void print_kmers ()
{

    ElType this_alpha;
    char kmer_out[KMER_LENGTH+1];
    FILE *f = fopen("kmers_rank0.log", "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    assert (kmer_proc_buf.size() == kmer_cnt_proc_buf.size());
    for (size_t it=0; it<kmer_proc_buf.size(); it++)
    {

         kmer_t kmer_1 = kmer_proc_buf[it];
         for (size_t l=0; l<KMER_LENGTH; l++) {
              this_alpha = kmerel (kmer_1, l);
              kmer_out[l] = el_to_char(this_alpha);
         }
         kmer_out[KMER_LENGTH] = '\0';
         fprintf(f, "kmer id: %lu, kmer: %s, count: %d\n", kmer_1, kmer_out, kmer_cnt_proc_buf[it]);
    }

    fclose(f);
}
*/

