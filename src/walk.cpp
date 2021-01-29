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

extern int rank, size;
extern int coverage;
extern int num_threads;
extern int num_contigs;
int num_contigs_c=0;


MnodeInfo retrieve_mn_sinfo(BasePairVector &succ_ext, BasePairVector &key)
{
    assert(key.size() == MN_LENGTH);
    BasePairVector next_ext;
    kmer_t next_mn_key=0;
    //MnodeInfo next_mn_info;

    if (succ_ext.size() >= MN_LENGTH) {
        next_mn_key = succ_ext.extract_succ(MN_LENGTH);

        uint64_t start = succ_ext.size()-MN_LENGTH;
        //for (kmer_t i=start; i<succ_ext.size(); i++)
        //     next_mn_key = mnmer_shift(next_mn_key, succ_ext[i]);

        if (start) {
            next_ext=key;
            next_ext.mnode_extract_succ(succ_ext, start);
            //for (kmer_t i=0; i<start; i++)
            //     next_ext.push_back(succ_ext[i]);
        }
        else
            next_ext=key;
    }
    else {
        uint64_t remainder = MN_LENGTH-succ_ext.size();
        next_mn_key = key.extract_succ(remainder);
        next_mn_key = ((next_mn_key << (succ_ext.size()*2)) | succ_ext.extract(succ_ext.size())) & (kmer_t)MN_MASK;

        //for (kmer_t i=succ_ext.size(); i<key.size(); i++) {
        //     next_mn_key = mnmer_shift(next_mn_key, key[i]);
        //}
        //for (kmer_t i=0; i<succ_ext.size(); i++)
        //     next_mn_key = mnmer_shift(next_mn_key, succ_ext[i]);

        next_ext.populate(key.extract(succ_ext.size()), succ_ext.size());
        //for (kmer_t i=0; i<succ_ext.size(); i++)
        //     next_ext.push_back(key[i]);
    }

    assert(next_ext.size() == succ_ext.size());
    //next_mn_info.search_mn = next_mn_key;
    //next_mn_info.search_ext = next_ext;

    //return next_mn_info;
    return MnodeInfo{next_mn_key, next_ext};

}

MnodeInfo retrieve_mn_pinfo(BasePairVector &pred_ext, BasePairVector &key)
{
    assert(key.size() == MN_LENGTH);
    BasePairVector next_ext;
    kmer_t next_mn_key=0;
    //MnodeInfo next_mn_info;

    if (pred_ext.size() >= MN_LENGTH) {
        
        next_mn_key = pred_ext.extract(MN_LENGTH); 
        //for (int i=0; i<MN_LENGTH; i++)
        //     next_mn_key = mnmer_shift(next_mn_key, pred_ext[i]);

        size_t remainder = pred_ext.size()-MN_LENGTH;
        if (remainder) {
            next_ext.mnode_extract_pred(pred_ext, remainder);
            //for(size_t i=MN_LENGTH; i<pred_ext.size(); i++)
            //    next_ext.push_back(pred_ext[i]);
            next_ext.append(key);
        }
        else
            next_ext=key;
            //next_ext.append(key);
    }
    else {
        size_t remainder = MN_LENGTH - pred_ext.size();

        next_mn_key = pred_ext.extract(pred_ext.size());
        next_mn_key = ((next_mn_key << (remainder*2)) | key.extract(remainder)) & (kmer_t)MN_MASK;

        //for (size_t i=0; i<pred_ext.size(); i++)
        //     next_mn_key = mnmer_shift(next_mn_key, pred_ext[i]);
        //for (size_t i=0; i<remainder; i++)
        //     next_mn_key = mnmer_shift(next_mn_key, key[i]);

        next_ext.populate(key.extract_succ(pred_ext.size()), pred_ext.size());        
        //for (size_t i=remainder; i<key.size(); i++)
        //    next_ext.push_back(key[i]);
    }

    assert(next_ext.size() == pred_ext.size());
    //next_mn_info.search_mn = next_mn_key;
    //next_mn_info.search_ext = next_ext;

    //return next_mn_info;
    return MnodeInfo{next_mn_key, next_ext};

}


//std::vector<std::pair<kmer_t,MacroNode>>::iterator find_mnode_exists

int  find_mnode_exists (kmer_t key, std::vector<std::pair<kmer_t,MacroNode>> &MN_map)
{
     std::vector<std::pair<kmer_t,MacroNode>>::iterator it;
     size_t pos=-1; 
   
     it = std::lower_bound(MN_map.begin(), MN_map.end(), key,
                      [] (const std::pair<kmer_t,MacroNode>& lhs, const kmer_t &rhs) {
                             return (lhs.first < rhs);
                      });

     //pos = it - MN_map.begin();
     pos = std::distance(MN_map.begin(), it);
     if (MN_map[pos].first != key) {
         char mn_str[MN_LENGTH+1];
         for (int k=0; k<MN_LENGTH; k++) {
              mn_str[k] = el_to_char(kmerel_mn(key, k));
         }
         mn_str[MN_LENGTH] = '\0';

         //printf("Error, key not found! rank: %d, pos: %d, MN_map[pos].first: %lu, key: %lu, key_str: %s, sizeof(key): %lu\n", 
         printf("Error, key not found! rank: %d, pos: %d, key_str: %s, sizeof(key): %lu\n", 
                                       rank, pos, mn_str, sizeof(key));
     }
     assert(MN_map[pos].first == key);
     
     return pos;

}

void print_output_contigs(std::vector<BasePairVector>& partial_contig_list,
                  const std::string& func)
{

    //printing partial contigs to file
    std::string output_filename(func + std::to_string(rank) + ".fa");
    FILE *fc = fopen(output_filename.c_str(), "w");
    if (fc == NULL)
    {
            printf("Error opening partial output contig file!\n");
            exit(1);
    }

    size_t contig_len_cutoff=400;

    for (size_t j=0; j<partial_contig_list.size(); j++) {
         if (partial_contig_list[j].size() > contig_len_cutoff) {
         std::string output_contig_name(">contig_" + std::to_string(__sync_add_and_fetch(&num_contigs, 1))
                   + "_l_" + std::to_string(partial_contig_list[j].size())
                   + "_r_" + std::to_string(rank));
              std::string out_contig;
              for (size_t k=0; k<partial_contig_list[j].size(); k++)
                   out_contig += el_to_char(partial_contig_list[j][k]);
              fprintf(fc, "%s\n", output_contig_name.c_str());
              fprintf(fc, "%s\n", out_contig.c_str());
        } 
    }
    fclose(fc);

}

void output(BasePairVector &partial_contig)
            //std::vector<BasePairVector> &local_contig_list) 
{
    num_contigs_c++;
    //std::string output_contig_name(">contig_" + std::to_string(num_contigs_c) + "_l_" + std::to_string(partial_contig.size())
    
    std::string output_contig_name(">contig_" + std::to_string(__sync_add_and_fetch(&num_contigs_c, 1)) 
                + "_l_" + std::to_string(partial_contig.size())
                + "_r_" + std::to_string(rank));
    /*
    ContigInfo ci;
    ci.contig_name = output_contig_name;
    ci.contig_data = partial_contig;
    */
    //local_contig_list.push_back(partial_contig);

    
    std::string out_contig;
    for (size_t i=0; i<partial_contig.size(); i++) {
         out_contig += el_to_char(partial_contig[i]);
    }
    
    std::string output_file("contigs_r" + std::to_string(rank) + "_t_" + std::to_string(omp_get_thread_num()) + ".fa");
    std::ofstream outfile;
    outfile.open(output_file, std::ios_base::app);
    outfile << output_contig_name;
    outfile << "\n";
    outfile << out_contig;
    outfile << "\n";
    

}

/* must set: 
 * a) ulimit -s unlimited
 * b) export OMP_STACKSIZE=128M
 * in order to increase the individual stack sizes of the threads
 * The size is proportional to the length of the largest contig
 */
void walk(BasePairVector &partial_contig, int freq, 
        int offset_in_prefix, int prefix_id, const MacroNode &mn,
        std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map,
        std::vector<BasePairVector> &local_contig_list)
{
    assert(prefix_id >=0 && prefix_id < mn.prefixes.size());
    size_t cur_sz = partial_contig.size();
#ifdef WALK_DEBUG 
    if (cur_sz == 500000)
    {
        printf("rank: %d fails with partial_contig size > 400\n", rank);
        std::string out_contig;
        for (size_t t=0; t<partial_contig.size(); t++) {
             out_contig += el_to_char(partial_contig[t]);
        }
        std::ofstream outfile;
        //outfile.open("debug_400k_contig.dat", std::ios_base::out);
        outfile.open("debug_500k_contig.dat", std::ios_base::out);
        outfile << out_contig;
        outfile << "\n";
    }
#endif
    assert(cur_sz < 500000);
    int freq_remaining = freq;

  for(int i=0, off=0; i<mn.prefix_begin_info[prefix_id].num_wires;
      off+=mn.wiring_info[mn.prefix_begin_info[prefix_id].prefix_pos+i].count, i++) 
  {
    int id = mn.wiring_info[mn.prefix_begin_info[prefix_id].prefix_pos+i].suffix_id;
    int sz = mn.wiring_info[mn.prefix_begin_info[prefix_id].prefix_pos+i].count;
    int off_in_suffix = mn.wiring_info[mn.prefix_begin_info[prefix_id].prefix_pos+i].offset_in_suffix;
   
#ifdef WALK_DEBUG 
    fprintf(pFile,"prefix_id: %d, id: %d, sz: %d, off_in_suffix: %d, offset_in_prefix: %d, freq_remaining: %d, freq: %d\n",
            prefix_id, id, sz, off_in_suffix, offset_in_prefix, freq_remaining, freq);
#endif
    //if(off + sz <= offset_in_prefix || off > offset_in_prefix + freq) {
    if(off + sz <= offset_in_prefix || off > offset_in_prefix + freq_remaining) {
      continue;
    }
    int off_in_wire = offset_in_prefix <= off ? 0 : offset_in_prefix - off;
    int next_off = off_in_suffix + off_in_wire;
    int freq_in_wire = std::min(freq_remaining, (sz - off_in_wire));
    partial_contig.append(mn.suffixes[id]);

    if(mn.suffixes_terminal[id]) {
#ifdef WALK_DEBUG 
      fprintf(pFile,"printing output contig\n");
#endif
      //output(partial_contig);
      //output(partial_contig, local_contig_list);
      local_contig_list.push_back(partial_contig);
    } else {
      BasePairVector succ_ext = mn.suffixes[id];
      BasePairVector key = mn.k_1_mer;
      MnodeInfo mn_info = retrieve_mn_sinfo(succ_ext, key);

      int pos = find_mnode_exists(mn_info.search_mn, global_MN_map);
#ifdef WALK_DEBUG 
      if (global_MN_map[pos].first != mn_info.search_mn) {
          printf("Error, key not found! rank: %d, pos: %d, MN_map[pos].first: %lu, key: %lu\n", rank, pos, 
                 global_MN_map[pos].first, mn_info.search_mn);
          
          std::string pred_key, s_ext;
          char new_key[MN_LENGTH+1];
          for (size_t k=0; k<key.size(); k++) {
               pred_key += el_to_char(key[k]);
          }
          for (size_t k=0; k<succ_ext.size(); k++) {
               s_ext += el_to_char(succ_ext[k]);
          }
          for (int k=0; k<MN_LENGTH; k++)
               new_key[k] = el_to_char(kmerel_mn(mn_info.search_mn, k));
          new_key[MN_LENGTH] = '\0';

          printf("pred_key: %s, succ_ext: %s, new_key: %s\n", pred_key.c_str(), s_ext.c_str(), new_key);

          assert(global_MN_map[pos].first == mn_info.search_mn);
      }
#endif
      MacroNode &next_mn = global_MN_map[pos].second; // lookup the next macro node
      std::vector<BasePairVector>::iterator viter = 
          std::find(next_mn.prefixes.begin(), next_mn.prefixes.end(), BasePairVector(mn_info.search_ext));
      if (viter == next_mn.prefixes.end()) {
          printf("Error, prefix not found\n");
          MPI_Finalize();
          exit(1);
      }

      int next_prefix_id = std::distance(next_mn.prefixes.begin(), viter); // find the prefix id this suffix id connects to

#ifdef WALK_DEBUG 
      std::string tmp_contig, tmp_key;
            for (size_t t=0; t<partial_contig.size(); t++)
                 tmp_contig += el_to_char(partial_contig[t]);

            for (size_t l=0; l<next_mn.k_1_mer.size(); l++)
                 tmp_key += el_to_char(next_mn.k_1_mer[l]);

      fprintf(pFile,"walk: partial_contig: %s, freq: %d, next_off: %d, next_pref_id: %d, off: %d, node: %s\n",
              tmp_contig.c_str(), freq_in_wire, next_off, next_prefix_id, off, tmp_key.c_str());

      if (freq_in_wire <= 0) {
         std::string tmp_key;
         for (size_t l=0; l<next_mn.k_1_mer.size(); l++)
                 tmp_key += el_to_char(next_mn.k_1_mer[l]);
         printf("Rank: %d, thread: %d, freq: %d, next_off: %d, next_pref_id: %d, off: %d, node: %s\n",
         rank, omp_get_thread_num(), freq_in_wire, next_off, next_prefix_id, off, tmp_key.c_str());
      }
      //printf("rank: %d, thread: %d, size of local_contig_list: %lu\n", rank, omp_get_thread_num(), local_contig_list.size()); 
#endif
      assert(freq_in_wire > 0);
      walk(partial_contig, freq_in_wire, next_off, next_prefix_id, next_mn, global_MN_map, local_contig_list);
      //walk(partial_contig, freq_in_wire, next_off, next_prefix_id, next_mn, global_MN_map);
    }

#ifdef WALK_DEBUG 
    std::string tmp_key;
    for (size_t l=0; l<mn.k_1_mer.size(); l++)
         tmp_key += el_to_char(mn.k_1_mer[l]);

    fprintf(pFile,"decrementing freq_in_wire, mn: %s, i: %d, off: %d, num_wires: %d, count: %d\n", tmp_key.c_str(), i, off,
            mn.prefix_begin_info[prefix_id].num_wires, mn.wiring_info[mn.prefix_begin_info[prefix_id].prefix_pos+i].count);
#endif

    freq_remaining -= freq_in_wire;
    //printf("rank: %d, partial_contig size: %lu, cur_sz: %lu\n", rank, partial_contig.size(), cur_sz);
    if (partial_contig.size()!= cur_sz)
        partial_contig.resize(cur_sz);
  }

}

void traverse_pakgraph (std::vector<std::pair<kmer_t,MacroNode>> &global_MN_map,
                        std::vector<BeginMN> &list_of_begin_kmers,
                        std::vector<BasePairVector> &partial_contig_list)
{

//#ifdef WALK

    double mn9 = MPI_Wtime();
   
    std::vector< std::vector<BasePairVector> > global_contig_list(num_threads);
#pragma omp parallel
{
    std::vector<BasePairVector> local_contig_list;

#pragma omp for schedule(guided)
    for (size_t i=0; i<list_of_begin_kmers.size(); i++) {
      
        kmer_t key = list_of_begin_kmers[i].node;

#ifdef WALK_DEBUG 
        if (key == 3986664039332806047)
        {
        FILE * pFile;
#endif
        //pFile = fopen ("bkmer_debug.txt","a");
        int terminal_prefix_id = list_of_begin_kmers[i].terminal_prefix_id;
        int pos = find_mnode_exists(key, global_MN_map);
        if (pos<0) {
            printf("Error, the key is absent\n");
            //printf("Error, the key: %lu is absent\n", key);
            MPI_Finalize();
            exit(1);
        }
        else 
        {
            MacroNode &mn_val = global_MN_map[pos].second;
            int freq = mn_val.prefix_count[terminal_prefix_id].second;
            BasePairVector partial_contig;
            for(size_t t=0; t<mn_val.prefixes[terminal_prefix_id].size(); t++)
                partial_contig.push_back(mn_val.prefixes[terminal_prefix_id][t]);
            //for(kmer_t t=0; t<key.size(); t++)
            //    partial_contig.push_back(key[t]);
            partial_contig.append(mn_val.k_1_mer);
            //assert(partial_contig.size() == (mn_val.prefixes[terminal_prefix_id].size()+key.size()));
            //assert(key == mn_val.k_1_mer);
#ifdef WALK_DEBUG 
            std::string tmp_contig, tmp_key;
            for (size_t t=0; t<partial_contig.size(); t++)
                 tmp_contig += el_to_char(partial_contig[t]);

            for (size_t l=0; l<mn_val.k_1_mer.size(); l++)
                 tmp_key += el_to_char(mn_val.k_1_mer[l]);

            fprintf(pFile,"Begin kmer: partial_contig: %s, freq: %d, terminal_pref: %d, node: %s\n",
                    tmp_contig.c_str(), freq, terminal_prefix_id, tmp_key.c_str());
#endif
            walk(partial_contig, freq, 0, terminal_prefix_id, mn_val, global_MN_map, local_contig_list);
            //walk(partial_contig, freq, 0, terminal_prefix_id, mn_val, global_MN_map);
        }

#ifdef WALK_DEBUG 
        }
#endif
    } // end of for loop

    global_contig_list[omp_get_thread_num()] = local_contig_list;
    local_contig_list.clear();

} // end of parallel region 
 
    double mn10 = MPI_Wtime ();
    double walk_time=0.0, global_walk_time=0.0;
    walk_time = mn10 - mn9;
    
    if (rank==0)
        fprintf(stderr, "Completed Walk phase\n");
        
    MPI_Reduce(&walk_time, &global_walk_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for contig generation (walk time) across all procs (secs): %f \n", 
                            (double)global_walk_time/(double)size);

    MPI_Barrier(MPI_COMM_WORLD);

    global_MN_map.clear();
    global_MN_map.shrink_to_fit();

    double mn13 = MPI_Wtime ();
    //printing contigs to file
    print_output_contigs(partial_contig_list, std::string("contigs_par_r"));
    partial_contig_list.clear();
    
    std::string output_filename("contigs_out_r" + std::to_string(rank) + ".fa");
    FILE *fc = fopen(output_filename.c_str(), "w");
    if (fc == NULL)
    {
            printf("Error opening output contig file!\n");
            exit(1);
    }

    size_t contig_len_cutoff=400;

    for (int i=0; i<num_threads; i++) {
         for (size_t j=0; j<global_contig_list[i].size(); j++) {
              if (global_contig_list[i][j].size() > contig_len_cutoff) {
              std::string output_contig_name(">contig_" + std::to_string(__sync_add_and_fetch(&num_contigs, 1))
                   + "_l_" + std::to_string(global_contig_list[i][j].size())
                   + "_r_" + std::to_string(rank));
              std::string out_contig;
              for (size_t k=0; k<global_contig_list[i][j].size(); k++)
                   out_contig += el_to_char(global_contig_list[i][j][k]);
              fprintf(fc, "%s\n", output_contig_name.c_str());
              fprintf(fc, "%s\n", out_contig.c_str());
              }
         }
    }
    fclose(fc);
    double mn14 = MPI_Wtime ();
    double cprint_time=0.0, global_cprint_time=0.0;
    cprint_time = mn14 - mn13;
    
    MPI_Reduce(&cprint_time, &global_cprint_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for printing all contigs (secs): %f \n", 
                            (double)global_cprint_time/(double)size);

    for (int i=0; i<num_threads; i++)
         global_contig_list[i].clear();

    list_of_begin_kmers.clear();

//#endif //end of WALK

}

