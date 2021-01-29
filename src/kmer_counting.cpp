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
extern int num_threads;
extern int num_batch_transfers;
uint64_t num_recalculate_lmer, global_num_recalculate_lmer;

extern std::vector<lmer_t> lmer_frequency;
extern std::vector<lmer_t> global_lmer_frequency;

extern std::vector<KmerPairs> kmer_proc_buf;


void sort_recv_buffer(std::vector<KmerPairs>& kmer_recv_buf, 
                      std::vector<int>& rcounts_kmer, 
                      std::vector<int>& rdisp_kmer)
{
     size_t rsize=0;
     std::vector<int> offsets(rdisp_kmer);
     offsets.push_back(rcounts_kmer[size-1]+rdisp_kmer[size-1]);
     for (int t=0; t<rcounts_kmer.size(); t++) rsize += rcounts_kmer[t];
     assert(kmer_recv_buf.size() == rsize);
     assert(offsets.size() == (size+1));

     while(offsets.size()>2) {
            //assert(offsets.back() == kmer_recv_buf.size());
            assert(offsets.front() == 0);
            std::vector<int> new_offsets;
            int x = 0;
            while(x+2 < offsets.size()) {
                    // mergesort (offsets[x],offsets[x+1]) and (offsets[x+1],offsets[x+2])
                    std::inplace_merge(kmer_recv_buf.begin()+offsets[x]
                                 ,kmer_recv_buf.begin()+offsets[x+1]
                                 ,kmer_recv_buf.begin()+offsets[x+2] // this *might* be at the end
                                 ,[](const auto& i, const auto& j) {return i.seq < j.seq;} 
                                 );
                    // now they are sorted, we just put offsets[x] and offsets[x+2] into the new offsets.
                    // offsets[x+1] is not relevant any more
                    new_offsets.push_back(offsets[x]);
                    new_offsets.push_back(offsets[x+2]);
                    x += 2;
            }
            // if the number of offsets was odd, there might be a dangling offset
            // which we must remember to include in the new_offsets
            if(x+2==offsets.size()) {
                    new_offsets.push_back(offsets[x+1]);
            }
            // assert(new_offsets.front() == 0);
            //assert(new_offsets.back() == kmer_recv_buf.size());
            offsets.swap(new_offsets);

    }
    offsets.clear();
    offsets.shrink_to_fit();
}  


/*void SortAndAggregate(std::vector<kmer_t>& arr, std::vector<int>& count)
{

     std::vector<int> add;
     std::vector<size_t> indices(count.size());
     std::iota(indices.begin(), indices.end(), 0);

     std::sort(indices.begin(), indices.end(), Comp(arr));
     sort(arr.begin(), arr.end());

     assert(indices.size() == arr.size());
     std::vector<int> temp(count.size());
     for (int it = 0; it < (int)indices.size(); it++){
             temp[it] = count[ indices[it] ];
     }

     count = temp;
     temp.clear();
     indices.clear();

     int aggr=count[0];
     kmer_t prev=arr[0];
     for(int i = 1; i < (int)(arr.size()); i++)
     {
          if (arr[i] == prev)
              aggr += count[i];
           else {
                add.push_back(aggr);
                aggr = count[i];
                prev = arr[i];
           }
     }
     add.push_back(aggr);

     arr.erase( unique( arr.begin(), arr.end() ), arr.end() );
     assert(arr.size() == add.size());
     count = add;

}
*/

void SortAndAggregate(std::vector<KmerPairs>& arr)
{

    std::vector<KmerPairs>::iterator low,up, it;
    //auto cmp = [](const KmerPairs& i, const KmerPairs& j) { return i.seq < j.seq; };
    //size_t new_size = std::set<KmerPairs, decltype(cmp)>(arr.begin(), arr.end(), cmp).size();
    std::vector<KmerPairs> new_arr;
    //new_arr.reserve(new_size);

/*    
#ifdef MEMOPT
    auto cmp = [](const KmerPairs& i, const KmerPairs& j) { return i.seq < j.seq; };
    size_t new_size = std::set<KmerPairs, decltype(cmp)>(arr.begin(), arr.end(), cmp).size();
    std::vector<KmerPairs> new_arr (new_size);
    size_t pos=0;
#else
    std::vector<KmerPairs> new_arr;
#endif
*/ 

    for( it = arr.begin(); it != arr.end(); )
          {
              kmer_t key = (*it).seq;
              low=std::lower_bound (arr.begin(), arr.end(), key,
                                   [] (const KmerPairs& lhs, kmer_t rhs) {
                             return (lhs.seq < rhs);
                             });
 
              up= std::upper_bound (arr.begin(), arr.end(), key,
                                   [] (kmer_t rhs, const KmerPairs& lhs) {
                             return (rhs < lhs.seq);
                             });

              //size_t count = (up-full_vector.begin()) - (low-full_vector.begin());
              //if (count >= 1) {
              int sum=0;
              for (auto itr=low; itr!= up; itr++) {
                        sum += (*itr).k_count;
              }
/*
#ifdef MEMOPT
              new_arr[pos].seq= key;
              new_arr[pos].k_count=sum;
              pos++;
*/
              new_arr.push_back(KmerPairs{key, sum});

              it = up;
          }

     arr=new_arr;
     new_arr.clear();
     new_arr.shrink_to_fit();
}



void transfer_kmers (std::vector<int>& scounts_kmer, 
                     std::vector<KmerPairs> &kmer_send_buf) 
{

    int ssize=0, rsize=0;// disp=0;
    std::vector<int> rcounts_kmer (size,0);
    std::vector<int> rdisp_kmer (size,0);
    std::vector<int> sdisp_kmer (size,0);

    for (int t=0; t<size; t++) ssize += scounts_kmer[t];

    sdisp_kmer[0] = 0;
    for (int i=1; i<size; i++) sdisp_kmer[i] = scounts_kmer[i-1] + sdisp_kmer[i-1];

     //create contiguous derived data type
     MPI_Datatype rowtype;
     MPI_Type_contiguous(sizeof(KmerPairs), MPI_BYTE, &rowtype);
     MPI_Type_commit(&rowtype);

    MPI_Barrier(MPI_COMM_WORLD);

    double t7 = MPI_Wtime ();
    MPI_Alltoall (scounts_kmer.data(), 1, MPI_INT, rcounts_kmer.data(), 1, MPI_INT, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    double t8 = MPI_Wtime ();
    alltoall_time += (t8 - t7);

    for (int t=0; t<size; t++) rsize += rcounts_kmer[t];
    rdisp_kmer[0] = 0;
    for (int i=1; i<size; i++) rdisp_kmer[i] = rcounts_kmer[i-1] + rdisp_kmer[i-1];

    std::vector<KmerPairs> kmer_recv_buf (rsize);

    MPI_Barrier(MPI_COMM_WORLD);

    double t9 = MPI_Wtime ();
     MPI_Alltoallv(kmer_send_buf.data(), scounts_kmer.data(), sdisp_kmer.data(), rowtype,
                   kmer_recv_buf.data(), rcounts_kmer.data(), rdisp_kmer.data(), rowtype, 
                   MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    double t10 = MPI_Wtime ();
    alltoallv_time += (t10 - t9);

    kmer_send_buf.clear();
    kmer_send_buf.shrink_to_fit();

    double t21 = MPI_Wtime ();

    // sort the recv buffer
    double tm1 = MPI_Wtime ();
    //sort_recv_buffer(kmer_recv_buf, rcounts_kmer, rdisp_kmer);
    std::sort(kmer_recv_buf.begin(), kmer_recv_buf.end(),
              [](const auto& i, const auto& j) {return i.seq < j.seq;}
             );

    double tm2 = MPI_Wtime ();
    unpack_rbuf_sort += (tm2 - tm1);

    // Insert and sort
    double tm3 = MPI_Wtime ();
    if (kmer_proc_buf.size())
    {
       size_t offset = kmer_proc_buf.size();
       kmer_proc_buf.insert(kmer_proc_buf.end(), kmer_recv_buf.begin(), kmer_recv_buf.end());

       std::inplace_merge(kmer_proc_buf.begin(),
                  kmer_proc_buf.begin()+offset,
                  kmer_proc_buf.end(), // this *might* be at the end
                  [](const auto& i, const auto& j) { return i.seq < j.seq; }
                 );
    }
    else
      kmer_proc_buf.insert(kmer_proc_buf.end(), kmer_recv_buf.begin(), kmer_recv_buf.end());
       
    double tm4 = MPI_Wtime ();
    unpack_rbuf_insert += (tm4 - tm3);

    kmer_recv_buf.clear();
    kmer_recv_buf.shrink_to_fit();

    double tm5 = MPI_Wtime ();
    SortAndAggregate (kmer_proc_buf);
    double tm6 = MPI_Wtime ();
    unpack_rbuf_acc += (tm6 - tm5);


    double t22 = MPI_Wtime ();
    unpack_rbuf_time += (t22-t21);
    num_batch_transfers++;
  
    // free datatype
    MPI_Type_free(&rowtype);
    
    rcounts_kmer.clear();
    rdisp_kmer.clear();
    sdisp_kmer.clear();
     
 
}


void recalculate_min_lmer (kmer_t kmer_in, lmer_t *m_lmer, lmer_t *m_lmer_freq, int *m_pos)
{
    lmer_t min_lmer=0, tmp_lmer=0;
    lmer_t min_lmer_freq=0, tmp_lmer_freq=0;
    int min_pos=0, k=0;

    for (k=0; ((KMER_LENGTH-1) - k) >= (LMER_LENGTH-1); k++) {
        lmer_t lmer_out=0;
        for(int j=k; j<LMER_LENGTH+k; j++) {
            lmer_out = kmer_to_lmer (kmer_in, j, lmer_out);
        }

        tmp_lmer = lmer_out;
        tmp_lmer_freq = global_lmer_frequency[tmp_lmer];

        if (k == 0) {
            min_lmer = tmp_lmer;
            min_lmer_freq = tmp_lmer_freq;
            min_pos = 0;
        }
        else {
           if (tmp_lmer_freq < min_lmer_freq) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
               min_pos = k;
           }
        }
    }
    assert (k == (KMER_LENGTH-LMER_LENGTH+1));

    *m_lmer = min_lmer;
    *m_lmer_freq = min_lmer_freq;
    *m_pos = min_pos;
}


void Sliding_window_l (const char *ptr, size_t length) {

#ifdef LMER_DEBUG2
      ElType this_alpha;
      char kmer_out[lmer_len+1];
#endif 

  size_t p=0;
  //lmer_frequency.reserve(pow(4, LMER_LENGTH));

  /*find start of a read*/
  for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

  lmer_t kmer = 0;

  while(p<length) {
    assert(ptr[p]=='>'); /*this will be true*/

    /*skip till newline*/
    for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
    p++; /*skip the newline*/

    if(p+LMER_LENGTH > length) break; /*too short a read*/
    kmer = 0;
    int i;
    for(i=0; ptr[p]!='\n' && i<LMER_LENGTH-1; i++) {
      kmer = lmer_shift(kmer, char_to_el(ptr[p++]));
      //kmer = kmer_cons(kmer, i, char_to_el(ptr[p++]));
    }

    while(p<length && ptr[p]!='\n') {
      kmer = lmer_shift(kmer, char_to_el(ptr[p++]));

#ifdef LMER_DEBUG2
      printf("lmer: %lu,");
      for(int j=0; j<LMER_LENGTH; j++) {
          this_alpha = lmerel (kmer, j);
          lmer_out[j] = el_to_char(this_alpha);
      }
      lmer_out[lmer_len] = '\0';
      printf(" lmer_out: %s \n", lmer_out);
#endif

      lmer_frequency[kmer]++;

    }
    p++; /*skip the newline*/
  }

}

void Sliding_window (const char *ptr, size_t length, int *n_kmers, 
                     std::vector<std::vector<kmer_t>> &partial_kmer_counts)
{
#ifdef LMER_DEBUG
    FILE *fp_d;
    char debug_file_name[25];
    char proc_id[3];

    sprintf(proc_id, "%d", rank); 
    strcpy(debug_file_name,"debug_p");
    strcpy(&debug_file_name[strlen(debug_file_name)],proc_id);
    strcpy(&debug_file_name[strlen(debug_file_name)],".log");
    fp_d = fopen (debug_file_name, "w");
    
    /*check to see if it opened okay */
    if (fp_d == NULL)
    {
		printf ("Error opening proc %d 's dump file \n", rank);
		exit (0);
    }
#endif

  size_t p=0;
  std::vector<int> scounts_kmer (size,0);
  std::vector<KmerPairs> kmer_send_buf;
  size_t kpos=0;

  //char kmer_out[WINDW_SIZE_PLUS];

  /*find start of a read*/
  for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

  int num_kmers=*n_kmers;
  kmer_t kmer = 0; 
  lmer_t lmer_out = 0, min_lmer=0, tmp_lmer=0;
  uint64_t min_lmer_freq=0, tmp_lmer_freq=0;
  int min_pos=0, tmp_pos=0;

  while(p<length) {
    assert(ptr[p]=='>'); /*this will be true*/

    /*skip till newline*/
    for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
    p++; /*skip the newline*/

    if(p+KMER_LENGTH > length) break; /*too short a read*/

    kmer=0, lmer_out=0;
    min_lmer=0, min_lmer_freq=0;
    tmp_lmer=0, tmp_lmer_freq=0;
    min_pos=0, tmp_pos=0;
    int i;

    //double TM1 = MPI_Wtime();
    for(i=0; ptr[p]!='\n' && i<KMER_LENGTH-1; i++) {
      kmer = kmer_shift(kmer, char_to_el(ptr[p]));

      if (i<LMER_LENGTH-1) 
          lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));
      else {
           lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));

           tmp_lmer = lmer_out;
           tmp_lmer_freq = global_lmer_frequency[tmp_lmer];
           tmp_pos = i-(LMER_LENGTH-1);

           if (i == LMER_LENGTH-1) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
           }
           else {
                if (tmp_lmer_freq < min_lmer_freq) {
                    min_lmer = tmp_lmer;
                    min_lmer_freq = tmp_lmer_freq;
                    min_pos = tmp_pos;
                }
           }
      }
      p++;
    }
    //double TM2 = MPI_Wtime();
    //kmer_shift_time += (TM2-TM1);

    while(p<length && ptr[p]!='\n') {
      //double TM3 = MPI_Wtime();
      kmer = kmer_shift(kmer, char_to_el(ptr[p]));
      lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));
      uint64_t lmer_out_freq = global_lmer_frequency[lmer_out];
   
      //double TM5 = MPI_Wtime();
      if (min_pos < 0) {
          recalculate_min_lmer(kmer, &min_lmer, &min_lmer_freq, &min_pos);
          num_recalculate_lmer++;
      }
      //double TM6 = MPI_Wtime();
      //kmer_recal_time += (TM6-TM5);

      if (lmer_out_freq < min_lmer_freq) {
          min_lmer = lmer_out;
          min_lmer_freq = lmer_out_freq;
          min_pos = KMER_LENGTH-LMER_LENGTH;
          //min_lmer_freq = global_lmer_frequency[lmer_out];
      }
      p++;
      min_pos--;
      //double TM4 = MPI_Wtime();
      //kmer_shift_time += (TM4-TM3);

#ifdef LMER_DEBUG
      fprintf(fp_d, "kmer: %lu, lmer: %lu, lmer_freq: %lu, min_lmer: %lu, min_freq: %lu\n", kmer,lmer_out,lmer_out_freq,min_lmer,min_lmer_freq);
#endif
#if 0
      /* Murmerhash method */
      /*
      calculate owner using MurmurHash and update corresponding scounts_kmer
      got = tmp_k_map.find(kmer);
      if ( got == tmp_k_map.end() ) 
          scounts_kmer[(int) MurmurHash64A ((char*)&kmer, sizeof(kmer), SEED, size)]++;

      tmp_k_map[kmer]++;
      */

#endif
      double T1 = MPI_Wtime();
#ifdef NOOPT
      kmers_per_proc[retrieve_proc_id(min_lmer)].push_back(kmer);
#endif
      partial_kmer_counts[retrieve_proc_id(min_lmer)].push_back(kmer);
      double T2 = MPI_Wtime();
      vec_insert_time += (T2-T1);

      num_kmers++;

      if (num_kmers > MAX_KMER_COUNT){
          double t25 = MPI_Wtime ();
          // initiate collective communication to pass k-mers and their respective counts to rightful owners
          // calculate global owner of each k-mer and populate the k-mer and count to their respective 'p'th local buffers
          // if global position of a k-mer in my k_map != me, delete k-mer from my k_map
          // iterate through k-mers recieved after collective communication ends, and add k-mers to my k_map 
          // reset num_kmers count to 0.
          
          double T5 = MPI_Wtime();
#ifdef NOOPT
          for (int t=0; t<size; t++)
          {
		  int counter=1;
                  sort(kmers_per_proc[t].begin(), kmers_per_proc[t].end());
                  kmer_t prev=kmers_per_proc[t][0];
		  for(int i = 1; i < (int)(kmers_per_proc[t].size()); i++)
		  {
			 if (kmers_per_proc[t][i] == prev) {
			     counter++;
                         } else {
			     kmer_cnt_tmp_buf[t].push_back(counter);
			     counter=1;
                             prev=kmers_per_proc[t][i];
			 } 
		     
		  }
                  kmer_cnt_tmp_buf[t].push_back(counter);

                  kmers_per_proc[t].erase( unique( kmers_per_proc[t].begin(), kmers_per_proc[t].end() ), kmers_per_proc[t].end() );
                  assert(kmers_per_proc[t].size() == kmer_cnt_tmp_buf[t].size());
                  scounts_kmer[t] = kmer_cnt_tmp_buf[t].size(); 
                  //scounts_kmer[t] = std::accumulate(kmer_cnt_tmp_buf[t].begin(), kmer_cnt_tmp_buf[t].end(), 0); 
          }
#else
          for (int t=0; t<size; t++)
          {
		  int counter=1;
                  kpos=0;
                  sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                  kmer_t prev=partial_kmer_counts[t][0];
		  for(int i = 1; i < (int)(partial_kmer_counts[t].size()); i++)
		  {
			 if (partial_kmer_counts[t][i] == prev) {
			     counter++;
                         } else {
                             kmer_send_buf.push_back(KmerPairs{prev, counter});
                             //kmers_per_proc[t].push_back(KmerPairs{prev, counter});
			     kpos++;
                             counter=1;
                             prev=partial_kmer_counts[t][i];
			 }
		     
		  }
                  //kmers_per_proc[t].push_back(KmerPairs{prev, counter});
                  kmer_send_buf.push_back(KmerPairs{prev, counter});
                  kpos++;

                  //kmers_per_proc[t].erase( unique( kmers_per_proc[t].begin(), kmers_per_proc[t].end() ), kmers_per_proc[t].end() );
                  //assert(kmers_per_proc[t].size() == kmer_cnt_tmp_buf[t].size());
                  //scounts_kmer[t] = kmers_per_proc[t].size();
                  scounts_kmer[t] = kpos;
                  partial_kmer_counts[t].clear();
                  partial_kmer_counts[t].shrink_to_fit();
          }
          /*
          for (int t=0; t<size; t++)
          {
                  sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                  std::vector<kmer_t>::iterator low,up;
		  for(auto it=partial_kmer_counts[t].begin(); it != partial_kmer_counts[t].end(); )
		  {
                      low=std::lower_bound (partial_kmer_counts[t].begin(), partial_kmer_counts[t].end(), *it);
                      up= std::upper_bound (partial_kmer_counts[t].begin(), partial_kmer_counts[t].end(), *it);
                      int ncount = std::count (low, up, *it);
                      kmers_per_proc[t].push_back(KmerPairs{*it, ncount});
                      it = up;
		  }
                  scounts_kmer[t] = kmers_per_proc[t].size(); 
          }
          */
#endif
          double T6 = MPI_Wtime();
          tmap_insert_time += (T6-T5);

          /*double T5 = MPI_Wtime();
          for(std::vector<kmer_t>::iterator it = kmer_tmp_buf.begin(); it != kmer_tmp_buf.end(); ++it) {
              if (tmp_k_map[*it]++ == 0)
                  scounts_kmer[(int) MurmurHash64A ((char*)&(*it), sizeof(*it), SEED, size)]++;
          }
          double T6 = MPI_Wtime();
          tmap_insert_time += (T6-T5);
          */

          //clear the partial counts
          /*
          for(int k=0; k<size; k++)
              partial_kmer_counts[k].clear();         
          */

#ifdef NOOPT
          transfer_kmers (scounts_kmer, kmers_per_proc, kmer_cnt_tmp_buf);
#else
          transfer_kmers (scounts_kmer, kmer_send_buf);
#endif
          num_kmers = 0;

#ifdef NOOPT
          for (int t=0; t<size; t++) {
               kmers_per_proc[t].clear();
               kmer_cnt_tmp_buf[t].clear();
          }
#endif
          //memset(scounts_kmer, 0, size*sizeof(*scounts_kmer));
          scounts_kmer.clear(); 
          double t26 = MPI_Wtime ();
          sl_win_time_int += (t26-t25);
       } // end of if condition
    } // end of while loop
#ifdef LMER_DEBUG
    fprintf(fp_d, "----------------------------\n");
#endif
    p++; /*skip the newline*/
  }

#ifdef LMER_DEBUG
          fclose (fp_d);
#endif

  *n_kmers = num_kmers;
  scounts_kmer.shrink_to_fit();

}

void process_remaining_kmers(
                     std::vector<std::vector<kmer_t>> &partial_kmer_counts) 
{

     std::vector<int> scounts_kmer (size,0);
     std::vector<KmerPairs> kmer_send_buf;
     size_t kpos=0;

      double T7 = MPI_Wtime();
          for (int t=0; t<size; t++)
          {
             if(partial_kmer_counts[t].size())
             {
		  int counter=1;
                  kpos=0;
                  sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                  kmer_t prev=partial_kmer_counts[t][0];
		  for(int i = 1; i < (int)(partial_kmer_counts[t].size()); i++)
		  {
			 if (partial_kmer_counts[t][i] == prev) {
			     counter++;
                         } else {
                             kmer_send_buf.push_back(KmerPairs{prev, counter});
			     kpos++;
			     counter=1;
                             prev=partial_kmer_counts[t][i];
			 }
		     
		  }
                  kmer_send_buf.push_back(KmerPairs{prev, counter});
                  kpos++;

                  //kmers_per_proc[t].erase( unique( kmers_per_proc[t].begin(), kmers_per_proc[t].end() ), kmers_per_proc[t].end() );
                  //assert(kmers_per_proc[t].size() == kmer_cnt_tmp_buf[t].size());
                  //scounts_kmer[t] = kmers_per_proc[t].size(); 
                  scounts_kmer[t] = kpos;
                  partial_kmer_counts[t].clear();
                  partial_kmer_counts[t].shrink_to_fit();
              }
           }
      /*
      for(int t=0; t<size; t++)
      {
           if(partial_kmer_counts[t].size())
           {
                  sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                  std::vector<kmer_t>::iterator low,up;
		  for(auto it=partial_kmer_counts[t].begin(); it != partial_kmer_counts[t].end(); )
		  {
                      low=std::lower_bound (partial_kmer_counts[t].begin(), partial_kmer_counts[t].end(), *it);
                      up= std::upper_bound (partial_kmer_counts[t].begin(), partial_kmer_counts[t].end(), *it);
                      int ncount = std::count (low, up, *it);
                      kmers_per_proc[t].push_back(KmerPairs{*it, ncount});
                      it = up;
		  }
                  scounts_kmer[t] = kmers_per_proc[t].size(); 

             
             // if (kmers_per_proc[t].size() != kmer_cnt_tmp_buf[t].size())
             // {
             //     fprintf(stderr, "counters not matching, for rank: %d, t: %d!! kmers_per_proc size: %lu, kmer_cnt_tmp_buf: %lu\n",
             //             rank, t, kmers_per_proc[t].size(), kmer_cnt_tmp_buf[t].size());

             //     
             //     char pros_id[3];
             //     char outfile_name[25];

             //     sprintf(pros_id, "%d", rank); 
             //     strcpy(outfile_name,"tmp_kmers_");
             //     strcpy(&outfile_name[strlen(outfile_name)],pros_id);
             //     strcpy(&outfile_name[strlen(outfile_name)],".log");
             //     FILE *ft = fopen(outfile_name, "w");
             //     if (ft == NULL)
             //     {
             //        printf("Error opening tmp file!\n");
             //        exit(1);
             //     }
             //     
             //     for (size_t m=0; m<kmers_per_proc[t].size(); m++)
             //          fprintf(ft, "m: %lu, kmer: %lu\n", m, kmers_per_proc[t][m]);
             //     fprintf(ft, "----------- end of kmers ------------\n");
             //     for (size_t m=0; m<kmer_cnt_tmp_buf[t].size(); m++)
             //          fprintf(ft, "m: %lu, cnt: %d\n", m, kmer_cnt_tmp_buf[t][m]);

             //     fclose(ft);

             // }
             // assert(kmers_per_proc[t].size() == kmer_cnt_tmp_buf[t].size());
             // scounts_kmer[t] = kmer_cnt_tmp_buf[t].size();
               
           }
          }*/
          double T8 = MPI_Wtime();
          tmap_insert_time += (T8-T7);

          //clear the partial counts
          for(int k=0; k<size; k++)
              partial_kmer_counts[k].clear();         
 
          transfer_kmers (scounts_kmer, kmer_send_buf);
         
          scounts_kmer.clear();
          scounts_kmer.shrink_to_fit();

}

void print_kmer_count_timers()
{

    MPI_Reduce(&alltoall_time, &global_alltoall_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for Alltoall across all procs (secs): %f \n", 
                            (double)global_alltoall_time/(double)size);

    MPI_Reduce(&pack_sbuf_time, &global_pack_sbuf_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for pack_sbuf_time across all procs (secs): %f \n", 
                            (double)global_pack_sbuf_time/(double)size);

    MPI_Reduce(&alltoallv_time, &global_alltoallv_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for AlltoallV across all procs (secs): %f \n", 
                            (double)global_alltoallv_time/(double)size);

    MPI_Reduce(&unpack_rbuf_time, &global_unpack_rbuf_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_time/(double)size);

    MPI_Reduce(&unpack_rbuf_sort, &global_unpack_rbuf_sort, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time:sort across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_sort/(double)size);

    MPI_Reduce(&unpack_rbuf_insert, &global_unpack_rbuf_insert, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time:insert across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_insert/(double)size);

    MPI_Reduce(&unpack_rbuf_acc, &global_unpack_rbuf_acc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time:acc across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_acc/(double)size);

    MPI_Reduce(&sl_win_time_int, &global_sl_win_time_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for partial Sliding Window across all procs (secs): %f \n", 
                            (double)global_sl_win_time_int/(double)size);

    MPI_Reduce(&sl_win_time, &global_sl_win_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for Sliding Window across all procs (secs): %f \n", 
                            (double)global_sl_win_time/(double)size);

    MPI_Reduce(&sl_lmer_freq, &global_sl_lmer_freq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for lmer Sliding Window across all procs (secs): %f \n", 
                            (double)global_sl_lmer_freq/(double)size);

    MPI_Reduce(&vec_insert_time, &global_vec_insert_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for vector insert across all procs (secs): %f \n", 
                            (double)global_vec_insert_time/(double)size);

    MPI_Reduce(&tmap_insert_time, &global_tmap_insert_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for temp_map insert across all procs (secs): %f \n", 
                            (double)global_tmap_insert_time/(double)size);

    MPI_Reduce(&num_recalculate_lmer, &global_num_recalculate_lmer, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average number of re-calculations of min l-mer across all procs (secs): %f \n", 
                            (double)global_num_recalculate_lmer/(double)size);

    //MPI_Reduce(&kmer_recal_time, &global_kmer_recal_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //if (rank == 0) printf ("Average time for recalculating min-mers across all procs (secs): %f \n", 
    //                        (double)global_kmer_recal_time/(double)size);

    //MPI_Reduce(&kmer_shift_time, &global_kmer_shift_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //if (rank == 0) printf ("Average time for k-mer construction across all procs (secs): %f \n", 
    //                        (double)global_kmer_shift_time/(double)size);

    MPI_Reduce(&kmer_count_time, &global_kmer_count_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for performing k-mer counting across all procs (secs): %f \n", 
                            (double)global_kmer_count_time/(double)size);

    uint64_t all_proc_kmer_count = 0;
    uint64_t tmp_kmer_count = kmer_proc_buf.size();
    //assert(kmer_proc_buf.size() == kmer_cnt_proc_buf.size());

    MPI_Reduce(&tmp_kmer_count, &all_proc_kmer_count, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) { 
        printf("Total distinct k-mer entries across all proc's: %lu\n", all_proc_kmer_count); 
        printf("Number of batch iterations: %d \n", num_batch_transfers);
    }

}

void free_kmer_count_buffers()
{
   
    kmer_proc_buf.clear();
    kmer_proc_buf.shrink_to_fit();

}

void perform_kmer_counting (const char *read_data, size_t length)
{

    double start_t = MPI_Wtime ();

    //calculate frequencies of all l-mers in the Read dataset
    double time_l1 = MPI_Wtime ();
    Sliding_window_l(read_data, length);
    double time_l2 = MPI_Wtime ();
    sl_lmer_freq = time_l2 - time_l1;

    // Perform Allreduce to obtain global lmer counts across all proc's
    int num_lmers = pow(4, LMER_LENGTH);
    //global_lmer_frequency.reserve(num_lmers);

    MPI_Allreduce (lmer_frequency.data(), global_lmer_frequency.data(), num_lmers, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef NOOPT
    std::vector< std::vector<kmer_t> > kmers_per_proc; //, std::vector<kmer_t>);
    std::vector< std::vector<int> > kmer_cnt_tmp_buf; //, std::vector<int>);
    for (int i=0; i<size; i++) {
         kmers_per_proc.push_back(std::vector<kmer_t> ());
         kmer_cnt_tmp_buf.push_back(std::vector<int> ());
    }
#else
    std::vector< std::vector<kmer_t> > partial_kmer_counts(size);
#endif
    int num_kmers = 0;

    double t23 = MPI_Wtime ();

    if (rank==0)
        fprintf(stderr, "Starting sliding window\n");
        
#ifdef NOOPT
    Sliding_window (rdata.read_data, rdata.read_data_size, &num_kmers, kmers_per_proc, kmer_cnt_tmp_buf);
#else
    Sliding_window (read_data, length, &num_kmers, partial_kmer_counts);
#endif

    double t24 = MPI_Wtime ();
    sl_win_time = t24 - t23;

    if (rank==0)
        fprintf(stderr, "Completed sliding window\n");
        
    // initiate communication for the residual kmers
    //printf("rank: %d, num_kmers: %d\n", rank, num_kmers);
    if (num_kmers) {
         
          if (rank==0) fprintf(stderr, "remaining sliding window, num_kmers: %d\n", num_kmers);
#ifdef NOOPT
          process_remaining_kmers (kmers_per_proc, kmer_cnt_tmp_buf);
#else
          process_remaining_kmers (partial_kmer_counts);
#endif
    }

#ifdef NOOPT
    kmers_per_proc.clear();
    kmers_per_proc.shrink_to_fit();
    kmer_cnt_tmp_buf.clear();
    kmer_cnt_tmp_buf.shrink_to_fit();
#else
    for (int i=0; i<size; i++)
         partial_kmer_counts[i].clear();
    partial_kmer_counts.shrink_to_fit(); 
#endif

    lmer_frequency.clear();
    global_lmer_frequency.clear();

    double end_t = MPI_Wtime ();

    kmer_count_time = end_t - start_t;

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0)
        fprintf(stderr, "Completed kmer counting\n");

    print_kmer_count_timers();
        
}


