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

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

extern double kmer_count_time, global_kmer_count_time;
extern double alltoall_time, global_alltoall_time;
extern double alltoallv_time, global_alltoallv_time;
extern double partial_contig_time, global_partial_contig_time;
extern double preproc_time, global_preproc_time;
extern double gather_time, global_gather_time;
extern double barrier_time, global_barrier_time;
extern double sort_time, global_sort_time;

extern double contigs_time, global_contigs_time;
extern double io_time, cleanup_time, vector_time;
extern double count_time, global_count_time;
extern double read_input_time, global_read_input_time;
extern double pack_sbuf_time, global_pack_sbuf_time;
extern double unpack_rbuf_time, global_unpack_rbuf_time;
extern double unpack_rbuf_sort, global_unpack_rbuf_sort;
extern double unpack_rbuf_insert, global_unpack_rbuf_insert;
extern double unpack_rbuf_acc, global_unpack_rbuf_acc;
extern double sl_win_time, global_sl_win_time;
extern double sl_win_time_int, global_sl_win_time_int;
extern double sl_lmer_freq, global_sl_lmer_freq;

extern double vec_insert_time, global_vec_insert_time;
extern double tmap_insert_time, global_tmap_insert_time;

//P2 Timers
extern double p2alltoall_time, p2global_alltoall_time;
extern double p2alltoallv_time, p2global_alltoallv_time;

