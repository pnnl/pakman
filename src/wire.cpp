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

void initiate_mnode_wiring(std::vector<std::pair<kmer_t,MacroNode>> &MN_map)
{

//#ifdef DO_WIRE
    double mn5 = MPI_Wtime ();
    /* Iterate through the macro_nodes and initiate the wiring */

#pragma omp parallel
{
    size_t chunk = (size_t)(MN_map.size()/num_threads);
    size_t lo = chunk*omp_get_thread_num();
    size_t hi = std::min(chunk*(omp_get_thread_num()+1), MN_map.size());
    if (omp_get_thread_num() == 0)
        lo = 0;

    if (omp_get_thread_num() == num_threads-1)
        hi = MN_map.size();

    for (size_t i=lo; i<hi; i++)
    {
         MacroNode &mn = MN_map[i].second;
         mn.prefixes.push_back(BasePairVector{});
         mn.prefixes_terminal.push_back(true);
         mn.prefix_count.push_back(std::make_pair(-1,-1));
         mn.suffixes.push_back(BasePairVector{});
         mn.suffixes_terminal.push_back(true);
         mn.suffix_count.push_back(std::make_pair(-1,-1));

         int num_p = mn.prefixes.size();
         int num_s = mn.suffixes.size();

         mn.wiring_info.resize(num_p+num_s+1);
         mn.prefix_begin_info.resize(num_p);

         mn.setup_wiring();
         /* Check can be commented out if for performance tests */
         //check(mn);
    }
}

    double mn6 = MPI_Wtime ();
    double wire_time=0.0, global_wire_time=0.0;
    wire_time = mn6 - mn5;
    
    MPI_Reduce(&wire_time, &global_wire_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for initial wiring of MN nodes across all procs (secs): %f \n\n", 
                            (double)global_wire_time/(double)size);

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0)
        fprintf(stderr, "Completed wiring of Macro nodes\n");
        
//#endif // end of DO_WIRE

}

