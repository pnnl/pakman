# PaKman: A Scalable Algorithm for Generating Genomic Contigs on Distributed Memory Machines

We address the problem of performing large-scale genome assemblies on a distributed memory parallel computer.
PaKman presents a solution for the two most time-consuming phases in the full genome assembly pipeline, namely, k-mer Counting and
Contig Generation. A key aspect of our algorithm is its graph data structure (PaK-Graph), which comprises fat nodes (or what we call
"macro-nodes") that reduce the communication burden during contig generation.


#Dependencies:
----------------
PaKman has the following dependencies:
* MPI library (preferably MPI-3 compatible)
* C++14 (or greater) compliant compiler


Note: 
At this time, PaKman requires the input files to be in the FASTA format and utilizes single-end reads to generate the contigs. 
We are working to extend the functionality to  perform scaffolding accepting paired-end reads as input. 


#Build: 
----------------
Please specify the k-mer length at the time of building:

For e.g: make ksize=32

If ksize is not specified at the time of build, PaKman will build with default size of k=32.

Note:
a) At present, PaKman supports k-mer lengths of k<=64.
b) Pass/update -DLMER_LENGTH in the Makefile to update the LMER size. Default value is set to 8 (recommended).


#Execute:
----------------
mpiexec -np $procs $BINARY -f $INPUT -b $MAX_BATCH_SIZE -r $AVG_READ_LEN -c $COVERAGE -t $BUCKET_CUTOFF -n $MERGE_CUTOFF
For example:
mpiexec -np 4 ./pkmer -f ~/string-graph/inputs/Ecoli_reads_80x.fasta -b 100000000 -r 100 -c 80 -t 21 -n 100000

#Mandatory input arguments to PaKman:
------------------------------
1. -f <INPUT>: input reads file in .fasta format
2. -b <MAX_BATCH_SIZE>: number of k-mers to include in a batch for collective communication. Default value set to 100M (100,000,000). If memory is short, consider reducing to 50M.
3. -r <AVG_READ_LENGTH>: average length of the short reads.
4. -c <COVERAGE>: coverage of the input reads dataset
5. -t <BUCKET_CUTOFF>: number of buckets to consider while determining the cutoff from the k-mer frequency histogram. Default value set to 21.
6. -n <MERGE_CUTOFF>: number of nodes at which the iterative phase of merging macro-nodes is concluded, Default value set to 100,000.

#Contact:
---------
Please contact the following for any queries or support:
-Priyanka Ghosh, PNNL (priyanka dot ghosh at pnnl dot gov)


#CITE:
---------
1) Ghosh, Priyanka, Sriram Krishnamoorthy, and Ananth Kalyanaraman. "PaKman: A Scalable Algorithm for Generating Genomic Contigs on Distributed Memory Machines." IEEE Transactions on Parallel and Distributed Systems (TPDS) vol. 32, no. 5, pp. 1191-1209, 2021. DOI: 10.1109/TPDS.2020.3043241. 

2) Ghosh, Priyanka, Sriram Krishnamoorthy, and Ananth Kalyanaraman. "PaKman: Scalable Assembly of Large Genomes on Distributed Memory Machines." In 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS), pp. 578-589. IEEE, 2019.

