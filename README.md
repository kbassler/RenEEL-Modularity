# RenEEL-Modularity
Reduced Network Extremal Ensemble Learning (RenEEL) Algorithm for Modularity Maximization 

If you use this code, please cite the paper:
J. Guo, P. Singh, K. E. Bassler, Scientific Reports 9, 14234 (2019).

This is an implementation of the Reduced Network Extremal Ensemble Learning (RenEEL) scheme for community detection in complex networks. The example network used here for illustration is the Zachary’s Karate Club network (W. W. Zachary, Journal of Anthropological Research. 33 (4): 452–473 (1977)).

For comments/questions and reporting any bugs that you encounter in the program please contact Kevin E. Bassler (bassler@uh.edu)

Usage: 
To use the RenEEL scheme for maximizing modularity (Q), follow the steps below:

1. Prepare data

1.1 Use an edgelist file including 2 columns separated by space or tab with no header. The network should be unweighted, undirected. Self-loops will be ignored.\
(See example input file karate.txt)\
1.2 Use bash script (work.sh) to generate the three files required by the program. 

Example:
	sh work.sh karate.txt 


2. Compile code

compile main.c, help.c and rg.c with required libraries (math).

Example:
	gcc-9 main.c help.c rg.c  -lm

This will generate file a.out.

3. Run a.out with 3 arguments.

argument 1: Positive Integer, parameter for Randomized Greedy  (usually 2)\
argument 2: Positive Integer, maximum and initial ensemble size of partitions used in RenEEL\
argument 3: Positive Integer, ensemble size of partitions of the reduced network for iteration part in RenEEL\
(seed of random number will be generated using system time)

Example:
	./a.out 2 10 5 

4. Collect results

file 1: partition.txt\ 
file 2: result.txt, a copy will also be printed to stdout\
file 3. records.txt each run will generate a single line for record log.
