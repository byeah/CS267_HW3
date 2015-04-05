#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>

#include "packingDNAseq.h"
#include "kmer_hash.h"


int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
    //printf("Thread %d of %d: hello UPC world\n",MYTHREAD, THREADS);

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	///////////////////////////////////////////
	// Your code for input file reading here //
	///////////////////////////////////////////
	char *input_UFX_name = argv[1];
	int64_t nKmers = getNumKmersInUFX(input_UFX_name);
	int64_t startP = MYTHREAD*(nKmers / THREADS);
	int64_t endP = (MYTHREAD+1)*(nKmers / THREADS);
	if (endP > nKmers) endP = nKmers;
	total_chars_to_read = (endP - startP) * LINE_SIZE;
	
    working_buffer = (unsigned char*) malloc(total_chars_to_read * sizeof(unsigned char));
    FILE *inputFile = fopen(input_UFX_name, "r");
    fseek(inputFile,startP*LINE_SIZE,SEEK_SET);
    int64_t cur_chars_read = fread(working_buffer, sizeof(unsigned char),total_chars_to_read , inputFile);
    fclose(inputFile);
    
	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
	upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	upc_barrier;
	traversalTime += gettime();

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}