#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>

#include "packingDNAseq.h"
#include "kmer_hash.h"
#include "kmer_hash_upc.h"

//Shared data
shared unsigned char*shared pointers[THREADS];

//Private data
int64_t nKmers, total;
unsigned char *buffer;

void read_data(char *input_UFX_name)
{
	total = getNumKmersInUFX(input_UFX_name);
	int64_t perThread = (total + THREADS - 1) / THREADS;
	int64_t startP = MYTHREAD * perThread;
	int64_t endP = (MYTHREAD + 1) * perThread;
	if (endP > total) endP = total;
	nKmers = endP - startP;
	int64_t total_chars_to_read = nKmers * LINE_SIZE;
	
    buffer = malloc(total_chars_to_read * sizeof(unsigned char));
    FILE *inputFile = fopen(input_UFX_name, "r");
    fseek(inputFile,startP*LINE_SIZE,SEEK_SET);
    int64_t cur_chars_read = fread(buffer, sizeof(unsigned char), total_chars_to_read, inputFile);
    fclose(inputFile);

    if (MYTHREAD == 0)
        printf("Finish Reading\n");
}

int64_t getHashVal(int64_t size,const unsigned char* kmer)
{
    char packedKmer[KMER_PACKED_LENGTH];
    packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
    int64_t hashval = hashkmer(size, (char*) packedKmer);
    return hashval;
}

void prepare_localdata(shared int64_t* pos,shared int64_t* length)
{
    int64_t N = nKmers * LINE_SIZE;

    pointers[MYTHREAD] = upc_alloc((N + 1) * sizeof(unsigned char));
    unsigned char* working_buffer = (unsigned char*) pointers[MYTHREAD];
    int *id = malloc(nKmers * sizeof(int));
    int64_t *cnt = malloc(THREADS * sizeof(int64_t));

    memset(cnt, 0, sizeof(int64_t)*THREADS);
    
    int64_t ptr = 0, i = 0;
    
    while (ptr < N)
    {
        int64_t hVal = getHashVal(total, buffer + ptr);
        id[i] = hVal % THREADS;
        cnt[hVal % THREADS]++;
        i++;
        ptr += LINE_SIZE;
    }
    //for(int k=0;k<THREADS;k++)
        //printf("%d to %d:  %lld\n",MYTHREAD,k,cnt[k]);
        
    for (int k = 0; k < THREADS; k++)
        length[MYTHREAD * THREADS + k] = cnt[k];

    for (int k = 1; k < THREADS; k++)
        cnt[k] += cnt[k-1];
        
    ptr = 0, i = 0;
    while (ptr < N)
    {
        int64_t offset = (--cnt[id[i]]) * LINE_SIZE;
        memcpy(working_buffer + offset, buffer + ptr, LINE_SIZE * sizeof(char));
        i++;
        ptr += LINE_SIZE;
    }

    for (int k = 0; k < THREADS; k++)
        pos[MYTHREAD * THREADS + k] = cnt[k];
    
    upc_barrier;

    free(buffer);
    free(id);
    free(cnt);
}


void shuffle_data()
{
    shared int64_t* pos = upc_all_alloc(THREADS, THREADS * sizeof(int64_t));
    shared int64_t* length = upc_all_alloc(THREADS, THREADS * sizeof(int64_t));

    prepare_localdata(pos, length);
    
    //Pull remote data
    nKmers = 0;
    
    for (int i = 0; i < THREADS; i++)
        nKmers += length[i * THREADS + MYTHREAD];

    printf("%d: nKmers: %lld\n", MYTHREAD, nKmers);
    
    unsigned char* ptr = malloc((nKmers * LINE_SIZE + 1) * sizeof(char));
    buffer = ptr;

    for (int i = 0; i < THREADS; i++){
        shared unsigned char* pointer = pointers[i];
        upc_memget(ptr, pointer + (pos[i * THREADS + MYTHREAD] * LINE_SIZE * THREADS), length[i * THREADS + MYTHREAD] * LINE_SIZE * sizeof(char));
        ptr += length[i*THREADS+MYTHREAD] * LINE_SIZE;
        //printf("%d from %d at %d len %d\n",MYTHREAD,i,pos[i*THREADS+MYTHREAD],length[i*THREADS+MYTHREAD]);
    }
    /*printf("%d\n",sizeof(char));
    char s[LINE_SIZE+1];
    shared unsigned char* pointer = pointers[1];
    upc_memget(s,pointer+LINE_SIZE*THREADS,LINE_SIZE);
    s[LINE_SIZE]=0;
    printf("%s",s);*/
    
    /*char file[5]="out00";
    file[3]=(char)(MYTHREAD+48);file[4]=0;
    printf("%s\n",file);
    FILE *fout = fopen(file,"w");
    buffer[nKmers*LINE_SIZE]=0;
    fprintf(fout,"%s",buffer);
    fclose(fout);*/
}

hash_table_upc_t *hashtable;
start_kmer_upc_t *startKmersList = NULL;

shared hash_table_upc_t hashtables[THREADS];
hash_table_upc_t* local_tables;

void preprocessing()
{
    /* Create a hash table */
    hashtable = create_hash_table_upc(nKmers);
    /* Process the working_buffer and store the k-mers in the hash table */
    /* Expected format: KMER LR ,i.e. first k characters that represent the kmer, then a tab and then two chatacers, one for the left (backward) extension and one for the right (forward) extension */
    int64_t N = nKmers*LINE_SIZE, ptr = 0;
    while (ptr < N) {
      /* working_buffer[ptr] is the start of the current k-mer                */
      /* so current left extension is at working_buffer[ptr+KMER_LENGTH+1]    */
      /* and current right extension is at working_buffer[ptr+KMER_LENGTH+2]  */
    
      char left_ext = (char) buffer[ptr + KMER_LENGTH+1];
      char right_ext = (char) buffer[ptr + KMER_LENGTH+2];
      /* Add k-mer to hash table */
      add_kmer_upc(hashtable, &buffer[ptr], left_ext, right_ext);
      /* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
      if (left_ext == 'F') {
         addKmerToStartList_upc(&hashtable->heap, &startKmersList);
      }
    
      /* Move to the next k-mer in the input working_buffer */
      ptr += LINE_SIZE;
    }
    //!!!@!@!@!@!@!hashtables[MYTHREAD] = *hashtable;
    upc_memput(hashtables+MYTHREAD,hashtable,sizeof(hash_table_upc_t));
    //printf("%lld\n",(*hashtable).size);
    //printf("%lld\n",hashtables[MYTHREAD].size);
    upc_barrier;
    local_tables = (hash_table_upc_t*)malloc(THREADS*sizeof(hash_table_upc_t));
    for(int i=0;i<THREADS;i++)
        upc_memget(local_tables+i,hashtables+i,sizeof(hash_table_upc_t));
        //local_tables[i]=hashtables[i];
}

char cur_contig[MAXIMUM_CONTIG_SIZE];
kmer_upc_t tmp;

kmer_upc_t* lookup_kmer_upc(const unsigned char *kmer)
{
    int64_t hashval = getHashVal(total,kmer);
    return lookup_kmer_helper_upc(local_tables+hashval%THREADS,kmer,&tmp);
}

void traversal()
{
    char file[20];
    sprintf(file,"pgen%d.out",MYTHREAD);
    FILE* fout = fopen(file, "w");
    char unpackedKmer[KMER_LENGTH+1];
    /* Pick start nodes from the startKmersList */
    
    start_kmer_upc_t *curStartNode = startKmersList;
    
    while (curStartNode != NULL ) {
      /* Need to unpack the seed first */
      kmer_upc_t* cur_kmer_ptr = curStartNode->kmerPtr;
      unpackSequence((unsigned char*) cur_kmer_ptr->kmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);
      /* Initialize current contig with the seed content */
      memcpy(cur_contig ,unpackedKmer, KMER_LENGTH * sizeof(char));
      int posInContig = KMER_LENGTH;
      char right_ext = cur_kmer_ptr->r_ext;
    
      /* Keep adding bases while not finding a terminal node */
      while (right_ext != 'F') {
         cur_contig[posInContig] = right_ext;
         posInContig++;
         /* At position cur_contig[posInContig-KMER_LENGTH] starts the last k-mer in the current contig */
         cur_kmer_ptr = lookup_kmer_upc((const unsigned char *) &cur_contig[posInContig-KMER_LENGTH]);
         //printf("!!Fin look\n");
         //cur_kmer_ptr = lookup_kmer_upc(hashtable, (const unsigned char *) &cur_contig[posInContig-KMER_LENGTH],&tmp);
         if (cur_kmer_ptr == NULL){
             printf("NULL\n");
             break;
         }
         right_ext = cur_kmer_ptr->r_ext;
      }
    
      /* Print the contig since we have found the corresponding terminal node */
      cur_contig[posInContig] = '\0';
      fprintf(fout,"%s\n", cur_contig);
      //contigID++;
      //totBases += strlen(cur_contig);
      /* Move to the next start node in the list */
      curStartNode = curStartNode->next;
    }
    
    fclose(fout);
}

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
	init_LookupTable();

	read_data(argv[1]);
    
	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
	shuffle_data();
	preprocessing();
	//printf("Finish preprocessing\n");
	upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	traversal();
	
	upc_barrier;
	traversalTime += gettime();

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD == 0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}