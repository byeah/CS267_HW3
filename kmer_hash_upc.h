#ifndef KMER_HASH_UPC_H
#define KMER_HASH_UPC_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include <upc.h>
#include "contig_generation.h"

#define LOAD_FACTOR_UPC 2

typedef struct kmer_upc_t kmer_upc_t;
struct kmer_upc_t{
   char kmer[KMER_PACKED_LENGTH];
   char l_ext;
   char r_ext;
   int64_t next;
};

typedef shared kmer_upc_t* kmer_refer;

typedef struct memory_heap_upc_t memory_heap_upc_t;
struct memory_heap_upc_t {
   int64_t posInHeap;
   kmer_refer heap;
   kmer_upc_t *local_heap;
};


typedef struct hash_table_upc_t hash_table_upc_t;
struct hash_table_upc_t {
   int64_t size;           // Size of the hash table
   shared int64_t *table;			// Entries of the hash table are pointers to buckets
   int64_t *local_table;			// Entries of the hash table are pointers to buckets
   memory_heap_upc_t heap;
};


/* Start k-mer data structure */
typedef struct start_kmer_upc_t start_kmer_upc_t;
struct start_kmer_upc_t{
   kmer_upc_t *kmerPtr;
   start_kmer_upc_t *next;
};

/* Creates a hash table and (pre)allocates memory for the memory heap */
hash_table_upc_t * create_hash_table_upc(int64_t nEntries)//, memory_heap_upc_t *memory_heap)
{
   hash_table_upc_t *result;
   int64_t n_buckets = (int64_t)(nEntries * LOAD_FACTOR_UPC);

   result = (hash_table_upc_t*) malloc(sizeof(hash_table_upc_t));
   result->size = n_buckets;
   result->table = upc_alloc(n_buckets*sizeof(int64_t));
   result->local_table = (int64_t*) result->table;
   for(int i=0;i<result->size;i++)
        result->local_table[i]=-1;
    
   
   if (result->table == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", n_buckets, sizeof(bucket_t));
      exit(1);
   }
   memory_heap_upc_t *memory_heap = &result->heap;
   memory_heap->heap = upc_alloc(nEntries * sizeof(kmer_upc_t));
   memory_heap->local_heap = (kmer_upc_t*)memory_heap->heap;
   if (memory_heap->heap == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
      exit(1);
   }
   memory_heap->posInHeap = 0;
   
   return result;
}


/* Looks up a kmer in the hash table and returns a pointer to that entry */
kmer_upc_t* lookup_kmer_helper_upc(hash_table_upc_t* hashtable,const unsigned char *kmer,kmer_upc_t* tmp)
{
   char packedKmer[KMER_PACKED_LENGTH];
   //char tmp[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(hashtable->size, (char*) packedKmer);
   //kmer_refer result;
   //int64_t cur;
   int64_t pos = hashtable->table[hashval*THREADS];
   //result = hashtable->heap.heap+pos*THREADS;
   //for (; result!=NULL; ) {
   for(;pos>=0;){
      //*tmp = hashtable->heap.heap[pos*THREADS];
      upc_memget(tmp,hashtable->heap.heap+pos*THREADS,sizeof(kmer_upc_t));
      //upc_memget(tmp,hashtable->heap.heap[pos*THREADS],KMER_PACKED_LENGTH * sizeof(char));
      if ( memcmp(packedKmer, tmp->kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
         return tmp;
      }
      pos = tmp->next;
   }
   return NULL;
}

/* Adds a kmer and its extensions in the hash table (note that a memory heap should be preallocated. ) */
int add_kmer_upc(hash_table_upc_t *hashtable, const unsigned char *kmer, char left_ext, char right_ext)
{
   /* Pack a k-mer sequence appropriately */
   memory_heap_upc_t *memory_heap = &hashtable->heap;
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(hashtable->size, (char*) packedKmer);
   int64_t pos = memory_heap->posInHeap;
   //printf("%lld\n",hashval);
   
   /* Add the contents to the appropriate kmer struct in the heap */
   memcpy((memory_heap->local_heap[pos]).kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
   (memory_heap->local_heap[pos]).l_ext = left_ext;
   (memory_heap->local_heap[pos]).r_ext = right_ext;
   
   /* Fix the next pointer to point to the appropriate kmer struct */
   (memory_heap->local_heap[pos]).next = hashtable->local_table[hashval];
   /* Fix the head pointer of the appropriate bucket to point to the current kmer */
   hashtable->local_table[hashval] = pos;
   
   /* Increase the heap pointer */
   memory_heap->posInHeap++;
   
   return 0;
}


/* Adds a k-mer in the start list by using the memory heap (the k-mer was "just added" in the memory heap at position posInHeap - 1) */
void addKmerToStartList_upc(memory_heap_upc_t *memory_heap, start_kmer_upc_t **startKmersList)
{
   start_kmer_upc_t *new_entry;
   kmer_upc_t *ptrToKmer;
   
   int64_t prevPosInHeap = memory_heap->posInHeap - 1;
   ptrToKmer = &(memory_heap->local_heap[prevPosInHeap]);
   new_entry = (start_kmer_upc_t*) malloc(sizeof(start_kmer_upc_t));
   new_entry->next = (*startKmersList);
   new_entry->kmerPtr = ptrToKmer;
   (*startKmersList) = new_entry;
}

/* Deallocation functions */
int dealloc_heap_upc(memory_heap_upc_t *memory_heap)
{
   upc_free(memory_heap->heap);
   return 0;
}

int dealloc_hashtable_upc(hash_table_upc_t *hashtable)
{
   upc_free(hashtable->table);
   dealloc_heap_upc(&hashtable->heap);
   return 0;
}


#endif // KMER_HASH_H
