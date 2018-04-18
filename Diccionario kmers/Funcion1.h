#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#define MAX_ID_LEN 900
#define MAX_POS_LEN 90000
#define MAX_KMER_LEN 9000000
#define MAX_KMER_SIZE 9000
#define MAX_SEQ_LEN 90000000
#define KMEROS 16



struct Sequence {
    char id[MAX_ID_LEN];
    char *seq;
    int len;
};

struct Dictionary{
    char **kmer;
    int pos[MAX_KMER_SIZE][KMEROS];
    int numPos[KMEROS];
    int numFilas;
};

struct Kmer{
    char seq[MAX_KMER_SIZE];
    int pos;
};


int loadSeq(FILE *f, struct Sequence *S, int flag);
int findKmers(struct Sequence *S, int K, struct Kmer *kmer);
void bubbleSortKmer(struct Kmer *kmers, int total);
int ListPositionsForKmer(struct Kmer *kmers, int totalkmers, struct Dictionary *dic);
int busquedaBinaria(char *kmerdic1, struct Dictionary dic2);
void diccionarioGlobal(FILE* fout, struct Dictionary dic1, struct Dictionary dic2, double tot);
void terror( char * msg);
int funcionQ(const void *kmer1, const void *kmer2);



