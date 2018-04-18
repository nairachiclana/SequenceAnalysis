// sfunc.h: header file for sequence library
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//User definitions
#define MAXSID 2000000
#define MAXSEQSIZE 90000000000000

//Structure definitions
struct Sequence{
	
	char *id;
	char *seq;
	char *seqRev;
	int len;

};

// Function prototype defintions
int loadSeq(FILE *f, struct Sequence *S, int flag);
void terror(const char * msg);
int kmerIndex(struct Sequence,int , double K);
void kmerIndex2Word(int index, int K, char *);
void printKmers(int , double *, double *, char *, char*, int,double len, double time,   FILE *fout);
void computeKmers(struct Sequence ,double K, double *);
void computeKmersRev(struct Sequence ,double K, double *);
void seqToNum(struct Sequence );
void anadeSeq(struct Sequence , int);
int kmerIndexRev(struct Sequence ,int , double);
void seqRevToNum(struct Sequence );

