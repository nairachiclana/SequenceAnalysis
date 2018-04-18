
#include "sfunc.h"
#include <math.h>
#include <time.h>

int main(int ac, char** av){

   FILE *f;
   FILE *fout;
   
   double **freq, **freqRev, K,tot;
   struct Sequence *S;
   char **kmer, **kmerRev;
   K=atoi(av[2]);
   char* fileName=av[3];
   tot= pow(4,K);
   int Long;
   double len; 
   clock_t t1, t2;
   t1 = clock();



   if (ac!=4) terror ("Missing arguments. Argument order: kmers file.IN K file.OUT");

   if((f=fopen(av[1],"rt"))==NULL) terror("Could not open input sequence file");
    fseek( f, 0L, SEEK_END );
    Long = ftell( f );
    rewind(f);

    printf("open \n");
    if((S = (struct Sequence *) malloc(MAXSID * sizeof(S)))==NULL) terror("Not enough memory to allocate sequences");

    int num_secuencias=loadSeq(f, S, 1);
   
   if((freq = (double**) malloc(num_secuencias*sizeof(double*)))==NULL) terror("Not enough memory to allocate frequencies");
   if((freqRev = (double**) malloc(num_secuencias*sizeof(double*)))==NULL) terror("Not enough memory to allocate reverse frequencies");
   if((kmer = (char**) malloc(Long*sizeof(char*)))==NULL) terror("Not enough memory to allocate kmeros");
   if((kmerRev = (char**) malloc(Long*sizeof(char*)))==NULL) terror("Not enough memory to allocate reverse kmeros");

  for(int i=0; i<num_secuencias; i++) {
      //printf("num secuencia %d\n ", i);
      anadeSeq(S[i], S[i].len); 
        
   
      if((freq[i]= calloc(tot,sizeof(double)))==NULL) terror("Not enough memory 1 ");
      if((freqRev[i]= calloc(tot,sizeof(double)))==NULL) terror("Not enough memory 1 ");
    

      seqToNum(S[i]); 
      seqRevToNum(S[i]); 
     
      computeKmers(S[i], K,freq[i]); 
      computeKmersRev(S[i], K,freqRev[i]);
       
      if((kmer[i]= (char*) malloc(MAXSID*sizeof(char)))==NULL) terror ("Not enough memory2"); //imprimir kameros y espacio de memoria que necesitan
      if((kmerRev[i]= (char*) malloc(MAXSID*sizeof(char)))==NULL) terror ("Not enough memory2"); 
      char *fileName =av[3];

    free(S[i].id);
    free(S[i].seq);
}

    t2 =clock();
    double tt=(double)(t2-t1)/CLOCKS_PER_SEC;

   fout= fopen(fileName,"wt");
   if (fout == NULL) terror("can't open file");

   for(int i=0; i<num_secuencias;i++){
      //printf(fout, "Secuencia %d\n", i );
     printKmers(tot, freq[i], freqRev[i], kmer[i],  kmerRev[i], K,S[i].len, tt, fout);
     
    free(freq[i]);
    free(freqRev[i]);
    free(kmer[i]);
    free(kmerRev[i]);

   }

  fclose(fout);
  fclose(f);
  free(S);

   
   return 0;
}