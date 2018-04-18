#include <ctype.h>
#include "Funcion1.h"
#include <math.h>


int main(int ac, char** argv){


    FILE *filein1, *filein2, *fileout;
    struct Sequence S1;
    struct Sequence S2;

    struct Kmer *kmers1;
    struct Kmer *kmers2;

    struct Dictionary dic1;
    struct Dictionary dic2; 

    int  K, tot;



    if (ac!=5) terror("USE:  fileIN1 fileIN2 K fileOUT \n");

    fprintf(stderr, "entrada bien\n");

    K=atoi(argv[3]);
    tot = pow(4,K);

  
  if((dic1.kmer= (char **) malloc(tot*sizeof(char*)))==NULL) terror("Not enough memory for dictionary");
  if((dic2.kmer= (char **) malloc(tot*sizeof(char*)))==NULL) terror("Not enough memory for dictionary");

    
    if((S1.seq = (char*) malloc(MAX_SEQ_LEN*sizeof(char)))==NULL)  terror("Not enough memory for sequence");
    if((S2.seq = (char*) malloc(MAX_SEQ_LEN*sizeof(char)))==NULL)  terror("Not enough memory for sequence");


    if((kmers1 = (struct Kmer*) malloc(MAX_KMER_LEN*sizeof(struct Kmer)))==NULL) terror("Not enough memory for kmers");
    if((kmers2 = (struct Kmer*) malloc(MAX_KMER_LEN*sizeof(struct Kmer)))==NULL) terror("Not enough memory for kmers");

     fprintf(stderr, "Memorias\n");

    if ((filein1= fopen(argv[1],"r"))==NULL) terror("Error opening input file");
    if ((filein2= fopen(argv[1],"r"))==NULL) terror("Error opening input file");

    printf("\nFileINS abiertos \n");
  

    loadSeq(filein1, &S1,0);//read sequence from input file
    fclose(filein1);
    loadSeq(filein2, &S2,0);
    fclose(filein2);
    printf("\n loads\n");
  

    int totalKmers1=findKmers(&S1,K, kmers1);//write kmers found in output file
    int totalKmers2=findKmers(&S2,K, kmers2);
     printf("\n findKmers\n");
   

    qsort((void*)kmers1, totalKmers1, sizeof(struct Kmer), funcionQ);
    //bubbleSortKmer(kmers1, totalKmers1);
    //put all same kmers together
    //bubbleSortKmer(kmers2, totalKmers2); 
    //printf("\n bubble2\n");
    qsort(kmers2, totalKmers2, sizeof(struct Kmer), funcionQ);
    printf("\n QSORTS 2\n"); 

    for(int i=0; i<20; i++) {
        for(int j=0;j<50; j++){
            printf("%c", kmers1[i].seq[j]);
        }
        printf("\n");
    }
    
    for(int i=0; i<20; i++) {
        printf("Kmers2:\n %s\n", kmers2[i].seq);
    }
 
    int totalDifKmers1 = ListPositionsForKmer(kmers1,totalKmers1,&dic1); //for each kmer, list all of its positions
    int totalDifKmers2 = ListPositionsForKmer(kmers2,totalKmers2,&dic2); 


 
    if((fileout = fopen(argv[4],"w"))==NULL) terror("Error opening output file");
    diccionarioGlobal(fileout, dic1, dic2, tot);


    printf("Done\n");

    fclose(fileout);
    free(dic1.kmer);
    free(dic2.kmer);
    free(S1.seq);
    free(S2.seq);
    free(kmers1);
    free(kmers2);


    return 0;
}

