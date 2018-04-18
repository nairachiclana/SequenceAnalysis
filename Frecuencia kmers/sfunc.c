#include "sfunc.h"
#include <math.h>
#include <ctype.h>

int loadSeq(FILE *f, struct Sequence *S, int flag){ 
char c;
  int cont=0;
  
  if (flag==0) while ((c=fgetc(f))!='>'); //encuentra primera secuencia

  while( !feof(f)) { 
    while((c=fgetc(f))!='>' && !feof(f)) {//esta dentro de una secuencia
      if((S[cont].id= (char*) malloc(MAXSID*sizeof(char)))==NULL) terror("Not enough memory to allocate the sequence's id");
        if((S[cont].seq= (char*) malloc(MAXSID*sizeof(char)))==NULL) terror("Not enough memory to allocate sequence");
        if((S[cont].seqRev= (char*) malloc(MAXSID*sizeof(char)))==NULL) terror("Not enough memory to allocate sequence");


      //pasar espacios hasta ID
      while((c=(char)fgetc(f))== ' ');
        
      //leer ID hata espacio o salto linea
      int k=0, len=0;
      while(k<MAXSID && c!='\n' && c!= ' ') {
        S[cont].id[k++]=c;
        c=(char)fgetc(f);
      }

      c=(char)fgetc(f);
    
      while(c!='>' && !feof(f)) {
        c=toupper(c);
        if (c>='A' && c<='Z') S[cont].seq[len++] = c;
        c=(char)fgetc(f);
          
      } 
      
      S[cont].len=len;
      cont++;
    } 
  }
  return cont;
}

void terror(const char * msg){
   printf("%s\n", msg);
   exit(-1);
}

void seqToNum(struct Sequence S){
  char c ;
  int i;
  for (i=0;i<S.len;i++){
    switch(toupper(S.seq[i])){
      case 'A':S.seq[i]=0;
      break;
      case 'C':S.seq[i]=1;
      break;
      case'G':S.seq[i]=2;
      break;
      case'T':S.seq[i]=3;
      break;
      default: S.seq[i]=9;
    }
  }
}

void seqRevToNum(struct Sequence S){
  char c ;
  int i;
  for (i=0;i<S.len;i++){
    switch(toupper(S.seqRev[i])){
      case 'A':S.seqRev[i]=0;
      break;
      case 'C':S.seqRev[i]=1;
      break;
      case'G':S.seqRev[i]=2;
      break;
      case'T':S.seqRev[i]=3;
      break;
      default: S.seqRev[i]=9;
    }
  }
}

void anadeSeq(struct Sequence S, int len){


  for(int i=0; i<len ;i++){
    char c= toupper(S.seq[i]);
    
    if( c=='A') { S.seqRev[len-i+1]='T';}
    else if ( c=='T') { S.seqRev[len-i+1]='A';}
    else if ( c=='C') { S.seqRev[len-i+1]='G';}
    else if ( c=='G') { S.seqRev[len-i+1]='C';}
  }
}


int kmerIndex(struct Sequence S,int j, double K){  
   double l=0;
   int i=0;
   double total=0; 

   for (i=j;i<=(j+K-1);i++) { 
      total=total +(S.seq[i] * pow(4,K-l-1));
      l++;
   }
  return total;  
}

int kmerIndexRev(struct Sequence S,int j, double K){  
   double l=0;
   int i=0;
   double total=0; 

   for (i=j;i<=(j+K-1);i++) { 
    total=total +(S.seqRev[i] * pow(4,K-l-1)); 
    l++;
   }
  return total;  
}

void computeKmers(struct Sequence S,double K, double *freq){
  int i=0;
  int value=0;
  for (i=0;i<S.len-K+1;i++){
      value= kmerIndex(S,i,K);
      if (value!=-1){ 
         freq[value]++; 
      }        
  }
}

void computeKmersRev(struct Sequence S,double K, double *freq){
  int i=0;
  int value=0;
  for (i=0;i<S.len-K+1;i++){
      value= kmerIndexRev(S,i,K);
      if (value!=-1){ 
         freq[value]++; 
      }        
  }
}


void kmerIndex2Word(int index, int K, char *kmer){ 
   int resto;
   int cociente = index; 
   char alph[]={'A','C','G','T'};
   int i=0;
      for(i=0;i<K;i++){
         while(cociente>=4){ // 4 is the alphabet size
            resto=cociente%4;
            cociente=cociente/4;
            kmer[K-i-1]=alph[resto];
            i++;
         }
          kmer[K-i-1]=alph[cociente];
      }
}

void printKmers(int tot, double *freq, double*freqRev, char *kmer, char *kmerRev,int K,double len, double time,  FILE *fout){
   
int i;
 fprintf(fout,"Frecuencias para K = %d  en la cadena normal\n"  ,K);
   
   for (i=0;i<tot;i++){
      kmerIndex2Word(i,K, kmer);
      fprintf(fout,"Frecuencia de %s : %f  \n"  ,kmer, freq[i]/len);
      printf("%f\n",freq[i]/len);
      
    }

  fprintf(fout,"Frecuencias para K = %d  en la cadena reversa\n"  ,K);
     for (i=0;i<tot;i++){
      kmerIndex2Word(i,K, kmerRev);
      fprintf(fout,"Frecuencia de %s : %f  \n"  ,kmerRev, freqRev[i]/len);
      printf("%f\n",freqRev[i]/len);
      
      }

   fprintf(fout,"Longitud sec : %f \n"  , len);
   fprintf(stderr, "Tiempo de ejecucion = %f --> para K= %d\n", time, K);
   fprintf(fout,"Tiempo de ejecución = %f --> PARA K = %d\n ", time, K);
}




