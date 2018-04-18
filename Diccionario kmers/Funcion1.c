#include <ctype.h>
#include "Funcion1.h"

int loadSeq(FILE *f, struct Sequence *S, int flag){ 
    char c;
 
    if (flag==0) while ((c=fgetc(f))!='>'); 
    while((c=(char)fgetc(f))== ' '); 
    int k=0, len=0;
    //leer id
    while (k<MAX_ID_LEN && c!='\n' && c!=' ') {c=(char)fgetc(f);}
    c=(char)fgetc(f);
    while (c!='>' && !feof(f)) {
        c=toupper(c);
        if (c>='A' && c<='Z') { 
           len++;
            S->seq[len]=c; 
        }
        c=(char)fgetc(f);
    }

    S->len=len;
    return len;
}


int findKmers(struct Sequence *S, int K, struct Kmer *kmer){

    int i,j;
    for(i=0;i<(S->len-K+1);i++) {
            for(j=i;j<=(i+K-1);j++) 
                kmer[i].seq[j]=S->seq[j];
            kmer[i].pos=i;     
    }
    return i; 
}

int funcionQ(const void *kmer1, const void *kmer2) {
    return strncmp(((struct Kmer*)kmer1)->seq,((struct Kmer*)kmer2)->seq, 2);
}


int KmersUnicos(struct Kmer *kmers, int totalkmers, struct Dictionary *dic){
    int i;
    int m=0, n=0, j=0; 

    if(totalkmers > 0) { //El primero nunca está repetido
        strcpy(dic->kmer[0],kmers[0].seq);
        dic->pos[n][m]=kmers[0].pos;
        n++;
        dic->numPos[j]=1; //en pos j ( kmero0) solo ha salido una vez
    }
    
    for(i=1;i<totalkmers;i++){ //ya ordenados con qsort

        if(strcmp(kmers[i-1].seq,kmers[i].seq)!=0){ 
            j++;
            dic->numPos[j]=1;
            m++;
            n=0;
            strcpy(dic->kmer[m],kmers[i].seq); 
            dic->pos[n][m] = kmers[i].pos;
            n++;
            //pos[posiciones kmero i][kmero i del que miro posiciones]
        }

        else{   //kmers are equal
            dic->numPos[j]++;
            dic->pos[n][m]=kmers[i].pos;
            n++;
        }
    }
     
    return m;//kmeros diferentes
}

int busquedaBinaria(char *kmerdic1, struct Dictionary dic2) {
    
    int inf=0;
    int sup=16;
    int middle; 
    int i=0;
    while (inf<=sup) {
        middle = (inf+sup)/2;
        if (strcmp(kmerdic1,dic2.kmer[i])==0) { return i;}
        
        else if(strcmp(dic2.kmer[i],kmerdic1)<0) {
            sup=middle;
            middle=(inf+sup)/2;
        }

        else if(strcmp(dic2.kmer[i],kmerdic1)>0) {
            inf=middle;
            middle=(inf+sup)/2;
        }
        i++;
    }
    return -1;
}


void diccionarioGlobal(FILE* fout, struct Dictionary dic1, struct Dictionary dic2, double tot) {
 
    int posk;
    int numPos1; 
    int numPos2;

    for(int i=0; i<tot; i++) {  
        numPos1=dic1.numPos[i];
        posk = busquedaBinaria(dic1.kmer[i],dic2); 
        numPos2 =dic2.numPos[posk];
        
        if (posk !=-1) { //coincidencia

           for(int k=0; k<numPos1;k++ ){ //NUM DE POSICIONES QUE TIENE KMER[i]
                for(int j =0; j <numPos2;j++ ) { //NUM DE POSICIONES DE KMER[POSK] 
                    fprintf(fout, "%s: (%i,%i), ", dic1.kmer[i], dic1.pos[k][i], dic2.pos[j][posk]);
            }
        }
        fprintf(fout, "\n");
        }
    }
}

void terror( char * msg){
    printf("%s\n", msg);
    exit(-1);
}








