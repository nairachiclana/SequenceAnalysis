
from Bio import SeqIO

f1='/Users/nairachiclana/Desktop/2.Diccionario/fastas/Myco02.fasta'
f2='/Users/nairachiclana/Desktop/2.Diccionario/fastas/Myco05.fasta'
k= int(input('Enter your k length: '))

def parseFasta(file):
    fasta_sequences = SeqIO.parse(open(file), 'fasta')
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
        #print(sequence + '\n')
        return(sequence)

seq1=parseFasta(f1)
seq2=parseFasta(f2)

def makeDictionary(k,seq):
        dic = {} #empty dictionary kmer:pos(i)
        kmer=([seq[i:i+k] for i in range(0, len(seq)-k+1, 1)])
        for i in range(1,len(kmer)):
            key = kmer[i]
            dic.setdefault(key, [])
            dic[key].append(i)
        return(dic)

dic1=makeDictionary(k,seq1)
dic2=makeDictionary(k,seq2)

def writeDictionary(dict, fo, sep):
    with open(fo, "a") as f: #appending mode
        f.writelines('{}:{}'.format(k, v) for k, v in dict.items())
        f.write('\n')
        f.write(i + " " + sep.join([str(x) for x in dict[i]]) + "\n")

def ComputeHits(d1,d2):
    matches=list()
    for key in d1:
        list1=d1.get(key) #posiciones para cada kmer (key)
        list2=d2.get(key)
        s1=set(list1) #posiciones no repetidas
        s2=set(list2)
        match =s1.intersection(s2) #para cada kmero, pos dic1 y dic2
        matches.extend(match)
    return(matches) #para todos los kmeros

matchesTotal=ComputeHits(dic1,dic2)
print("Hits:",matchesTotal)
print("Total Hits:",len(matchesTotal))





