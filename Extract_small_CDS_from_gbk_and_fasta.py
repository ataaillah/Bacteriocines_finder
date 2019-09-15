#! /usr/bin/env python
#-------------------------------------------------
#README : 
#ce script extrait les sequences CDS codant pour des peptide entre 10 et 80 aa en prenant
#en arguement 1: un fichier .gbk(output de BactgeneSHOW par ex) et en argument 2:
#le fichier .fasta correspondant
#-------------------------------------------------
import time
tps1=time.clock()
import re
import sys
filename_1=sys.argv[1]
filename_2=sys.argv[2]
#if len(sys.argv) == 1 :
#	sys.exit("ERREUR: il faut exactement un argument !")
filin=open(filename_1,'r')
CDS_regex=re.compile("^ *(CDS) *(complement)?\(?([0-9]+)\.\.([0-9]+)\)?")
seq_regex=re.compile("[AGTC]{,70}")
CDS_list=[]
pep_CDS_list=[]
num_cds=0
num_pep_cds=0
sequence=""
fasta_width=80
tps2=0
tps1=0
#les fonctions utilisee:
def find_CDS(i):
	if CDS_regex.search(i):
		res_CDS_list=CDS_regex.search(i)
		CDS_list.append(res_CDS_list.group(1))
		CDS_list.append(int(res_CDS_list.group(3)))
		CDS_list.append(int(res_CDS_list.group(4)))
	return
def find_pep_CDS(i):
	size_cds=0
	size_cds=CDS_list[i+2]-CDS_list[i+1]
	if  30 < size_cds <= 240:
		global num_pep_cds
		num_pep_cds+=1
		pep_CDS_list.append(CDS_list[i])
		pep_CDS_list.append(CDS_list[i+1])
		pep_CDS_list.append(CDS_list[i+2])
	return
def extract_sequence(i):
	if seq_regex.search(i):
		seq_line=seq_regex.search(i)
		seq=seq_line.group(0)
		global sequence
		sequence+=seq
	return
def extract_pep_CDS(i):
	fasta_seq=""
	fasta_seq=sequence[pep_CDS_list[i+1]-1:pep_CDS_list[i+2]-1]
	filout.write("\n".join([fasta_seq[i:i+fasta_width] for i in xrange(0,len(fasta_seq),fasta_width)])+ "\n")
	return
# mettre dans une liste les CDS, debut, fin:
for i in filin:
	find_CDS(i)
for i in xrange(0,len(CDS_list),3):
	find_pep_CDS(i)	
filin.close()
filin=open(filename_2,'r')
for i in filin:
	extract_sequence(i)
filout=open('pep_CDS_'+filename_2,'w')
for i in xrange(0,len(pep_CDS_list),3):
	filout.write("> CDS  "+str(pep_CDS_list[i+1])+".."+str(pep_CDS_list[i+2]) + "\n")
	extract_pep_CDS(i)
tps2=time.clock()
print "le temps (en seconds)  d'execution de ce script est:"
print (tps2-tps1)
print "le nombre de CDS correspondants au peptides entre 10 et 80 aa est:"
print num_pep_cds
