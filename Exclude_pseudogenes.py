#! /usr/bin/env python
#-------------------------------------------------
#README : 
#ce script extrait les query d'un output BLAST ayant 0 hits found, il prend comme premier argument le output
#de BLAST et comme deuxieme argument le input de BLAST pour produire un fichier contenant les query ayant
#0 hits found et leurs sequences
#-------------------------------------------------
import re
import sys
sequence=""
debut_fin_regex=re.compile("# Query: pep([0-9]+)\.{0,2}\|?([0-9]+)")
prot_seq_regex=re.compile("^[ARNDCEQGHILKMFPSTWYV]+")
filename_1=sys.argv[1]
debut_fin_list=[]
filename_2=sys.argv[2]
list_sequence=[]
fasta_width=80
fasta_seq=""

# extraire et concatener les sequences proteiques correspondantes au peptides:
filin_2=open(filename_2, 'r')
for i in filin_2:
	if prot_seq_regex.search(i):
		res_prot_seq_regex=prot_seq_regex.search(i)
		seq=res_prot_seq_regex.group(0)
		sequence+=seq
#transformer la chaine de caractere en liste de caracteres
for i in sequence:
	list_sequence.append(i)

#recuperer les debut fin de peptides:
filin_1=open(filename_1, 'r')
lignes=filin_1.readlines()
for i in xrange(len(lignes)):
	if debut_fin_regex.search(lignes[i]):
		res_debut_fin_regex=debut_fin_regex.search(lignes[i])
		debut_fin_list.append(int(res_debut_fin_regex.group(1)))
		debut_fin_list.append(int(res_debut_fin_regex.group(2)))
filin_1.close()
#recuperer les sequences proteiques ayant 0 hits found
filout=open(filename_1+'_+NON_pseudogene','w')
filin_1=open(filename_1, 'r')
lignes=filin_1.readlines()
for i in xrange(len(lignes)):
	prot_size=0
	if debut_fin_regex.search(lignes[i]) and lignes[i+2][2:] == "0 hits found\n":
		res_debut_fin_regex=debut_fin_regex.search(lignes[i])
		prot_size=(int(res_debut_fin_regex.group(2))-int(res_debut_fin_regex.group(1))-2)/3
		list_prot=list_sequence[0:prot_size]
		
		fasta_seq="".join(list_prot)
		
		del list_sequence[0:prot_size]
		filout.write(">"+str(res_debut_fin_regex.group(1))+".."+str(res_debut_fin_regex.group(2))+ " NON_pseudogene" + "\n")
		filout.write("\n".join([fasta_seq[i:i+fasta_width] for i in xrange(0,len(fasta_seq),fasta_width)])+ "\n")
	elif debut_fin_regex.search(lignes[i]) and lignes[i+2][2:8] == "Fields": 
		res_debut_fin_regex=debut_fin_regex.search(lignes[i])
		prot_size=(int(res_debut_fin_regex.group(2))-int(res_debut_fin_regex.group(1))-2)/3
		del list_sequence[0:prot_size]
		
filin_2.close()
print debut_fin_list
print sequence


filin_2.close()
filout.close()
