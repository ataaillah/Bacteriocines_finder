#! /usr/bin/env python
#-------------------------------------------------
#README : 
#ce script extrait les sequences codantes pour des peptide entre 10 et 60 aa en
#prenant en arguement 1: un fichier .gbk(output de BactgeneSHOW par ex) et en
#argument 2: le fichier .fasta correspondant au genome d'interet
#-------------------------------------------------
import time
tps1=time.clock()
import re
import sys
import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
filename_1=sys.argv[1]
filename_2=sys.argv[2]

filin=open(filename_1,'r')
CDS_plus_regex=re.compile("^ *(CDS) *([0-9]+)\.\.([0-9]+)")
CDS_minus_regex=re.compile("^ *(CDS) *(complement)\(([0-9]+)\.\.([0-9]+)\)")
seq_regex=re.compile("([AGTC]{,70})")
CDS_plus_list=[]
CDS_minus_list=[]
pep_CDS_plus_list=[]
pep_CDS_minus_list=[]
num_cds=0
size_cds_plus=0
size_cds_minus=0
num_pep_cds_plus=0
num_pep_cds_minus=0
sequence=""
fasta_width=80
tps2=0
tps1=0

# mettre dans une liste les CDS, debut, fin:
for i in filin:
	if CDS_plus_regex.search(i):
		res_CDS_plus_list=CDS_plus_regex.search(i)
		CDS_plus_list.append(res_CDS_plus_list.group(1))
		CDS_plus_list.append(int(res_CDS_plus_list.group(2)))
		CDS_plus_list.append(int(res_CDS_plus_list.group(3)))
	if CDS_minus_regex.search(i):
		res_CDS_minus_list=CDS_minus_regex.search(i)
		CDS_minus_list.append(res_CDS_minus_list.group(1))
		CDS_minus_list.append(int(res_CDS_minus_list.group(3)))
		CDS_minus_list.append(int(res_CDS_minus_list.group(4)))
for i in xrange(0,len(CDS_plus_list),3):
	size_cds_plus=CDS_plus_list[i+2]-CDS_plus_list[i+1]
	if  30 < size_cds_plus <= 180:
		num_pep_cds_plus+=1
		pep_CDS_plus_list.append(CDS_plus_list[i])
		pep_CDS_plus_list.append(CDS_plus_list[i+1])
		pep_CDS_plus_list.append(CDS_plus_list[i+2])
for i in xrange(0,len(CDS_minus_list),3):
	size_cds_minus=CDS_minus_list[i+2]-CDS_minus_list[i+1]
	if  30 < size_cds_minus <= 180:
		num_pep_cds_minus+=1
		pep_CDS_minus_list.append(CDS_minus_list[i])
		pep_CDS_minus_list.append(CDS_minus_list[i+1])
		pep_CDS_minus_list.append(CDS_minus_list[i+2])	
filin.close()
filin=open(filename_2,'r')
for i in filin:
	if seq_regex.search(i):
		seq_line=seq_regex.search(i)
		seq=seq_line.group(0)
		sequence+=seq
filout=open('pep_CDS_translated_'+filename_2,'w')
for i in xrange(0,len(pep_CDS_plus_list),3):
	filout.write(">pep"+str(pep_CDS_plus_list[i+1])+"|"+str(pep_CDS_plus_list[i+2])+"|"+str(i)+ " pep_CDS_translated"+"\n")
	fasta_seq=""
	fasta_seq=sequence[pep_CDS_plus_list[i+1]-1:pep_CDS_plus_list[i+2]-1]
	fasta_seq_DNA=Seq(fasta_seq, IUPAC.unambiguous_dna)
	fasta_seq_DNA_translated=fasta_seq_DNA.translate()
	filout.write("\n".join([str(fasta_seq_DNA_translated)[i:i+fasta_width] for i in xrange(0,len(fasta_seq_DNA_translated),fasta_width)])+ "\n")
for i in xrange(0,len(pep_CDS_minus_list),3):
	filout.write(">pep"+str(pep_CDS_minus_list[i+1])+"|"+str(pep_CDS_minus_list[i+2])+"|"+str(i)+" pep_CDS_mins_reverse_complement_translated" +"\n")
	fasta_seq=""
	fasta_seq=sequence[pep_CDS_minus_list[i+1]:pep_CDS_minus_list[i+2]]
	fasta_seq_DNA=Seq(fasta_seq, IUPAC.unambiguous_dna)
	fasta_seq_DNA_rev_compl=fasta_seq_DNA.reverse_complement()
	fasta_seq_DNA_rev_compl_translated=fasta_seq_DNA_rev_compl.translate()
	filout.write("\n".join([str(fasta_seq_DNA_rev_compl_translated)[i:i+fasta_width] for i in xrange(0,len(str(fasta_seq_DNA_rev_compl_translated)),fasta_width)])+ "\n")
	print "sequence"
	print fasta_seq


############################################################################################################################################
# README:
# pour extraire et trier selon la taille les sequences proteiques a partir d'un fichier fasta

filename_3=sys.argv[3]
pseudo_out=sys.argv[4]
E_value=sys.argv[5]


	
filin=open(filename_3,'r')

prot_id_regex=re.compile("location=(complement)?\(?([0-9]+)\.\.([0-9]+)\)?")
prot_seq_regex=re.compile("^[ARNDCEQGHILKMFPSTWYV]+")
prot_id_list=[]
list_sequence=[]
num_prot=0
seq=""
sequence=""
fasta_seq=""
fasta_width=80
#recuperer les coordonnes l'ID et nom des proteines:

for i in filin:
	if prot_id_regex.search(i):
		res_prot_id_regex=prot_id_regex.search(i)
		prot_id_list.append(int(res_prot_id_regex.group(2)))
		prot_id_list.append(int(res_prot_id_regex.group(3)))

# extraire et concatener les sequences proteiques correspondantes au genes:
		
		
	if prot_seq_regex.search(i):
		res_prot_seq_regex=prot_seq_regex.search(i)
		seq=res_prot_seq_regex.group(0)
		sequence+=seq
#transformer la chaine de caractere en liste de caracteres
for i in sequence:
	list_sequence.append(i)

print sequence
print len(sequence)
print list_sequence
print len(list_sequence)
print prot_id_list
#recuperer les sequences proteiques superieurs a 100 aa
filout=open(filename_3+'_+q_100aa','w')
for i in xrange(0,len(prot_id_list),2):
	prot_size=((prot_id_list[i+1]-prot_id_list[i])-2)/3
	if  prot_size >= 100:
		list_prot=list_sequence[0:((prot_id_list[i+1]-prot_id_list[i])-2)/3]
		fasta_seq="".join(list_prot)
		print fasta_seq
		del list_sequence[0:((prot_id_list[i+1]-prot_id_list[i])-2)/3]
		filout.write(">"+str(prot_id_list[i])+".."+str(prot_id_list[i+1])+ " prot_plus_que_100aa" + "\n")
		filout.write("\n".join([fasta_seq[i:i+fasta_width] for i in xrange(0,len(fasta_seq),fasta_width)])+ "\n")
	else:
		seq_prot=list_sequence[0:((prot_id_list[i+1]-prot_id_list[i])-2)/3]
		del list_sequence[0:((prot_id_list[i+1]-prot_id_list[i])-2)/3]

filout.close()
filin.close()
#######################################################################################################################
if pseudo_out=="-pseudo_out":
	import os
	os.system("makeblastdb -in "+filename_3+" -parse_seqids -dbtype prot -title "+filename_3)
	output_blast="BLASTP_predict_pep_vs_annotated_prot_+q100aa"+filename_2+E_value+".txt"
	input_blast="pep_CDS_translated_"+filename_2
	os.system("blastp -query "+input_blast+" -db "+filename_3+" -task blastp -outfmt 7 -evalue "+E_value+" >  "+output_blast)

#######################################################################################################################


#-------------------------------------------------
#README : 
#ce script extrait les query d'un output BLAST ayant 0 hits found, il prend comme premier argument le output de BLAST et comme deuxieme argument le input de BLAST pour produire un fichier contenant les query ayant 0 hits found et leurs sequences
#-------------------------------------------------

	sequence=""
	debut_fin_regex=re.compile("# Query: pep([0-9]+)\.{0,2}\|?([0-9]+)")
	prot_seq_regex=re.compile("^[ARNDCEQGHILKMFPSTWYV]+")
	filename_1=output_blast
	debut_fin_list=[]
	filename_2=input_blast
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
	output_blast_non_pseudogene=filename_1+'_+NON_pseudogene'
	filout=open(output_blast_non_pseudogene,'w')
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
######################################################################################################
intergenic_ORF=sys.argv[6]
if intergenic_ORF == "-intergenic_ORF":
#-------------------------------------------------
#README : 
# extrait les sequences codantes predites par BactgeneSHOW, qui font parties de l'espace intergenique du genome d'interet .
#il prend en premier argument, le fichier .gbk du genome d'interet et en deuxieme argument le fichier contenant les genes
#a verifier s'il sont intergeniques, il prend en 3 eme argument le nombre de paire de base tolere a l'interieur des region
#intergenique.
#-------------------------------------------------


	filename_1=sys.argv[7]
	filename_2=output_blast_non_pseudogene
	X=sys.argv[8]

	gene_plus_minus_regex=re.compile("^ *gene *(complement)?\(?([0-9]+)\.\.([0-9]+)\)?")
	my_regex=re.compile("^>([0-9]+)\.\.([0-9]+)")
	my_list_intergenic=[]
	my_list=[]
	list_sequence=[]
	gene_list_plus_minus=[]




	filin=open(filename_1,'r')
	for i in filin:
		if gene_plus_minus_regex.search(i):
			res_gene_plus_minus=gene_plus_minus_regex.search(i)
			gene_list_plus_minus.append(int(res_gene_plus_minus.group(2)))
			gene_list_plus_minus.append(int(res_gene_plus_minus.group(3)))
	

	filin.close()

	filin=open(filename_2,'r')
	for i in filin:
		if my_regex.search(i):
			res_my_regex=my_regex.search(i)
			my_list.append(int(res_my_regex.group(1)))
			my_list.append(int(res_my_regex.group(2)))
		
	filin.close()

	for i in xrange(0,len(my_list),2):
		for j in xrange(1,len(gene_list_plus_minus),2):
			if (gene_list_plus_minus[j]-int(X))<= my_list[i] <= (gene_list_plus_minus[j+1]+int(X)) and (gene_list_plus_minus[j]-int(X))<= my_list[i+1] <=(gene_list_plus_minus[j+1]+int(X)):
				my_list_intergenic.append(my_list[i])
				my_list_intergenic.append(my_list[i+1])
		



	filout=open('intergenic_CDS_'+filename_2,'w')
	filin=open(filename_2,'r')
	lignes=filin.readlines()
	for i in xrange(0,len(lignes),2):
		for j in xrange(0,len(my_list_intergenic),2):
			if lignes[i].find(str(my_list_intergenic[j]))!= -1:
				filout.write(lignes[i])
				filout.write(lignes[i+1])
	filin.close()
	filout.close()

