
############################################################################################################################################
# README:
# pour extraire et trier selon la taille les sequences proteiques a partir d'un fichier
# fasta,il prend comme argument le fichier multifasta contenant les proteines d'un genome
# YOU NEED BLASTP IN A STANDALONE FORMAT INSTALLED IN YOU CAMPUTER!!!!
# E-value: The Expect value (E) is a parameter that describes the number of hits one can "expect" to see by chance when searching a database of a particular size. 
# commande format: ./Extract_prot_fasta_pseudogenes_out.py <multifasta.fa> -pseudo_out <E-value>
filename_3=sys.argv[1]
pseudo_out=sys.argv[2]
E_value=sys.argv[3]


	
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
