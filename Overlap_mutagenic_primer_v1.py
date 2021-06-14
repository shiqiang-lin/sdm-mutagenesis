#usage:
#python3.7 Overlap_mutagenic_primer_v1.py DNAsequence.txt mutaa.txt species
#
#species can be one of the following Ec, Yeast, Insect,
#Ce, Dm, Human, Mouse, Rat, Pig, Pp, At, Streptomyces,
#Zm, Tobacco, Sc, Cg.

"""
This script searches possible mutagenc primers for overlap extension PCR.
The resulting overlap primers are sorted according to the theoretical Tm value of
the "overlap" sequence during the overlapping stage. Also, the length, GC content
and Tm value of each primer are provided to help designing the gene primers and
set annealing temperature for the PCR run. We hope these listed primers are helpful
for preparing mutagenic primers with overlap extension PCR.
In addition, the Overlap_mutagenic_primer_v1.py and the CodonPreference.txt should be
in the same folder for proper running of the program.
"""

import sys
import os,shutil
import string
from operator import itemgetter
from Bio.Seq import Seq
from Bio.SeqUtils import GC,MeltingTemp

#get current path
current_path = os.getcwd()
print("Current path is %s." % current_path)
print("\n")

#read gene file
DNA_file = open(sys.argv[1])
DNA_file_sequence = DNA_file.read()
DNA_file.close()

#filter char and print gene sequence
DNA_sequence_list = []
for i in DNA_file_sequence:
    if i.isalpha():
        i = i.upper()
        DNA_sequence_list.append(i)
gene_sequence_str = ''.join(DNA_sequence_list)

print("The gene sequence is")   
for j in range(0,len(gene_sequence_str),100):
    DNA_string_100_per_line = gene_sequence_str[j:j+100]
    print(DNA_string_100_per_line)
print("The lenght of gene sequence is %s bps." % len(gene_sequence_str))
print("\n")

#see if gene length is enough for running
if len(gene_sequence_str) < 399:
    print("The gene sequence is less than 399 and has not enough length for the script to run.\n")
    print("In this situation, the whole plasmid PCR mutagenesis might be utilized.\n")
    print("Another way is to include the DNA sequences upstream and downstream of gene and recalculate\n")
    print("the position of mutagenic amino acids. Thank you!\n")
    sys.exit()
    
    
"""
get file of codon preference table and extract species of interest to species_codon_preference_table_list[]
each item is a str with length 11 like this 'GGG\tG\t0.21\n', can be extracted according to the 5th char
"""

species = sys.argv[3]
species_str = species + ","
codon_preference_table_list = []

codon_file = open("CodonPreference.txt")
codon_file_lines = codon_file.readlines()
codon_file.close()

for line in codon_file_lines:
    if species_str in line:
        species_line_index = codon_file_lines.index(line)
        species_codon_preference_table_list = codon_file_lines[species_line_index+2:species_line_index+66]

        print("The species selected is %s" % codon_file_lines[species_line_index])
        print("The codon preference table is as follows")
        for i in species_codon_preference_table_list:
            i = i.strip()
            print(i)
            
if len(species_codon_preference_table_list) == 0:
    print("Species not found. Please rerun with appropriate species. Thank you!")
    sys.exit()
    

#get amino acid mutations
amino_acid_mutation_file = open(sys.argv[2])
amino_acid_mutation_file_lines = amino_acid_mutation_file.readlines()
amino_acid_mutation_file.close()

for amino_acid_mutation in amino_acid_mutation_file_lines:
    amino_acid_mutation = amino_acid_mutation.strip()
    amino_acid_mutation_line_length = len(amino_acid_mutation)
    orignal_amino_acid = amino_acid_mutation[0]
    new_amino_acid = amino_acid_mutation[amino_acid_mutation_line_length-1]
    amino_acid_position = int(str(amino_acid_mutation[1:amino_acid_mutation_line_length-1]))

    #check if the position is too near to the ends of protein sequence
    """
    Three conditions need to be satisfied, including
    1) amino_acid_position < 66
    2) amino_acid_position > len(gene)/3-65
    3) len(gene) > 400
    """
    if amino_acid_position < 66 or amino_acid_position > int(len(gene_sequence_str)/3-65):
        print("Amino acid mutation file contains a line with mutation position within the first or last 65.")
        print(amino_acid_mutation)
        print("Mutagenic primers can be directly designed without overlap extension PCR.\n")
        print("Similar message may occur if there is still other such positions.\n")
        print("Please delete them all from the txt file of amino acid mutation and run script again. Thank you!\n")
        sys.exit()
    

    else:
        print("The orginal %s in position %i will be changed to %s." % (orignal_amino_acid, \
                                                                    amino_acid_position, \
                                                                    new_amino_acid))
    

    #get the codon with the highest preference value
    #    format: GGG\tG\t0.21\n
    potential_codons_list = []
    for codon_line in species_codon_preference_table_list:
        codon_line = codon_line.strip()
        if new_amino_acid == codon_line[4]:
            potential_codons_list.append(codon_line)

    mutagnenic_codon_str = potential_codons_list[0]
    if len(potential_codons_list) == 1:
        print("The potential codon is as following")
        print(mutagnenic_codon_str)       
  
    else:
        print("The potential codons are as follows")
        for i in range(len(potential_codons_list)):
            print(potential_codons_list[i])
            if (i >=1) and (float(potential_codons_list[i][6:10]) > float(potential_codons_list[0][6:10])):
                mutagnenic_codon_str = potential_codons_list[i]

    print("The codon used is %s" % mutagnenic_codon_str)
    mutagenic_codon_used_str = mutagnenic_codon_str[0:3]


    """
    get the tetrad primers
    """
    p = amino_acid_position   # int type
    primer_trimers_list = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    overlap_for_left_str = gene_sequence_str[3*p-12-i:3*p-3]
                    overlap_for_right_str = gene_sequence_str[3*p:3*p+15+j]
                    overlap5_Seq = (overlap_for_left_str + \
                                    mutagenic_codon_used_str+ \
                                    overlap_for_right_str)

                    overlap_rev_left_str = gene_sequence_str[3*p-18-l:3*p-3]
                    overlap_rev_right_str = gene_sequence_str[3*p:3*p+9+k]
                    overlap3_Seq = Seq(overlap_rev_left_str + mutagenic_codon_used_str + overlap_rev_right_str).reverse_complement()

                    overlap_Seq = Seq(gene_sequence_str[3*p-12-i:3*p-3] + mutagenic_codon_used_str + gene_sequence_str[3*p:3*p+9+k])

                    Tm_overlap5 = round(MeltingTemp.Tm_NN(overlap5_Seq),2)
                    Tm_overlap3 = round(MeltingTemp.Tm_NN(overlap3_Seq),2)
                    Tm_overlap = round(MeltingTemp.Tm_NN(overlap_Seq),2)

                    primer_trimer_list = []
                    primer_trimer_list.append(overlap5_Seq)
                    primer_trimer_list.append(overlap3_Seq)
                    primer_trimer_list.append(Tm_overlap)
                    
                    primer_trimers_list.append(primer_trimer_list)
                    primer_trimer_list = []
                                
 
    
        
    #sort primer_sixers_list by K value first then by the first element.
    primer_trimers_list_sorted = sorted(primer_trimers_list, key = itemgetter(2,0))

    #print the original primer list
    print("The number of original primers is %s." %len(primer_trimers_list))
    print("The original primers are shown in the following")
    for i in range(len(primer_trimers_list)):
        for j in range(3):
            print(primer_trimers_list[i][j])
        
    #print the sorted primer list
    print("The number of sorted primers is %s." %len(primer_trimers_list_sorted))
    print("The sorted primers are shown in the following")
    for i in range(len(primer_trimers_list_sorted)):
        for j in range(3):
            print(primer_trimers_list_sorted[i][j])
    
           
    """
    This section print each combination of tetrad primers to txt file, plus necessary parameters of each primer,
    including name, sequence, length, GC, Tm and K value.
    """
    
     
    file_name_str = amino_acid_mutation + "$" + ".txt"
    primer_file = open(file_name_str,'w')
    print("primer".ljust(20,' '),\
          "sequence".ljust(40,' '),\
          "length".ljust(10,' '),\
          "GC".ljust(10,' '),\
          "Tm".ljust(10,' '),\
          "Overlap_Tm".ljust(15,' '),\
          sep="",\
          file=primer_file
          )

    for i in range(len(primer_trimers_list_sorted)):
        for j in range(2):
            print_primer_sequence_str = primer_trimers_list_sorted[i][j]
            
            if j == 0:
                print("\n",file=primer_file)
                S0_print_primer_line_str = amino_acid_mutation + "_overlap_5"
            else:
                S0_print_primer_line_str = amino_acid_mutation + "_overlap_3"

            S1_print_primer_line_str = str(print_primer_sequence_str)
            S2_print_primer_line_str = str(len(S1_print_primer_line_str))
            S3_print_primer_line_str = str(round(GC(S1_print_primer_line_str),2))
            S4_print_primer_line_str = str(round(MeltingTemp.Tm_NN(S1_print_primer_line_str),2))
            S5_print_primer_line_str = str(primer_trimers_list_sorted[i][2])

            print(S0_print_primer_line_str.ljust(20,' '),\
                  S1_print_primer_line_str.ljust(40,' '),\
                  S2_print_primer_line_str.ljust(10,' '),\
                  S3_print_primer_line_str.ljust(10,' '),\
                  S4_print_primer_line_str.ljust(10,' '),\
                  S5_print_primer_line_str.ljust(15,' '),\
                  sep="",\
                  file=primer_file
                  )
            
    primer_file.close()
    

gene_file = sys.argv[1]
folder_name = gene_file.split(".")[0] + "_primers"
dirs = os.listdir(current_path)

if(folder_name not in dirs):
    os.system('mkdir temp_foldername')
    os.system('mv *$.txt temp_foldername')
    os.rename("temp_foldername", folder_name)

else:
    shutil.rmtree(folder_name)
    os.system('mkdir temp_foldername')
    os.system('mv *$.txt temp_foldername')
    os.rename("temp_foldername", folder_name)









