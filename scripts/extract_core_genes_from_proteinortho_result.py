import os
import csv
import random
from Bio import SeqIO



#please note that Bakta was used to generate both the .faa files used with proteinortho
#and the gbff files located in input_path_bakta  

# Define your INPUT #
name_of_pr_ortho_file = "name_of_pr_ortho_file.proteinortho.tsv"
input_path = "/input_path_to_proteinortho_folder/"
input_path_bakta = "/input_path_to_bakta_annotated_sequences_gbff/"
z = #add a number  # number of genomes/sequences used in the proteinortho 


output_path_RAW = "/path_to_some_remote_folder_on_your_pc/coregenes_locustags/"
final_output_path = "/path_to_some_remote_folder_on_your_pc/coregenes/"
alignment_path = "/path_to_some_remote_folder_on_your_pc/coregenes_alignment_prank/"
cleaned_alignment_path = "//path_to_some_remote_folder_on_your_pc/coregenes_alignment_prank_cleaned/"

for folder in [output_path_RAW, final_output_path, alignment_path, cleaned_alignment_path]:
    os.makedirs(folder, exist_ok=True)
    

# parse genome/sequence list from the header #
with open(os.path.join(input_path, name_of_pr_ortho_file)) as f:
    proteinortho_results = csv.reader(f, delimiter='\t')
    header = next(proteinortho_results)
    ICElist = [name.replace(".faa", "") for name in header[3:3+z]]
    

# Extract core genes locus tag #
with open(os.path.join(input_path, name_of_pr_ortho_file)) as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)  # skip header
    for row in reader:
        if int(row[0]) == z and int(row[1]) == z:
            locus_tags = row[3:3+z]
            with open(os.path.join(output_path_RAW, row[3] + '_singlegene.txt'), 'w') as out1:
                out1.write(','.join(locus_tags))

print("Core gene locus tags extracted")


# Extract core gene sequences from the gbff files"

num = 1
for filename in os.listdir(output_path_RAW):
    gene_path = os.path.join(output_path_RAW, filename)
    with open(gene_path, 'r') as f:
        coregenes = f.read().strip().split(',')

    gene_output_file = os.path.join(final_output_path, filename + "_coregenes.fasta")
    with open(gene_output_file, "w") as ofile:
        print(f"🔍 Extracting gene {num} from genomes...")
        num += 1
        for element in ICElist:
            ofile.write(f">{element}\n")
            gbk_file = os.path.join(input_path_prokka, f"{element}.gbff")
            found = False

            for record in SeqIO.parse(gbk_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        locustag = feature.qualifiers.get("locus_tag", [""])[0]
                        if locustag in coregenes:
                            DNAseq = feature.extract(record.seq)
                            ofile.write(str(DNAseq) + "\n")
                            found = True
                            break
                if found:
                    break

print("Core gene sequences extracted.")



# align with prank all the core genes #
for filename in os.listdir(final_output_path):
   
    if filename.endswith(".fasta"):
        name=filename.replace('.fasta', '')
        command_prank='prank -d={}{} -o=alignment_path/al_{}.fasta'.format(final_ouput_path,filename,name)
        os.system(command_prank)


# remove gaps #
working_directory='alignment_path/'
for filename in os.listdir(working_directory):
    name=filename.replace('.fasta', '')
    print("cleaning    " + name)
    command='goalign clean sites -i {}{} > cleaned_alignment_path/{}.fasta'.format(working_directory,filename,name)
    os.system(command)


# concatenate the alignments 
command='goalign concat -i cleaned_alignment_path/*.fasta -o PRANK_clean_concatenated_core_aln.fasta'
os.system(command)


##this alignment can be used now to build a tree, we used fasttree