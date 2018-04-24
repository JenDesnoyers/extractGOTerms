#Hello, my name is Jennifer Densoyers. This is my script for BINF6410 Assignment #1. The goal of the program is to find GO Terms for differentially expressed genes in Species X after treatment. 
#However, Species X doesn't have any assigned GO Terms, so a closely related species (Sorghum bicolor), which has GO Terms, will be used.   

#The following code line takes the 'speciesx_genes.fa' file, which is a fasta file, and makes each fasta entry one line. The reults are then writen into the 'speciesx_genes_oneline.fa' file. 
#It is important for the file to be one line per entry in order to use the grep command later on. Grep looks at one line at a time, so if an entry is split between multiple lines grep may miss
#entries or return more than there actually are (grep would treat a sequence spread across 4 lines as 4 separate sequences).

sed 's/\(^>.*$\)/#\1#/' speciesx_genes.fa | tr -d "\n" | sed 's/$/#/' | tr '#' "\n" | sed '/^$/d' > speciesx_genes_oneline.fa 

#Again, the following line takes the 'Sorghum_bicolor.Sorbi1.25.cds.all.fa' file and makes the file one line per entry, and then writes the results to the file called 
#'Sorghum_bicolor.Sorbi1.25.cds.all.oneline.fa'.

sed 's/\(^>.*$\)/#\1#/' Sorghum_bicolor.Sorbi1.25.cds.all.fa | tr -d "\n" | sed 's/$/#/' | tr '#' "\n" | sed '/^$/d' > Sorghum_bicolor.Sorbi1.25.cds.all.oneline.fa 

#The next code line uses grep, which searches for a provided pattern, in order to find the sequences (from the 'speciesx_genes_oneline.fa') that match the differentially expressed gene IDs
#(from the 'dif_expressed.txt'). The 'pattern' that grep searches for is the differentially expressed gene IDs in the 'dif_expressed.txt' file, and searches for this pattern in the
#'speciesx_genes_oneline.fa' file.The result from this command are then written to the file called 'dif_expressed_speciesx_seq.txt'. -A1 makes grep return the entire line that matches the pattern. Without this parameter, grep would just return ~ and not the sequence that matches the ID
#in the 'dif_expressed.txt' file. -f lets grep know that it will be looking at files.  

grep -A1 -f dif_expressed.txt speciesx_genes_oneline.fa > dif_expressed_speciesx_seq.txt

#The next code line loads the BLAST module onto Sharcnet so thatBLAST searches can be performed. 

module load blast/2.2.28+

#The next code line creates a database using the 'Sorghum_bicolor.Sorbi1.25.cds.all.oneline.fa' so that the Species X differentially expressed gene sequences can be BLASTed against the Sorghum
#sequences in order to find the closest match for each differentially expressed gene. The matched Sorghum sequences will then be used to find GO Terms for Sorghum, which should also apply to 
#Species X because the two species are closely related.

makeblastdb -in Sorghum_bicolor.Sorbi1.25.cds.all.oneline.fa -dbtype nucl

#The next code line performs a BLAST search. It BLASTs the sequences in the 'dif_expressed_speciesx_seq.txt' file against the database created above in order to find Sorghum sequenes that are similar 
#to the Species X sequences. The -outfmt 6 paramter returns the result in tabular form so that command like cut, which depend on delimiters such as tab, can be used. 
#The '| sort -k1,1 -k12,12nr -k11,11 | sort -u -k1,1 -- merge' portion retrieves only the top hit for each BLAST search. The results are then written into the 'tophitsBLAST.txt' file. 

blastn -db Sorghum_bicolor.Sorbi1.25.cds.all.oneline.fa -query dif_expressed_speciesx_seq.txt -outfmt 6 | sort -k1,1 -k12,12nr -k11,11 | sort -u -k1,1 -- merge > tophitsBLAST.txt  

#The next code line cuts the 'tophitsBLAST.txt' file so that only the 2nd column remains (the column with Sorghum IDs) and then writes the cut results to the file named 
#'tophitsBLAST_onlySorghumGeneIDs.txt' so that the original 'tophitsBLAST.txt' file isn't altered. This line of code is necessary because the GO_analysis.R program will not work with the file 
#'tophitsBLAST.txt' file, it only requires the Sorghum IDs from column #2.

cut -f2 tophitsBLAST.txt > tophitsBLAST_onlySorghumGeneIds.txt

#The next code lines load the require modules in order to use R in Unix.
module unload r intel gcc
module load r

#The next code line runs the script named GO_analysis_no_download.R, which will identify GO Terms that are highly represented by the Sorghum genes that relate to the differentially expressed 
#Species X genes. The results from this GO analysis will be written to the file named 'GOTerms_Sorghum_RelatedToSpeciesX.txt'. 

Rscript GO_analysis_no_download.R  

