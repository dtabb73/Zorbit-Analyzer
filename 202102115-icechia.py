import numpy as np
import pandas as pd
from io import StringIO
import re
import csv
from csv import reader, writer
import sys
import os
import glob
import fnmatch
from os import path
import matplotlib
from matplotlib import pyplot as plt

print("You are using Zorbit Analyzer v0.1")
directory_path = input("Please enter the path to the directory of your files. All files should be in the same location: ") #Asks users for path
os.chdir(directory_path)
x = input('Input your Interproscan output gff3 file(s):') #Asks users for gff3 input
if "*" in x: #Handles the case of *.gff3
    gff3_input = glob.glob("*.gff3")
else: 
    y = re.sub('[|; ]', ', ', x) #Substitutes possible gff3 file delimeters with commas
    gff3_input = re.split(', ', y) #Splits gff3 input into a list
for i in gff3_input:
    if os.path.exists(i): #Checks existence of gff3 file
        pass
    else:
        print("There does not seem to be a file by that name. Please check your path/filename and try again")
        sys.exit()
fasta_input = input('Input your fasta file:') #Asks users for fasta input file
if os.path.exists(fasta_input): #Checks existence of fasta input file
    pass
else:
    print("There does not seem to be a file by that name. Please check your path/filename and try again")
    sys.exit()
if fnmatch.fnmatch(fasta_input, '*fastq*'):
    print("Zorbit Analyzer is not specifically constructed to handle fastq files but will try. If errors convert to fasta format")
ortho_input = input ('Input your ProteinOrtho output file:') #Asks users for ProteinOrtho input
if os.path.exists(ortho_input): #Checks existence of ProteinOrtho input
    pass
else:
    print("There does not seem to be a file by that name. Please check your path/filename and try again")
    sys.exit()
ortho_input_file_name = input ('Input your ProteinOrtho input file name (faa). Leave blank if unknown though will run slower:') #Asks users for ProteinOrtho output file
while True: 
    file_to_write = input('Input your desired ZorbitAnalyzer output file name: ') #Asks users for output file
    if file_to_write != '': #Checks to see if user entered a file name
        break
    else:
        print("You did not enter an output file name") #Repeatedly asks for output file name if not given
        continue
Choice = ['yes', 'y', 'no', 'n']
flag = True
while flag is True:
    exclusion_flag = input("Would you like to exclude sequences that do not have either Interproscan or ProteinOrtho hits? (Yes/No) ").lower()
    for i in Choice:
        if exclusion_flag.startswith(i):
            flag = False
            break
        else: 
            continue
if exclusion_flag.startswith('y'):
    exclusion_flag = 1
else:
    exclusion_flag = 0
print("Analyzing files") #Lets user know input portion has completed

pdortho = pd.read_csv(ortho_input, "/t", engine="python") #Creates ProteinOrtho pd
test_file = 'test.txt'
test2_file = 'test2.txt'
test3_file = 'test3.txt'

#Testing open/closing files
def try_file(input_file): #Defining function that creates/opens user output file and truncates it before closing it
    try: 
        open(input_file, 'w+').close()
    except IOError:
        print("Unable to open output file")
try_file('file_to_write.txt') #Creates/opens output file and truncates it before closing it
try_file('test.txt') #Creates/opens test file and truncates it before closing it
try_file('gff3_file_to_write.txt') #Creates/opens gff3 output file and truncates it before closing it
try_file('gff3_statsfile_to_write.txt') #Creates/opens gff3 output file and truncates it before closing i
try_file('fasta_file_to_write.txt') #Creates/opens fasta output file and truncates it before closing it
try_file('ortho_file_to_write.txt')  #Creates/opens ProteinOrtho output file and truncates it before closing it
try_file('ortho_file_to_write2.txt') #Creates/opens a second ProteinOrtho output file and truncates it before closing it
try_file('zorbit_statistics.txt') #Creates/opens a statistics file and truncates it before closing it

#Defining variables for later use
fasta_file_to_write = 'fasta_file_to_write.txt' #Defining the interim fasta file to write
gff3_file_to_write = 'gff3_file_to_write.txt' #Defining the interim gff3 file to write
gff3_statsfile_to_write = 'gff3_statsfile_to_write.txt'
ortho_file_to_write = 'ortho_file_to_write.txt' #Defining the interim Protein Ortho file to write
zorbit_statistics = 'zorbit_statistics.txt' #Defining the Zorbit Statistics variable
string_to_remove1 = '##' #Removes header and gene introduction lines
string_to_remove2 = 'polypeptide' #Removes redundant polypeptide line
string_to_remove3 = 'MobiDBLite' #Removes results from MobiDBLite database
string_to_end = '##FASTA' #Sets end of file as the start of the fasta/code part of gff3 files

#fasta
fasta_file = None
fastq_file = None
fasta_type = "amino_acid"
fastq_start_character = '@'
fasta_start_character = '>' #Setting start character for fasta information line
fastq_third_line_character ='+'
fna_type = "fna"
if fna_type in fasta_input:
    fasta_type = "nucleotide"
with open(fasta_input, 'r') as fasta: #Opening fasta input file to read
    for line in fasta: #reading lines in fasta file
        if line.startswith(fasta_start_character): #Altering lines with > but not sequence lines
            fasta_file = fasta_input
            break
        elif line.startswith(fastq_start_character): #Altering lines with @ but not sequence lines (for fastq)
            fastq_file = fasta_input
            fasta_type = "nucleotide"
            break
        else:
            print("The fasta input file does not seem to have typical fasta or fastq format")
            sys.exit()
if fasta_file is not None: #Checking to see if fasta input was fasta file (should not be empty)
    print("Working on fasta file")
    with open(fasta_input, 'r') as fasta: #Opening fasta input file to read
        with open(fasta_file_to_write, 'a') as f: #Opens the output file to append
            for line in fasta: #reading lines in fasta file
                if line.startswith(fasta_start_character): #Altering lines with > but not sequence lines
                    fasta_nostart = re.sub('>', '\n', line) #Removing > symbol and replacing with carriage return from each occurrence
                    fasta_nospace = ', '.join(fasta_nostart.rsplit('\n',1)) #Removes carriage return (before aa or na code) and replaces with comma
                    fasta_csv = ', '.join(fasta_nospace.split(' ',1)) #Removes first space (after Trinity output name) and replaces with comma
                    f.write(fasta_csv) #Writes output to file
                else:
                    if not line.isspace(): #Will not write blank lines
                        sequence_no_carriage = re.sub('\n', '', line) #Removes carriage return from before the sequence data
                        sequence_no_line_break = re.sub('\r', '', sequence_no_carriage) #Removes line break from before the sequence data
                        f.write(sequence_no_line_break) #Writes the sequence line without line breaks or carriage returns
                    else:
                        continue
elif fastq_file is not None: #Checking to see if fasta input was fastq file (should not be empty)
    print("Working on fastq file")
    with open(fasta_input, 'r', encoding="latin-1") as fasta: #Opening fasta input file to read
        with open(fasta_file_to_write, 'a', encoding="latin-1") as f: #Opens the output file to append
            for i, line in enumerate(fasta): #reading lines in fasta file
                if i == 0: # Dealing with first line differently (no line break)
                    fasta_nostart = re.sub('@', '', line) #Removing @ symbol from each occurrence and replaces with nothing
                    fasta_nospace = ', '.join(fasta_nostart.rsplit('\n',1)) #Removes carriage return (before aa or na code) and replaces with comma
                    fasta_csv = ', '.join(fasta_nospace.split(' ',1)) #Removes first space (after Trinity output name) and replaces with comma
                    f.write(fasta_csv) #Writes output to file
                elif line.startswith(fastq_start_character): #Altering lines with @ but not sequence lines (for fastq)
                    fasta_nostart = re.sub('@', '\n', line) #Removing @ symbol from each occurrence and replaces with carriage return
                    fasta_nospace = ', '.join(fasta_nostart.rsplit('\n',1)) #Removes carriage return (before aa or na code) and replaces with comma
                    fasta_csv = ', '.join(fasta_nospace.split(' ',1)) #Removes first space (after Trinity output name) and replaces with comma
                    f.write(fasta_csv) #Writes output to file
                elif i % 4 == 1: #Writing line 2/4 (sequence file) to output file
                    sequence_no_carriage = re.sub('\n', '', line) #Removes carriage return from before the sequence data
                    sequence_no_line_break = re.sub('\r', '', sequence_no_carriage) #Removes line break from before the sequence data
                    f.write(sequence_no_line_break) #Writes the sequence line without line breaks or carriage returns
                else:
                    pass
else:
    print("The input file does not seem to be in typical fasta or fastq format. Please check and try again") #Ending if atypical fasta/fastq format
    sys.exit()

for i in gff3_input: #Cleaning up gff3 file prior to conversion to dataframe
    with open(i, 'r') as stack:
        with open(gff3_file_to_write, 'a') as f:
               for line in stack:
                if string_to_end in line: #Closes file at the start of the sequence data without including
                 f.close()
                 break
                elif string_to_remove1 in line: #Removing header and gene introduction lines (if present)
                    continue
                elif string_to_remove2 in line: #Removing polypeptide line (if present)
                    continue
                elif string_to_remove3 in line: #Removing MobiDBLite database (if present)
                    continue
                else:
                     f.write(line)    
for i in gff3_input: #Saving unedited gff3 input into file for statistics purposes later
    with open(i, 'r') as stack:
        with open(gff3_statsfile_to_write, 'a') as f:
               for line in stack:
                    if string_to_end in line: #Closes file at the start of the sequence data without including
                        f.close()
                        break
                    elif string_to_remove1 in line: #Removing header and gene introduction lines (if present)
                        continue
                    else:
                        f.write(line) 
fasta_column_names = ['SeqID', 'Information', 'Sequence'] #Defining the list of fasta column names to pass to the dataframe
fastapd = pd.read_csv(fasta_file_to_write, names=fasta_column_names, engine = "python", header=None) #Creating a Pandas dataframe from the fasta output csv
SeqID_list = fastapd["SeqID"].tolist() #Saving contents of the SeqID column to a list
fasta_row_number = len(fastapd) #Counting the number of rows in the fasta dataframe for the statistics output
with open(zorbit_statistics, 'a') as f:
    f.write("The number of sequences in the fasta is " + str(fasta_row_number) + "\n")


#Start orthopd
print("Working on ProteinOrtho dataframe")
orthopd = pd.read_csv(ortho_input, sep='\t', engine="python", na_values="*") #Creates a Pandas dataframe from ProteinOrtho input csv
ortho_column_names = list(orthopd.columns)
#Defining the SeqID column
if ortho_input_file_name != "":
    orthopd.columns = ["SeqID" if col.startswith(ortho_input_file_name) else col for col in orthopd.columns] #Renaming the fasta input column in ProteinOrtho dataframe to SeqID to match other dataframes
else: pass
#Attempting to identify which column corresponds to the input fasta
fasta_input_split = fasta_input.split('.', 1)[0] #Trying to delete file handle from the fasta input file in case there was .fasta versus .faa, etc
orthopd_pruned = orthopd.drop(columns=['# Species', 'Genes', 'Alg.-Conn.']) #Creating a new dataframe without the first three columns which will always have data in each row in order to id longest column
if orthopd.columns.astype(str).str.contains("SeqID").any(): #Checking to see if fasta input file name is in the ProteinOrtho column name list
    print("Found fasta Sequence ID column in ProteinOrtho file")
else:
    print("Trying to find fasta file in ProteinOrtho file through other means")
    orthopd.columns = ["SeqID" if col.startswith(fasta_input_split) else col for col in orthopd.columns] #Using the input fasta file name as a guess for the faa file name
    if orthopd.columns.astype(str).str.contains("SeqID").any(): #Breaks loops if the column name has been found/replaced
        print("Found fasta Sequence ID column in ProteinOrtho file")       
    else: 
        print("Attempting another way of identifying fasta file column. This may take some time")
        orthopd_fasta_column_name = orthopd_pruned.count().idxmax() #Finding column with the least number of NaN which is likely the input fasta
        for l in SeqID_list: #Searching to see if any values from the fastapd SeqID column (l) are in the putative SeqID ProteinOrtho column
            if orthopd[orthopd_fasta_column_name].astype(str).str.contains(l).any(): 
                orthopd.rename(columns=lambda x: x.replace(orthopd_fasta_column_name, "SeqID"), inplace=True) #Renaming the ProteinOrtho column with fasta sequence names as SeqID
                break
            else:
                print("Final method to identify fasta file column. This may take hours")
                orthopd = orthopd.drop(orthopd[(orthopd['Genes'] == 1)].index) #Gets rid of rows with just a single gene found in order to speed up full frame search
                for l in SeqID_list: #Searching to see if any values from the fastapd SeqID column (l) are in the ProteinOrtho dataframe
                    for i in orthopd.columns:
                        if orthopd[i].astype(str).str.contains(l).any():
                            orthopd.rename(columns=lambda x: x.replace(i, "SeqID"), inplace=True) #Renaming the ProteinOrtho column with fasta sequence names as SeqID
                            break
orthopd = orthopd.drop(orthopd[(orthopd['SeqID'].isna())].index)#Removing SeqID rows with NaN
#Splitting the duplicated entries in the SeqID column and making new rows with a SeqID member on each but with same data otherwise
def pir2(df, c): #Defining function to split the SeqID column at each comma and place one of each split value onto a new, otherwise duplicated row
    colc = df[c].astype(str).str.split(',')
    clst = colc.values.astype(object).tolist()
    lens = [len(l) for l in clst]
    j = df.columns.get_loc(c)
    v = df.values
    n, m = v.shape
    r = np.arange(n).repeat(lens)
    return pd.DataFrame(
        np.column_stack([v[r, 0:j], np.concatenate(clst), v[r, j+1:]]),
        columns=orthopd.columns
    )

orthopd3 = pir2(orthopd, "SeqID")  #Running column split function on the SeqID column on orthopd
print("Beginning data analysis on the ProteinOrtho dataframe")

#Graph Algebraic Connectivity
orthopd_algconn_nozero = orthopd3[orthopd3['Alg.-Conn.'] != 0] #Removing zero and one counts in orthopd for graph
orthopd_algconn_noone = orthopd_algconn_nozero[orthopd_algconn_nozero['Alg.-Conn.'] != 1] #Getting the count of each Alg.Conn in the gff3 dataframe
orthopd_algconn_noone['Alg.-Conn.'].plot.hist(grid=True, bins=100, 
                   color='#607c8e')
plt.title('Distribution of Algebraic Connectivity without Unity')
plt.xlabel('Degree of Connectivity')
plt.ylabel('Number of Genes with Degree of Connectivity')
plt.tight_layout()
plt.savefig("ProteinOrtho_AlgConn_graph_noone.png")#Saving graph to file
plt.clf()
orthopd_algconn_nozero['Alg.-Conn.'].plot.hist(grid=True, bins=100, 
                   color='#607c8e')
plt.title('Distribution of Algebraic Connectivity')
plt.xlabel('Degree of Connectivity')
plt.ylabel('Number of Genes with Degree of Connectivity')
plt.tight_layout()
plt.savefig("ProteinOrtho_AlgConn_graph.png")#Saving graph to file
plt.clf()

#Graph Gene Counts
orthopd_gene_count_values = orthopd3['Genes'].value_counts() #Getting the count of each database in the gff3 dataframe
orthopd_gene_count_values.plot(kind='bar') #Graphing the database counts
plt.title('Graph of Gene Counts')
plt.xlabel('Number of Shared transcripts')
plt.ylabel('Number of Genes with same frequency')
plt.tight_layout()
plt.savefig("ProteinOrtho_gene_graph.png")#Saving graph to file
plt.clf()


#Start gff3pd
print("Working on gff3 dataframe")
gff3pd_column_names = ['SeqID', 'Database', 'Match type', 'Start', 'Stop', 'Score', 'Strand', 'Phase', 'Match information'] #Renaming static gff3 columns
statsgff3pd = pd.read_csv(gff3_statsfile_to_write, sep='\t', names=gff3pd_column_names, header=None, engine="python") #Creating a dataframe for gff3 stats
gff3pd_original_row_number = len(statsgff3pd) #Counting the number of rows in the original gff3pd dataframe for the statistics output
with open(zorbit_statistics, 'a') as f: #Writing the number of rows in the original gff3pd dataframe to the statistics output
    f.write("The number of sequences in the original gff3 file is " + str(gff3pd_original_row_number) + "\n")
gff3pd = pd.read_csv(gff3_file_to_write, sep='\t', names=gff3pd_column_names, header=None, engine = "python") #Creating a Pandas dataframe from the gff3 output csv
gff3pd_row_number = len(gff3pd) #Counting the number of rows in the final gff3 file dataframe for the statistics output
gff3pd_max_score = gff3pd['Score'].max() #Finding maximum value in Score column of gff3 dataframe
gff3pd_without_null = gff3pd[gff3pd['Score'] !=  "."] #Finding minimum value in Score column of gff3 dataframe
gff3pd_without_null_or_zero = gff3pd_without_null[gff3pd_without_null['Score'] != 0.0]
gff3pd_min_score = gff3pd_without_null_or_zero['Score'].min() 
statsgff3pd_without_null = statsgff3pd[statsgff3pd['Score'] !=  "."]
statsgff3pd_max_score = statsgff3pd_without_null['Score'].max()
with open(zorbit_statistics, 'a') as f:
    f.write("The number of sequences in the gff3 file after removal of MobiDBLite and duplicates is " + str(gff3pd_row_number) + "\n") #Adding cleaned gff3 stastitics to file
    f.write("The range of quality scores for the gff3 file range from " + str(gff3pd_min_score) + " to " + str(gff3pd_max_score) + "\n")#Adding range of scores to statistics file
    f.write("The maximum quality score for the original gff3 file is " + str(statsgff3pd_max_score) + "\n")

#Graph database distribution
gff3pd_database_count_values = gff3pd['Database'].value_counts() #Getting the count of each database in the gff3 dataframe
gff3pd_database_count_values.plot(kind='bar') #Graphing the database counts
plt.title('Distribution of Database hits')
plt.xlabel('Database name')
plt.ylabel('Number of Database hits')
plt.tight_layout()
plt.savefig("Gff3_database_graph.png")#Saving graph to file
plt.clf()

#Preparing dataframes for merging
print("Preparing dataframes for merge")
gff3pd['SeqID'] = gff3pd['SeqID'].astype(str) #Setting column type as string
orthopd3['SeqID'] = orthopd3['SeqID'].astype(str) #Setting column type as string
fastapd['SeqID'] = fastapd['SeqID'].astype(str) #Setting column type as string

#Dealing with fna versus faa
protein_flag = 0
if fasta_type == "nucleotide": #Checking to see if the fasta_type is nucleotide
    gff3pd_split = gff3pd['SeqID'].str.rsplit('_', n=2, expand=True) #Removing the extra two numbers after the fasta SeqID to allow match
    gff3pd['SeqID'] = gff3pd_split[0] #Setting the gff3 SeqID column as the split column
    orthopd_split = orthopd3['SeqID'].str.rsplit('_', n=2, expand=True) #Removing the extra two numbers after the fasta SeqID to allow match
    orthopd['SeqID'] = orthopd_split[0] #Setting the ProteinOrtho SeqID column as the split column
else:
    #Pulling out reading frame information
    protein_flag = 1
    gff3pd['SeqID2'] = gff3pd['SeqID']
    gff3pd_split = gff3pd['SeqID2'].str.rsplit('_', n=1, expand=True) #Removing the extra number after the fasta SeqID 
    gff3pd['SeqID2'] = gff3pd_split[0] #Setting the gff3 SeqID column as the split column
    gff3pd_split = gff3pd['SeqID2'].str.rsplit('_', n=1, expand=True) #Splitting the frame number out
    gff3pd['SeqID2'] = gff3pd_split[0] #Setting the gff3 SeqID column
    gff3pd['Reading_Frame'] = gff3pd_split[1] #Setting the gff3 Frame column
    gff3pd = gff3pd.drop(['SeqID2'], axis=1)
    orthopd3['SeqID2'] = orthopd3['SeqID']
    orthopd_split = orthopd3['SeqID2'].str.rsplit('_', n=1, expand=True) #Removing the extra two numbers after the fasta SeqID to allow match
    orthopd3['SeqID2'] = orthopd_split[0] #Setting the ProteinOrtho SeqID column as the split column
    orthopd_split = orthopd3['SeqID2'].str.rsplit('_', n=1, expand=True) #Splitting the frame number out
    orthopd3['SeqID2'] = orthopd_split[0] #Setting the orthopd SeqID column
    orthopd3['Reading_Frame'] = orthopd_split[1] #Setting the gff3 Frame column
    orthopd = orthopd3.drop(['SeqID2'], axis=1)
 

#Merging
print("Combining dataframes") 
gff3_ortho_merge = pd.merge(orthopd, gff3pd, how='outer', on=['SeqID']) #Merging the ProteinOrtho and interproscan dataframes
all_merge = pd.merge(gff3_ortho_merge, fastapd, how='outer', on=['SeqID']) #Merging the fasta dataframe with the combined ProteinOrtho/Interproscan dataframes



#Adding marks to merged dataframe to make fasta
all_merge['SeqID'] = all_merge['SeqID'].apply(lambda x: f'>{x}') #Placing > at the beginning of each new line and a tab at the end of SeqID
all_merge['Sequence'] = all_merge['Sequence'].apply(lambda x: f'\n{x}') #Placing a new line before the Sequence data
all_merge = all_merge[ ['SeqID'] + [ col for col in all_merge.columns if col != 'SeqID' ] ] #Moving SeqID to the far left of the dataframe
all_merge = all_merge[ [ col for col in all_merge.columns if col != 'Sequence' ] + ['Sequence'] ] #Moving Sequence to the far right of the dataframe

#Statistics on the merged dataframe
all_merge_both = all_merge.drop(all_merge[((all_merge['Database'].isna()) | (all_merge['Genes'] == 1))].index)
all_merge_neither = all_merge.drop(all_merge[((all_merge['Database'].notna()) | (all_merge['Genes'] !=1))].index)
all_merge_just_ortho = all_merge.drop(all_merge[((all_merge['Database'].notna()) | (all_merge['Genes'] == 1))].index)
all_merge_just_inter = all_merge.drop(all_merge[((all_merge['Database'].isna()) | (all_merge['Genes'] !=1))].index)


all_merge_all = len(pd.unique(all_merge['SeqID'])) #Calculating the number of unique sequences
all_merge_both = len(pd.unique(all_merge_both['SeqID'])) #Calculating unique sequences with both interproscan and proteinortho hits
all_merge_neither = len(pd.unique(all_merge_neither['SeqID'])) #Calculating unique sequences without interproscan or proteinortho hits
all_merge_just_ortho = len(pd.unique(all_merge_just_ortho['SeqID'])) #Calculating unique sequences with proteinortho but not interproscan hits
all_merge_just_inter = len(pd.unique(all_merge_just_inter['SeqID'])) #Calculating unique sequences with interproscan but no proteinortho hits

#Writing merged dataframe statistics to file
with open(zorbit_statistics, 'a') as f:
    f.write("The total number of unique sequences is: " + str(all_merge_all) + "\n")
    f.write("The number of unique sequences with both ProteinOrtho and Interproscan hits is: " + str(all_merge_both) + "\n")
    f.write("The number of unique sequences with neither ProteinOrtho nor Interproscan hits is: " + str(all_merge_neither) + "\n")
    f.write("The number of unique sequences with only ProteinOrtho and not Interproscan hits is: " + str(all_merge_just_ortho) + "\n")
    f.write("The number of unique sequences with only Interproscan and not ProteinOrtho hits is: " + str(all_merge_just_inter) + "\n")

#Excluding based on user preference
if exclusion_flag == 1:
    all_merge = all_merge.drop(all_merge[(all_merge['Genes'] == 1 & all_merge['Database'].isna())].index)
    print("Dropped rows")
    #Statistics on the merged dataframe
    all_merge_both = all_merge.drop(all_merge[((all_merge['Database'].isna()) | (all_merge['Genes'] == 1))].index)
    all_merge_neither = all_merge.drop(all_merge[((all_merge['Database'].notna()) | (all_merge['Genes'] !=1))].index)
    all_merge_just_ortho = all_merge.drop(all_merge[((all_merge['Database'].notna()) | (all_merge['Genes'] == 1))].index)
    all_merge_just_inter = all_merge.drop(all_merge[((all_merge['Database'].isna()) | (all_merge['Genes'] !=1))].index)
    all_merge_all = len(pd.unique(all_merge['SeqID'])) #Calculating the number of unique sequences
    all_merge_both = len(pd.unique(all_merge_both['SeqID'])) #Calculating unique sequences with both interproscan and proteinortho hits
    all_merge_neither = len(pd.unique(all_merge_neither['SeqID'])) #Calculating unique sequences without interproscan or proteinortho hits
    all_merge_just_ortho = len(pd.unique(all_merge_just_ortho['SeqID'])) #Calculating unique sequences with proteinortho but not interproscan hits
    all_merge_just_inter = len(pd.unique(all_merge_just_inter['SeqID'])) #Calculating unique sequences with interproscan but no proteinortho hits
    #Writing merged dataframe statistics to file
    with open(zorbit_statistics, 'a') as f:
        f.write("The total number of unique sequences (after dropping) is: " + str(all_merge_all) + "\n")
        f.write("The number of unique sequences with both ProteinOrtho and Interproscan hits (after dropping) is: " + str(all_merge_both) + "\n")
        f.write("The number of unique sequences with neither ProteinOrtho nor Interproscan hits (after dropping) is: " + str(all_merge_neither) + "\n")
        f.write("The number of unique sequences with only ProteinOrtho and not Interproscan hits (after dropping) is: " + str(all_merge_just_ortho) + "\n")
        f.write("The number of unique sequences with only Interproscan and not ProteinOrtho hits (after dropping) is: " + str(all_merge_just_inter) + "\n")
else:
    pass

if protein_flag == 1: #Taking care of aa sequence fastas
    #Plotting Frame data as a histogram
    all_merge = all_merge.drop(['Reading_Frame_y'], axis=1)
    all_merge = all_merge.rename(columns = {"Reading_Frame_x": "Reading_Frame"})
    all_merge['Reading_Frame'] = all_merge['Reading_Frame'].astype(str) #Setting column type as string
    all_merge_frame_counts = all_merge['Reading_Frame'].value_counts() #Getting the count of Frames in the all_merge dataframe
    all_merge_frame_counts.plot(kind='bar') #Graphing the frame counts
    plt.title('Distribution of Frame hits')
    plt.xlabel('Frame number')
    plt.ylabel('Number of Frame hits')
    plt.tight_layout()
    plt.savefig("Reading_frame.png")#Saving graph to file
    plt.clf()
    #Calculating number of duplicated, multi-frame sequences there are
    all_merge_duplicates = all_merge
    all_merge_split = all_merge['SeqID'].str.rsplit('_', n=2, expand=True) #Removing the extra two numbers after the fasta SeqID to allow match
    all_merge_duplicates['SeqID'] = all_merge_split[0] #Setting the ProteinOrtho SeqID column as the split column
    all_merge_two_columns = all_merge_duplicates[['SeqID', 'Reading_Frame']]
    all_merge_two_columns.to_csv(test2_file, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE, escapechar="\\")
    all_merge_unique_sequence_count = len(pd.unique(all_merge_two_columns['SeqID']))
    quotient = all_merge_unique_sequence_count/all_merge_all
    percentage = quotient*100
    with open(zorbit_statistics, 'a') as f:
            f.write("The number of unique sequences is: " + str(all_merge_unique_sequence_count) + "\n")
            f.write("The percentage of unique sequences is: " +  str(percentage) + "%" + "\n")


#Writing to file and cleaning up
all_merge.to_csv(test_file, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE, escapechar="\\") #Printing final dataframe to csv file with tab delimeters and no index or header
with open('test.txt', 'r') as infile, \
    open(file_to_write, 'w') as outfile:
    data = infile.read()
    data = data.replace("\\", "")
    outfile.write(data)

    

os.remove("test.txt") #Cleaning up files
os.remove("gff3_statsfile_to_write.txt")
os.remove("fasta_file_to_write.txt")
os.remove("gff3_file_to_write.txt")
os.remove("file_to_write.txt")
os.remove("ortho_file_to_write.txt")
os.remove("ortho_file_to_write2.txt")
print("Zorbit Analysis Complete")
input("Press any key to exit program")
