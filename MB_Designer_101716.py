#!/anaconda/bin/python3
import csv
import sys
import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
import subprocess
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
print("\n"*5)
undscr = "->"*30
print(undscr)
print("\n"+"WARNING: Previous files will be overwritten!  Save them in a "+"\n"+"different location than the current file, or rename them to "+"\n"+"ensure they are not misused (e.g. use probes from a different target)."+"\n")
print(undscr)

# request ss-count file path and name
filename = input('Enter a file name: ')

#use the path of input to save all files
userpath = os.path.dirname(filename)

#check if the input file is txt and if Us are present not Ts
with open(filename) as infile:
    if 'T' in infile.read():
        print('Wrong type of file (not txt) OR you may have used DNA parameters to fold your target RNA!')
        sys.exit('Try again using RNA parameters!')

#clean-up the file remove whitespaces from left side
with open(filename) as infile, open(userpath+'/usefile2.txt', 'w') as outfile:
    for line in infile:
        line = line.lstrip(" ")
        if not line:
            continue
        outfile.write(line.upper())
        
#extract the number of suboptimal structures
with open (userpath+'/usefile2.txt', 'r') as f:
    allval = f.readlines()
    so = int(allval[0])

# convert to csv file format  
with open(userpath+'/usefile2.txt', "r") as in_txt:
    next(in_txt)
    reader = csv.reader(in_txt, delimiter = " ")
    with open(userpath+'/tocsv2.csv','w') as csv_file:
        writer = csv.writer(csv_file, delimiter = ",", lineterminator = '\n')
        writer.writerows(reader)
    
# request the desired length of the probe sequence
while True:
        try:
                probe=int(input("Enter the length of probe; a number between 14 and 28: "))
        except:
                print('You must type a number between 14 and 28, try again:')
                continue

        else:
                if probe <14 or probe >28:
                        print('The value you entered is incorrect!')
                        sys.exit('Try again!')
                else:
                        probe = probe
                        break
                    
#extract information as lists and total number of bases for the target
with open(userpath+'/tocsv2.csv', 'r') as f:
    bl = [[],[],[]]
    reader = csv.reader(f)
    for row in reader:
        for col in range(3):
            bl[col].append(row[col])
sscount = bl[1]
position = bl[0]
max_base = bl[0][-1]
bases = bl[2]
seq = ''.join(bases)
size = len(seq)

# split sequence into fragments of probe size
def itersplit_into_x_chunks(argum, size, chunksize): # we assume here that x is an int and > 0
    for pos in range(0, size-chunksize):
        yield argum[pos:pos+chunksize]

#list of target sequence split into fragments of probe length as list of each probe
size = len(seq)
result = list(itersplit_into_x_chunks(seq, size, probe))

basecomplement = str.maketrans({'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'})
basesl = []

#RNA complement
def reverse_complement(seq):
    return seq.translate(basecomplement)[::-1]

#reverse complement of each probe sequence
for i in result:
     i = reverse_complement(i)
     basesl.append(i)

#calculate melting temperature for each probe
Tml = []
for i in basesl:
     Tmx = mt.Tm_NN(i, dnac1 = 50000, dnac2 = 50000, Na = 100, nn_table = mt.RNA_NN1, saltcorr = 1)
     Tml.append(int(Tmx))

result_bases = list(itersplit_into_x_chunks(bases, size, probe)) #list of lists of each base for each probe

#base number as j and list of these numbers as jl, list of percent of Gs and Cs as perl
j = 0
perl = []
jl = []

for i in result_bases:
    j += 1
    cs = i.count('C')
    gs = i.count('G')
    per = int((cs+gs)/probe*100)
    perl.append(per)
    jl.append(j)

#calculate sscount for each probe from values of individual bases
size2=len(sscount)
result2 = list(itersplit_into_x_chunks(sscount, size2, probe))

#calulate the overall sscount for each probe as fraction (1 for fully single stranded)
sumsl = []
for i in result2:
    i = list(map(int, i))
    sums = sum(i)/(probe*so)
    sumsl.append(sums)

#put together all data as indicated in header
with open(userpath+'/allcsv.csv', 'w') as csv_file:
    writer = csv.writer(csv_file, lineterminator='\n')
    writer.writerow(["Base number", "%GC", "sscount", "Probe sequence", "Tm"])
    rows = zip(jl,perl,sumsl,basesl,Tml)
    for row in rows:
        writer.writerow(row)

# consider only certain regions of target? maybe use only CDS
df = pd.read_csv(userpath+'/allcsv.csv')
tg_start = int(input("If a specific region within the target is needed, please enter the number of start base, or 1: "))
tg_end = int(input("  and the number of end base or max number of bases " + str(max_base) + ": "))
tg_diff = tg_end - tg_start
assert(tg_start>0 and tg_end>0), "The base numbers should be positive and larger then zero!"

#consider only probes within the requested region of target RNA, but check that the program can finish
if tg_diff < 0:
    print("The number of end base cannot be smaller than the number of the start base!")
    sys.exit('Try again!')
elif tg_diff < probe:
    print("You have to enter a region with a size larger than the probe size!")
    sys.exit('Try again!')
#sort descending by sscount = larger sscount more accessible target region
elif tg_start >=1 and tg_end <= int(max_base):
    slice2 = df[tg_start:tg_end]
    for row in slice2:
        slice2 = slice2.sort_values(by='sscount', ascending=False)
        slice2.to_csv(userpath+'/Full_List_sorted.csv', index=False)

#filter sorted probes to keep only probes with desired GC percentage
with open(userpath+'/Full_List_sorted.csv','r') as fin:
    next(fin)
    reader = csv.reader(fin)
    with open(userpath+'/filter1.csv', 'w') as fout:
        writer = csv.writer(fout, lineterminator='\n')
        for row in reader:
            if int(row[1]) >30 and int(row[1])<56:
                writer.writerow(row)

#combine in csv file is this needed?
with open(userpath+'/combined_file.csv', 'w') as outcsv:
    writer = csv.writer((outcsv), lineterminator = '\n')
    #writer.writerow(["Base number", "%GC", "sscount", "Probe sequence", "Tm"])
    with open(userpath+'/filter1.csv', 'r', newline='') as incsv:
        reader = csv.reader(incsv)
        writer.writerows(row for row in reader)

#sort by sscount?
flistsort = pd.read_csv(userpath+'/Full_List_sorted.csv', sep=',', usecols=[0,1,2,3,4])
    
#new file with only sequences of probes for calculating free energies using oligoscreen
flistsort2 = pd.read_csv(userpath+'/Full_List_sorted.csv', usecols=[3])
flistsort2.to_csv(userpath+'/Sorted_list.csv', index=False)

with open(userpath+"/Sorted_list.csv",'r') as f:
        with open(userpath+"/updated_test.lis",'w') as f1:
                next(f) # skip header line
                for line in f:
                        f1.write(line)
subprocess.check_output(["/Users/irinacatrina/Desktop/RNAstructure/exe/oligoscreen", userpath + "/updated_test.lis", userpath+ "/updated_test.csv"])
read_oligosc = pd.read_csv(userpath+'/updated_test.csv', delimiter = '\t', usecols=[1,2])
read_oligosc.to_csv(userpath+'/updated_test2.csv', index=False)

test3 = pd.concat([flistsort, read_oligosc], axis=1)
test3.to_csv(userpath+'/final_out.csv', index=False)

#keep only probes that meet the energy requirements
with open(userpath+'/final_out.csv', 'r') as fin, open(userpath+'/sortedout.csv', 'w') as fout:
    reader = csv.reader(fin)
    writer = csv.writer(fout, lineterminator='\n')
    next(fin, None)
    writer.writerow(["Base number", "%GC", "sscount", "Probe sequence", "Tm", "DGbimol", "DGunimol"])
    for row in reader:
        if float(row[6]) >-2 and float(row[5])>-5:
            writer.writerow(row)

#determine the total number of probes that meet the eg criteria for the selected target (region or full)
with open(userpath+'/sortedout.csv', 'r') as flistsort3:
    reader = csv.reader(flistsort3)
    row_no = sum(1 for row in reader)-1
    print("Maximum number of possible probes is: "+str(row_no)+"\n")

try:
    no_pb = int(input ('How many probes do you want to save? Enter a number between 2 and 50: '))

except:
    print('You must type a number between 2 and 50!')
    sys.exit('Try again!')

else:
    if no_pb > int(tg_end-tg_start):
        print("This number is too large! You cannot enter a number larger then "+ str(tg_end-tg_start)+ " !")
        sys.exit('Try again!')
    elif row_no==0:
            print("No probes meet the criteria for the selected region, please expand the search region.")
            sys.exit('Try again!')
    elif no_pb > row_no and row_no > 0:
        print("Only "+str(row_no)+" meet the criteria.  Instead of "+ str(no_pb)+", " + str(row_no)+ " probe(s) will be considered")
        if row_no > 1:
            input1 = open(userpath+'/sortedout.csv', 'r').read().split('\n')
            output = open(userpath + '/eg_sorted1.csv', 'w')
            output.write('\n'.join(input1))
            output.close()
        elif row_no==1:
            input1 = open(userpath+'/sortedout.csv', 'r')
            output.write(input1)
            output.close()
        
        with open(userpath+'/eg_sorted1.csv','r') as f, open(userpath+"/eg_sorted2.csv",'w') as f1:
            next(f) # skip header line
            for line in f:
                f1.write(line)
    elif no_pb >0 and no_pb <= 50:
        input1 = open(userpath+'/sortedout.csv', 'r').read().split('\n')
        outputData = input1[:no_pb+1]
        output = open(userpath + '/eg_sorted1.csv', 'w')
        output.write('\n'.join(outputData))
        output.close()
        with open(userpath+'/eg_sorted1.csv','r') as f, open(userpath+"/eg_sorted2.csv",'w') as f1:
            next(f) # skip header line
            for line in f:
                f1.write(line)

read1 = pd.read_csv(userpath+'/eg_sorted2.csv', delimiter = ',', usecols=[3])
read1.to_csv(userpath+'/forblast.csv', index=False)

#write the fasta file containing the final sequences for blast
with open(userpath+'/forblast.csv', 'r') as file, open(userpath+"/blast_picks.fasta",'w') as f1:
    for line in file:
        f1.write('>'+'\n'+line)


blastm = input('Do you want to perform blast separately? (highly recommended! as direct qblast may take 2h for 30 probes) y/n: ')

if blastm == 'n':
    # what organism should be used for qblast
    org_blastno = int(input("Enter the organism for BLAST: 1 for D. melanogaster, 2 for H. sapiens and 3 for M. musculus: "))
    if org_blastno == 2:
        org_blast = ('Homo sapiens')
    elif org_blastno == 3:
        org_blast = ('Mus musculus')
    elif org_blastno == 1:
        org_blast = ('Drosophila melanogaster')
    else:
        print("You have to enter 1 or 2 or 3!")
        sys.exit('Try again!')
        #perform qblast fpr selected probes
        #print ('The organism for blast is: '+ org_blast)
    #record2 = SeqIO.parse(userpath+'/blast_picks.fasta', format='fasta')
    handle = open(userpath+"/blast_picks.fasta","r")
    save_file = open(userpath+"/my_blast.xml", "a")
    for record2 in SeqIO.parse(handle, "fasta"):
        result_handle = NCBIWWW.qblast("blastn", "refseq_rna", record2, entrez_query = "\"" + org_blast + "\""+'[ORGANISM]', nucl_reward=1, nucl_penalty = -3)
            #   , format_type = "XML"
            #    , gapcosts = '5 2'
        save_file.write(result_handle.read())
#    save_file = (userpath+"/my_blast.xml")
     
    handle.close()
    save_file.close()
    results3 = open(userpath+"/my_blast.xml", 'r')
    records= NCBIXML.parse(results3)
#handle.close()
elif blastm == 'y':
    print ('Use file blast_picks.fasta to perform blast with refseq-rna database, and desired organism')
    save_file = input('Enter path and file name for saved blast XML file: ')
    results3 = open(save_file, 'r')
    records= NCBIXML.parse(results3)

#userpath+"/my_blast.xml"
pick = 0
query1 = []
with open (userpath+'/blast_results.csv', 'w') as output:
    writer = csv.writer(output)
    writer.writerow(["Pick#", "Positives", "Gaps"])
    for record in records:
        pick +=1
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                query1.append(hsp.query)
                assert hsp.positives <= probe, "A different probe length was used for blast to obtain the current XML file!"
                if hsp.gaps == 0:
                    if hsp.frame == (1, -1):
                        if hsp.positives < probe:
                            writer.writerow([pick, hsp.positives, hsp.gaps])
#print(pick)
assert pick == no_pb, "The number of queries does not match the number of probes, wrong XML file?" #check if the correct file was used for blast

with open (userpath+'/blast_picks.fasta', 'r') as filebst:
    linebst = filebst.readlines()
    str1 = linebst[1].rstrip()
    str2 = str1.replace('U','T')

#print(query1[0])
#print(str2)
assert query1[0] == str2, "The first query is not the same as the first sequence in the blast_picks.fasta file" #check if the correct file was used for blast

    
records.close()
results3.close()

df = pd.read_csv(userpath+'/blast_results.csv', sep = ",", index_col = None, engine = 'python')
grouped = df.groupby('Pick#').first().reset_index()

a = pd.DataFrame(grouped)
a.to_csv(userpath+"/top_mb_picks.csv", index=False)

#extract the information for final output from other files
with open(userpath+'/top_mb_picks.csv') as f:
    bl = [[],[],[]]
    reader = csv.reader(f)
    next(f)
    for row in reader:
        for col in range(2):
            bl[col].append(row[col])
pick_no =bl[0]
positives =bl[1]

with open (userpath+'/eg_sorted2.csv') as f3:
    dl = [[]]
    reader = csv.reader(f3)
    for row in reader:
        for col in range (1):
            dl[col].append(row[col])
base_no = dl[0]

with open(userpath+'/forblast.csv', 'r') as f1:
    cl =[[]]
    reader = csv.reader(f1)
    for row in reader:
        for col in range(1):
            cl[col].append(row[col])
probe_seq = cl[0]

with open (userpath+'/blast_results_picks.csv', 'w') as output2:
    writer = csv.writer(output2)
    writer.writerow(["Pick#", "Base Number","Positives","Probe Sequence"])
    rows = zip(pick_no, base_no, positives, probe_seq)
    for row in rows:
        writer.writerow(row)

df = pd.read_csv(userpath+'/blast_results_picks.csv')
for row in df:
    df = df.sort_values(by="Positives", ascending=True)
    df.to_csv(userpath+'/Picks_Sorted.out', index=False)

i=0
mb_pick = pd.read_csv(userpath+'/Picks_Sorted.out', sep=',', usecols=[1,3])
mb_pick.to_csv(userpath+'/mb_picks.csv', index=False, header = False)

print("\n"*2+"Number of suboptimal structures in the selected ss-count file: "  + filename + '=  ' + str(so)+"\n"*2)

#design the stem for the molecular beacon
with open(userpath+'/mb_picks.csv') as ff:
    tw =[[],[]]
    reader = csv.reader(ff)
    for row in reader:
        for col in range(2):
            tw[col].append(row[col])
    bs_ps =tw[0]
    fseq =tw[1]
    seq_slc = []
    for item in fseq:
        i += 1
        seql = list(item)
        cseq0 = seql[0].translate(basecomplement)
        cseq1 = seql[1].translate(basecomplement)
        cseq2 = seql[2].translate(basecomplement)
        stem1 = ['C', 'G', 'U', 'C', 'G']
        seq_slc = []
        for p in range(0,probe-5):
            seq_slc.append(seql[p:p+5])
        if stem1 in seq_slc:
            stem = 'GCACG'
        elif cseq0 == seql[-1] and cseq1 == seql[-2] and cseq2 == seql[-3]:
            stem = 'CG'
        elif cseq0 == seql[-1] and cseq1 == seql[-2]:
            stem = 'CGC'
        elif seql[0] == 'A' and seql[-1] == 'U':
            stem = 'CGCC'
        elif seql[0] == 'U' and seql[-1] == 'A':
            stem = 'CGCC'
        elif seql[0] == 'C' and seql[-1] == cseq0 or seql[0] == 'G' and seql[-1] == cseq0:
            stem = 'CGAC'
        else:
            stem = 'CGACG'
                               
        stemr = reverse_complement(stem)
        aseq = (str(i) + " molecular beacon sequence at base number " + bs_ps[i-1] + ' is:  '+ stem + item + stemr)
        with open(userpath+'/Seq'+str(i)+'.seq', 'w') as fiseq:
            fiseq.write(';'+'\n'+ str(i) + ' molecular beacon' + '\n' + stem.strip() + item.strip() + stemr.strip()+'1')

        subprocess.check_output(["/Users/irinacatrina/Desktop/RNAstructure/exe/fold", userpath+ "/Seq"+str(i)+".seq" , userpath+ "/Seq"+str(i)+ ".ct"])
        subprocess.check_output(["/Users/irinacatrina/Desktop/RNAstructure/exe/draw", userpath+ "/Seq"+str(i)+".ct", userpath+ "/Seq"+str(i)+ ".svg", '--svg', '-n', '1'])
        print(aseq+"\n")
        with open (userpath + '/Final_molecular_beacons.csv', 'a') as outputf:
            outputf.write(aseq+'\n')

#remove results that are highly structured        
for j in range(1, no_pb+1):     
    with open (userpath+"/Seq"+str(j)+".ct", 'r') as gin:
        linesa = gin.readlines()
        egdraw = float(linesa[0][16:20])
        if egdraw < -4.0 or egdraw > -2.5:
            os.remove(userpath+"/Seq"+str(j)+".seq")
            os.remove(userpath+"/Seq"+str(j)+".ct")
            os.remove(userpath+"/Seq"+str(j)+".svg")

print("\n"+"This information can be also be found in the file Final_molecular_beacons.csv"+"\n")
print("\n"+"Check the structure for the selected probes using your favorite browser by opening the corresponding SVG files!")
print("\n"+"If no SVG file are found, increase the number of probes and/or target region!")
#remove intermediate files
os.remove(userpath+'/usefile2.txt')
os.remove(userpath+'/tocsv2.csv')
os.remove(userpath+'/allcsv.csv')  
os.remove(userpath+'/updated_test.lis')
os.remove(userpath+'/Sorted_list.csv')
os.remove(userpath+'/Full_List_sorted.csv')
os.remove(userpath+'/filter1.csv')
os.remove(userpath+'/updated_test2.csv')
os.remove(userpath+'/sortedout.csv')
os.remove(userpath+'/final_out.csv')
os.remove(userpath+'/forblast.csv')
os.remove(userpath+'/eg_sorted1.csv')
os.remove(userpath+'/top_mb_picks.csv')
os.remove(userpath+'/mb_picks.csv')
os.remove(userpath+'/updated_test.csv')
os.remove(userpath+'/blast_results.csv')
os.remove(userpath+'/blast_results_picks.csv')
