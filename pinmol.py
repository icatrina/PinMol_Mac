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

def readSScountFile(filename): #sscount file for molecular beacon design
    userpath = os.path.dirname(filename) #use the path of input to save all files
    with open(filename) as infile: #check if the input file is txt and if Us are present not Ts
        if 'T' in infile.read():
            print('Wrong type of file (not txt) OR you may have used DNA parameters to fold your target RNA!')
            sys.exit('Try again using RNA parameters!')
    baselines = []
    with open(filename) as infile: #clean-up the file remove whitespaces from left side
        for line in infile:
            line = line.lstrip(" ")
            if not line:
                continue
            baselines.append(line.upper())
    so = int(baselines[0]) #number of suboptimal structures
    return (so, userpath, baselines)

def seqTarget(f): #sequence of target & sscount for each probe as fraction (1 for fully single stranded)
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
    return(sscount, position, max_base, bases, seq, size)

def itersplit_into_x_chunks(argum, size, chunksize): #split sequence in chunks of probe size
    for pos in range(0, size-chunksize):
        yield argum[pos:pos+chunksize]

def reverse_complement(seq): #generate RNA complement
    return seq.translate(basecomplement)[::-1]

def probeLength(probe): #input desired probe length; limited to range[14,28)]
    if probe <14 or probe >28:
        print('The value you entered is incorrect!')
        sys.exit('Try again!')
    else:
        probe = probe
    return (probe)

def seqProbes(mb_seq, mb_size, mb_sscount, probe):
    result = list(itersplit_into_x_chunks(mb_seq, mb_size, probe))
    basesl = []
    for i in result:
        i = reverse_complement(i)
        basesl.append(i)

    Tml = []
    for i in basesl:
        Tmx = mt.Tm_NN(i, dnac1 = 50000, dnac2 = 50000, Na = 100, nn_table = mt.RNA_NN1, saltcorr = 1)
        Tml.append(int(Tmx))
    result_bases = list(itersplit_into_x_chunks(mb_bases, mb_size, probe)) #list of lists of each base for each probe
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
    size2=len(mb_sscount)
    result2 = list(itersplit_into_x_chunks(mb_sscount, size2, probe))
    sumsl = []
    for i in result2:
        i = list(map(int, i))
        sums = sum(i)/(probe*mb_so)
        sumsl.append(sums)
    return (jl, perl, sumsl, basesl, Tml) #put together all data as indicated in header

def regionTarget(tg_start, tg_end): #if only a region of the target needs to be considered
    tg_diff = tg_end - tg_start
    assert(tg_start>0 and tg_end>0), "The base numbers should be positive and larger then zero!"
    if tg_diff < 0: #consider only probes within the requested region of target RNA, but check that the program can finish
        print("The number of end base cannot be smaller than the number of the start base!")
        sys.exit('Try again!')
    elif tg_diff < probe:
        print("You have to enter a region with a size larger than the probe size!")
        sys.exit('Try again!')
    elif tg_start >=1 and tg_end <= int(mb_max_base):
        df = pd.read_csv(mb_userpath+'/allcsv.csv')
        slice2 = df[tg_start:tg_end]
        for row in slice2:
            slice2 = slice2.sort_values(by='sscount', ascending=False) #sort descending by sscount = larger sscount more accessible target region
        return(slice2)

def numberProbes(no_pb):  #how many probes should be retained; limited to range [2, 50]
    if no_pb > int(tg_end-tg_start):
        print("This number is too large! You cannot enter a number larger then "+ str(tg_end-tg_start)+ " !")
        sys.exit('Try again!')
    elif row_no==0:
            print("No probes meet the criteria for the selected region, please expand the search region.")
            sys.exit('Try again!')
    elif no_pb > row_no and row_no > 0:
        print("Only "+str(row_no)+" meet the criteria.  Instead of "+ str(no_pb)+", " + str(row_no)+ " probe(s) will be considered")
        if row_no > 1:
            no_pb = row_no
            input1 = open(mb_userpath+'/sortedoutby3.csv', 'r').read().split('\n')
            output = open(mb_userpath +'/eg_sorted1.csv', 'w')
            output.write('\n'.join(input1))
            output.close()
        elif row_no==1:
            no_pb = row_no
            input1 = open(mb_userpath+'/sortedoutby3.csv', 'r')
            output.write(input1)
            output.close()

        with open(mb_userpath+'/eg_sorted1.csv','r') as f, open(mb_userpath+'/eg_sorted2.csv','w') as f1:
            next(f) # skip header line
            for line in f:
                f1.write(line)
    elif no_pb >0 and no_pb <= 50:
        no_pbf = no_pb
        input1 = open(mb_userpath+'/sortedoutby3.csv', 'r').read().split('\n')
        outputData = input1[:no_pb+1]
        output = open(mb_userpath +'/eg_sorted1.csv', 'w')
        output.write('\n'.join(outputData))
        output.close()
        with open(mb_userpath+'/eg_sorted1.csv','r') as f, open(mb_userpath+'/eg_sorted2.csv','w') as f1:
            next(f) # skip header line
            for line in f:
                f1.write(line)
    return(no_pb, no_pbf)

def blastChoice(blastm): #perform blast separately? HIGHLY RECOMMENDED!!!
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
        handle = open(mb_userpath+'/blast_picks.fasta','r')
        save_file = open(mb_userpath+'/my_blast.xml', 'a')
        for record2 in SeqIO.parse(handle, "fasta"):
            result_handle = NCBIWWW.qblast("blastn", "refseq_rna", record2, entrez_query = "\"" + org_blast + "\""+'[ORGANISM]', nucl_reward=1, nucl_penalty = -3)
                #   , format_type = "XML"
                #    , gapcosts = '5 2'
            save_file.write(result_handle.read())         
        handle.close()
        save_file.close()
        results3 = open(mb_userpath+'/my_blast.xml', 'r')
        #records= NCBIXML.parse(results3)
    elif blastm == 'y':
        print ("\n"*2+'Please use the file blast_picks.fasta to perform blast with refseq-rna database, and desired organism.'+'\n'+' For targets other than mRNAs make sure you use the Nucleotide collection (nr/nt) instead!')
        save_file = input('Enter path and file name for saved blast XML file: ')
        results3 = open(save_file, 'r')
    records= NCBIXML.parse(results3)

    pick = 0
    query1 = []
    with open (mb_userpath+'/blast_results.csv', 'w') as output:
        writer = csv.writer(output)
        writer.writerow(["Pick#", "Positives", "Gaps"])
        for record in records:
            pick +=1
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    query1.append(hsp.query)
                    assert hsp.positives <= probe, "A different probe length was used for blast to obtain the current XML file!"
                    if hsp.positives < probe and hsp.frame == (1, -1) and hsp.gaps == 0:
                        writer.writerow([pick, hsp.positives, hsp.gaps])
    
    with open (mb_userpath+'/blast_picks.fasta', 'r') as filebst:
        linebst = filebst.readlines()
        str1 = linebst[1].rstrip()
        str2 = str1.replace('U','T')
        str3 = linebst[-1].rstrip()
        str4 = str3.replace('U','T')
    records.close()
    results3.close()
    return (pick, query1, str2, str4)
    
def stemDesign(): #design the stem of the final probes
    i = 0
    with open(mb_userpath+'/mb_picks.csv') as ff:
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
            #print(sl)
            seql = list(item)
            seq_slc0 = seql[0:probe+1]
            seq_slc = list(itersplit_into_x_chunks(seq_slc0, probe+1, 4))
            cseq0 = seql[0].translate(basecomplement)
            cseq1 = seql[1].translate(basecomplement)
            cseq2 = seql[2].translate(basecomplement)
            stem13 = ['U', 'U', 'G', 'C']
            stem14 = ['G', 'U', 'G', 'U']
            stem15 = ['G', 'C', 'G', 'G']
            sls = 'CU'
            slq = 'UC'
            sl1 = list(sls)
            sl2 = list(slq)
            slco1 = []
            slco2 = []
            slco3 = []
            for slt in sl1:
                for slz in sl2:
                    X = slt
                    Y = slz
                    stem1 = [X, 'U', Y, 'G']
                    stem2 = ['G', X, 'U', Y]
                    stem3 = ['G', 'A', 'G', X]
                    stem4 = [X, 'G', 'A', 'G']
                    stem5 = ['G', 'U', X, 'G']
                    stem6 = [X, 'G', 'U', Y]
                    stem7 = ['G', 'A', X, 'G']
                    stem8 = [X, 'G', 'A', Y]
                    stem9 = [X, 'U', 'G', Y]
                    stem10 = ['G', X, 'U', 'G']
                    stem11 = [X, 'A', 'G', Y]
                    stem12 = ['G', X, 'A', 'G']
                    
                    slco1.append(stem1)
                    slco1.append(stem2)
                    slco1.append(stem3)
                    slco1.append(stem4)
                    slco2.append(stem5)
                    slco2.append(stem6)
                    slco2.append(stem7)
                    slco2.append(stem8)
                    slco3.append(stem9)
                    slco3.append(stem10)
                    slco3.append(stem11)
                    slco3.append(stem12)
                slcol1 = list(slco1)
                slcol2 = list(slco2)
                slcol3 = list(slco3)
                #for p in range(0,probe-3):
                #print(seq_slc)
            if  cseq0 == seql[-1] and cseq1 == seql[-2] and cseq2 == seql[-3]:
                stem = 'GG'
                
            elif cseq0 == seql[-1] and cseq1 == seql[-2]:
                stem = 'CCG'
                
            elif seql[0] == 'U' and seql[-1] == 'A' and stem15 not in seq_slc:
                stem = 'GCCG'
                
            elif seql[0] == 'U' and seql[-1] == 'A' and stem15 in seq_slc:
                stem = 'CCGG'
                
            elif seql[0] == 'A' and seql[-1] == 'U':
                stem = 'CGCC'
            
            elif seql[0] == 'C' and seql[-1] == cseq0 or seql[0] == 'G' and seql[-1] == cseq0:
                stem = 'CGAG'
            
            elif seql[0] == 'U' and seql[-1] == 'G' and stem13 not in seq_slc:
                stem = 'CGCGA'
                
            elif seql[0] == 'G' and seql[-1] == 'U':
                stem = 'CGCGA'
                
            elif any(s in seq_slc for s in (slcol3)) and stem14 not in seq_slc:
                stem = 'GCACG'
            
            elif any(s in seq_slc for s in (slco1))and stem14 not in seq_slc:
                stem = 'CGACG'
            
            elif any(s in seq_slc for s in (slcol2)) and stem14 not in seq_slc:
                stem = 'GCAGC'
                
            else:
                stem = 'CGAGC'
            

            stemr = reverse_complement(stem)
            #print(slcol2)
            aseq = (str(i) + " MB sequence at base number " + bs_ps[i-1] + ' is:  '+ stem + item + stemr)
            print(aseq + '\n')
            with open (mb_userpath+'/Final_molecular_beacons.csv', 'a') as outputf:
                outputf.write(aseq+'\n') 
            with open(mb_userpath+'/Seq'+str(i)+'.seq', 'w') as fiseq:
                fiseq.write(';'+'\n'+ str(i) + ' at base # ' + bs_ps[i-1] + ' molecular beacon' + '\n' + stem.strip() + item.strip() + stemr.strip()+'1')  
        return(i)
    


if __name__ == "__main__":
    filename = input('Enter a file name: ') # request ss-count file path and name
    mb_so, mb_userpath, mb_baselines = readSScountFile(filename)

    while True: #probe length?
        try:
            probe=int(input("Enter the length of probe; a number between 14 and 28: "))
        except:
            print('You must type a number between 20 and 28, try again:')
            continue

        else:
            probe = probeLength(probe)
            break
    
    with open (mb_userpath+'/usetxtfile.txt', 'w') as outfile:
        for line in mb_baselines:
            outfile.write(line)

    with open(mb_userpath+'/usetxtfile.txt', 'r') as infile, open(mb_userpath+'/mb_tocsv.csv', 'w') as csv_file:
        next(infile)
        reader = csv.reader(infile, delimiter = ' ')
        writer = csv.writer(csv_file, delimiter = ",", lineterminator = '\n')
        writer.writerows(reader)
        
    os.remove(mb_userpath+'/usetxtfile.txt')

    with open(mb_userpath+'/mb_tocsv.csv', 'r') as f:
        mb_sscount, mb_position, mb_max_base, mb_bases, mb_seq, mb_size = seqTarget(f)
        
    os.remove(mb_userpath+'/mb_tocsv.csv')
        
    basecomplement = str.maketrans({'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'})
    mb_jl, mb_perl, mb_sumsl, mb_basesl, mb_Tml = seqProbes(mb_seq, mb_size, mb_sscount, probe)
    with open(mb_userpath+'/allcsv.csv', 'w') as csv_file:
        writer = csv.writer(csv_file, lineterminator='\n')
        writer.writerow(["Base number", "%GC", "sscount", "Probe sequence", "Tm"])
        rows = zip(mb_jl,mb_perl,mb_sumsl,mb_basesl,mb_Tml)
        for row in rows:
            writer.writerow(row)

    tg_start = int(input("If a specific region within the target is needed, please enter the number of start base, or 1: "))
    tg_end = int(input("  and the number of end base or max number of bases " + str(mb_max_base) + ": "))
    mb_slice2 = regionTarget(tg_start, tg_end)
    mb_slice2.to_csv(mb_userpath+'/Full_List_sorted.csv', index=False)
    
    os.remove(mb_userpath+'/allcsv.csv')
 
    with open(mb_userpath+'/Full_List_sorted.csv','r') as ctin,open(mb_userpath+'/ss_strcnt.csv', 'w') as ctout:
        next(ctin)
        reader2 = csv.reader(ctin)
        writer2 = csv.writer(ctout, lineterminator='\n')
        for row2 in reader2:
            if float(row2[2]) >0.5:
                writer2.writerow(row2)
                
    with open(mb_userpath+'/Full_List_sorted.csv','r') as fin,open(mb_userpath+'/filter1.csv', 'w') as fout:
        next(fin)
        reader = csv.reader(fin)
        writer = csv.writer(fout, lineterminator='\n')
        for row in reader:
            if int(row[1]) >30 and int(row[1])<56:
                writer.writerow(row)
        
                
    #sort by sscount?
    flistsort = pd.read_csv(mb_userpath+'/filter1.csv', sep=',', usecols=[0,1,2,3,4])
    #new file with only sequences of probes for calculating free energies using oligoscreen
    flistsort2 = pd.read_csv(mb_userpath+'/filter1.csv', usecols=[3])
    flistsort2.to_csv(mb_userpath+'/Sorted_list.csv', index=False)
    
    
    with open (mb_userpath+'/filter1.csv', 'r') as GC_probes, open(mb_userpath+'/Full_List_sorted.csv', 'r') as all_probes:     
        no_GCprobes = sum(1 for row in GC_probes)
        no_probes = sum(1 for row in all_probes)
        

    with open(mb_userpath+'/Sorted_list.csv','r') as f, open(mb_userpath+'/updated_test.lis','w') as f1:
        next(f) # skip header line
        for line in f:
            f1.write(line)
            
    os.remove(mb_userpath+'/Sorted_list.csv')

    subprocess.check_output(["/Users/irinacatrina/Desktop/RNAstructure/exe/oligoscreen", mb_userpath+'/updated_test.lis', mb_userpath+'/updated_test.csv'])
    read_oligosc = pd.read_csv(mb_userpath+'/updated_test.csv', delimiter = '\t', usecols=[1,2,3])
    read_oligosc.to_csv(mb_userpath+'/updated_test2.csv', index=False)

    os.remove(mb_userpath+'/updated_test.lis')
    os.remove(mb_userpath+'/updated_test.csv')

    test3 = pd.concat([flistsort, read_oligosc], axis=1)
    test3.to_csv(mb_userpath+'/final_out.csv', index=False)
    #keep only probes that meet the energy requirements
    with open(mb_userpath+'/final_out.csv', 'r') as fin, open(mb_userpath+'/sortedout.csv', 'w') as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout, lineterminator='\n')
        next(fin, None)
        writer.writerow(["Base number", "%GC", "sscount", "Probe sequence", "Tm", "DGbimol", "DGunimol", "DGduplex"])
        for row in reader:
            if float(row[6]) >-2.5 and float(row[5])>-7.5:
                writer.writerow(row)
                
    df = pd.read_csv(mb_userpath+'/sortedout.csv')
    for row in df:
        sort3 = df.sort_values(['sscount','DGunimol', 'DGbimol'], ascending=[False, False, False]) #sort descending by sscount = larger sscount more accessible target region
        sort3.to_csv(mb_userpath+'/sortedoutby3.csv', index=False)
        
    os.remove(mb_userpath+'/updated_test2.csv')
    os.remove(mb_userpath+'/final_out.csv')        
        
    #determine the total number of probes that meet the eg criteria for the selected target (region or full)
    with open(mb_userpath+'/sortedoutby3.csv', 'r') as flistsort3:
        reader = csv.reader(flistsort3)
        row_no = sum(1 for row in reader)-1
        print("Maximum number of possible probes is: "+str(row_no)+"\n")

    try: #how many probes?
        no_pb = int(input ('How many probes do you want to save? Enter a number between 2 and 50: '))

    except:
        print('You must type a number between 2 and 50!')
        sys.exit('Try again!')

    else:
        no_pb, no_pbf = numberProbes(no_pb)
        
    read1 = pd.read_csv(mb_userpath+'/eg_sorted2.csv', delimiter = ',', usecols=[3])
    read1.to_csv(mb_userpath+'/forblast.csv', index=False)
    
    #write the fasta file containing the final sequences for blast
    with open(mb_userpath+'/forblast.csv', 'r') as file, open(mb_userpath+'/blast_picks.fasta','w') as f1:
        for line in file:
            f1.write('>'+'\n'+line)
    
    blastm = input('Do you want to perform blast separately? (highly recommended! as direct qblast may take 2h for 30 probes) y/n: ')
    pick, query1, str2, str4 = blastChoice(blastm)
    assert pick == no_pb, "The number of queries does not match the number of probes, wrong XML file?" #check if the correct file was used for blast 
    assert query1[0] == str2, "The first query is not the same as the first sequence in the blast_picks.fasta file" #check if the correct file was used for blast    
    assert query1[-1] in str4, "The last query is not contained in the last sequence in the blast_picks.fasta file" #check if the correct file was used for blast
    
    df = pd.read_csv(mb_userpath+'/blast_results.csv', sep = ",", index_col = None, engine = 'python')
    df_grouped = df.groupby(['Pick#']).agg({'Positives':'max'})
    df_grouped = df_grouped.reset_index()
    df = pd.merge(df, df_grouped, how='left', on=['Pick#'])
    a = pd.DataFrame(df_grouped)
    a.to_csv(mb_userpath+'/top_mb_picks.csv', index=False)

    with open(mb_userpath+'/top_mb_picks.csv') as f: #extract the information for final output from other files
        bl = [[],[],[]]
        reader = csv.reader(f)
        next(f)
        for row in reader:
            for col in range(2):
                bl[col].append(row[col])
    pick_no =bl[0]
    positives =bl[1]

    with open (mb_userpath+'/eg_sorted2.csv') as f3:
        dl = [[],[],[]]
        reader = csv.reader(f3)
        for row in reader:
            for col in range (3):
                dl[col].append(row[col])
    base_no = dl[0]
    sscntl = dl[2]

    with open(mb_userpath+'/forblast.csv', 'r') as f1:
        cl =[[]]
        reader = csv.reader(f1)
        for row in reader:
            for col in range(1):
                cl[col].append(row[col])
    probe_seq = cl[0]
    
    #os.remove(mb_userpath+'/top_mb_picks.csv')    
    os.remove(mb_userpath+'/forblast.csv')
    os.remove(mb_userpath+'/eg_sorted2.csv')

    with open (mb_userpath+'/blast_results_picks.csv', 'w') as output2:
        writer = csv.writer(output2)
        writer.writerow(["Pick#", "Base Number","Positives","Probe Sequence", "ss-count fraction"])
        rows = zip(pick_no, base_no, positives, probe_seq, sscntl )
        for row in rows:
            writer.writerow(row)

    df = pd.read_csv(mb_userpath+'/blast_results_picks.csv')
    for row in df:
        df = df.sort_values(['Positives', 'Pick#'], ascending=[True, True])
        df.to_csv(mb_userpath+'/Picks_Sorted.csv', index=False)
        

    mb_pick = pd.read_csv(mb_userpath+'/Picks_Sorted.csv', sep=',', usecols=[1,3])
    mb_pick.to_csv(mb_userpath+'/mb_picks.csv', index=False, header = False)

    i = stemDesign() #design the stem for the molecular beacon
    for x in range(1, int(i)+1):
        subprocess.check_output(["/Users/irinacatrina/Desktop/RNAstructure/exe/fold", mb_userpath+"/Seq"+str(x)+".seq" , mb_userpath+"/Seq"+str(x)+ ".ct"])
        subprocess.check_output(["/Users/irinacatrina/Desktop/RNAstructure/exe/draw", mb_userpath+"/Seq"+str(x)+".ct", mb_userpath+"/Seq"+str(x)+ ".svg", '--svg', '-n', '1'])    
    for j in range(1, int(i)+1):     #remove results that are highly structured
        with open (mb_userpath+"/Seq"+str(j)+".ct", 'r') as gin:
            linesa = gin.readlines()
            egdraw = float(linesa[0][16:20])
            no_bs = int(linesa[0][3:5])
            paired = int(linesa[1][23:26])
            #print (egdraw, no_bs, paired)
            if egdraw < -7.2 or egdraw > -2.5:
                os.remove(mb_userpath+"/Seq"+str(j)+".svg")
            elif no_bs != paired:
                os.remove(mb_userpath+"/Seq"+str(j)+".svg")
        os.remove(mb_userpath+"/Seq"+str(j)+".seq")
        os.remove(mb_userpath+"/Seq"+str(j)+".ct")      
    
    with open(mb_userpath+'/ss_strcnt.csv', 'r') as count_ss:
        no_ss = sum(1 for row in count_ss)
    
    print((no_probes-1), no_GCprobes, row_no, no_ss)

    print("\n"+"This information can be also be found in the file Final_molecular_beacons.csv"+"\n")
    print("\n"+"Check the structure for the selected probes using your favorite browser by opening the corresponding SVG files!")
    print("\n"+"If no SVG files are found, increase the number of probes and/or target region!")
    #remove intermediate files
    
    #os.remove(mb_userpath+'/grk'+'_'+str(probe)+'_50_'+str(mb_so)+'so'+'/Full_List_sorted.csv')
    #os.remove(mb_userpath+'/grk'+'_'+str(probe)+'_50_'+str(mb_so)+'so'+'/filter1.csv')
    #os.remove(mb_userpath+'/grk'+'_'+str(probe)+'_50_'+str(mb_so)+'so'+'/sortedout.csv')
    os.remove(mb_userpath+'/sortedoutby3.csv')
    os.remove(mb_userpath+'/ss_strcnt.csv')
    os.remove(mb_userpath+'/mb_picks.csv')
    #os.remove(mb_userpath+'/blast_results.csv')
    #os.remove(mb_userpath+'/blast_results_picks.csv')
