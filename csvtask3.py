import csv
import os

#function for printing the 1-digit sequence
def printsq(sequence):
    val={'GLY':'G',
        'ALA':'A',
        'ARG':'R',
        'ASN':'N',
        'ASP':'D',
        'GLN':'Q',
        'GLU':'E',
        'HIS':'H',
        'ILE':'I',
        'LEU':'L',
        'LYS':'K',
        'CYS':'C',
        'MET':'M',
        'PHE':'F',
        'PRO':'P',
        'SER':'S',
        'THR':'T',
        'TRP':'W',
        'TYR':'Y',
        'VAL':'V'
        }
    return val[sequence]



file_path =r'C:\Users\Ishanee\Desktop\sumintprogs\1hh9.pdb'
with open(file_path, 'r') as file:
        s=''
        l=''
        lens=''
        cntx=0
        exp=''
        z=0
        brk='no'
        loc=''
        wb=''
        prot=''
        rval=''
        c=0
        chain=''
        res=''
        lst=''
        t=''
        modres=''
        HETLOC=[]
        spcgrp=''
        hetbrk='no'
        for line in file:
            #PDB id
            if line[0:6]=="HEADER":
                PDB_ID=line[62:66]

            #setting initial value for z to help in calculating chain break
            if line[0:4]=="ATOM" and z==0:
                 z=int((line[23:26]).strip())

            #for hetatm replaced by what
            if line.startswith('MODRES'):
                modres=modres+line[24:27]+'\n'
                HETLOC.append(line[19:22])

            #for space group
            if line[0:6]=="CRYST1":
                 spcgrp="SPACE GROUP:"+line[55:66].replace(" ","")
                 
            #uniprot id     elif
            if line[0:5]=="DBREF" and line[26:29]=="UNP":
                prot=prot+line[29:40].strip()+"("+line[7:13]+")"+"\n"+"\n"

            #to calculate chain length, chain composition and no. of chains elif
            if line[0:3]=="TER":
                cntx=cntx+1
                l=l+str(len(s))+'|'
                
                for i in "GARNDQEHILKCMFPSTWYV":
                    lens=lens+i+":"+str(s.count(i))+"|"
                lens=lens+"\n"
                s=''
                z=0

            #experiment type- X-ray or NMR
            elif line[0:6]=="REMARK" and line[12:27]=="EXPERIMENT TYPE":
                exp=line[44:100]
                exp=exp.strip()

            #for printing resolution and rvalue related to x-ray type experiments
            if line[0:6]=="REMARK" and line[11:21]=="RESOLUTION":
                    
                    res=line[11:21]+":"+line[22:31].strip()
                    #print(res)
            if line[0:6]=="REMARK" and line[13:20]=="R VALUE":
                    rval=("R-VALUE"+line[46:60]).strip()

            #for NMR and NMR related experiment details    
            elif exp=="NMR":
                #calculating no. of models by countimg end model "ENMDL" sequences
                if line[0:6]=="ENDMDL":
                    c=c+1

                #might have many subchains ( not completely correct)
                if line[0:4]=="ATOM" and line[21:22] not in chain:
                    chain=chain+line[21:22]
                    
            #for generating the 1-digit sequence by checking for C-alpha carbon
            if line[0:4]=="ATOM" and "CA" in line[13:16]:
                s=s+printsq(line[17:20])


            #for which hetatm
            elif t[0:4]=="ATOM" and line[0:6]=="HETATM" and line[23:26] in HETLOC:
                 lst=lst+line[17:20]+"\n"
                 hetbrk="yes"
                 #line[0:6]="ATOM  "

            #for calculating chain break, which chain and the location of the break
            elif line[0:4]=="ATOM": # or (line[0:6]=="HETATM" and t[0:4]=="ATOM"):
                temp=int((line[23:26]).strip())
                if z+1==temp:
                    z=z+1
                if z==temp:
                    continue
                else:
                    brk='yes'
                    loc=loc+str(z)+"-"+line[23:26]+"\n"
                    wb=wb+str(cntx+1)+'\n'
                    z=temp

            #for storing the previous line
            t=line
        #print(s)

        if exp=="NMR":
            deets="No. of models:"+str(c)
            index1=l.index("|")
            l=l[0:index1]
            index2=lens.index("\n")
            lens=lens[0:index2]
        else:
            deets=res+"\n"+rval+"\n"+spcgrp
        
        if brk=="no":
            loc="-"
            wb='-'
            hetbrk="no"



# Data to write
data = [ [PDB_ID, cntx, l, lens, brk, wb, loc, hetbrk, lst, modres, exp, deets, prot] ]

# CSV file path
file_path = r'C:\Users\Ishanee\Desktop\ops\tab.csv'

# Check if the file already exists
file_exists = os.path.exists(file_path)

# Open the CSV file for appending
with open(file_path, 'a', newline='\n') as csvfile:
    writer = csv.writer(csvfile, dialect='excel', quoting=csv.QUOTE_ALL)
    
    # Write the header only if the file is new or empty
    if not file_exists or os.path.getsize(file_path) == 0:
        writer.writerow(["PDB_ID", "No. of chains", "Chain Length", "Composition of Chain",
                         "Chain Break", "Which Chain Break", "Location of Break",
                         "Chain break due to Hetatm", "Which Hetatm", "Replaced", 
                         "Experiment", "Experiment Details", "UniProt ID"])
    
    # Write the data rows
    for row in data:
        writer.writerow(row)
