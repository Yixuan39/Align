#!/bin/python3
import os
import re
"""
This program select the file form "TGD_CDS" according to Dr. Conant's text file "TGD_DrGaDupl_1plusloss_90per"
"""

def num_sort(test_string):
    return list(map(int, re.findall(r'\d+', test_string)))[0]

def read(file = "../content_file/TGD_DrGaDupl_1plusloss_90per.txt"):
    """
    This function read the pillar numbers where there is at least one gene missing
    but including zebrafish, stickleback fish, and gar, according to Gavin's orthology
    inference output (TGD_DrGaDupl_1plusloss_90per.txt)
    """
    numberList=[]
    with open(file, "r") as f:
        for line in f: # same as readlines() then loop over, it's better because readlines() will read all lines and cache in memory
            number = line.split()[0] # first element as number, do not assume number length
            numberList.append(number)
    return numberList




def GarMatching2(path, files, GarNames = "../content_file/TGD_reextract_ances_genes.txt"):
    """
    This function catch the first gene of each dataset and match with correct pillar numbers
    Then, according to Gavin's Gar gene names file, return a dictionary of pillar number mathching
    gar sequence
    """
    oldGeneDic = {}
    GarGeneDic = {}
    for number in files:
        with open(path+str(number),"r") as f:
            for line in f:
                if line[0] == '>':
                    firstGene = line[1:-1] # "-1" is to delete the "\n" at the end
                    oldGeneDic[firstGene] = number
                    break
    #print(oldGeneDic)

    with open(GarNames, "r") as G:
        GarList = []
        for line in G:
            line = line[:-1] # bug "48" fixed by delete the "\n" at the end of the line
            line = line.split("\t")
            GarList.append(line)
        for i in oldGeneDic:
            for j in GarList:
                if i in j:
                    subdic2 = {str(oldGeneDic[i]):j[1]}
                    GarGeneDic.update(subdic2)
    return GarGeneDic

def copyGarSeq(name,path_in, path_out, dic, garSeq = "../content_file/lepisosteus_oculatus_dna.fas"):
    """
    This function add the gar sequence to a individual dataset
    """
    os.makedirs(path_out, exist_ok=True)
    #dic = GarMatching2(List)
    # try:
    with open(path_in+str(name),"r") as read:
        r = read.readlines()
        rowNum=0
        with open(garSeq, "r") as Gar:
            g = Gar.readlines()
            for line in g:
                if dic.get(name) in line:
                    r.append(g[rowNum][:-5]+"\n")
                    r.append(g[rowNum + 1])
                    r.append("\n")
                rowNum += 1
        with open(path_out + "/"+str(name),"w") as w:
            for line in r:
                w.write(line)


def addGar(List, path_in, path_out):
    """
    add all gar sequence separately to all 5588 dataset.
    """
    os.makedirs(path_out, exist_ok=True)
    dic = GarMatching2(files=List, path=path_in)
    for number in List:
        copyGarSeq(str(number), path_in, path_out, dic)



def select(List,direcetory):
    """
    This function selects the appropriate files from the whole data folder, according to
    Gavin's text file.
    """

    os.makedirs("../algndna_new/pir/"+direcetory, exist_ok=True)
    numbers = List
    for number in numbers:
        print(number)
        with open("../input_file/TGD_CDS_withGar/"+str(number),"r") as old:
            with open("../algndna_new/pir/"+direcetory+"/"+str(number), "w") as new:
                o = old.readlines()
                for i in o:
                    new.write(i)

def clean_seq(file, path_in = "algndna_new/output", path_out = "Aligned_files"):
    with open(path_in + "/"+ file,"r") as old:
        out = []
        for line in old.readlines():
            if ">" not in line:
                line = line[:-1]
            else:
                line = "\n" + line
            out.append(line)
        out[0]=out[0][1:] # Delete the first "\n"
        out = ''.join(out)
        out = out.split('\n')
        out = removeGap(out)
        with open(path_out + "/" + file.replace("fas", "fasta"),"w") as new:
            for line in out:
                new.write(line.upper() + '\n')

def listGap(file):
    # a list of index where there is at least one gap in a column
    out = []
    for line in file:
        if ">" not in line:
            index = 0
            for i in line:
                if i == "-":
                    out.append(index)
                index += 1
        out = list(set(out))
    return out

def removeGap(file):
    list = listGap(file)
    out = []
    for line in file:
        newline = ""
        if ">" in line:
            newline = line
        else:
            for i in range(0,len(line)):
                if i not in list:
                    newline += line[i]
        out.append(newline)
    return out

def nameDic():
    with open("TGD_DrGaDupl_1plusloss_90per.txt","r") as new:
        n=new.readlines()
        dic = {}
        for line in n:
            subdic1={}
            newline = line.split("\t")
            for i in newline[2:10]:
                subdic1.update({i:"01"})
            for j in newline[10:18]:
                subdic1.update({j:"02"})
            dic.update({newline[0]:subdic1})
    return dic

def matchName(List):
    dic = nameDic()
    number = List
    for i in number:
        try:
            with open("reformatFolder/Pillar"+str(i)+"_gapsClean.fas","r") as old:
                rename = []
                o = old.readlines()
                for line in o:
                    if ">" in line:
                        if "ENSLOCG" in line:
                            line = ">ENSLOCG01" + "\n"
                        else:
                            line = line[0:8] + dic[str(i)][line[1:19]] +"\n"
                    rename.append(line)
                with open("renamed/Pillar"+str(i)+"R.fasta","w") as new:
                    for line in rename:
                        new.write(line)
        except:
            print("no Pillar"+str(i))


def onlyLetter(List):
    out = []
    for i in List:
        out.append("".join(re.findall("[a-zA-Z]+", i)))
    return out

def check(path, n, outgroup='ENSLOCG'):
    files = os.listdir(path)
    out_files = []
    for file in files:
        lines = open(path + file, 'r').readlines()
        names = [line for line in lines if '>' in line if outgroup not in line]
        names = onlyLetter(names)
        n_paralog = {name:names.count(name) for name in names}
        if len([i for i in n_paralog.values() if i < 2]) <= n:
            out_files.append(file)
    return out_files


def main():
    """Description of main() - what does this function do?  Does it run a 
		program?  Does it execute test code?"""
    # 45 Full gene dataset
    input_folder = '../input_file/TGD_CDS/'
    List = os.listdir(input_folder) # All files from Gavin
    # List = ["Pillar"+i+"_CDS.fas" for i in read()]
    directory = "selectedFiles"
    List = check(input_folder, 0)
    List = [i.replace("R.fasta","_CDS.fas") for i in os.listdir('../input_file/IGC_R') if 'fasta' in i]
    print(len(List))
    print("adding Gar sequence")
    addGar(List,path_in=input_folder, path_out="../input_file/TGD_CDS_withGar")
    select(List,directory)


    os.makedirs("../algndna_new/pir/outputs", exist_ok=True)
    os.makedirs("../algndna_new/output", exist_ok=True)

    for file in List:
        print("aligning" + file)
        os.system("cd ../algndna_new/pir\nls \nt_coffee -other_pg seq_reformat -in "+directory+"/"
                      + str(file) + " -action +translate -output fasta_aln > ./outputs/" + file.replace('_CDS.fas', '')
                  + "_pro.fas\nls\nt_coffee ./outputs/" + file.replace('_CDS.fas', '') + "_pro.fas -output=pir\n")

        os.system("cd ../algndna_new\n ./algndna_new .pir/" + file.replace('_CDS.fas', '') + "_pro.pir output/"
                  + file.replace('_CDS.fas', '') + ".fas -c:pir/" + directory + "/" + file + " -s\n")

# Clean Gaps
    path_in = "algndna_new/output"
    List = os.listdir(path_in)

    output_folder = 'Aligned_files'
    os.makedirs(output_folder, exist_ok=True)
    for file in List:
        try:
            clean_seq(file, path_out=output_folder)
        except:
            print(file)

if __name__ == '__main__':
    main()



