import os
import re
"""
This program select the file form "TGD_CDS" according to Dr. Conant's text file "TGD_DrGaDupl_1plusloss_90per"
"""

def num_sort(test_string):
    return list(map(int, re.findall(r'\d+', test_string)))[0]

def read():
    """
    This function read the pillar numbers where there is at least one gene missing
    but including zebrafish, stickleback fish, and gar, according to Gavin's orthology
    inference output (TGD_DrGaDupl_1plusloss_90per.txt)
    """
    numberList=[]
    with open("TGD_DrGaDupl_1plusloss_90per.txt", "r") as f:
        for line in f: # same as readlines() then loop over, it's better because readlines() will read all lines and cache in memory
            number = line.split()[0] # first element as number, do not assume number length
            numberList.append(number)
    return numberList




def GarMatching2(List):
    """
    This function catch the first gene of each dataset and match with correct pillar numbers
    Then, according to Gavin's Gar gene names file, return a dictionary of pillar number mathching
    gar sequence
    """
    oldGeneDic = {}
    GarGeneDic = {}
    for number in List:
        with open("TGD_CDS/"+str(number),"r") as f:
            for line in f:
                if line[0] == '>':
                    firstGene = line[1:-1] # "-1" is to delete the "\n" at the end
                    oldGeneDic[firstGene] = number
                    break
    #print(oldGeneDic)

    with open("TGD_reextract_ances_genes.txt", "r") as G:
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






def copyFile(name, directory, dic):
    """
    This function add the gar sequence to a individual dataset
    """
    os.makedirs(directory, exist_ok=True)
    #dic = GarMatching2(List)
    # try:
    with open("TGD_CDS/"+str(name),"r") as read:
        r = read.readlines()
        rowNum=0
        with open("lepisosteus_oculatus_dna.fas", "r") as Gar:
            g = Gar.readlines()
            for line in g:
                if dic.get(name) in line:
                    r.append(g[rowNum][:-5]+"\n")
                    r.append(g[rowNum + 1])
                    r.append("\n")
                rowNum += 1
        with open(directory + "/"+str(name),"w") as w: # "w+" creates the file if it doesn't exist
            for line in r:
                w.write(line)
    # except:
    #     print("no")
    #     pass


def addGar(List):
    """
    add all gar sequence separately to all 5588 dataset.
    """
    try:
        os.mkdir("TGD_CDS_withGar")
    except:
        pass
    dic = GarMatching2(List=List)
    for number in List:
        copyFile(str(number), "TGD_CDS_withGar",dic)



def select(List,direcetory):
    """
    This function selects the appropriate files from the whole data folder, according to
    Gavin's text file.
    """

    os.makedirs("algndna_new/pir/"+direcetory, exist_ok=True)
    numbers = List
    for number in numbers:
        print(number)
        with open("TGD_CDS_withGar/"+str(number),"r") as old:
            with open("algndna_new/pir/"+direcetory+"/"+str(number), "w") as new:
                o = old.readlines()
                for i in o:
                    new.write(i)





def main():
    """Description of main() - what does this function do?  Does it run a 
		program?  Does it execute test code?"""
    # 45 Full gene dataset
    input_folder = './TGD_CDS'
    List = os.listdir(input_folder)
    # List = ["Pillar"+i+"_CDS.fas" for i in read()]
    print(List)
    directory = "selectedFiles"

    print("adding Gar sequence")
    addGar(List)
    select(List,directory)

    os.makedirs("algndna_new/pir/outputs", exist_ok=True)
    os.makedirs("algndna_new/output", exist_ok=True)

    for file in List:
        print("aligning" + file)
        os.system("cd algndna_new/pir/ \nt_coffee -other_pg seq_reformat -in "+directory+"/"
                      + str(file) + " -action +translate -output fasta_aln > ./outputs/" + file.replace('_CDS.fas', '')
                  + "_pro.fas\nls\nt_coffee ./outputs/" + file.replace('_CDS.fas', '') + "_pro.fas -output=pir\n")
        print("cd algndna_new\n ./algndna_new pir/" + file.replace('_CDS.fas', '') + "_pro.pir output/"
                  + file.replace('_CDS.fas', '') + ".fas -c:pir/" + directory + "/" + file + " -s\n")
        os.system("cd algndna_new\n ./algndna_new pir/" + file.replace('_CDS.fas', '') + "_pro.pir output/"
                  + file.replace('_CDS.fas', '') + ".fas -c:pir/" + directory + "/" + file + " -s\n")


if __name__ == '__main__':
    main()



