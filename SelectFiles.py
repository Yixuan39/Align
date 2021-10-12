import os

"""
This program select the file form "TGD_CDS" according to Dr. Conant's text file "TGD_DrGaDupl_1plusloss_90per"
"""

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




def GarMatching2():
    """
    This function catch the first gene of each dataset and match with correct pillar numbers
    Then, according to Gavin's Gar gene names file, return a dictionary of pillar number mathching
    gar sequence
    """
    oldGeneDic = {}
    GarGeneDic = {}
    for number in range(0,5589):
        with open("TGD_CDS/Pillar"+str(number)+"_CDS.fas","r") as f:
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






def copyFile(number, directory, dic):
    """
    This function add the gar sequence to a individual dataset
    """
    try:
        os.mkdir(directory)
    except:
        pass
    #dic = GarMatching2()
    try:
        with open("TGD_CDS/Pillar"+str(number)+"_CDS.fas","r") as read:
            r = read.readlines()
            rowNum=0
            with open("lepisosteus_oculatus_dna.fas", "r") as Gar:
                g = Gar.readlines()
                for line in g:
                    if dic.get(str(number)) in line:
                        r.append(g[rowNum][:-5]+"\n")
                        r.append(g[rowNum + 1])
                        r.append("\n")
                    rowNum += 1
            with open(directory + "/Pillar"+str(number)+"_CDS.fas","w+") as w: # "w+" creates the file if it doesn't exist
                for line in r:
                   w.write(line)
    except:
        pass


def addGar():
    """
    add all gar sequence separately to all 5588 dataset.
    """
    try:
        os.mkdir("TGD_CDS_withGar")
    except:
        pass
    dic = GarMatching2()
    for number in range(0, 5589):
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
        with open("TGD_CDS_withGar/Pillar"+str(number)+"_CDS.fas","r") as old:
            with open("algndna_new/pir/"+direcetory+"/Pillar"+str(number)+"_CDS.fas", "w") as new:
                o = old.readlines()
                for i in o:
                    new.write(i)




def createScriptsIn1(List,directory):
    """
    Create folders and scripts of commands needed for t-coffee and algndna_new software
    """
    try:
        os.mkdir("algndna_new/pir")
    except:
        pass
    try:
        os.mkdir("algndna_new/pir/outputs")
    except:
        pass
    try:
        os.mkdir("algndna_new/output")
    except:
        pass
    number = List
    with open("algndna_new/pir/tcoffeeCommand.sh","w") as f:
        f.write("#!/bin/bash\n\n")
        for i in number:
            # step 1
            f.write("t_coffee -other_pg seq_reformat -in "+directory+"/Pillar"
                      + str(i) + "_CDS.fas -action +translate -output fasta_aln > outputs/Pillar" + str(i) + "_pro.fas\n")
            # step 2
            f.write("t_coffee outputs/Pillar" + str(i) + "_pro.fas -output=pir\n")
    with open("algndna_new/algndnaCommand.sh","w") as f:
        f.write("#!/bin/bash\n\n")
        for i in number:
            f.write(" ./algndna_new pir/Pillar"+str(i)+"_pro.pir output/Pillar"+str(i)+".fas -c:pir/"+directory+"/Pillar"+str(i)+"_CDS.fas -s\n")
            #new.write(" ./algndna_new new/Pillar" + str(i) + "_pro.pir output/Pillar" + str(i) + ".fas -c:selectedFiles/Pillar" + str(i) + "_CDS.fas -s\n")









def main():
    """Description of main() - what does this function do?  Does it run a 
		program?  Does it execute test code?"""
    # 45 Full gene dataset
    List45 = ['211', '214', '222', '223', '337', '479', '521', '526', '561', '735', '755', '852', '1050',
              '1053', '1215', '2129', '2158', '2210', '2214', '2321', '2358', '2371', '2382', '2861',
              '3278', '3295', '3309', '3337', '3346', '3347', '3390', '3994', '4025', '4031', '4063',
              '4268', '4287', '4494', '4553', '4570', '4932', '5153', '5233', '5316', '5550']
    # 194 missing value dataset
    List = List45 + read()
    # List = [*range(0,5589)]
    print(List)
    directory = "selectedFiles"
    print("step1")
    createScriptsIn1(List,directory)
    print("step2")
    #addGar()
    print("step3")
    select(List,directory)



if __name__ == '__main__':
    main()



