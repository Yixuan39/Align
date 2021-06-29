import os



"""
This program select the file form "TGD_CDS" according to Dr. Conant's text file "TGD_DrGaDupl_1plusloss_90per"
"""

def read():
    numberList=[]
    with open("TGD_DrGaDupl_1plusloss_90per.txt", "r") as f:
        for line in f: # same as readlines() then loop over, it's better because readlines() will read all lines and cache in memory
            number = line.split()[0] # first element as number, do not assume number length
            numberList.append(number)
    return numberList




def GarMatching():
    dictionary = {}
    with open("TGD_reextract_ances_genes.txt","r") as old:
        o = old.readlines()
        for line in o:
            newline = line.split("\t")
            for string in newline[1:]:
                if "ENSLOCG" in string:
                    dic = {newline[0]:string}
                    dictionary.update(dic)
    return dictionary

def GarMatching2():
    numbers = range(0,5589)
    oldGeneDic = {}
    GarGeneDic = {}
    for number in numbers:
        with open("TGD_CDS/Pillar"+str(number)+"_CDS.fas","r") as old:
            firstGene = old.readlines()
            firstGene = firstGene[0][1:-1] # "-1" is to delete the "\n" at the end
            subdic = {firstGene:number}
            oldGeneDic.update(subdic)
    #print(oldGeneDic)

    with open("TGD_reextract_ances_genes.txt", "r") as G:
        GarList = []
        Gar = G.readlines()
        for line in Gar:
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
            with open(directory + "/Pillar"+str(number)+"_CDS.fas","w") as w:
                for line in r:
                   w.write(line)
    except:
        pass


def addGar():
    try:
        os.mkdir("TGD_CDS_withGar")
    except:
        pass
    ranges = range(0, 5589)
    numbers = []
    for i in ranges:
        i = str(i)
        numbers.append(i)
    dic = GarMatching2()
    for number in numbers:
        copyFile(number, "TGD_CDS_withGar",dic)
        
def select(List,direcetory):
    try:
        os.mkdir("algndna_new/pir/"+direcetory)
    except:
        pass
    numbers = List
    for number in numbers:
        with open("TGD_CDS_withGar/Pillar"+str(number)+"_CDS.fas","r") as old:
            with open("algndna_new/pir/"+direcetory+"/Pillar"+str(number)+"_CDS.fas", "w") as new:
                o = old.readlines()
                for i in o:
                    new.write(i)




def createScriptsIn1(List,directory):
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


"""def jointList():
    newList = []
    dic = GarMatching2()
    GarList = []
    for i in dic:
        GarList.append(i)
    OldList = read()
    for i in GarList:
        if i in OldList:
            newList.append(i)
    with open("newList.txt","w") as new:
        for number in newList:
            new.write(number+"\n")


def newList():
    newList = []
    with open("newList.txt","r") as o:
        old=o.readlines()
        for line in old:
            newList.append(line[:-1])
    return newList"""






def main():
    """Description of main() - what does this function do?  Does it run a 
		program?  Does it execute test code?"""
    List45 = ['211', '214', '222', '223', '337', '479', '521', '526', '561', '735', '755', '852', '1050',
              '1053', '1215', '2129', '2158', '2210', '2214', '2321', '2358', '2371', '2382', '2861',
              '3278', '3295', '3309', '3337', '3346', '3347', '3390', '3994', '4025', '4031', '4063',
              '4268', '4287', '4494', '4553', '4570', '4932', '5153', '5233', '5316', '5550']
    List = read()
    directory = "selectedFiles"

    createScriptsIn1(List,directory)
    addGar()
    select(List,directory)



if __name__ == '__main__':
    main()




"""
dictionary = {}
    with open("TGD_DrGaDupl_1plusloss_90per.txt","r") as o:
        Fulldic = {}
        Old = o.readlines()
        FullList = []
        for line in Old:
            line = line.split("\t")
            for i in line:
                if i == "NONE":
                    line.remove(i)
                FullList.append(line)

    with open("TGD_reextract_ances_genes.txt","r") as G:
        GarList = []
        Gar = G.readlines()
        for line in Gar:
            line = line.split("\t")
            for i in line:
                if i == "NONE":
                    line.remove(i)
                if "\n" in i:
                    newi=i[:-1]
                    line.remove(i)
                    line.append(newi)
            GarList.append(line)

    for i in FullList:
        print(i)
        #for j in GarList:

            #print(j[1])"""