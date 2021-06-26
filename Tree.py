
from SelectFiles import *
import os

def qualify(List):
    number = List
    for i in number:
        with open("renamed/Pillar"+str(i)+"R.fasta","r") as old:
            outgroup = 0
            otherspecies = 0
            o = old.readlines()
            nameList = []
            for line in o:
                if ">" in line:
                    if "ENSLOCG" in line:
                        outgroup += 1
                    if line[1:8] not in nameList:
                        nameList.append(line[1:8])
                    else:
                        otherspecies += 1
            if outgroup != 1:
                print(i)
            if otherspecies < 2:
                print(i)

# all data are qualified


def Diff(li1, li2):
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif


def nameDic():

    with open("TGD_DrGaDupl_1plusloss_90per.txt","r") as new:
        n=new.readlines()
        dic = {}
        for line in n:
            newList = []
            newline = line.split("\t")
            for i in newline[2:18]:
                if i != "NONE":
                    newList.append(i[0:7])
            new = set([x for x in newList if newList.count(x) > 1])
            dic.update({newline[0]:list(new)})
    return dic




def buildTree(directory,number):
    dic = nameDic()
    gene = ["ENSTRUG","ENSTNIG","ENSGACG","ENSXMAG","ENSORLG","ENSONIG","ENSDARG","ENSAMXG"]
    node = "*"

    fakegene = []
    newgene = []
    for i in gene:
        if i in dic[number]:
            fakegene.append("gene")
            newgene.append(i)
        else:
            fakegene.append("")
    treeStru = "(((((" + fakegene[0] + "," + fakegene[1] + ")"+node+"," + fakegene[2] + ")"+node+",((" + fakegene[3] + "," + fakegene[4] + ")"+node+"," + \
               fakegene[5] + ")"+node+")"+node+",(" + fakegene[6] + "," + fakegene[7] + ")"+node+")"+node+",ENSLOCG)"+node
    j = 0
    while j < 1:
        if "(," in treeStru:
            treeStru = treeStru.replace("(,", "(")
            j = 0
        elif ",)*" in treeStru:
            treeStru = treeStru.replace(",)*", ")*")
            j = 0
        elif "(,)" in treeStru:
            treeStru = treeStru.replace("(,)", "()")
            j = 0
        elif "()*" in treeStru:
            treeStru = treeStru.replace("()*", "")
            j = 0
        elif "(gene)*" in treeStru:
            treeStru = treeStru.replace("(gene)*", "gene")
            j = 0
        elif "((gene,gene)*)*" in treeStru:
            treeStru = treeStru.replace("((gene,gene)*)*", "(gene,gene)*")  # special case
            j = 0
        else:
            j = 1
    #print(treeStru)
    for k in newgene:
        treeStru=treeStru.replace("gene",k,1)
    nodeNum = treeStru.count("*") - 1
    for l in range(0,nodeNum+1):
        treeStru=treeStru.replace("*","N"+str(nodeNum-l),1)
    #print(newgene)
    print(treeStru)
    with open(directory+"/"+str(number)+"/Pillar"+str(number)+".newick","w") as new:
        new.write(treeStru)



#        while "(,)" in treeStru:
#            treeStru=treeStru.replace("(,)","()")
#        while "(," in treeStru:
#            treeStru=treeStru.replace("(,","(")
#        while ",)N" in treeStru:
#            treeStru=treeStru.replace(",)N",")N")


#(((((ENSTRUG,ENSTNIG)N7,ENSGACG)N6,((ENSXMAG,ENSORLG)N5,ENSONIG)N4)N3,(ENSDARG,ENSAMXG)N2)N1,ENSLOCG)N0
#(()N1,ENSLOCG)N0
#(((((ENSTRUG,ENSTNIG)N7,ENSGACG)N6,((ENSXMAG,ENSORLG)N5,ENSONIG)N4)N3,(ENSDARG,ENSAMXG)N2)N1,ENSLOCG)N0



def main():
    """Description of main() - what does this function do?  Does it run a
		program?  Does it execute test code?"""

    List45 = ['211', '214', '222', '223', '337', '479', '521', '526', '561', '735', '755', '852', '1050',
              '1053', '1215', '2129', '2158', '2210', '2214', '2321', '2358', '2371', '2382', '2861',
              '3278', '3295', '3309', '3337', '3346', '3347', '3390', '3994', '4025', '4031', '4063',
              '4268', '4287', '4494', '4553', '4570', '4932', '5153', '5233', '5316', '5550']
    List = read()
    #nameDic()
    directory = "TGD_CDS_new"
    numbers = List
    for number in numbers:
        file = os.path.exists("renamed/Pillar"+str(number)+"R.fasta")
        if file:
            try:
                os.makedirs(directory+"/" + str(number),exist_ok = True)
                with open("renamed/Pillar"+str(number)+"R.fasta","r") as old:
                    with open(directory+"/"+str(number)+"/Pillar"+str(number)+".fasta","w") as new:
                        for line in old:
                            new.write(line)
                buildTree(directory,number)
            except:
                print("no Pillar" + str(number))
    #MGList = []
    #for i in numbers:
    #    MGList.append("MG"+i)
    #print(MGList)


if __name__ == '__main__':
    main()



"""
ENSTRUG = "ENSTRUG"
    ENSTNIG = "ENSTNIG"
    ENSGACG = "ENSGACG"
    ENSXMAG = "ENSXMAG"
    ENSORLG = "ENSORLG"
    ENSONIG = "ENSONIG"
    ENSDARG = "ENSDARG"
    ENSAMXG = "ENSAMXG"
"""