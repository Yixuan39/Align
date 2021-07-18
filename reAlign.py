import os
from SelectFiles import *



def oneLine(List):
    # make each sequence into one line
    numbers = List
    try:
        os.mkdir("reformatFolder")
    except:
        pass
    for number in numbers:
        try:
            with open("algndna_new/output/Pillar"+str(number)+".fas","r") as old:
                out = []
                for line in old.readlines():
                    if ">" not in line:
                        line = line[:-1]
                    else:
                        line = "\n" + line
                    out.append(line)
                out[0]=out[0][1:] # Delete the first "\n"
                with open("reformatFolder/Pillar"+str(number)+"_oneLine.fas","w") as new:
                    for line in out:
                        new.write(line)
        except:
            pass

def listGap(number):
    # a list of index where there is at least one gap in a column
    try:
        with open("reformatFolder/Pillar" + str(number) + "_oneLine.fas", "r") as old:
            out = []
            for line in old.readlines():
                if ">" not in line:
                    index = 0
                    for i in line:
                        if i == "-":
                            out.append(index)
                        index += 1
                out = list(set(out))
        return out
    except:
        pass



def removeGap(number):
    list = listGap(number)
    try:
        with open("reformatFolder/Pillar" + str(number) + "_oneLine.fas", "r") as old:
            out = []
            for line in old.readlines():
                newline = ""
                if ">" in line:
                    newline = line
                else:
                    for i in range(0,len(line)):
                        if i not in list:
                            newline += line[i]
                out.append(newline)
            with open("reformatFolder/Pillar" + str(number) + "_gapsClean.fas", "w") as new:
                for i in out:
                    new.write(i)
    except:
        pass

def removeGapAll(List):
    oneLine(List)
    numbers = List
    for number in numbers:
        removeGap(number)

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
    try:
        os.mkdir("renamed")
    except:
        pass
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








def check(List):
    number = List
    #number = range(0,5589)
    output = []
    unQ = []
    for i in number:
        with open("algndna_new/output/Pillar"+str(i)+".fas","r") as o:
            paralog1 = 0
            paralog2 = 0
            outgroup = 0
            old = o.readlines()
            for line in old:
                if ">" in line:
                    if 'ENSGACG' in line:
                        paralog1 += 1
                    elif 'ENSDARG' in line:
                        paralog2 += 1
                    elif 'ENSLOCG' in line:
                        outgroup += 1
            if paralog1 == 2:
                if paralog2 == 2:
                    if outgroup == 1:
                        output.append(i)
    for i in output:
        if i not in List:
            unQ.append(i)
    return unQ




def main():
    """Description of main() - what does this function do?  Does it run a 
		program?  Does it execute test code?"""
    List45 = ['211', '214', '222', '223', '337', '479', '521', '526', '561', '735', '755', '852', '1050',
              '1053', '1215', '2129', '2158', '2210', '2214', '2321', '2358', '2371', '2382', '2861',
              '3278', '3295', '3309', '3337', '3346', '3347', '3390', '3994', '4025', '4031', '4063',
              '4268', '4287', '4494', '4553', '4570', '4932', '5153', '5233', '5316', '5550']
    List = read()



    #if len(check(List)) == 0:
    #    print("all data have 2 paralog1, 2 paralog2, and an outgroup")
    oneLine(List)
    removeGapAll(List)
    matchName(List)



if __name__ == '__main__':
    main()


