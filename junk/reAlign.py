import os
from SelectFiles import *



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
                new.write(line + '\n')

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
        with open("algndna_new/output/" + i,"r") as o:
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
    List45 = ['210','211', '214', '222', '223', '337', '479', '521', '526', '561', '735', '755', '852', '1050',
              '1053', '1215', '2129', '2158', '2210', '2214', '2321', '2358', '2371', '2382', '2861',
              '3278', '3295', '3309', '3337', '3346', '3347', '3390', '3994', '4025', '4031', '4063',
              '4268', '4287', '4494', '4553', '4570', '4932', '5153', '5233', '5316', '5550']
    List = read() + List45

    path_in = "algndna_new/output"
    List = os.listdir(path_in)

    if len(check(List)) == 0:
        print("all data have 2 paralog1, 2 paralog2, and an outgroup")
    else:
        print("*")
    output_folder = 'Aligned_files'
    os.makedirs(output_folder, exist_ok=True)
    for file in List:
        try:
            clean_seq(file, path_out=output_folder)
        except:
            print(file)




if __name__ == '__main__':
    main()


