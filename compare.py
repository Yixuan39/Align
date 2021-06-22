import filecmp
from SelectFiles import *

def check1():
    print("compare input selectedFiles folders old vs new")
    for i in newList():
        try:
            if filecmp.cmp("compare/selectedFiles_1/Pillar"+str(i)+"_CDS.fas",
                           "compare/selectedFiles_2/Pillar"+str(i)+"_CDS.fas") == False:
                print(i)
        except:
            print(str(i)+"excluded")
            #pass

def checkCDS():
    for i in range(0,5589):
        with open("TGD_CDS/Pillar"+str(i)+"_CDS.fas","r") as a:
            with open("TGD_CDS_withGar/Pillar"+str(i)+"_CDS.fas","r") as b:
                old = a.readlines()
                new = b.readlines()
                if len(old) == len(new[:-2]):
                    print(old)
                    #print(new)
                    for n in range(0,len(old)):
                        if old[n] != new[:-2][n]:
                            print("Pillar"+str(i)+"line"+str(n)+"not same.")
                else:
                    print(len(old))
                    print(len(new))
                break




def main():
    """Description of main() - what does this function do?  Does it run a 
		program?  Does it execute test code?"""
    checkCDS()
    """print("compare tcoffee outcomes")
    for i in range(0, 5589):

        try:
            if filecmp.cmp("algndna_newold/new/Pillar" + str(i) + "_pro.pir",
                    "algndna_new/pir/Pillar" + str(i) + "_pro.pir") == False:
                print(i)
        except:
            #print(str(i)+"excluded")
            pass
    print("compare compare input file")
    for i in range(0, 5589):
        try:
            if filecmp.cmp("TGD_CDS_withGar_Oldversion/Pillar" + str(i) + "_CDS.fas",
                    "TGD_CDS_withGar/Pillar" + str(i) + "_CDS.fas") == False:
                print(i)
        except:
            print(str(i)+ "not in")"""


if __name__ == '__main__':
    main()
