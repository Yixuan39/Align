from SelectFiles import *
import os

if __name__ == '__main__':
    input_folder = 'TGD_CDS_withGar/'
    output_folder = 'MAFFT_result/'
    os.makedirs(output_folder, exist_ok=True)
    files = os.listdir(input_folder)
    for file in files:
        mafft_cline = 'mafft ' + input_folder + file + ' > ' + output_folder + file
        os.system(mafft_cline)

# Clean Gaps
    path_in = output_folder
    path_out = 'mafft_Aligned_files'
    List = os.listdir(path_in)

    os.makedirs(path_out, exist_ok=True)
    for file in List:
        try:
            clean_seq(file,path_in=path_in, path_out=path_out)
        except:
            print(file)