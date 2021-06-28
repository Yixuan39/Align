# IGC 2 record

## Introduce the data

The data used in this project is downloaded from Dr. Gavin Conant’s website (http://wgd.statgen.ncsu.edu). The folder of the data is called “TGD_CDS” which contain gene of fishes. 

## Process of aligning the input files

1. Add Gar sequences

   The original gene sequence file, TGD_CDS, does not include the Gar (ENSLOCG) sequence. So attach the gar data to each of the 5588 files. The algorithm for this part is for each data file, we catch the name of the first gene (first line), then match the gene name with the dictionary Dr. Conant gave us. In the file named “TGD_reextract_ances_genes.txt”, each line represent a gene, and in each gene, gene names are tab separated (the first number in each line is the row number). 

2. Select according to Gavin’s instruction

   We select the fish data file according to Gavin’s list.  The number at the beginning of each line is the pillar number of the fish data. The selected data in the list were insured to have 2 paralogs of stickleback and 2 paralogs of zebrafish, with a 90% confident orthology inference.


3. Align the sequence using t-coffee and algndna_new software

   The alignment procedure has 3 steps:

   1. Use t-coffee to translate each DNA data to amino acid data

      ```bash
      t_coffee -other_pg seq_reformat -in selectedFiles/Pillar2148_CDS.fas -action +translate -output fasta_aln > outputs/Pillar2148_pro.fas
      ```

   2. Use t-coffee to convert each data into “.pir ” format. (For Gavin’s software to recognize)

      ```bash
       t_coffee outputs/Pillar2148_pro.fas -output=pir
      ```

   3. Use algndna_new to get the aligned sequence, where the input file include both the file in pir format and the original data (with gar) in fas format

      ``` bash
      ./algndna_new new/Pillar2148_pro.pir output/Pillar2148.fas -c:selectedFiles/Pillar2148_CDS.fas -s
      ```

4. Clean the gaps

   1. Make each sequence into one line
   2. Create an index list. Within one file, for each gene sequence, if there is a gap in the sequence line, record the index. Remove the duplicated item in the list, then we get a list of index whenever there is a gap in a column of the sequence matrix.
   3. Delete the item by index

5. Change the names to “01” and “02” according to Gavin’s file

   Refer to the picture of No.2, for each line of the file, there are 16 gene positions. If there is a missing gene, a “NONE” will appear in the corresponding position. Then, according to Gavin’s orthology inference outcome, the first 8 position will be renamed as paralog 1, and the others will be renamed as paralog 2.

6. Produce tree

   The algorithm I produce the tree for each gene is to set up a tree with full gene first, then delete the gene name where it does not have 2 paralogs in the gene sequence file (except for the outgroup). When there is a missing gene, in the “read in” tree, the “gene” will be changed to “”. After this, I clean up the tree where it appears any pattern that violate the rule of a newick tree.

   **Original tree structure**

   ```
   (((((ENSTRUG,ENSTNIG)N7,ENSGACG)N6,((ENSXMAG,ENSORLG)N5,ENSONIG)N4)N3,(ENSDARG,ENSAMXG)N2)N1,ENSLOCG)N0
   ```

   **for read in**

   ```python
   # "gene" represents a gene name, "*" represents a node name
   (((((gene,gene)*,gene)*,((gene,gene)*,gene)*)*,(gene,gene)*)*,ENSLOCG)*
   (((((,gene)*,gene)*,((gene,)*,gene)*)*,(gene,gene)*)*,ENSLOCG)* # Example
   ```

   Algorithm 1

   ```python
   i = 0
   while i < 1:
     	if "(," in treeStru:
       	treeStru=treeStru.replace("(,","(")
       	i = 0
     	elif ",)*" in treeStru:
       	treeStru=treeStru.replace(",)*",")*")
       	i = 0
     	elif "(,)" in treeStru:
       	treeStru=treeStru.replace("(,)","()")
       	i = 0
     	elif "()*" in treeStru:
       	treeStru=treeStru.replace("()*","")
       	i = 0
     	elif "(gene)*" in treeStru:
       	treeStru=treeStru.replace("(gene)*","gene")
       	i = 0
     	elif "((gene,gene)*)*" in treeStru:
       	treeStru=treeStru.replace("((gene,gene)*)*","(gene,gene)*") # special case
       	i = 0
     	else:
       	i = 1
   ```

   Algorithm 2

   ```python
   i = 0
   while i < 10:
       treeStru=treeStru.replace("(,","(")
       treeStru=treeStru.replace(",)*",")*")
       treeStru=treeStru.replace("(,)","()")
       treeStru=treeStru.replace("()*","")
       treeStru=treeStru.replace("(gene)*","gene")
       treeStru=treeStru.replace("((gene,gene)*)*","(gene,gene)*") # special case
       i += 1
   ```

   

7. Run through cluster

8. Optimization

## insight of the outcomes

1. Tau parameters

2. Calculate the tau proportion

   Formula $\frac{branch[1->2]+branch[2->1]}{branch[1->2]+branch[2->1]+branch[mut]}$ , we exclude (ENSLOCG, N0), (N0, N1)
   
3. Paralog exchange

   Using species: ENSDARG, ENSGACG, ENSLOCG (outgroup)


