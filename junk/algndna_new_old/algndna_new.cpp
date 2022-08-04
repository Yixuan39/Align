//Copyright 1999-2002 Gavin Conant

#include <iostream>
#include <fstream>
#include "string.h"
#include "gen_dna_funcs.h"
#include "gen_code.h"
#include "read_seq.h"
#include "write_seq.h"
#include "algndna_funcs.h"
#define DEBUG


using namespace::std;

int main(int argc, char **argv)
{
	int i,j,k, numtaxa, numchars, exon, DNAsize, dummyt, exonpos, cds_in_file, codon[3], num_DNA_dummy, DNA_id;
	char dump[5],  proname[81], outname[81], **pro_ids, DNAname[100],  param_file[100], read_file[300], inchar, nuc_file[300];

	Data_info current_info;
	String_for_list an_id;
	List<String_for_list> list_gene_ids;
	ifstream paramin;
	BOOL end_set=FALSE, backcheck_ok=TRUE, all_single_exon=FALSE, single_nuc_file=FALSE;
	Gene_struct *curtaxa, new_taxainfo;
	Sequence_dataset *current_data, *DNA_seq, *final_sequence, *new_protein_sequences;
	Read_Sequence *get_seq;
	Read_FASTA *get_DNA;
	Write_Sequence *write_seq;
	Exchange current_exchange;
	Genetic_code *current_code;

	if (argc<4) {
		cerr<<"Syntax: algndna <protein alignment file>  <output file> <Protein ID1> <Protein ID2> ... -f:<DNA location file>/(-s) (-c:DNA file)\n";
	}

	else {
      
      strcpy(proname, argv[1]);
      strcpy(outname, argv[2]);  
      
      //Reads in the list of possible gene ids
      while(list_gene_ids.get_list_length()<argc-3 && argv[list_gene_ids.get_list_length()+3][0]!='-') {
		  strcpy(an_id.string, argv[list_gene_ids.get_list_length()+3]);
		  list_gene_ids.add_to_list(an_id);
	  }
      
  
      pro_ids=new char* [list_gene_ids.get_list_length()];
      list_gene_ids.return_to_start();

      for (i=0; i<list_gene_ids.get_list_length(); i++) {
		  pro_ids[i]=new char[20];
		  strcpy(pro_ids[i], list_gene_ids.get_current_item()->string);
	  }

      current_info.set_gene_ids(list_gene_ids.get_list_length(), pro_ids);
      
    
      for (i=0; i<list_gene_ids.get_list_length(); i++)
		  delete[] pro_ids[i];
      delete pro_ids;

   
      parse_args(argc, argv, param_file, all_single_exon, single_nuc_file, nuc_file, &current_exchange);
      current_code=create_genetic_code(&current_exchange);

	  if (param_file[0] != '\0') {
		current_info.set_global_paramfile(param_file);
		get_config_file_data(&current_info);
	  }

	  switch (guess_dataformat(proname, strlen(proname))) {
		case NEXUS:
		  get_seq=new Read_Nexus();
		  break;
		case PIR:
		  get_seq=new Read_PIR();
		  break;
		case PHYLIP:
		  get_seq=new Read_Phylip_interleave();
		  break;
		case FASTA:
		  get_seq=new Read_FASTA();
		  break;
	  }

     
      switch(guess_dataformat(outname, strlen(outname))) {
		case NEXUS:
		  write_seq=new Write_Nexus(outname, NUCLEIC);
		  break;
		case PIR:
		  write_seq=new Write_PIR(outname, NUCLEIC);
		  break;
		case PHYLIP:
		  write_seq=new Write_Phylip_interleave(outname, NUCLEIC);
		  break;
		case FASTA:
		  write_seq=new Write_FASTA(outname, NUCLEIC);
		  break;
	  }

      
      
      current_data=get_seq->get_dataset(numtaxa, numchars, proname, TRUE);
      if (current_data==0 || write_seq->fail()==TRUE) {
		  delete write_seq;
		  return(-1);
	  }

      final_sequence = new Sequence_dataset(numtaxa, 3*numchars);
      
      read_file[299]='\0';
      
		if (single_nuc_file == TRUE) {
			get_DNA=new Read_FASTA;
			DNA_seq=get_DNA->get_dataset(numtaxa, num_DNA_dummy, nuc_file, NUCLEIC);
			delete get_DNA;
		}
		
	
      
      //Does alignment
      for(i=0; i<numtaxa; i++) {
		  (*final_sequence)[i].Assign_name((*current_data)[i].Sequence_name());  
	  
		  if (single_nuc_file == FALSE) {
			  DNA_id=0;
			  get_next_dna_seq((*current_data)[i].Sequence_name(), &current_info, DNAname, TRUE);
			  

			  cout<<"Reading DNA sequence from |"<<DNAname<<"|"<<endl;
			  get_DNA=new Read_FASTA;
			  DNA_seq=get_DNA->get_dataset(dummyt, DNAsize, DNAname, FALSE);
			  delete get_DNA;
			  cout<<"Got :"<<(*DNA_seq)[0].Sequence_size()<<" bases\n";
			  if (DNA_seq==0)
				return(-1);

			  cout<<"Have config file: "<<current_info.check_for_config_file()<<" Type is filenames: "<<current_info.check_for_filenames()<<" Simple option: "<<all_single_exon<<endl;
			  if (all_single_exon==TRUE) {
					cout<<"Assuming all DNA sequence files are already spliced\n";
				  new_taxainfo.numexons=1;
				  new_taxainfo.starts[0]=1;
				  new_taxainfo.ends[0]=(*DNA_seq)[0].Sequence_size();
				  current_info.set_taxa_info((*current_data)[i].Sequence_name(), new_taxainfo.numexons, 
				       new_taxainfo.starts, new_taxainfo.ends);
	    
			  }
			  else {
				  if(current_info.check_for_config_file()==FALSE || 
					 (current_info.check_for_filenames()==TRUE && current_info.check_for_exons()==FALSE)) {
					  //We try to read from the genbank records if we have no config file
	      	 	     
					  get_next_dna_seq((*current_data)[i].Sequence_name(), &current_info, param_file, FALSE);
					  paramin.clear();
					  paramin.open(param_file);
					  cout<<"Current Param file: "<<param_file<<endl;
		
					  cds_in_file=0;
		
					  paramin.getline(read_file, 299);
		
					  while (paramin.eof()!=TRUE && word_match(read_file, 299, "FEATURES", 8)!=1)
						  paramin.getline(read_file, 299);
		
					  while (paramin.eof()!=TRUE) {
						  paramin.getline(read_file, 299);
						  if (word_match(read_file, 299, "  CDS  ", 7)==1) {
							  cds_in_file++;
							  get_coding_region(read_file, paramin, &current_info, &new_taxainfo);
						  }
					  }
					  paramin.close();
		
					  if (cds_in_file>1) {
						  paramin.clear();
						  paramin.open(param_file);
		    
						  paramin.getline(read_file, 299);
		    
						  while (word_match(read_file, 299, "FEATURES", 8)!=1)
							  paramin.getline(read_file, 299);
		    
		    
						  if (find_correct_cds((*current_data)[i].Sequence_name(), &current_info, paramin)!=TRUE) {
							  cerr<<"Could not find any identifer in file "<<param_file<<". Aborting\n";
							  return(-1);
						  }
						  paramin.close();
					  }
					  else
						  current_info.set_taxa_info((*current_data)[i].Sequence_name(), new_taxainfo.numexons, 
													 new_taxainfo.starts, new_taxainfo.ends);
				  } //end reading genbank file
			  }	    
			  curtaxa= current_info.get_taxa_info((*current_data)[i].Sequence_name());
	    

#if defined (DEBUG)  
			  cout<<"Name: "<<(*current_data)[i].Sequence_name()<<" "<<" info file: "<<param_file
			  <<" Num exons: "<<curtaxa->numexons<<" Size: "<<(*current_data)[i].Sequence_size()<<endl<<flush;
			  for(j=0; j<curtaxa->numexons; j++)
				  cout<<"Start: "<<curtaxa->starts[j]<<"  End: "<<curtaxa->ends[j]<<" ";
			  cout<<endl<<endl;
#endif	 
		  }
		  else {
			  DNA_id=-1;
			  for(j=0; j<DNA_seq->Num_sequences(); j++) {
				  if (strcmp((*final_sequence)[i].Sequence_name(), (*DNA_seq)[j].Sequence_name()) == 0) DNA_id=j;
			  }
			  
			  if (DNA_id == -1) {
				  cerr<<"Cannot find sequence "<< (*final_sequence)[i].Sequence_name()<<" in the provided nucleotide file...Aborting\n";
				  return(-1);
			  }
			  
			  new_taxainfo.numexons=1;
			  new_taxainfo.starts[0]=1;
			  new_taxainfo.ends[0]=(*DNA_seq)[DNA_id].Sequence_size();
			  current_info.set_taxa_info((*current_data)[i].Sequence_name(), new_taxainfo.numexons, 
										 new_taxainfo.starts, new_taxainfo.ends);
			  curtaxa= current_info.get_taxa_info((*current_data)[i].Sequence_name());
			  
#if defined (DEBUG)  
			  cout<<"Name: "<<(*current_data)[i].Sequence_name()<<" "<<" ID #"<<DNA_id
			  <<" Num exons: "<<curtaxa->numexons<<" Size: "<<(*current_data)[i].Sequence_size()<<endl<<flush;
			  for(j=0; j<curtaxa->numexons; j++)
				  cout<<"Start: "<<curtaxa->starts[j]<<"  End: "<<curtaxa->ends[j]<<" ";
			  cout<<endl<<endl;
#endif	
			  
		  }

		  exon=0;
		  exonpos=0;

		  //Making the new sequence
		  for (j=0; j<numchars; j++){
			  for (k=0; k<3; k++) {
				  if ((exonpos+curtaxa->starts[exon])>curtaxa->ends[exon]) {
					  ++exon;
					  exonpos=0;
				  }
	
				  if ((*current_data)[i][j]==readchar_to_aa('-')) {
					  (*final_sequence)[i].Assign_site(3*j+k, readchar_to_base('-'));
					  //cout<<"-";
				  }
				  else  { 
					  if ((curtaxa->starts[exon]+exonpos-1) < (*DNA_seq)[DNA_id].Sequence_size()) {
						  (*final_sequence)[i].Assign_site(3*j+k, (*DNA_seq)[DNA_id][curtaxa->starts[exon]+exonpos-1]);
		  // cout<<num_to_base((*DNA_seq)[0][curtaxa->starts[exon]+exonpos-1]);
					  }
					  else
						  (*final_sequence)[i].Assign_site(3*j+k, readchar_to_base('N'));
					  exonpos++;
				  }
			  }  
		  } //End for to create final sequence
	  
	  // cout<<endl;
	  paramin.close();
	}  //End for-looping through taxa


        new_protein_sequences=new Sequence_dataset(numtaxa, numchars, PROTEIN);

	cout<<"Testing sequence alignment of "<<numtaxa<<" taxa and "<<numchars<<" characters\n";
	for(i=0; i<numtaxa; i++)
	  {
	    (*new_protein_sequences)[i].Assign_name((*final_sequence)[i].Sequence_name());
	    
	    for(j=0; j<numchars; j++)
	      {
		if (num_to_base((*final_sequence)[i][j*3])!='-')
		  {
		    codon[0]=(*final_sequence)[i][j*3];
		    codon[1]=(*final_sequence)[i][j*3+1];
		    codon[2]=(*final_sequence)[i][j*3+2];

		    (*new_protein_sequences)[i].Assign_site(j, readchar_to_aa(current_code->return_AA(codon)));
		}
	      else
		(*new_protein_sequences)[i].Assign_site(j, readchar_to_aa('-'));
	    }
	 
	  }
      for(i=0; i<numtaxa; i++) 
		if (!((*new_protein_sequences)[i]==(*current_data)[i])) 
		  backcheck_ok=FALSE;
      
     
	write_seq->write_dataset(numtaxa, numchars*3, final_sequence);
		if (single_nuc_file == TRUE) delete DNA_seq;
	
	if (backcheck_ok == FALSE)
	  cout<<"ERROR: Aligned DNA sequence doesn't match orignal protein sequence\n";
      delete final_sequence;
    
    } 
  
  
} //End main











