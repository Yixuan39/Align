//Copyright 1999-2002 Gavin Conant

#include <fstream>
#include <string.h>
#include <iostream>
#include "algndna_funcs.h"

using namespace::std;


Gene_struct& Gene_struct::operator=(Gene_struct & assign_from)
{
  int i;

  strcpy(name, assign_from.name);
  numexons=assign_from.numexons;

  for(i=0; i<40; i++)
    {
      starts[i]=assign_from.starts[i];
      ends[i]=assign_from.ends[i];
    }
  return(*this);
}  //End Gene_struct::operator= 



int Gene_struct::operator==(Gene_struct test)
{
  if (strcmp(name, test.name)==0)
    return(0);
  else
    return(-1);
}


Taxa_File& Taxa_File::operator=(Taxa_File & assign_from)
{
  id_len=assign_from.id_len;
  strcpy(file_prefix, assign_from.file_prefix);
  strcpy(Taxa_name, assign_from.Taxa_name);
  strcpy(id, assign_from.id);
  return(*this);
} 


int Taxa_File::operator==(Taxa_File test)
{
  if (strcmp(Taxa_name, test.Taxa_name)==0)
    return(0);
  else
    return(-1);
}



String_for_list& String_for_list::operator=(String_for_list & assign_from)
{
  strcpy(string, assign_from.string);
  return(*this);
}

int String_for_list::operator==(String_for_list test)
{
  if (strcmp(string, test.string)==0)
    return(0);
  else
    return(-1);
}




Data_info::Data_info()
{
  have_filenames=FALSE;
  config_file=FALSE;

  num_gene_ids=0;

  seq_filenames=new List<Taxa_File>();
  taxa_data = new List<Gene_struct>();
}


BOOL Data_info::check_for_filenames()
{
  return(have_filenames);
}

BOOL Data_info::check_for_exons()
{
  return(have_exons);
}


BOOL Data_info::check_for_config_file()
{
  return(config_file);
}
  
char * Data_info::get_global_paramfile()
{
  return(global_param_file);
}


char * Data_info::get_gene_id_num(int num)
{
  if (num<num_gene_ids && num>=0)
    return(gene_ids[num]);
  else
    {
      cerr<<num<<" is an invalid gene identification number\n";
      return(0);
    }
}  //End Data_info::get_gene_id


int Data_info::get_num_gene_ids()
{
  return(num_gene_ids);
}


Gene_struct* Data_info::get_taxa_info(char *Taxaname)
{
  strcpy(dummy.item()->name, Taxaname);

  requested_taxa=taxa_data->find_data(dummy.item());

  if (requested_taxa==0)
    {  
      cerr<<"Couldn't find taxa "<<Taxaname<<endl;
      return(0);
    }
  else
      return(requested_taxa->item());
}  //End Data_info::get_taxa_info
  


void Data_info::set_taxa_filename(char *Taxa_name, char *filename, char *id)
{

  strcpy(assigner.Taxa_name, Taxa_name);
  strcpy(assigner.file_prefix, filename);
  strcpy(assigner.id, id);
  assigner.id_len=strlen(id);

  seq_filenames->add_to_list(assigner);
}



void Data_info::set_taxa_info(char name[], long numexons, long starts[], long ends[])
{
  int i;

  strcpy(dummy.item()->name, name);
  dummy.item()->numexons=numexons;

  for (i=0; i<numexons; i++)
    {
      dummy.item()->starts[i]=starts[i];
      dummy.item()->ends[i]=ends[i];
    }
  
  taxa_data->add_to_list((*dummy.item()));
}



void Data_info::set_gene_ids(int num, char **id)
{
  int i;
  num_gene_ids=num;

  gene_ids=new char*[num_gene_ids];

  for (i=0; i<num_gene_ids; i++)
    {
      gene_ids[i]=new char[20];
      strcpy(gene_ids[i], id[i]);
    }
}  //End Data_info::set_gene_id



void Data_info::set_have_filenames(BOOL set)
{
  have_filenames=set;
}
void Data_info::set_have_exons(BOOL set)
{
  have_exons=set;
}


void Data_info::set_global_paramfile(char *name)
{
  strcpy(global_param_file, name);
  config_file=TRUE;
}


void Data_info::taxa_file_name(char *taxa_name, char *filename)
{
  strcpy(dummy_file.item()->Taxa_name, taxa_name);

  file_search=seq_filenames->find_data(dummy_file.item());

  strcpy(filename, file_search->item()->file_prefix);
}



Data_info::~Data_info()
{
  int i;

  if (num_gene_ids!=0)
    {
      for (i=0; i<num_gene_ids; i++)
	delete[] gene_ids[i];
      delete[] gene_ids;
    }
}



void parse_args(int argc, char **argv, char * param_file, BOOL &single_exon, BOOL &single_nuc_file, char *nuc_file, Exchange *curr_exchange)
{
	int i=0, j;

	curr_exchange->set_genetic_code(UNIVERSAL);
	param_file[0] = '\0';
	nuc_file[0]='\0';
	single_nuc_file=FALSE;

	//Get rid of non-optional args
	while((i< argc) && (argv[i][0] != '-')) i++;
	

  
	while (i<argc) {
		switch(argv[i][1]) {
			case 'f':
			case 'F':
				j=3;
				while(argv[i][j]!='\0' && (j<strlen(argv[i]))) {
					param_file[j-3]=argv[i][j];
					j++;
				}
				param_file[j-3]='\0';
				i++;
				break;
			case 's':
			case 'S':
				single_exon=TRUE;
				i++;
				break;
			case 'c':
			case 'C':
				single_nuc_file=TRUE;
				j=3;
				while(argv[i][j]!='\0' && (j<strlen(argv[i]))) {
					nuc_file[j-3]=argv[i][j];
					j++;
				}
				nuc_file[j-3]='\0';
				i++;
				break;
			case 'v':
			case 'V':
				switch (argv[i][3]) {
					case 'u':
					case 'U':
							curr_exchange->set_genetic_code(UNIVERSAL);
					break;
					case 'v':
					case 'V':
						curr_exchange->set_genetic_code(VERT_MITO);
						break;
					case 'y':
					case 'Y':
						curr_exchange->set_genetic_code(YEAST_MITO);
						break;
					case 'm':
					case 'M':
						if ((argv[i][4] == 'm') || (argv[i][4] == 'M'))
							curr_exchange->set_genetic_code(MOLD_MITO);
						else
							curr_exchange->set_genetic_code(MYCOPLASMA);
						break;
						break;
					case 'i':
					case 'I':
						curr_exchange->set_genetic_code(INVERT_MITO);
						break;
					case 'c':
					case 'C':
						curr_exchange->set_genetic_code(CILIATE_NUC);
						break;
					case 'e':
					case 'E':
						curr_exchange->set_genetic_code(ECHINO_MITO);
						break;
				}
			i++;
			break;     
		}
	}
}





void get_config_file_data(Data_info *info)
{
  BOOL first=TRUE;
  long nexons, estarts[40], eends[40];
  char config_type[10], taxaname[100], taxafilename[100], taxaid[20], dump[5];
  ifstream paramin;



  paramin.open(info->get_global_paramfile());
 
  paramin>>config_type;
	  
  if(strcmp("Filenames", config_type)==0)
    {
      info->set_have_filenames(TRUE);
      info->set_have_exons(FALSE);
     
      paramin>>taxaname>>taxafilename>>taxaid;
      info->set_taxa_filename(taxaname, taxafilename, taxaid);
      

      while(!paramin.eof())
	{
	  paramin>>taxaname>>taxafilename>>taxaid;
	  info->set_taxa_filename(taxaname, taxafilename, taxaid);
	}
      cout<<"Read filenames from config file "<<info->get_global_paramfile()<<endl; 
    }
  else
  {
    if(strcmp("Exons+Filenames", config_type)==0) {
      cout<<"Reading filenames and exons\n";
      info->set_have_filenames(TRUE);
      info->set_have_exons(TRUE);
     
      while(!paramin.eof())
	{
	 
	  paramin>>taxaname>>taxafilename;
	 
	  
	  cout<<"Reading exons for "<<taxaname<<"|"<<taxafilename<<endl;
	  nexons=0;
	  
	  paramin.get(dump[0]);
	  
	  while(dump[0]!='\n' && !(paramin.eof()) )
	    {
	      paramin>>estarts[nexons]>>eends[nexons];
	      nexons++;
	      paramin.get(dump[0]);
	    }
	  
	  info->set_taxa_info(taxaname, nexons, estarts, eends);
	  info->set_taxa_filename(taxaname, taxafilename, "");
	}
      cout<<"Read exons from config file "<<info->get_global_paramfile()<<endl;
    }
    else {
      info->set_have_filenames(FALSE);
      info->set_have_exons(TRUE);
     
      while(!paramin.eof())
	{
	  if (!first) 
	    paramin>>taxaname;
	  else {
	    first=FALSE;
	    strcpy(taxaname, config_type);
	  }
	  
	  cout<<"Reading exons for "<<taxaname<<"|"<<endl;
	  nexons=0;
	  
	  paramin.get(dump[0]);
	  
	  while(dump[0]!='\n' && !(paramin.eof()) )
	    {
	      paramin>>estarts[nexons]>>eends[nexons];
	      nexons++;
	      paramin.get(dump[0]);
	    }
	  
	  info->set_taxa_info(taxaname, nexons, estarts, eends);
	}
      cout<<"Read exons from config file "<<info->get_global_paramfile()<<endl;
    }
  }
    paramin.close();
}




void get_next_dna_seq (char *Taxa_name, Data_info *cinfo, char *filename, BOOL dna_file)
{
  int i;
  ifstream testfile;
 
  if (cinfo->check_for_filenames()==FALSE)
    {
      strcpy(filename, Taxa_name);

      if (dna_file==TRUE)
	strcat(filename, ".fas");
      else
	strcat(filename, ".all");
      
      testfile.open(filename);

      cout<<"Checking filename "<<filename<<endl;

      i=0;
      while (testfile.fail() && i<cinfo->get_num_gene_ids())
	{
	  strcpy(filename, Taxa_name);
	  strcat(filename, "_");
	  strcat(filename, cinfo->get_gene_id_num(i));
	  
	  
	  if (dna_file==TRUE)
	    strcat(filename, ".fas");
	  else
	    strcat(filename, ".all");
	  
	  testfile.open(filename);
	  i++;
	}      
      if (!testfile.fail())
	testfile.close();
      else
	cerr<<"Could not find file for taxa "<<Taxa_name<<endl;
    }


  else
    {
      cinfo->taxa_file_name(Taxa_name, filename);
      if (dna_file==TRUE)
	strcat(filename, ".fas");
      else
	strcat(filename, ".all");  
    }
  
}




void update_word(char *old_word, char newchar, int word_length)
{
  int m;
  
  for(m=0; m<word_length; m++)
    old_word[m]=old_word[m+1];

  old_word[word_length-1]=newchar;
}



void get_coding_region(char cds_line[],  ifstream &fin, Data_info *cinfo,  Gene_struct *curtaxa)
{
  int point=0;
 

  while (cds_line[point]!='S')
    point++;
  while (!((cds_line[point]>=48 && cds_line[point]<=57)  || (cds_line[point]== 'j' || cds_line[point] == 'J')))
    point++;
 
  if (cds_line[point]== 'j' || cds_line[point] == 'J')
    {
      curtaxa->numexons=0;
      point+=4;
    
      while (cds_line[point]!=')')
	{
	  point++;   
	  if (cds_line[point]=='\n' || cds_line[point]=='\0')
	    {
	      fin.getline(cds_line, 299);
	      point=0;
	      while (!(cds_line[point]>=48 && cds_line[point]<=57))
		point++;
	    }
	
	  if (cds_line[point]=='\n' || cds_line[point]=='\0')
	    fin.getline(cds_line, 299);

	  curtaxa->starts[curtaxa->numexons]=cds_line[point]-48;
	  point++;

	  while (cds_line[point] != '.')
	    {
	      curtaxa->starts[curtaxa->numexons]*=10;
	      curtaxa->starts[curtaxa->numexons]+=cds_line[point++]-48;
	    }

	  point+=2;
	  curtaxa->ends[curtaxa->numexons]=cds_line[point++]-48;
		
	  cout<<"Join ends :"<<cds_line[point]<<cds_line[point+1]<<cds_line[point+2]<<cds_line[point+3]<<endl;

	  while (cds_line[point] != ',' && cds_line[point] != ')')
	    {
	      curtaxa->ends[curtaxa->numexons]*=10;
	      curtaxa->ends[curtaxa->numexons]+=cds_line[point++]-48;
	    }	
	  curtaxa->numexons++;
	}
    }
	      
	      
  else
    { 
      curtaxa->numexons=1;
      
      while (!(cds_line[point]>=48 && cds_line[point]<=57))
	point++;
		    
		  
      curtaxa->starts[0]=cds_line[point++]-48;
	
      while (cds_line[point] != '.')
	{
	  curtaxa->starts[0]*=10;
	  curtaxa->starts[0]+=cds_line[point++]-48;
	}
		  
      point+=2;


      curtaxa->ends[0]=cds_line[point++]-48;
	
      while (cds_line[point] != '\n' && cds_line[point] != ' ' && cds_line[point] !='\0')
	{
	  curtaxa->ends[0]*=10;
	  curtaxa->ends[0]+=cds_line[point++]-48;
	}
    }
  
}  //End get_coding_region





BOOL find_an_id_in_file(Data_info *cinfo, char *first_line)
{
  int i;
  char read_file[300];
  BOOL found_in_line=FALSE;

  i=0;
  while(i<cinfo->get_num_gene_ids() && found_in_line==FALSE)
    { 
      if (word_match(first_line, 299, cinfo->get_gene_id_num(i), strlen(cinfo->get_gene_id_num(i)))==1)
	found_in_line=TRUE;
      i++;
    }
  return(found_in_line);
}




BOOL find_correct_cds(char *Taxaname, Data_info *cinfo, ifstream &fin)
{
  char read_file[300];
  BOOL found=FALSE;
  Gene_struct curtaxa;
 

  fin.getline(read_file, 299);
  while (word_match(read_file, 299, "  CDS  ", 7)!=1 && found==FALSE && !fin.eof()) 
    {
      if (word_match(read_file, 299, "GENE", 4)==1 ||
	  word_match(read_file, 299, "PRODUCT", 7)==1 ||
	  word_match(read_file, 299, "NOTE", 4)==1)
	found=find_an_id_in_file(cinfo, read_file);
	fin.getline(read_file, 299);
    } 
  

 while (word_match(read_file, 299, "  CDS  ", 7)!=1 && !fin.eof()) 
   fin.getline(read_file, 299);

 if (!fin.eof())
   get_coding_region(read_file, fin, cinfo, &curtaxa);


  if (found!=TRUE)
    {
      while (word_match(read_file, 299, "TRANSLATION", 11)!=1 && found==FALSE  && !fin.eof())
	{
	  if (word_match(read_file, 299, "GENE", 4)==1 || 
	      word_match(read_file, 299, "PRODUCT", 7)==1 ||
	      word_match(read_file, 299, "NOTE", 4)==1)
	    found=find_an_id_in_file(cinfo, read_file);
	  if (found==FALSE)
	    fin.getline(read_file, 299);
	
	}
    }
 
  if (found!=TRUE && !fin.eof())  //If the data isn't in the first coding block, recurse to another
    found=find_correct_cds(Taxaname, cinfo, fin);
    
  cinfo->set_taxa_info(Taxaname, curtaxa.numexons, curtaxa.starts, curtaxa.ends);
  return(found);
}
  



