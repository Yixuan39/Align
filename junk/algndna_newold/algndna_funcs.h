//Copyright 1999-2002 Gavin Conant

#include "gen_dna_funcs.h"
#include "read_seq.h"
#include "exchange.h"
#include "linked_list.cpp"


class Gene_struct {
 public:
  char name[700];
  long numexons, starts[40], ends[40];
  Gene_struct& operator=(Gene_struct & assign_from); 
  int operator==(Gene_struct test);
};


class Taxa_File {
 public:
  int id_len;
  char file_prefix[90], Taxa_name[100], id[20];
  Taxa_File& operator=(Taxa_File & assign_from); 
  int operator==(Taxa_File test);
};


class String_for_list {
 public:
  char string[20];
  String_for_list& operator=(String_for_list & assign_from);
  int operator==(String_for_list test);
};


class Data_info
{
 public:
  Data_info();
  BOOL check_for_filenames();
  BOOL check_for_exons();
  BOOL check_for_config_file();
  char * get_global_paramfile();
  char * get_gene_id_num(int num);
  int get_num_gene_ids();
  Gene_struct* get_taxa_info(char *Taxaname);

  void set_taxa_filename(char *Taxa_name, char *filename, char *id);
  void set_taxa_info(char name[], long numexons, long starts[], long ends[]);

  void set_gene_ids(int num, char **id);
  void set_have_filenames(BOOL set);
  void set_have_exons(BOOL set);
  void set_global_paramfile(char *name);
  void taxa_file_name(char *taxa_name, char *filename);
  ~Data_info();
 private:
  int num_gene_ids;
  char global_param_file[100], **gene_ids;
  BOOL have_filenames, config_file, have_exons;

  Taxa_File assigner;
  List_Element<Gene_struct> dummy, *requested_taxa;
  List_Element<Taxa_File> dummy_file, *file_search;
  List<Taxa_File> *seq_filenames;
  List<Gene_struct> *taxa_data;
  
};


void parse_args(int argc, char **argv, char * param_file, BOOL &single_exon, BOOL &single_nuc_file, char *nuc_file, Exchange *curr_exchange);
void get_config_file_data(Data_info *info);
void update_word(char *old_word, char newchar, int word_length);
void get_next_dna_seq(char *Taxa_name, Data_info *cinfo, char *filename, BOOL dna_file);
void get_coding_region(char *cds_line, ifstream &fin, Data_info *cinfo, Gene_struct *curtaxa);
BOOL find_an_id_in_file(Data_info *cinfo, char *first_line);
BOOL find_correct_cds(char *Taxaname, Data_info *cinfo, ifstream &fin);

