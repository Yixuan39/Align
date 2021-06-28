# Master Makefile for Sequence Utilities for Genome analysis
# Type make to build
include make.inc

O = o
SRC_DIR = ./

LIBRARY_DIR = -L./ 
INCLUDE = -I. 
CFLAGS = -g                 		                     # compiler switches to be applied to every module
OPTIM_SPEED = -O3             	  	                     # switches that give speed priority over size
OPTIM_SIZE = -O1              	  	                     # switches that give size priority over speed
MATH_LIB = -lm                  	                     # math library 
DNAFUNCS_LIB= -ldnafuncs

OPTIONS = $(CFLAGS) $(INCLUDE)


all:  algndna_new 

DNALIBRARY = libdnafuncs.a

DNALIBRARY_OBJS =		\
	read_seq.$(O)	\
	exchange.$(O) \
	gen_code.$(O) \
	gen_dna_funcs.$(O) \
	score_matrix.$(O) \
	write_seq.$(O) 

libdnafuncs : $(DNALIBRARY) 

$(DNALIBRARY): $(DNALIBRARY_OBJS) 
	ar  cr $(DNALIBRARY) $(DNALIBRARY_OBJS)
	ranlib $(DNALIBRARY)


ALGNDNA_NEW_OBJS = algndna_new.$(O)	\
			algndna_funcs.$(O) 


algndna_new: $(ALGNDNA_NEW_OBJS) $(DNALIBRARY)
	$(CC) $(LINUX_BUILD)  $(LIBRARY_DIR) \
	 -o algndna_new $(OPTIONS) $(ALGNDNA_NEW_OBJS) $(DNAFUNCS_LIB)


%.o: %.cpp
	$(CC)  $(OPTIONS) $(OPTIM_SPEED) -c $<


%.o: %.c
	$(CC)  $(OPTIONS) $(OPTIM_SPEED) -c $<





