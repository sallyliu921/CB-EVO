#---------------------------------------------------------------------------
#   Copyright 2002 by Paul D. Kundarewich, Michael Hutton, Jonathan Rose    
#   and the University of Toronto. 											
#   Use is permitted, provided that this attribution is retained  			
#   and no part of the code is re-distributed or included in any commercial	
#   product except by written agreement with the above parties.             
#                                                                           
#   For more information, contact us directly:                              
#  	   Paul D. Kundarewich (paul.kundarewich@utoronto.ca)					
#      Jonathan Rose  (jayar@eecg.toronto.edu)                              
#      Mike Hutton  (mdhutton@cs.toronto.edu, mdhutton@eecg.toronto.edu)    
#      Department of Electrical and Computer Engineering                    
#      University of Toronto, 10 King's College Rd.,                        
#      Toronto, Ontario, CANADA M5S 1A4                                     
#      Phone: (416) 978-6992  Fax: (416) 971-2286                           
#---------------------------------------------------------------------------


#/bin/sh

SHELL	= sh
CC 		= gcc -g
YACC	= bison
LEX		= flex
CP  	= cp
MV		= /bin/mv

CIRC	=	.
#PARTITION = ../hmetis-1.5-linux
#PARTITION = ../hmetis-1.5-sun4u-USparc

OBJ = main.o options.o util.o edges_and_nodes.o cluster.o sequential_level.o circuit.o circ_control.o symbol_table.o graph_constructor.o parse_blif.o token_blif.o lut.o graph_medic.o cycle_breaker.o drawer.o node_partitioner.o matrix.o delay_leveler.o degree_info.o statistic_reporter.o wirelength_character.o rand.o rnum.o
SRC	= main.cpp options.cpp util.cpp lut.cpp edges_and_nodes.cpp cluster.cpp sequential_level.cpp circuit.cpp circ_control.cpp symbol_table.cpp graph_constructor.cpp parse_blif.cpp token_blif.cpp graph_medic.cpp cycle_breaker.cpp drawer.cpp  node_partitioner.cpp matrix.cpp delay_leveler.cpp  degree_info.cpp statistic_reporter.cpp wirelength_character.cpp rand.cpp rnum.cpp
HDR	= circ.h output.h util.h lut.h options.h edges_and_nodes.h cluster.h sequential_level.h circuit.h circ_control.h symbol_table.h graph_constructor.h parse_blif.h graph_medic.h cycler_breaker.h drawer.h matrix.h node_partitioner.h delay_leveler.h degree_info.h statistic_reporter.h wirelength_character.h rand.h circ_version.h rnum.h

# The -I and -L options are directory search options
# The -I option says search this directory for include files
# The -L option says search this directory for files to link
# The -lm option includes the math libraries
# The -ll includes the lex libraries

INCLUDE	= -I$(CIRC)

CFLAGS = $(INCLUDE) -Wall -pedantic --std=c++11 #-m32 
LDFLAGS = -L$(PARTITION) -L. -lm

LEXFLAGS = 
YACCFLAGS = -d

LIBS	= -lstdc++ -lm #-lhmetis

EXE = ccirc

$(EXE): $(OBJ)
	$(CC) -c $(CFLAGS) "-DCURRENT_DATE=\"`date`\"" circ_version.cpp -o circ_version.o 
	$(CC) $(CFLAGS) $(OBJ) circ_version.o -o $(EXE) $(LDFLAGS) $(LIBS)

# the $@ and $< are automatic variables
# $^ is a list of all the dependencies of the variable
# $@ is the object file name
# $< is the source file name which should be 
#	the 1st dependency after the :
# 

%.o : %.cpp
	$(CC) -c $(CFLAGS) -o $@ $<

parse_blif.o: parse_blif.cpp token_blif.l parse_blif.y token_blif.cpp
	$(CC) -c $(CFLAGS) -o $@ $<

token_blif.o: token_blif.cpp parse_blif.y token_blif.l
	$(CC) -c $(CFLAGS) token_blif.cpp -o token_blif.o

parse_blif.cpp: parse_blif.y token_blif.l
	$(YACC) $(YACCFLAGS) parse_blif.y
	$(MV) parse_blif.tab.c parse_blif.cpp
	$(MV) parse_blif.tab.h parse_blif.h

token_blif.cpp: token_blif.l parse_blif.h
	$(LEX) $(LEXFLAGS) token_blif.l
	$(MV) lex.yy.c token_blif.cpp

circ_version.o : circ_version.h circ_version.cpp
	$(CC) -c $(CFLAGS) -o circ_version.o circ_version.cpp

clean:
	$(RM) $(OBJ)
	$(RM) parse_blif.cpp
	$(RM) parse_blif.h
	$(RM) token_blif.cpp
	$(RM) circ_version.o
	$(RM) $(EXE)
