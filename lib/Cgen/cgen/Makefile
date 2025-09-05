#!/bin/sh
#
#--------------------------------------------------------------------------*
# Copyright 2002 by Paul D. Kundarewich, Michael Hutton, Jonathan Rose     *
# and the University of Toronto. 										   *
# Use is permitted, provided that this attribution is retained  		   *
# and no part of the code is re-distributed or included in any commercial  *
# product except by written agreement with the above parties.              *
#                                                                          *
# For more information, contact us directly:                               *
#	 Paul D. Kundarewich (paul.kundarewich@utoronto.ca)					   *
#    Jonathan Rose  (jayar@eecg.toronto.edu)                               *
#    Mike Hutton  (mdhutton@cs.toronto.edu, mdhutton@eecg.toronto.edu)     *
#    Department of Electrical and Computer Engineering                     *
#    University of Toronto, 10 King's College Rd.,                         *
#    Toronto, Ontario, CANADA M5S 1A4                                      *
#    Phone: (416) 978-6992  Fax: (416) 971-2286                            *
#--------------------------------------------------------------------------*/


#
# Makefile for gen package.
#

GEN	= .

#		 gen_version.cpp gen_version.o
SRC		= main.cpp gen_control.cpp gen_options.cpp rand.cpp level_node.cpp \
			spec.cpp graph.cpp edges_and_nodes.cpp sequential_level.cpp move.cpp \
			matrix3D.cpp degree.cpp node_splitter.cpp simple_edge_assignment.cpp \
			post_process.cpp gen_drawer.cpp final_improvement.cpp util.cpp matrix.cpp fp.cpp \
			graph_medic.cpp
OBJ		= main.o gen_control.o gen_options.o rand.o level_node.o matrix3D.o \
			spec.o graph.o edges_and_nodes.o sequential_level.o move.o degree.o node_splitter.o \
			simple_edge_assignment.o post_process.o gen_drawer.o final_improvement.o util.o matrix.o fp.o \
			graph_medic.o
HDR		= gen.h gen_options.h gen_control.h matrix.h \
			types.h rand.h util.h output.h level_node.h spec.h graph.h edges_and_nodes.h \
			sequential_level.h move.h gen_version.h matrix3D.h degree.h node_splitter.h \
			simple_edge_assignment.h post_process.h gen_drawer.h final_improvement.h fp.h \
			graph_medic.h

PRE_COMP= 


EXE	= cgen


SHELL	= sh
CC	= gcc
YACC	= bison
LEX	= flex
RM	= /bin/rm -f
CP	= cp
MV	= mv

# Determine architecture.  On a msdos system you should make an "arch.bat"
# which echos "i486" or "i686" or something similar and make an "uname.bat" 
#that echos DOS

ARCH	= $(shell arch)
OS		= $(shell uname)

ifeq ($(OS),SunOS) 
	OS = solaris
else 
	ifeq ($(OS),Linux) 
	else
		ifeq ($(OS),DOS)
		else
			ARCH = unknown_arch
		endif
	endif
endif


INCLUDE = -I$(GEN) 
#CFLAGS	= $(INCLUDE) -D$(ARCH) -D$(OS) -DGEN -pedantic -Wall -g
#LDFLAGS	= -L. -g
CFLAGS	= $(INCLUDE) -D$(ARCH) -D$(OS) -DGEN -pedantic -Wall -O3
LDFLAGS	= -L.
LIBS = -lstdc++ -lm




#----------------------------------------------------------------------
# the $@ and $< are automatic variables
# $^ is a list of all the dependencies of the variable
# $@ is the object file name
# $< is the source file name which should be 
#	the 1st dependency after the :
# 


#all: 	$(ARCH) $(TARGET)

$(EXE): $(OBJ)
	$(CC) -c $(CFLAGS) "-DCURRENT_DATE=\"`date`\"" gen_version.cpp -o gen_version.o 
	$(CC) $(OBJ) gen_version.o -o $(EXE) $(LDFLAGS) $(LIBS) 


%.o : %.cpp
	$(CC) -c $(CFLAGS) -o $@ $< 

gen_version.o : gen_version.h gen_version.cpp
	$(CC) -c $(CFLAGS) -o gen_version.o gen_version.cpp

clean:
	$(RM) $(OBJ)
	$(RM) gen_version.o
	$(RM) $(EXE)

solaris:

Linux:

DOS:
	@echo "I have not attempted to compile the code with djgpp but it should work."
	@echo "If you want a windows version I have included a visual c++ makefile in this directory."

unknown_arch:
	@echo "Unknown architecture $(ARCH), assuming sun4."


# END 
