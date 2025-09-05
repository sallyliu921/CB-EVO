
# CCirc

Clustered And Iterative Synthetic Circuit Generation

Designed by: `FPGA CAD Software From the Toronto FPGA Research Group`

## Introduction

This archive contains CCirc a software tool to 
characterize physical properties of circuits and CGen a software
tool to generate synthetic circuit netlists in BLIF or VHDL format
ready to be place and routed.

There is a detailed manual covering everything from how to compile the programs
to how to use them in user_manual.pdf in the doc/ directory.

To quickly see CCirc characterize a sample circuit from the MCNC 
benchmark suite type:

cd bin
ccirc diffeq.blif --partitions 8

This will produce a diffeq.stats which contains the 
characterization results.

To quickly see CGen generate a synthetic circuit from 
a sample characterization file type:

cd bin
cgen diffeq.stats --vhdl_output


This will produce a diffeq_clone.vhd file that is ready to be placed and
routed.


If you liked the software and/or have a suggestion for its improvement send me and email.

cheers,

-- Paul Kundarewich, March 16, 2003

paul.kundarewich@utoronto.ca
http://www.eecg.toronto.edu/~jayar/software/Cgen/Cgen.html

==============================================================================
Contents of the archive:

Readme.txt:   This file.
user_manual.pdf: PDF version of the manual.

./doc		A user manual for the CCirc and CGen.


./ccirc

  *.cpp, *.h, *.l, *.y:  Source code for CCirc.


  diffeq.blif:  A sample netlist file from the MCNC bencmark set.  The logic
            block contains 4-input LUTs and FFs.

  Makefile

./cgen

  *.cpp, *.h:  Source code for CGen.


  Makefile

./bin		

	ccirc:	CCirc's executable.

	cgen:	CGen's executable.

  	diffeq.blif:  A sample tech-mapped MCNC bencmark circuit.

  	diffeq.stats:  A sample .stats file that can be fed into cgen to generate
				   a synthetic circuit.

./debug	

	ccirc-debug:	CCirc's executable compiled with debugging information.
	cgen-debug:		CGen's executable compiled with debugging information.

==============================================================================

## How to compile?

Please refer to the documentation: [user-manual](doc/user_manual.pdf)

Step 1: Make sure your version of gcc is at least 2.95.

Step 2: Download the hMetis library and save it somewhere.

Step 2: In the file Makefile set PARTITION to point to where the hMetis library is located.

Step 4: type make.

Step 5: Look at your new executable ccirc.