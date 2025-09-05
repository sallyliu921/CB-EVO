/*--------------------------------------------------------------------------*
 * Copyright 2002 by Paul D. Kundarewich, Michael Hutton, Jonathan Rose     *
 * and the University of Toronto. 											*
 * Use is permitted, provided that this attribution is retained  			*
 * and no part of the code is re-distributed or included in any commercial	*
 * product except by written agreement with the above parties.              *
 *                                                                          *
 * For more information, contact us directly:                               *
 *	  Paul D. Kundarewich (paul.kundarewich@utoronto.ca)					*
 *    Jonathan Rose  (jayar@eecg.toronto.edu)                               *
 *    Mike Hutton  (mdhutton@cs.toronto.edu, mdhutton@eecg.toronto.edu)     *
 *    Department of Electrical and Computer Engineering                     *
 *    University of Toronto, 10 King's College Rd.,                         *
 *    Toronto, Ontario, CANADA M5S 1A4                                      *
 *    Phone: (416) 978-6992  Fax: (416) 971-2286                            *
 *--------------------------------------------------------------------------*/



 /*
 *   Basic parser for BLIF.
 *   records it in generic graph-like data-structures for use by
 *   circuit analysis package.
 *
 *  P. Kundarewich, Sept 2000  
 *
 *  RESTRICTIONS:
 *
 *  - many BLIF things not supported, see "not_supported" and token_blif.l.
 *
 *  - for BLIF parsing, if opt.k is not set, we assume that the netlist 
 *    is not tech mapped and no analysis of LUTs is being done, so we 
 *    don't store them.  A warning is generated.
 *  
 */

/* ---------------------------------------------------------------------- */
/* ------------------------ PREAMBLE DECLARATIONS ----------------------- */
/* ---------------------------------------------------------------------- */

%{	

#ifdef VISUAL_C
#define MSDOS		// fixes a bug with flex and bison for win32
#endif

#include "lut.h"		    /* typedefs for function storage     */
#include "circuit.h"
#include "graph_constructor.h"
//#include <string>
#include <string.h>

extern char *yytext;
/* Default stack-size in BISON is 10000, I needed 20000 to read clma,    
 * and 30000 for spla.  May have to go bigger still for other circuits. 
 * This is the maximum that the parse stack is allowed to grow to.  It 
 * should be some small constant (like 4) times the number of .names  
 * and .latch specs                                                  
 *
 * Increased to 300000 for gen circuits.
 */
#define YYMAXDEPTH 300000

static string      	g_the_model_name; 
static string		g_latchin, g_latchout, g_latchclk;	
static LUT *		g_current_lut;
static string		g_entry;	    	// text of the current cube 
static VALUE_TYPE	g_value;	    	// value of the current cube
extern string		g_mytext;
static bool			g_got_graph = false;
		long		g_linenum = 1;    	/* in input file, for error messages */

static VARIABLE_STACK_TYPE * 	g_variable_name_stack;
GRAPH_CONSTRUCTOR * 			g_graph_constructor;	

int 	yylex();
int 	yyerror(char * error_text);
static void _init(string & model_name);
static void _cleanup();

extern OPTIONS * g_options;

#define YYDEBUG 1

%}	

/* ---------------------------------------------------------------------- */
/* ------------------------------ TOKEN LIST ---------------------------- */
/* ---------------------------------------------------------------------- */

%token MODEL_TOKEN
%token BLIFEND_TOKEN
%token INPUTS_TOKEN
%token OUTPUTS_TOKEN
%token CLOCK_TOKEN
%token NOTSUPPORTED_TOKEN
%token NAMES_TOKEN
%token LATCH_TOKEN


/* Generic tokens */
%token NEWLINE_TOKEN
%token STRING_TOKEN
%token ILLEGAL_TOKEN

/* ---------------------------------------------------------------------- */
/* ---------------------------- GRAMMAR BEGINS -------------------------- */
/* ---------------------------------------------------------------------- */

%%

/* Right now, the only type of circuit is a blif circuit.  */

circuit		: blif_circuit 	
			{ 
				_cleanup(); 
				return 0; 
			};
/* ---------------------------------------------------------------------- */
/* ------------------------------ BLIF GRAMMAR -------------------------- */
/* ---------------------------------------------------------------------- */

/*  
 *  # This is an example blif model, with latches.
 *  .model test
 *  # Preamble
 *  .inputs  i0 clk
 *  .inputs  i1 i2 
 *  .outputs out
 *  # Latches
 *  .latch a a1 re clk 2
 *  .latch b b1 re clk 2
 *  # Logic
 *  .names i0 i1 a
 *  11 1
 *  00 1
 *  .names il i2 b
 *  11 1
 *  00 1
 *  .names out a1 b1
 *  10 1
 *  .end
 *
 *  Subcircuits and gates are not implemented.
 *
 */

/*
	.model <decl-model-name>
	.inputs <decl-input-list>
	.outputs <decl-input-list>
	.clock <decl-clock-list>
	<command>
	.
	.
	.`
	<command>
	.end
*/
blif_circuit	: model_stmt	
				;

model_stmt 	: MODEL_TOKEN model_name newline model_spec BLIFEND_TOKEN
			;


model_name	: STRING_TOKEN 			{ _init(g_mytext); };

/*  We require that all .inputs statementes come in the preamble; not
 *  sure if BLIF grammar does too.  There can be any number of them, 
 *  and \<newline> at the end of the line implies a continuation (ignored
 *  by lexical analysis).  Also require that all .inputs are done before
 *  get to the .outputs, which follow the same rules.
 */

model_spec	: input_stmt output_stmt stmts 
			;


input_stmt	: INPUTS_TOKEN input_symbols newline 
			; 

input_symbols	: input_name input_symbols 
				| eps 
				;

input_name		: STRING_TOKEN 	{ g_graph_constructor->new_external_port(g_mytext, PORT::PI); };


output_stmt		: OUTPUTS_TOKEN output_symbols newline 
				;

output_symbols	: output_name output_symbols 
				| eps 
				;

output_name		: STRING_TOKEN 		{ g_graph_constructor->new_external_port(g_mytext, PORT::PO); };

/*  Remaining statements are either latches (.latch) or logic (.names) 
 *  or something unsupported or ignored.  Latches and locic can be 
 *  interleaved.  Note that we are counting on the newlines, in places,
 *  to resolve grammatical ambiguity --- blif is not easily parsable 
 *  if lexical analysis strips the newlines as whitespace; in fact, if
 *  we handle inomplete specifications I think it is ambiguous without them.
 */

stmts	: stmt stmts
		| eps
		;

stmt	: clock_spec 
		| names_spec 
		| latch_spec 
		| not_supported 
		;

// different from hutton
clock_spec	: CLOCK_TOKEN name_symbols { Warning("Ignoring .clocks stmt.  Assigning global clock."); };

not_supported	: NOTSUPPORTED_TOKEN { 
						Warning("Unsupported construct, will try to ignore it."); };

/*  To parse a logic element, we stack the arguments (because we don't know
 *  how big it is until we get the last one).  Then parse the table entries,
 *  putting them in a cubelist_t.  Once the table is in, _new_logic() will
 *  add the node to the graph, and instantiate it with its ROM.
 *  the ROM portion is for the XNF format
 */

names_spec	: NAMES_TOKEN name_symbols newline truth_table {
		    g_graph_constructor->new_combination_block(g_variable_name_stack, 
														g_current_lut);
			/* store it away */
		};

name_symbols	: name_symbol name_symbols  
				| eps
				;

name_symbol	: STRING_TOKEN { 
			debugif(DBLIF, "Got symbol '" << g_mytext << "' in logic");
		    // Just stack the symbol for now, resolve after we
			// have all of them, so we know who is the output.
		    g_variable_name_stack->push_back(g_mytext); 
		}
		;

truth_table	: table_line truth_table
		| eps
		;

table_line	: entry value newline {
		    // add the new cube and value to the current list.
			// the variable stack size is the number of inputs variables
		    g_graph_constructor->new_truth_table_entry(g_entry, g_value, 
									g_variable_name_stack->size() - 1,
									g_current_lut);
		} ;

entry		: STRING_TOKEN { 
		    // Store the cube until we have the cover

		    g_entry = g_mytext;
		    debugif(DBLIF, "got table entry '" << g_entry << "'");
		};

value		: got_value 
			| no_value 
			;

got_value	: STRING_TOKEN { 
			g_value = g_graph_constructor->new_value(g_mytext, g_current_lut);
			debugif(DBLIF, "got table value '" << g_value << "'");
		};

no_value	: eps {
		    debugif(DBLIF, "Function is a constant '" <<  g_entry << "'");
		    g_value = g_graph_constructor->new_value(g_entry, g_current_lut);
		    g_entry = string("");
		};

/* Latches are easy, but we require the optional rising-edge and timing   
 * stuff.  It is always there anyway, so I didn't take the time to make  
 * a grammar which could ignore them.  To make them optional again, need 
 * a couple different types of endings, using the newline as the end of  
 * the recursion.                                                        
 */

latch_spec	: LATCH_TOKEN latch_stuff 
		;

latch_stuff	: latch_in latch_out latch_typ latch_clk latch_n newline { 
		    g_graph_constructor->new_flip_flop(g_latchin, g_latchout, g_latchclk); 
		};

latch_in	: STRING_TOKEN 		{ g_latchin = g_mytext; };

latch_out	: STRING_TOKEN 		{ g_latchout = g_mytext; };

latch_typ 	: STRING_TOKEN 	// ignored:  don't care about type 
			;

latch_clk	: STRING_TOKEN 		{ 	g_latchclk = g_mytext; };

latch_n    	: STRING_TOKEN 	// ignored:  don't need this either
			;

/*  Newlines are *required* for the blif grammar.  I need them to notice 
 *  things like the end of a truth-table, because I have made the 
 *  end-value optional.  The epsilon production is just convenient.
 *  (Usually blif lines that extend past the end of line use \ to indicate 
 *	a continuing line)
 */


newline		: NEWLINE_TOKEN 
			| NEWLINE_TOKEN newline
			; 


eps		: /* */  
		;

%%	

/* ---------------------------------------------------------------------- */
/* -------------------------- SEMANTIC ACTIONS -------------------------- */
/* ---------------------------------------------------------------------- */


#include "circ.h"
#include "circuit.h"
#include "options.h"
#include "graph_constructor.h"

extern CIRCUIT * g_parsed_graph;
/*  
 *  Error routine for use by bison.  Note that it has to be called yyerror,  
 *  even though it is static global.
 */
int
yyerror(char * error_text) 
{
	Error("Parse error, line " << g_linenum << " of input: " << error_text);

    if (strncmp(yytext,"\n",1)) {
        Error("Offending token would be '" << yytext << "' or before.");
    }
    Fail("Cannot continue.  Goodbye.");

	return 0;
}


/*  
 *  Initialization -- setup global variables, and allocate the graph.
 */
static void
_init
(
	string & model_name
)
{
    g_graph_constructor = new GRAPH_CONSTRUCTOR(g_options);
	g_variable_name_stack = new VARIABLE_STACK_TYPE;
	g_got_graph = true;
	g_current_lut = 0;
}


/*  
 *  Free some stuff and set the global (shared with, but owned by the 
 *  caller) variable parsed_graph to our new digraph.
 */
static void
_cleanup()
{
	
    debugif(DCODE, "Parse completed.  Calling cleanup");
    g_variable_name_stack->clear();

	g_graph_constructor->delete_unusable_nodes();
    g_parsed_graph = g_graph_constructor->get_constructed_graph();

    if (! g_got_graph) 
	{
		g_parsed_graph = 0;
    }

}
