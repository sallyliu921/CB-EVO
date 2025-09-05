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


%{
#include "parse_blif.h"
#include <string>
using namespace std;
#include <assert.h>
#include <iostream>
static bool		g_in_comment;	    /* to keep track of comment newlines */
extern long		g_linenum;
string			g_mytext;	    	// a forward reference to yytext
%}
STRING_TOKEN_DEF 	[\[\]a-zA-Z0-9,=!@$%^&*_-]+
%option noyywrap
%option never-interactive
%option nounput
%%
[ \t]+			; /* skip whitespace */ 
\\\n			{ g_linenum += 1; } /* Skip these fake newlines, put in by sis or other tools*/
\n				{ 
					g_linenum += 1; 
					if (g_in_comment) {
						/* Don't mention newline after a comment */
						g_in_comment = false;
					} else {
						/* Tell the parser about it */
						return NEWLINE_TOKEN;
					}
				}
".model"		{ return MODEL_TOKEN; }
".end"			{ return BLIFEND_TOKEN; 	}
".inputs"		{ return INPUTS_TOKEN; 	}
".outputs"		{ return OUTPUTS_TOKEN;	}
".clock"		{ return CLOCK_TOKEN;	}
".names"		{ return NAMES_TOKEN;		}
".latch"		{ return LATCH_TOKEN;		}
\.gate.*$		{ 	g_in_comment = true;
					return NOTSUPPORTED_TOKEN;
				}
\.subckt.*$		{ 	g_in_comment = true;
					return NOTSUPPORTED_TOKEN;
				}
\.				{ 	g_mytext = string(yytext);
					cerr << "'.' not allowed. Find and delete any illegal periods around line " 
							<< g_linenum << endl << endl;
					exit(1);
				}
\:				{
					g_mytext = string(yytext);
					cerr << "':' not allowed. Find and delete any illegal : aournd line "
							<< g_linenum << endl << endl;
					exit(1);
				}

#.*$					{ g_in_comment=true; } /* skip comments */
{STRING_TOKEN_DEF}		{ g_mytext = string(yytext); return STRING_TOKEN; }
.						return ILLEGAL_TOKEN;

%%
