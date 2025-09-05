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



#include "util.h"

#include <strstream>
#include <string>
using namespace std;

#ifndef VISUAL_C
#include <sys/resource.h>
#include <sys/time.h>
#else
#include <time.h>
#endif

#include <assert.h>
 


string util_strip_file_extension
(
	const string & file_name
)
{
	string 	string_without_ending = ""; 
	char	ending_separator = '.';

	if (file_name.find(ending_separator) == file_name.npos)
	{
		// we don't have an extention in the string
		string_without_ending = file_name;
	}
	else
	{
		// we do have an extension in the string
		long position_after = file_name.find_last_of(ending_separator);
		assert(position_after >= 0);

		string_without_ending = file_name.substr(0, position_after);
	}

	return string_without_ending;
}

string util_strip_directory_name
(
	const string & file_name
)
{
	string 	string_without_directory = ""; 
	int 	last_directory_separator_position;
	char	directory_sepator = '/';


#ifdef VISUAL_C
	directory_sepator = '\\';
#endif
  
	if (file_name.find(directory_sepator) == file_name.npos)
	{
		// we don't have a directory separtor
		string_without_directory = file_name;
	}
	else
	{
		// we do have a directory separator
		last_directory_separator_position = file_name.find_last_of(directory_sepator);
	
		string_without_directory = 
			file_name.substr(last_directory_separator_position+1);
	}

	return string_without_directory;
}

//  Return the time since Jan 1, 1970, in clock_ticks.

long util_ticks()
{

#ifdef VISUAL_C
	time_t ticks;
    ticks = time(&ticks);
    return ticks;
#else
	struct timeval starttime;	
    struct timezone tz_dummy;
    clock_t ticks;
    gettimeofday(&starttime, &tz_dummy);
    ticks = starttime.tv_sec;
    return ticks;

#endif

}

/*
 *  Return the date and time, as a string.
 */
string util_time_string()
{

	//VISUAL_C might not have this.
    //Warn("Visual C can't do time strings, returning <notime>\n");
    //return "<notime>";
    //Warn("Unknown arch can't do time strings, returning <notime>");
    //return "<notime>";


    clock_t ticks;
    ticks = (clock_t) util_ticks();
	string time_string = ctime(&ticks);
    return time_string;
}


/*
 *  Return the cpu usage, in milliseconds.
 */
int util_cputime(void)
{
    int cputime = 0;
    
#ifndef VISUAL_C
    struct rusage rusage;
    (void) getrusage(RUSAGE_SELF, &rusage);
    cputime = (int) (rusage.ru_utime.tv_sec*1000 + 
	rusage.ru_utime.tv_usec/1000);

#else

	cout << "Warning: CPU time not implemented for Visual C\n";

#endif


    return cputime;
}


string  util_long_to_string(const long & number)
{
	strstream number_string;

	number_string << number << ends;

	return (number_string.str());
}
