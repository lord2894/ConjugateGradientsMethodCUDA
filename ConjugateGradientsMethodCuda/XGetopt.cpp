// XGetopt.cpp  Version 1.2
//
// Author:  Hans Dietrich
//          hdietrich2@hotmail.com
//
// Description:
//     XGetopt.cpp implements getopt(), a function to parse command lines.
//
// History
//     Version 1.2 - 2003 May 17
//     - Added Unicode stopport
//
//     Version 1.1 - 2002 March 10
//     - Added example to XGetopt.cpp module header 
//
// This software is released into the public domain.
// You are free to use it in any way you like.
//
// This software is provided "as is" with no expressed
// or implied warranty.  I accept no liability for any
// damage or loss of business that this software may cause.
//
///////////////////////////////////////////////////////////////////////////////
#include "XGetopt.h"

namespace ParsingOptions {
	char	*optarg;		// global argument pointer
	int		optind = 0; 	// global argv index

	int getopt(int argc, char *argv[], char *optstring)
	{
		static char *next = NULL;
		if (optind == 0)
			next = NULL;

		optarg = NULL;

		if (next == NULL || *next == '\0')
		{
			if (optind == 0)
				optind++;

			if (optind >= argc || argv[optind][0] != '-' || argv[optind][1] == '\0')
			{
				optarg = NULL;
				if (optind < argc)
					optarg = argv[optind];
				return EOF;
			}

			if (strcmp(argv[optind], "--") == 0)
			{
				optind++;
				optarg = NULL;
				if (optind < argc)
					optarg = argv[optind];
				return EOF;
			}

			next = argv[optind];
			next++;		// skip past -
			optind++;
		}

		char c = *next++;
		char *cp = strchr(optstring, c);

		if (cp == NULL || c == ':')
			return '?';

		cp++;
		if (*cp == ':')
		{
			if (*next != '\0')
			{
				optarg = next;
				next = NULL;
			}
			else if (optind < argc)
			{
				optarg = argv[optind];
				optind++;
			}
			else
			{
				return '?';
			}
		}

		return c;
	}
}