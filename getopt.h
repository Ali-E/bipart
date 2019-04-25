/*
Copyright 2008, 2009 Hamid Reza Chitsaz (chitsaz@cs.ucsd.edu)

    This file is part of pi-biRNA.

    pi-biRNA is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    pi-biRNA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with pi-biRNA; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

/*
	Desc: Command-line options class.

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow, 
		SFU, Computational Biology Lab

	Last Update: Nov 10, 2008
	Ver: 1.0
*/

#include <stdlib.h>
#include <string.h>


#define NEEDS_ARG 2
#define OPTIONAL_ARG 1
#define NO_ARG 0

#define FREE_ARG 27

class Option {
private:
	char shortForm;
	char *longForm;
	int argRequirement;
	char arg[1000];
	char description[1000];
public:
	Option(char s, char *l, int r, char *d)
	{
		shortForm = s;
		longForm = l;
		argRequirement = r;
		description[0] = 0;
		if(d) strcpy(description, d);
		arg[0] = 0;
	}
	char getShortForm(){ return shortForm;};
	char *getLongForm(){ return longForm;};
	int getArgRequirement(){ return argRequirement;};
	char *getArg() { return arg;};
	char *getDesc() { return description;};
	void setArg(char *a) {strcpy(arg, a);};
};

class GetOpt {
private:
	int argc;
	char **argv;
	void parseArgs();
	
	Option *options;
	int optionsNum;
	char helpstring[10000];


	Option *parsedOptions[1000];
	int parsedOptionsNum;
	int currentOp;
public:
	GetOpt(int a, char **v, Option *ops) {parsedOptionsNum = 0; currentOp = 0; argc = a; argv = v; options = ops; parseArgs();};
	~GetOpt(){for(int i=0; i < parsedOptionsNum; i++) free(parsedOptions[i]);};
	bool hasNext() {return currentOp < parsedOptionsNum;};
	Option *next(){return parsedOptions[currentOp++];}
	char *help();
};
