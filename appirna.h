/*
Copyright 2012 Hamid Reza Chitsaz (chitsaz@wayne.edu)

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
	Desc: Partition class drives the whole program.

	Author: Hamid Reza Chitsaz, Elmirasadat Forouzmand
		Wayne State University
		Algorithmic Biology Lab

	Last Update by Hamid Reza Chitsaz: Sep 19, 2012
	Ver: 1.0
*/

#ifndef APPART_H
#define APPART_H

#include <math.h>
#include <time.h>
#include <vector>

#include "ubpartitionfunction.h"
#include "getopt.h"
#include "config.h"

using namespace std;

Option OPTIONS[] = {
	Option('V', (char *)"version", NO_ARG, (char *)"prints the version"),
	Option('h', (char *)"help", NO_ARG, (char *)"shows this help"),
	Option('n', (char *)"NA", NEEDS_ARG, (char *)"=RNA (default) | DNA"),
	Option('e', (char *)"energy", NEEDS_ARG, (char *)"=|sig| (default) | sig"),
	Option('t', (char *)"tmin", NEEDS_ARG, (char *)"=Min Temperature (default 37)"),
	Option('i', (char *)"tinc", NEEDS_ARG, (char *)"=Temperature Inc (default 1)"),
	Option('T', (char *)"tmax", NEEDS_ARG, (char *)"=Max Temperature (default 37)"),
	Option('s', (char *)"suffix", NEEDS_ARG, (char *)"= suffix (default NULL)"),
	Option('d', (char *)"debug", NO_ARG, (char *)"debug"),
	Option('N', (char *)"sodium", NEEDS_ARG, (char *)"=Salt concentration"),
	Option('M', (char *)"magnesium", NEEDS_ARG, (char *)"=Mg concentration"),
	Option('p', (char *)"proc", NEEDS_ARG, (char *)"=number of CPUs"),
	Option('P', (char *)"polymer", NO_ARG, (char *)"polymer"),
	Option('q', (char *)"quiet", NO_ARG, (char *)"quiet mode"),
	Option('z', (char *)"zip", NO_ARG, (char *)"zip"),
	Option(2, (char *)"nodangle", NO_ARG, (char *)"no dangle energy (default)"),
	Option(10, (char *)"withisolates", NO_ARG, (char *)"allow isolate interacting base pairs"),
	Option(17, (char *)"store", NO_ARG, (char *)"store all tables in file"),
	Option(18, (char *)"retrieve", NO_ARG, (char *)"retrieve all tables from file"),
	Option(0, NULL, 0, NULL)
};

void version(const char* prog)
{
	printf("%s (%s) \n", prog, PACKAGE_STRING);
	puts("By Hamidreza Chitsaz, Elmirasadat Forouzmand");
	puts("Copyright (C) 2012");
	puts("Wayne State University");
	puts("Detroit, Michigan");
	exit(EXIT_SUCCESS);
}

class APPartition {
private:
	Alloc alloc;
	GetOpt *opts;
	int NA;
	double naConc;
	double mgConc;
	int polymer;
	int nodangle;
	double tMin;
	double tInc;
	double tMax;
	char *suffix;
	int zip;
	bool debug;
	bool quiet;
	int procNum;
	int files;
	char *filenames[20];
	bool no_isolate;
	int energyType;

	char dataDir[1000];
	char *prefix;
	bool store;
	bool retrieve;
	
	time_t end_time;

	vector<Sequence *> seq;
	Energy *energy;
	FILE *openfile(char *fn, int i, char *, char *);
	FILE *openfile(char *fn1, char *fn2, int i, int j, char *, char *);
	void sample(Table<double> *samp[2], unsigned int len);

public:
	APPartition(int argc, char** argv);
	~APPartition();
	void computeUpperbound();
	bool doInteraction(){return files > 1;};
};

#endif
