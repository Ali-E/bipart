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
	Desc: Partition class drives the whole program.

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow, 
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: Oct 25, 2018
	Ver: 3.4
*/

#ifndef PART_H
#define PART_H

#include <vector>
#include <math.h>
#include <time.h>

#include "partitionfunction.h"
#include "getopt.h"
#include "config.h"

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
	Option(3, (char *)"NotAA", NO_ARG, (char *)"do not compute AA partition function"),
	Option(4, (char *)"NotBB", NO_ARG, (char *)"do not compute BB partition function"),
	Option(5, (char *)"NotAB", NO_ARG, (char *)"do not compute AB partition function"),
	Option(10, (char *)"withisolates", NO_ARG, (char *)"allow isolate interacting base pairs"),
	Option(12, (char *)"A0", NEEDS_ARG, (char *)"=initial A concentration (if used, partition function will not be computed)"),
	Option(13, (char *)"B0", NEEDS_ARG, (char *)"=initial B concentration (if used, partition function will not be computed)"),
	Option(14, (char *)"ensemble", NO_ARG, (char *)"compute the ensemble energy (if used, partition function will not be computed)"),
	Option(15, (char *)"tm", NO_ARG, (char *)"compute AB melting temperature (if used, partition function will not be computed)"),
	Option(16, (char *)"prob", NO_ARG, (char *)"compute base pair probabilities"),
	Option(17, (char *)"store", NEEDS_ARG, (char *)"= P | Q | PQ, store all the probability and/or partition function tables in file"),
	Option(18, (char *)"retrieve", NEEDS_ARG, (char *)"= P | Q | PQ, retrieve all the probability and/or partition function tables from file"),
	Option(0, NULL, 0, NULL)
};

void version(const char* prog)
{
	printf("%s (%s) \n", prog, PACKAGE_STRING);
	puts("By Hamidreza Chitsaz");
	puts("Copyright (C) 2008-2018");
	puts("Colorado State University");
	puts("Fort Collins, Colorado");
	exit(EXIT_SUCCESS);
}

using namespace std;

class Partition {
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
	vector<char *> filenames;
	bool AA;
	bool BB;
	bool AB;
	bool no_isolate;
	int energyType;

	char dataDir[1000];
	char *prefix;
	double Na0, Nb0;
	bool ensemble;
	bool tm;
	bool prob;
	bool storeP, storeQ;
	bool retrieveP, retrieveQ;
	
	time_t end_time;

	vector<Sequence *> seq;
	int seq_num;
	Energy *energy;
	FILE *openfile(char *fn, int i, char *, char *);
	FILE *openfile(char *fn1, char *fn2, int i, int j, char *, char *);

public:
	Partition(int argc, char** argv);
	~Partition();
	void computeThermoValues();
	void computeConcentrations();
	void computedG();
	void computeTm();
	bool doPartition() {return seq_num != 0;}
	bool doConc() {return seq_num == 0 && (Na0 != 0.0 || Nb0 != 0.0);}
	bool doEns() {return ensemble;}
	bool doTm() {return tm;}
};

#endif
