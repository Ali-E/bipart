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
	Desc: OligoInteract class drives the whole program.

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow, 
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: April 25, 2009
	Ver: 1.2
*/

#ifndef INTERACT_H
#define INTERACT_H

#include <math.h>
#include <time.h>

#include "partitionfunction.h"
#include "jointprob.h"
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
	Option('m', (char *)"max", NEEDS_ARG, (char *)"=maximum candidate binding site pairs (default 1000)"),
	Option('w', (char *)"window", NEEDS_ARG, (char *)"=maximum sublength (default 25)"),
	Option('p', (char *)"proc", NEEDS_ARG, (char *)"=number of CPUs"),
	Option('P', (char *)"polymer", NO_ARG, (char *)"polymer"),
	Option('q', (char *)"quiet", NO_ARG, (char *)"quiet mode"),
	Option('z', (char *)"zip", NO_ARG, (char *)"zip"),
	Option(2, (char *)"nodangle", NO_ARG, (char *)"no dangle energy (default)"),
	Option(3, (char *)"nosites", NO_ARG, (char *)"do not dump candidate binding sites"),
	Option(10, (char *)"withisolates", NO_ARG, (char *)"allow isolate interacting base pairs"),
	Option(0, NULL, 0, NULL)
};

void version(const char* prog)
{
	printf("%s (%s)\n", prog, PACKAGE_STRING);
	puts("By Hamid Reza Chitsaz");
	puts("Copyright (C) 2008-2009");
	puts("Simon Fraser University");
	puts("Burnaby, BC Canada");
	exit(EXIT_SUCCESS);
}

class Interaction {
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
	int max_sublen;
	bool dump_sites;
	int sites_num;

	char dataDir[1000];
	char *prefix;
	
	Sequence *seq[100000];
	int seq_num;
	Energy *energy;
	void fname1(char **res, char *fn, int i);
	void fname2(char **res, char *fn1, char *fn2, int i, int j);

	int maxCandidates;
	
	void outputSites(RegionCollection *, FILE *, double);
	void filterHighProbIntRegions(RegionCollection *, PartitionFunction *[2], PartitionFunction *, JointProb *[2], int, int, int, int);
	void computeBestMatches(FILE *, int, RegionCollection *[2], PartitionFunction *[2], PartitionFunction *, JointProb *[2], int);
	void computeBestMatches(FILE *, RegionCollection *[2], PartitionFunction *[2], PartitionFunction *, JointProb *[2], int);
	double match(int sites_num, double *weights, int *assign);
public:
	Interaction(int argc, char** argv);
	~Interaction();
	void computeBindingSites();
};

#endif
