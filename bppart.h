/*
Copyright 2018 Hamid Reza Chitsaz (chitsaz@chitsazlab.org)

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
	Desc: BPPart class drives the whole program.

	Author: Hamid Reza Chitsaz and Ali Ebrahimpour Boroojeny
		Colorado State University
		Algorithmic Biology Lab

*/

#ifndef BPPART_H
#define BPPART_H

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>

#include "getopt.h"
#include "alloc.h"
#include "sequence.h"
#include "config.h"
#include "table.h"
#include "bpscore.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

long int total_tables_size;
int big_tables_num;


Option OPTIONS[] = {
	Option('V', (char *)"version", NO_ARG, (char *)"prints the version"),
	Option('h', (char *)"help", NO_ARG, (char *)"shows this help"),
	Option('s', (char *)"suffix", NEEDS_ARG, (char *)"= suffix (default NULL)"),
	Option('d', (char *)"debug", NO_ARG, (char *)"debug"),
	Option('p', (char *)"proc", NEEDS_ARG, (char *)"= number of CPUs"),
	Option('q', (char *)"quiet", NO_ARG, (char *)"quiet mode"),
	Option('A', (char *)"var2", NEEDS_ARG, (char *)"= var2"),
	Option('G', (char *)"var3", NEEDS_ARG, (char *)"= var3"),
	Option(0, NULL, 0, NULL)
};

void version(const char* prog)
{
	printf("%s (%s) \n", prog, PACKAGE_STRING);
	puts("By Hamidreza Chitsaz");
	puts("Copyright (C) 2018");
	puts("Colorado State University");
	puts("Fort Collins, Colorado");
	exit(EXIT_SUCCESS);
}


class BPPart {
private:
	Alloc alloc;
	GetOpt *opts;
	char *suffix;
	bool debug;
	bool quiet;
	int procNum;
	int files;
	vector<char *> filenames;
	FILE *logfile, *outfile;

	double var2;
	double var3;
	vector<Sequence *> seqs;
	Sequence *seq[2];
	int seq_num;
	unsigned char *sq1;
	unsigned char *sq2;
	int len1, len2, sublen1, sublen2;
	int current_seq;


	Table<double> *Q[2],*Qz[2], *QI, *QIa, *QIac, *QIs[2], *QIe, *QIaux[2], *QIm;
	char d2 = 2;
	BPScore scorer;

	FILE *openfile(char *fn, char *, char *, char *, char *);
	FILE *openfile(char *fn1, char *fn2, char *, char *, char *, char *);
	double score(int a, int b, double var2, double var3);
	double iscore(int a, int b, double var2, double var3);


public:
	BPPart(int argc, char** argv);
	~BPPart();
	void allocate();
	void compute();
	void release();
	bool more_pairs();
	void advance();
};

#endif
