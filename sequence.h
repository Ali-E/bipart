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

#ifndef SEQ_H
#define SEQ_H

#include <stdlib.h>
#include <stdio.h>
#include "alloc.h"

class Sequence {
private:
	Alloc alloc; 
	char *file_name;
	char *name;
	char *string;
	char short_filename[100];
	unsigned char *seq; /* [0-4] for [A,C,G,TU,N] */
	unsigned char *reverse_seq; /* [0-4] for [A,C,G,TU,N] */
	unsigned int len;
	bool loaded;
	unsigned char toNum(char c);
	int readSequence(FILE *);
	void checkArray(char** array, unsigned int* available, unsigned int used, unsigned int increment);
	void computeReverse();

public:
	Sequence(char *fn);
	Sequence(FILE*, char *fn);
	Sequence(char *, char *n, char *st, unsigned char *sq, unsigned int l);
	~Sequence();
	int load(FILE *f);
	int input(FILE* file);
	bool isLoaded() {return loaded;};
	char *getName();
	char *getFileName();
	char *getString();
	unsigned char *getSeq();
	unsigned char *getReverse();
	unsigned int getLen();
	Sequence *split();
};

#endif
