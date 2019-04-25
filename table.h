/*
Copyright 2008 Hamid Reza Chitsaz (chitsaz@wayne.edu)

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
	Desc: Table class stores dynamic programming values.

	Author: Hamid Reza Chitsaz
		Wayne State University, Algorithmic Biology Lab

	Last Update by Hamid Reza Chitsaz: Sep 14, 2012
	Ver: 2.6
*/

#ifndef TABLE_H
#define TABLE_H

#include <stdlib.h>
#include <stdio.h>

#include "alloc.h"
using namespace std;

extern long int total_tables_size;
extern int big_tables_num;

template<class T>
class Table {
private:
	Alloc alloc;
	T *Q;
	T IV;
	bool duplex;
	char dim;
	unsigned int len1, len2;
	unsigned int sublen1, sublen2;
	unsigned int len;
	unsigned int sublen;
	long unsigned int tableSize1, tableSize2, tableSize;
	bool reverse;
	long int allocate();
	inline long unsigned int addr_dup(unsigned int i1, unsigned int j1, unsigned int i2, unsigned int j2);
	inline long unsigned int singaddr2(unsigned int i, unsigned int j); 
	inline long unsigned int singaddr4(unsigned int i, unsigned int j, unsigned int d, unsigned int e); 
	void init(T iv);
	void init(T *);
	T *canvas() {return Q;}

public:
	Table<T>(unsigned int l1, unsigned int l2, T iv);
	Table<T>(unsigned int l1, unsigned int l2, unsigned int subl1, unsigned int subl2, T iv);
	Table<T>(unsigned int l, char *, T iv);
	Table<T>(unsigned int l, unsigned int subl, char *, T iv);
	Table<T>(Table *t);
	Table<T>();
	~Table();
	void Reverse();
	inline T element(int i, int j, int d, int e);
	inline T element(int i, int j);
	inline T *estar(int i, int j, int d, int e);
	inline T *estar(int i, int j);
	inline T & operator()(int i, int j, int d, int e);
	inline T & operator()(int i, int j);
	bool hasAny(T);
	bool hasOtherThan(T);
	void reset();
	T max(int *i, int *j);
	T max(int *i, int *j, int *d, int *e);
	void getParams(bool *_duplex, unsigned int *_len1, unsigned int *_len2, unsigned int *_sublen1, unsigned int *_sublen2, 
unsigned int *_len, unsigned int *_sublen, bool *_reverse, char *_dim);
	void store(FILE *);
	void retrieve(FILE *);
};

#include "table.C"
#endif
