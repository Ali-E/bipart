/*
Copyright 2012 Hamid Reza Chitsaz, Elmirasadat Forouzmand (chitsaz@wayne.edu)

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
	Desc: UBPartitionFunction class computes an upper bound for the partition function. It uses the Energy class to compute energies. This version accomodates a nonsaturated window sublength.

	Author: Hamid Reza Chitsaz, Elmirasadat Forouzmand
		Wayne State University 
		Algorithmic Biology Lab

	Last Update: Sep 14, 2012
	Ver: 1.0
*/

#ifndef UBPARFUNC_H
#define UBPARFUNC_H

#include "energy.h"
#include "sequence.h"
#include "table.h"
#include "alloc.h"
#include "partitionfunction.h"

class UBPartitionFunction {
private:	
	Alloc alloc;
	int nodangle;
	int zip;
	int NA;
	double naConc;
	double mgConc;
	int polymer;
	double tMin;
	double tInc;
	double tMax;
	char *suffix;
	bool no_isolate;
	int energyType;
	int procNum;
	bool DEBUG, quiet;
	ParFuncType parFuncType;


	Energy *energy;
	Sequence *sequence;
	Table<double> *perturbation[2]; // 0=positive base pair 1=negative base pair
	unsigned int sublen;

	UBPartitionFunction *parfunc[2];

	unsigned char *seqUp, *seqDown;

	Table<double> *Q[3], *Qm[3], *Qb; // color: 0=white 1=red 2=green


	Table<double> **auxQIa_dd[2][2], **auxQIa[2][2], **QIx[2][2], *QI, *QIb[2][2], *QIa, *QImm[3][3], *QImk[3][2], *QIkm[2][3], *QIkk[2][2], *QIa_nn[3][3], *QIa_nd[3][2],
		*QIa_dn[2][3], *QIa_dd[2][2], *QI_nn[3][3], *QI_nd[3][2], *QI_dn[2][3], *QI_dd[2][2], *QIa_r[2][2], *QIlr[2][2], *QI_r[2][2], 
		*QIhh[2][2], *QIhm[2][2],
		//first index for up, second index for down, white = 0, red = 1, green = 2
		*QIh[2], *QIy, *QIx_m[2][3], *QIx_k[2][2]; 
		//first index upOrDown (up = 0), second index color of the other side

	double RT;
	double Zz;

	double product(Table<double> *, Table<double> *, int, int, int, int);
	double product(Table<double> *, Table<double> **, int, int, int, int);
	double product(Table<double> **, Table<double> *, int, int, int, int);
	double product(Table<double> **, Table<double> **, int, int, int, int);
	double lsemiproduct(Table<double> *, Table<double> *, Table<double> *, int, int, int, int);
	double lsemiproduct(Table<double> **, Table<double> *, Table<double> *, int, int, int, int);
	double rsemiproduct(Table<double> *, Table<double> *, Table<double> *, int, int, int, int);
	double capproduct(Table<double> *, Table<double> *, Table<double> *, int, int, int, int, int);

	/*double sum(Table<double> *, Table<double> *, int, int, int, int);	
	double sum(Table<double> *, Table<double> **, int, int, int, int);
	double sum(Table<double> **, Table<double> *, int, int, int, int);
	double sum(Table<double> **, Table<double> **, int, int, int, int);
	double lsemisum(Table<double> *, Table<double> *, Table<double> *, int, int, int, int);
	double lsemisum(Table<double> **, Table<double> *, Table<double> *, int, int, int, int);
	double rsemisum(Table<double> *, Table<double> *, Table<double> *, int, int, int, int);
	double capsum(Table<double> *, Table<double> *, Table<double> *, int, int, int, int, int);*/

	double computeSingleParFunc();

	void computeSingleQbQm(int i, int j);
	void computeSingleQ(int i, int j);

	double QIx_(int i, int j, int i1, int j1, int i2, int j2) {return QIx_m[i][j]->element(i1, j1, i2, j2) + QIx_k[i][j]->element(i1, j1, i2, j2);}
	double computeDoubleParFunc();

	void computeBands(int i1, int j1, int i2, int j2);
	void computeXY(int i1, int j1, int i2, int j2);
	void computeAux(int i1, int j1, int i2, int j2);
	void computeImmkk(int i1, int j1, int i2, int j2);
	void computeIa_nndd(int i1, int j1, int i2, int j2);
	void computeI_nndd(int i1, int j1, int i2, int j2);
	void computeI(int i1, int j1, int i2, int j2);

	void selftest(int, int, int, int, int, int);

public:
	UBPartitionFunction(Energy *energy, Sequence *seq, Table<double> *pertrb[2], bool _debug = false, int _procNum = 16, bool _quiet = false, bool _no_isolate = true);	
	UBPartitionFunction(FILE *, Energy *energy, Sequence *seq, Table<double> *pertrb[2], bool _debug = false, int _procNum = 16, bool _quiet = false, bool _no_isolate = true);	
	UBPartitionFunction(Energy *energy, Sequence *seq, Table<double> *pertrb[2], unsigned int _sublen, bool _debug = false, int _procNum = 16, bool _quiet = false, bool _no_isolate = true);	
	UBPartitionFunction(Energy *energy, Table<double> *pertrb[2], UBPartitionFunction *p[2], bool _debug = false, int _procNum = 16, bool _quiet = false, bool _no_isolate = true);	
	UBPartitionFunction(FILE *, Energy *energy, Table<double> *pertrb[2], UBPartitionFunction *p[2], bool _debug = false, int _procNum = 16, bool _quiet = false, bool _no_isolate = true);	
	UBPartitionFunction(UBPartitionFunction *p);	
	~UBPartitionFunction();	

	double freeEnsembleEnergy(double = 1.0);
	double parFunc(double);

	UBPartitionFunction **getParFunc() {return parfunc;}
	ParFuncType getType() {return parFuncType;}
	Sequence *getSequence() {return sequence;}
	Table<double> *getPerturbation(int i) {return perturbation[i];}
	Energy *getEnergy() {return energy;}
	void getParams(bool *_debug, int *_procNum, bool *_quiet, bool *_no_isolate) {*_debug = DEBUG; *_procNum = procNum; *_quiet = quiet; *_no_isolate = no_isolate;}
	unsigned int getSubLen() {return sublen;}
	Table<double> *getQ(int i) {return Q[i];}
	Table<double> *getQm(int i) {return Qm[i];}
	Table<double> *getQb() {return Qb;}

	Table<double> *getQI() {return QI;}
	Table<double> *getQIa() {return QIa;}
	Table<double> *getQIb(int i, int j) {return QIb[i][j];}
	Table<double> *getQImm(int i, int j) {return QImm[i][j];}
	Table<double> *getQImk(int i, int j) {return QImk[i][j];}
	Table<double> *getQIkm(int i, int j) {return QIkm[i][j];}
	Table<double> *getQIkk(int i, int j) {return QIkk[i][j];}
	Table<double> *getQIa_nn(int i, int j) {return QIa_nn[i][j];}
	Table<double> *getQIa_nd(int i, int j) {return QIa_nd[i][j];}
	Table<double> *getQIa_dn(int i, int j) {return QIa_dn[i][j];}
	Table<double> *getQIa_dd(int i, int j) {return QIa_dd[i][j];}
	Table<double> *getQI_nn(int i, int j) {return QI_nn[i][j];}
	Table<double> *getQI_nd(int i, int j) {return QI_nd[i][j];}
	Table<double> *getQI_dn(int i, int j) {return QI_dn[i][j];}
	Table<double> *getQI_dd(int i, int j) {return QI_dd[i][j];}
	Table<double> *getQIa_r(int i, int j) {return QIa_r[i][j];}
	Table<double> *getQI_r(int i, int j) {return QI_r[i][j];}
	Table<double> *getQIlr(int i, int j) {return QIlr[i][j];}
	Table<double> *getQIhh(int i, int j) {return QIhh[i][j];}
	Table<double> *getQIhm(int i, int j) {return QIhm[i][j];}
	Table<double> *getQIx_m(int i, int j) {return QIx_m[i][j];}
	Table<double> *getQIx_k(int i, int j) {return QIx_k[i][j];}
	Table<double> *getQIh(int i) {return QIh[i];}
	Table<double> *getQIy() {return QIy;}
	Table<double> **getauxQIa(int i, int j) {return auxQIa[i][j];}
	Table<double> **getauxQIa_dd(int i, int j) {return auxQIa_dd[i][j];}
	Table<double> **getQIx(int i, int j) {return QIx[i][j];}

	void Reverse();
	void store(FILE *);
	void retrieve(FILE *);
};

#endif
