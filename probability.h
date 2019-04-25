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
	Desc: Probability class computes the base pair probabilities and also motif probabilities. It uses the 
PartitionFunction class to compute Q tables. 

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow, 
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: June 4, 2009
	Ver: 1.6
*/

#ifndef PROB_H
#define PROB_H

#include "energy.h"
#include "sequence.h"
#include "partitionfunction.h"
#include "table.h"
#include "alloc.h"

class Probability {
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

	bool finished;

	Energy *energy;

	PartitionFunction *parfunc;

	unsigned char *seqUp, *seqDown;
	Table<double> *bProb, *aProb[2];


	Table<double> *P[2][3], *Pm[2][3], *Pb[2]; // first index: the sequence, second index: color: 0=white 1=red 2=green

	Table<double> **auxPIa_dd[2][2], **auxPIa[2][2], **PIx[2][2], *PI, *PIb[2][2], *PIa, *PImm[3][3], *PImk[3][2], *PIkm[2][3], *PIkk[2][2], *PIa_nn[3][3], *PIa_nd[3][2],
		*PIa_dn[2][3], *PIa_dd[2][2], *PI_nn[3][3], *PI_nd[3][2], *PI_dn[2][3], *PI_dd[2][2], *PIa_r[2][2], *PIlr[2][2], *PI_r[2][2], 
		*PIhh[2][2], *PIhm[2][2],
		//first index for up, second index for down, white = 0, red = 1, green = 2
		*PIh[2], *PIy, *PIx_m[2][3], *PIx_k[2][2]; 
		//first index upOrDown (up = 0), second index color of the other side

	Table<double> *Q[2][3], *Qm[2][3], *Qb[2]; // first index: the sequence, second index: color: 0=white 1=red 2=green

	Table<double> **auxQIa_dd[2][2], **auxQIa[2][2], **QIx[2][2], *QI, *QIb[2][2], *QIa, *QImm[3][3], *QImk[3][2], *QIkm[2][3], *QIkk[2][2], *QIa_nn[3][3], *QIa_nd[3][2],
		*QIa_dn[2][3], *QIa_dd[2][2], *QI_nn[3][3], *QI_nd[3][2], *QI_dn[2][3], *QI_dd[2][2], *QIa_r[2][2], *QIlr[2][2], *QI_r[2][2], 
		*QIhh[2][2], *QIhm[2][2],
		//first index for up, second index for down, white = 0, red = 1, green = 2
		*QIh[2], *QIy, *QIx_m[2][3], *QIx_k[2][2]; 
		//first index upOrDown (up = 0), second index color of the other side

	double RT;

	void product(Table<double> *, Table<double> *, Table<double> *, Table<double> *, double, double, int, int, int, int);
	void product(Table<double> *, Table<double> **, Table<double> *, Table<double> **, double, double, int, int, int, int);
	void product(Table<double> **, Table<double> *, Table<double> **, Table<double> *, double, double, int, int, int, int);
	void product(Table<double> **, Table<double> **, Table<double> **, Table<double> **, double, double, int, int, int, int);
	void lsemiproduct(Table<double> *, Table<double> *, Table<double> *, Table<double> *, Table<double> *, Table<double> *, double, double, int, int, int, int);
	void lsemiproduct(Table<double> **, Table<double> *, Table<double> *, Table<double> **, Table<double> *, Table<double> *, double, double, int, int, int, int);
	void rsemiproduct(Table<double> *, Table<double> *, Table<double> *, Table<double> *, Table<double> *, Table<double> *, double, double, int, int, int, int);
	void capproduct(Table<double> *, Table<double> *, Table<double> *, Table<double> *, Table<double> *, Table<double> *, int, int, int, int, int);

	void processP(int i, int j, int);
	void processPbPm(int i, int j, int);

	double PIx_(int i, int j, int i1, int j1, int i2, int j2) {return PIx_m[i][j]->element(i1, j1, i2, j2) + PIx_k[i][j]->element(i1, j1, i2, j2);}

	unsigned int computeProb(unsigned int, time_t);

	void processBands(int i1, int j1, int i2, int j2);
	void processXY(int i1, int j1, int i2, int j2);
	void processAux(int i1, int j1, int i2, int j2);
	void processImmkk(int i1, int j1, int i2, int j2);
	void processIa_nndd(int i1, int j1, int i2, int j2);
	void processI_nndd(int i1, int j1, int i2, int j2);
	void processI(int i1, int j1, int i2, int j2);
	void processIa(int i1, int j1, int i2, int j2);
	void processIb(int i1, int j1, int i2, int j2);

	void Reverse(int);
	
	void buildTables();
	void buildEmptyTables();

	void selftest(int, int, int, int, int, int);
public:
	Probability(Energy *energy, PartitionFunction *pfunc, bool _debug = false, int _procNum = 16, bool _quiet = false, bool _no_isolate = true);	
	Probability(FILE *, FILE *, Energy *energy, PartitionFunction *pfunc, time_t, bool _debug = false, int _procNum = 16, bool _quiet = false, bool _no_isolate = true);	
	~Probability();
	void print(FILE *);
	void store(FILE *);
	void retrieve(FILE *);
	bool isFinished() {return finished;};
};

#endif
