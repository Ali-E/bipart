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
	Desc: PartitionFunction class computes the partition function for two single strands and also the double strand. It uses the Energy class to compute energies. This version accomodates a nonsaturated window sublength.

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow, 
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: June 5, 2009
	Ver: 3.8
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

#include "config.h"
#include "partitionfunction.h"
#include "getopt.h"
#include "sequence.h"


int inline max(int a, int b)
{
	int ret = (a > b) ? a : b;
	return ret;
}

PartitionFunction::PartitionFunction(PartitionFunction *p)
{
	double temp;
	energy = p->getEnergy();
	energy->getParams(&NA, &polymer, &naConc, &mgConc, &suffix, &temp, &zip, &nodangle, &energyType);
	RT = energy->getRT();

	p->getParams(&DEBUG, &procNum, &quiet, &no_isolate);

	parFuncType = p->getType();
	sequence = p->getSequence();
	sublen = p->getSubLen();
	parfunc[0] = parfunc[1] = NULL;
	QI = NULL;

	for(int j = 0; j < 3; j++)
	{
		Q[j] = new Table<double>(p->getQ(j));
		Qm[j] = new Table<double>(p->getQm(j));
	}
	Qb = new Table<double>(p->getQb());

	if(sequence->getLen() == sublen)
		Zz = Q[0]->element(0, sequence->getLen() - 1);
	else
		Zz = 1;
}

PartitionFunction::PartitionFunction(Energy *en, Sequence *seq, bool _debug, int _procNum, bool _quiet, bool _no_isolate)
{
	double temp;
	energy = en;
	energy->getParams(&NA, &polymer, &naConc, &mgConc, &suffix, &temp, &zip, &nodangle, &energyType);
	RT = energy->getRT();

	DEBUG = _debug;
	procNum = _procNum;
	quiet = _quiet;
	no_isolate = _no_isolate;

	parFuncType = SINGLE;
	sequence = seq;
	sublen = seq->getLen();
	parfunc[0] = parfunc[1] = NULL;
	QI = NULL;
	char d2 = 2;


	if(!quiet) printf("Allocating memory for single partition function of %s... \n", sequence->getName());
	for(int j = 0; j < 3; j++)
	{
		Q[j] = new Table<double>(sequence->getLen(), sublen, &d2, 0.0);
		Qm[j] = new Table<double>(sequence->getLen(), sublen, &d2, 0.0);
	}
	Qb = new Table<double>(sequence->getLen(), sublen, &d2, 0.0);

	if(!quiet) printf("Allocation and initialization finished. \n");
	computeSingleParFunc();
}

PartitionFunction::PartitionFunction(FILE *inp, Energy *en, Sequence *seq, bool _debug, int _procNum, bool _quiet, bool _no_isolate)
{
	double temp;
	energy = en;
	energy->getParams(&NA, &polymer, &naConc, &mgConc, &suffix, &temp, &zip, &nodangle, &energyType);
	RT = energy->getRT();

	DEBUG = _debug;
	procNum = _procNum;
	quiet = _quiet;
	no_isolate = _no_isolate;

	parFuncType = SINGLE;
	sequence = seq;
	sublen = seq->getLen();
	parfunc[0] = parfunc[1] = NULL;
	QI = NULL;
	char d2 = 2;


	for(int j = 0; j < 3; j++)
	{
		Q[j] = new Table<double>();
		Qm[j] = new Table<double>();
	}
	Qb = new Table<double>();

	retrieve(inp);
}

PartitionFunction::PartitionFunction(Energy *en, Sequence *seq, unsigned int _sublen, bool _debug, int _procNum, bool _quiet, bool _no_isolate)
{
	double temp;
	energy = en;
	energy->getParams(&NA, &polymer, &naConc, &mgConc, &suffix, &temp, &zip, &nodangle, &energyType);
	RT = energy->getRT();

	DEBUG = _debug;
	procNum = _procNum;
	quiet = _quiet;
	no_isolate = _no_isolate;

	parFuncType = SINGLE;
	sequence = seq;
	sublen = (_sublen < seq->getLen()) ? _sublen : seq->getLen();
	parfunc[0] = parfunc[1] = NULL;
	QI = NULL;
	char d2 = 2;


	if(!quiet) printf("Allocating memory for single partition function of %s... \n", sequence->getName());
	for(int j = 0; j < 3; j++)
	{
		Q[j] = new Table<double>(sequence->getLen(), sublen, &d2, 0.0);
		Qm[j] = new Table<double>(sequence->getLen(), sublen, &d2, 0.0);
	}
	Qb = new Table<double>(sequence->getLen(), sublen, &d2, 0.0);

	if(!quiet) printf("Allocation and initialization finished. \n");
	computeSingleParFunc();
}

PartitionFunction::PartitionFunction(Energy *en, PartitionFunction *par[2], bool _debug, int _procNum, bool _quiet, bool _no_isolate)
{
	double temp;
	energy = en;
	energy->getParams(&NA, &polymer, &naConc, &mgConc, &suffix, &temp, &zip, &nodangle, &energyType);
	RT = energy->getRT();

	DEBUG = _debug;
	procNum = _procNum;
	quiet = _quiet;
	no_isolate = _no_isolate;
	parFuncType = DUPLEX;
	sequence = NULL;
	Qb = NULL;
	for(int j = 0; j < 3; j++)
	{
		Q[j] = NULL;
		Qm[j] = NULL;
	}

	for(int i = 0; i < 2; i++)
		parfunc[i] = par[i];

	if(!quiet) printf("Allocating memory for double partition function of %s and %s... \n", parfunc[0]->getSequence()->getName(), parfunc[1]->getSequence()->getName());

	unsigned int l1 = parfunc[0]->getSequence()->getLen();
	unsigned int l2 = parfunc[1]->getSequence()->getLen();
	unsigned int subl1 = parfunc[0]->getSubLen();
	unsigned int subl2 = parfunc[1]->getSubLen();


	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(i+j) QIa_nn[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			if(i+j) QI_nn[i][j] = new Table<double>(l1, l2, subl1, subl2,0.0);
		}
	}

	for(int i = 0; i < 2; i++)
	{
		QIh[i] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QIx_m[i][2] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QImm[i][2] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QImm[2][i] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QImk[i][1] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QImk[2][i] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QIkm[1][i] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QIkm[i][2] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QIa_nd[i][1] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QIa_nd[2][i] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QI_nd[i][1] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QI_nd[2][i] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QIa_dn[1][i] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QIa_dn[i][2] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QI_dn[1][i] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		QI_dn[i][2] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		for(int j = 0; j < 2; j++)
		{
			QIb[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			QIhh[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			QIhm[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			QIx_m[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			QIx_k[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			if(i+j) QIkk[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			if(i+j) QIa_r[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			if(i+j) QI_r[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			if(i+j) QIlr[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			if(i+j) QIa_dd[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
			if(i+j) QI_dd[i][j] = new Table<double>(l1, l2, subl1, subl2, 0.0);
		}
	}
	QImm[2][2] = new Table<double>(l1, l2, subl1, subl2, 0.0);
	QIa_nd[1][0] = new Table<double>(l1, l2, subl1, subl2, 0.0);
	QI_nd[1][0] = new Table<double>(l1, l2, subl1, subl2, 0.0);
	QIa_dn[0][1] = new Table<double>(l1, l2, subl1, subl2, 0.0);
	QI_dn[0][1] = new Table<double>(l1, l2, subl1, subl2, 0.0);
	QIy = new Table<double>(l1, l2, subl1, subl2, 0.0);
	QI = new Table<double>(l1, l2, subl1, subl2, 0.0);
	QIa = new Table<double>(l1, l2, subl1, subl2, 0.0);

	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
		{
			auxQIa[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*5);
			*auxQIa[i][j] = QIa_nn[i][j];
			*(auxQIa[i][j] + 1) = QIa_nd[i][j];
			*(auxQIa[i][j] + 2) = QIa_dn[i][j];
			*(auxQIa[i][j] + 3) = QIa_dd[i][j];
			*(auxQIa[i][j] + 4) = NULL;
			auxQIa_dd[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*3);
			*auxQIa_dd[i][j] = QIa_dd[i][j];
			*(auxQIa_dd[i][j] + 1) = QIb[i][j];
			*(auxQIa_dd[i][j] + 2) = NULL;
			QIx[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*3);
			*QIx[i][j] = QIx_m[i][j];
			*(QIx[i][j] + 1) = QIx_k[i][j];
			*(QIx[i][j] + 2) = NULL;
		}
	if(!quiet) printf("Allocation of %ld MB of memory and its initialization finished. \n", total_tables_size / (1024*1024));
	computeDoubleParFunc();
}

PartitionFunction::PartitionFunction(FILE *inp, Energy *en, PartitionFunction *par[2], bool _debug, int _procNum, bool _quiet, bool _no_isolate)
{
	double temp;
	energy = en;
	energy->getParams(&NA, &polymer, &naConc, &mgConc, &suffix, &temp, &zip, &nodangle, &energyType);
	RT = energy->getRT();

	DEBUG = _debug;
	procNum = _procNum;
	quiet = _quiet;
	no_isolate = _no_isolate;
	parFuncType = DUPLEX;
	sequence = NULL;
	Qb = NULL;
	for(int j = 0; j < 3; j++)
	{
		Q[j] = NULL;
		Qm[j] = NULL;
	}

	for(int i = 0; i < 2; i++)
		parfunc[i] = par[i];

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(i+j) QIa_nn[i][j] = new Table<double>();
			if(i+j) QI_nn[i][j] = new Table<double>();
		}
	}

	for(int i = 0; i < 2; i++)
	{
		QIh[i] = new Table<double>();
		QIx_m[i][2] = new Table<double>();
		QImm[i][2] = new Table<double>();
		QImm[2][i] = new Table<double>();
		QImk[i][1] = new Table<double>();
		QImk[2][i] = new Table<double>();
		QIkm[1][i] = new Table<double>();
		QIkm[i][2] = new Table<double>();
		QIa_nd[i][1] = new Table<double>();
		QIa_nd[2][i] = new Table<double>();
		QI_nd[i][1] = new Table<double>();
		QI_nd[2][i] = new Table<double>();
		QIa_dn[1][i] = new Table<double>();
		QIa_dn[i][2] = new Table<double>();
		QI_dn[1][i] = new Table<double>();
		QI_dn[i][2] = new Table<double>();
		for(int j = 0; j < 2; j++)
		{
			QIb[i][j] = new Table<double>();
			QIhh[i][j] = new Table<double>();
			QIhm[i][j] = new Table<double>();
			QIx_m[i][j] = new Table<double>();
			QIx_k[i][j] = new Table<double>();
			if(i+j) QIkk[i][j] = new Table<double>();
			if(i+j) QIa_r[i][j] = new Table<double>();
			if(i+j) QI_r[i][j] = new Table<double>();
			if(i+j) QIlr[i][j] = new Table<double>();
			if(i+j) QIa_dd[i][j] = new Table<double>();
			if(i+j) QI_dd[i][j] = new Table<double>();
		}
	}
	QImm[2][2] = new Table<double>();
	QIa_nd[1][0] = new Table<double>();
	QI_nd[1][0] = new Table<double>();
	QIa_dn[0][1] = new Table<double>();
	QI_dn[0][1] = new Table<double>();
	QIy = new Table<double>();
	QI = new Table<double>();
	QIa = new Table<double>();

	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
		{
			auxQIa[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*5);
			*auxQIa[i][j] = QIa_nn[i][j];
			*(auxQIa[i][j] + 1) = QIa_nd[i][j];
			*(auxQIa[i][j] + 2) = QIa_dn[i][j];
			*(auxQIa[i][j] + 3) = QIa_dd[i][j];
			*(auxQIa[i][j] + 4) = NULL;
			auxQIa_dd[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*3);
			*auxQIa_dd[i][j] = QIa_dd[i][j];
			*(auxQIa_dd[i][j] + 1) = QIb[i][j];
			*(auxQIa_dd[i][j] + 2) = NULL;
			QIx[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*3);
			*QIx[i][j] = QIx_m[i][j];
			*(QIx[i][j] + 1) = QIx_k[i][j];
			*(QIx[i][j] + 2) = NULL;
		}
	retrieve(inp);
}

PartitionFunction::~PartitionFunction()
{
	if(parFuncType == SINGLE)
	{
		for(int j = 0; j < 3; j++)
		{
			delete Q[j];
			delete Qm[j];
		}
		delete Qb;
	} else
	{
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				if(i+j) delete QIa_nn[i][j];
				if(i+j) delete QI_nn[i][j];
			}
		}
		delete QImm[2][2];
		delete QIa_nd[1][0];
		delete QI_nd[1][0];
		delete QIa_dn[0][1];
		delete QI_dn[0][1];
		delete QIy;
		delete QI;
		delete QIa;
		for(int i = 0; i < 2; i++)
		{
			delete QIkm[1][i];
			delete QIkm[i][2];
			delete QImk[i][1];
			delete QImk[2][i];
			delete QIa_nd[i][1];
			delete QIa_nd[2][i];
			delete QI_nd[i][1];
			delete QI_nd[2][i];
			delete QIa_dn[1][i];
			delete QIa_dn[i][2];
			delete QI_dn[1][i];
			delete QI_dn[i][2];
			delete QIh[i];
			delete QIx_m[i][2];
			delete QImm[i][2];
			delete QImm[2][i];
			for(int j = 0; j < 2; j++)
			{
				free(auxQIa[i][j]);
				free(auxQIa_dd[i][j]);
				free(QIx[i][j]);
				delete QIb[i][j];
				delete QIhh[i][j];
				delete QIhm[i][j];
				if(i+j) delete QIkk[i][j];
				if(i+j) delete QIa_r[i][j];
				if(i+j) delete QI_r[i][j];
				if(i+j) delete QIlr[i][j];
				if(i+j) delete QIa_dd[i][j];
				if(i+j) delete QI_dd[i][j];
				delete QIx_m[i][j];
				delete QIx_k[i][j];
			}
		}
	}
}

double PartitionFunction::freeEnsembleEnergy(double h)
{
	return -RT*log(Zz*h);
}

double PartitionFunction::Z(double h)
{
	return Zz*h;
}

double PartitionFunction::intProb()
{
	return 0.0;
}

void PartitionFunction::Reverse()
{
	if(parFuncType == SINGLE)
	{
		for(int j = 0; j < 3; j++)
		{
			Q[j]->Reverse();
			Qm[j]->Reverse();
		}
		if(Qb) Qb->Reverse();
	}	
}

double PartitionFunction::computeSingleParFunc()
{
	if(parFuncType == DUPLEX)
		return 0.0;

	if(!quiet) printf("Computing single partition function for %s whose length is %d with sublength %d... \n", sequence->getName(), sequence->getLen(), sublen);

	unsigned int len = sequence->getLen();
	for(unsigned int l = 1; l <= sublen; l++)
	{
#pragma omp parallel for num_threads(procNum)
		for(int i = 0; (unsigned int)i <= len - l; i++)
		{
			
			int j = i + l - 1;
			if(DEBUG) printf("i:%d\tj:%d\n", i, j);
			computeSingleQbQm(i, j);
			computeSingleQ(i, j);
		}
		
		if(!quiet && !(l % 10)) printf("Length %d done.\n", l);		
	}

	if(!quiet) printf("Computing single partition function finished. \n");

	if(len == sublen)
		Zz = Q[0]->element(0, len - 1);
	else
		Zz = 1;
	return Zz;
}

void PartitionFunction::computeSingleQbQm(int i, int j)
{
	double *b = Qb->estar(i, j);
	unsigned char *seq = sequence->getSeq();
			
	if(j > i+3 && energy->base_pair(seq[i], seq[j]))
	{
		*b = exp(-energy->Eh(i, j, seq)/RT);
       		if (DEBUG) printf("Qb %lf %lf %lf %lf\n", energy->Eh(i,j,seq), RT, energy->Eh(i,j,seq)/RT, *b);
		for(int d = i+1; d <= j-5; d++)
			for(int e = d+4; e <= j-1; e++)
			{
				if (DEBUG) printf("Qb - b [%d,%d]\n", d, e);
				if (e == j-1 && d == i+1)
					*b += Qb->element(d, e)*exp(-energy->Es(i, j, seq)/RT);
					
				else 
				{
					*b += Qb->element(d, e)*exp(-energy->Ebi(i, j, d, e, seq)/RT);
					if(d > i+1)
						*b += Qb->element(d, e)*Qm[2]->element(i + 1, d - 1)*exp(-energy->Emulti(1, j-e-1, 1)/RT); 
				}
			}
	}

	double *mw = Qm[0]->estar(i, j);
	double *mr = Qm[1]->estar(i, j);
	double *mg = Qm[2]->estar(i, j);

	for(int d = i; d <= j - 4; d++)
		for(int e = d + 4; e <= j; e++)
		{
			if (d > i)
			{ 
				*mw += Qm[0]->element(i, d-1)*Qb->element(d, e); 
				*mr += Qm[1]->element(i, d-1)*Qb->element(d, e)*exp(-energy->Ekissing(0, j-e, 1)/RT);
				*mg += Qm[2]->element(i, d-1)*Qb->element(d, e)*exp(-energy->Emulti(0, j-e, 1)/RT);
			}					
			*mw += Qb->element(d, e); 
			*mr += Qb->element(d, e)*exp(-energy->Ekissing(0, j-i-(e-d), 1)/RT); 
			*mg += Qb->element(d, e)*exp(-energy->Emulti(0, j-i-(e-d), 1)/RT); 
		}
	
}

void PartitionFunction::computeSingleQ(int i, int j)
{
	*(Q[0]->estar(i, j)) = 1.0;
	*(Q[1]->estar(i, j)) = exp(-energy->Ekissing(0, j - i + 1, 0)/RT);
	*(Q[2]->estar(i, j)) = exp(-energy->Emulti(0, j - i + 1, 0)/RT);
	

	for(int d = i; d <= j-4; d++)
		for(int e = d+4; e <= j; e++)
		{
			double temp1 = 1.0, temp2 = 1.0, temp3 = 1.0;
			if (d > i)
			{
				temp1 = Q[0]->element(i, d-1);
				temp2 = Q[1]->element(i, d-1);
				temp3 = Q[2]->element(i, d-1);
			}
			*(Q[0]->estar(i, j)) += temp1*(Qb->element(d, e));
   			*(Q[1]->estar(i, j)) += temp2*(Qb->element(d, e))*exp(-energy->Ekissing(0, j-e, 1)/RT);
   			*(Q[2]->estar(i, j)) += temp3*(Qb->element(d, e))*exp(-energy->Emulti(0, j-e, 1)/RT);
		}
}


/* Hamid: The two RNAs interact 3'->5' and 5'->3'. But I assume that one RNA is reversed beforehand so the indices now interact in the same direction. */

double PartitionFunction::computeDoubleParFunc()
{
	Sequence *sequenceUp = parfunc[0]->getSequence();
	Sequence *sequenceDown = parfunc[1]->getSequence();

	seqUp = sequenceUp->getSeq();	
	seqDown = sequenceDown->getReverse();

	unsigned int lenUp = sequenceUp->getLen();
	unsigned int lenDown = sequenceDown->getLen();
	unsigned int sublenUp = parfunc[0]->getSubLen();
	unsigned int sublenDown = parfunc[1]->getSubLen();
	parfunc[1]->Reverse();	

	if(!quiet) printf("Computing double partition function for %s and %s whose lengths are %d and %d with sublengths %d and %d... \n", sequenceUp->getName(), sequenceDown->getName(), sequenceUp->getLen(), sequenceDown->getLen(), sublenUp, sublenDown);

	if(DEBUG) printf("tables: %d mem: %ld\n", big_tables_num, total_tables_size);
	for(unsigned int l1 = 1; l1 <= sublenUp; l1++)
	{
		for(unsigned int l2 = 1; l2 <= sublenDown; l2++)
		{
			if(((l1 - 1) * sublenDown + l2) % 10 == 0) if(!quiet) printf("Iteration %d out of %d done.\n", (l1 - 1) * sublenDown + l2, sublenUp*sublenDown);

/*			for(int i1 = 0; (unsigned int)i1 <= lenUp - l1; i1++)
				for(unsigned int i2 = 0; i2 <= lenDown - l2; i2++)*/
#pragma omp parallel for num_threads(procNum)
			for(int idx = 0; (unsigned int)idx < (lenUp - l1 + 1)*(lenDown - l2 + 1); idx++)
			{
				int i1 = idx / (lenDown - l2 + 1);
				int i2 = idx % (lenDown - l2 + 1);
				int j1 = i1 + l1 - 1;
				int j2 = i2 + l2 - 1;
					
				if(DEBUG) printf("%d\t%d\t%d\t%d\tbands\t", i1, j1, i2, j2);
				computeBands(i1, j1, i2, j2);
				if(DEBUG) printf("xy\t");
				computeXY(i1, j1, i2, j2);
				if(DEBUG) printf("aux\t");
				computeAux(i1, j1, i2, j2);
				if(DEBUG) printf("Ib, Ia, I\t");
				computeI(i1, j1, i2, j2);
				if(DEBUG) printf("mmkk\t");
				computeImmkk(i1, j1, i2, j2);
				if(DEBUG) printf("Ia_nndd\t");
				computeIa_nndd(i1, j1, i2, j2);
				if(DEBUG) printf("I_nndd\n");
				computeI_nndd(i1, j1, i2, j2);
			}
		}
/*		for(unsigned int l2 = 1; l2 <= l1; l2++)
		{
			for(int idx = 0; (unsigned int)idx < (lenUp - l1 + 1)*(lenDown - l2 + 1); idx++)
			{
				int i1 = idx / (lenDown - l2 + 1);
				int i2 = idx % (lenDown - l2 + 1);
				int j1 = i1 + l1 - 1;
				int j2 = i2 + l2 - 1;
				selftest(i1, j1, i2, j2, lenUp, lenDown);
			}
		} */
	}

	if(!quiet) printf("Computing double partition function finished. \n");

	parfunc[1]->Reverse();	
	
	if(sublenUp < lenUp || sublenDown < lenDown)
		Zz = 0.0;
	else
		Zz = QI->element(0, sequenceUp->getLen() - 1, 0, sequenceDown->getLen() - 1);
	return Zz;
}

double PartitionFunction::product(Table<double> *T1, Table<double> *T2, int i1, int j1, int i2, int j2)
{
	double result = 0.0;

	for(int d = i1; d <= j1-1; d++)
		for(int e = i2; e <= j2-1; e++)
			result += T1->element(i1, d, i2, e)*T2->element(d+1, j1, e+1, j2);

	return result;
}

double PartitionFunction::product(Table<double> **T1, Table<double> **T2, int i1, int j1, int i2, int j2)
{
	double result = 0.0;

	for(int d = i1; d <= j1-1; d++)
		for(int e = i2; e <= j2-1; e++)
		{
			double sum1 = 0.0;
			int  k = 0;
			while(T1[k]) 
			{
				sum1 += T1[k]->element(i1, d, i2, e);
				k++;
			}
			double sum2 = 0.0;
			k = 0;
			while(T2[k]) 
			{
				sum2 += T2[k]->element(d+1, j1, e+1, j2);
				k++;
			}
			result += sum1*sum2;
		}
	
	return result;
}

double PartitionFunction::product(Table<double> **T1, Table<double> *T2, int i1, int j1, int i2, int j2)
{
	double result = 0.0;

	for(int d = i1; d <= j1-1; d++)
		for(int e = i2; e <= j2-1; e++)
		{
			double sum1 = 0.0;
			int  k = 0;
			while(T1[k]) 
			{
				sum1 += T1[k]->element(i1, d, i2, e);
				k++;
			}
			result += sum1*T2->element(d+1, j1, e+1, j2);
		}

	return result;
}

double PartitionFunction::product(Table<double> *T1, Table<double> **T2, int i1, int j1, int i2, int j2)
{
	double result = 0.0;

	for(int d = i1; d <= j1-1; d++)
		for(int e = i2; e <= j2-1; e++)
		{
			double sum2 = 0.0;
			int k = 0;
			while(T2[k]) 
			{
				sum2 += T2[k]->element(d+1, j1, e+1, j2);
				k++;
			}
			result += T1->element(i1, d, i2, e)*sum2;
		}


	return result;
}

double PartitionFunction::rsemiproduct(Table<double> *Tu, Table<double> *Td, Table<double> *T2, int i1, int j1, int i2, int j2)
{
	double result = 0.0;

	for(int d = i1; d <= j1; d++)
		for(int e = i2; e <= j2; e++)
		{
			double temp1 = 1.0, temp2 = 1.0;
			if(d > i1) temp1 = Tu->element(i1, d-1);
			if(e > i2) temp2 = Td->element(i2, e-1);
			result += temp1*temp2*T2->element(d, j1, e, j2);
		}

	return result;
}

double PartitionFunction::lsemiproduct(Table<double> *T1, Table<double> *Tu, Table<double> *Td, int i1, int j1, int i2, int j2)
{
	double result = 0.0;
	for(int d = i1; d <= j1; d++){
		double temp1 = 1.0;
		if(d < j1) temp1 = Tu->element(d+1,j1);
		result += T1->element(i1, d,i2,j2)*temp1;
	}
	for(int e = i2; e < j2; e++){
		double temp1 = Td->element(e+1,j2);
		result += T1->element(i1,j1,i2,e)*temp1;
	}
	return result;
}

double PartitionFunction::lsemiproduct(Table<double> **T1, Table<double> *Tu, Table<double> *Td, int i1, int j1, int i2, int j2)
{
	double result = 0.0;
	for(int d = i1; d <= j1; d++){
		double sum1 = 0.0;
		int  k = 0;
		while(T1[k]) 
		{
			sum1 += T1[k]->element(i1, d, i2, j2);
			k++;
		}
		double temp1 = 1.0;
		if(d < j1) temp1 = Tu->element(d+1,j1);
		result += sum1*temp1;
	}
	for(int e = i2; e < j2; e++){
		double sum1 = 0.0;
		int  k = 0;
		while(T1[k]) 
		{
			sum1 += T1[k]->element(i1, j1, i2, e);
			k++;
		}
		result += sum1*Td->element(e+1,j2);
	}
	return result;
}

void PartitionFunction::computeBands(int i1, int j1, int i2, int j2)
{
	if(!energy->base_pair(seqUp[i1], seqDown[i2]))
		return;

	//QIh
	if(energy->base_pair(seqUp[j1], seqDown[j2]))
	{
		for(int i = 0; i < 2; i++)
		{
			double *qih = QIh[i]->estar(i1, j1, i2, j2);
			if(i1 == j1 && i2 == j2 && !no_isolate)
				*qih = 1.0;

			if(j1 == i1+1 && j2 == i2+1)
				if(i == 0)
					*qih += exp(-energy->Eintstack(0, i1, i2, seqUp, seqDown)/RT);
				else
					*qih += exp(-energy->Ekissingstack(i1, i2, seqUp, seqDown)/RT);
			else
				if(i == 0)
					*qih += exp(-energy->Eintbi(i1, i2, j1, j2, seqUp, seqDown)/RT);
				else
					*qih += exp(-energy->Ekissingbi(i1, i2, j1, j2, seqUp, seqDown)/RT);

			for(int d = i1+1; d < j1 && d < i1 + MAX_HYBRID_LEN; d++)
				for(int e = i2+1; e < j2 && e < i2 + MAX_HYBRID_LEN; e++)
				{
					double temp = QIh[i]->element(d, j1, e, j2);

					if(d == i1+1 && e == i2+1)
						if(i == 0)
							*qih += temp*exp(-energy->Eintstack(0, i1, i2, seqUp, seqDown)/RT);
						else
							*qih += temp*exp(-energy->Ekissingstack(i1, i2, seqUp, seqDown)/RT);
					else
						if(i == 0)
							*qih += temp*exp(-energy->Eintbi(i1, i2, d, e, seqUp, seqDown)/RT);
						else
							*qih += temp*exp(-energy->Ekissingbi(i1, i2, d, e, seqUp, seqDown)/RT);

				}
		}
	}

	//QIhh and QIhm
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
		{
			double *qihh = QIhh[i][j]->estar(i1, j1, i2, j2);
			double *qihm = QIhm[i][j]->estar(i1, j1, i2, j2);

			for(int d = i1; d <= j1; d++)
				for(int e = i2; e <= j2; e++)
				{
					double temp1 = 1.0, temp2 = 1.0, temp = 0.0, 
					temp3 = (i+j) ? exp(-energy->Ekissing(1, 0, 0)/RT) : exp(-energy->Eintstackpenalty()/RT), temp4;
					double empty = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, j1 - d, 0)/RT);
					if(d < j1 && e < j2)
					{
						temp1 = parfunc[0]->getQ(i)->element(d+1, j1);
						temp2 = parfunc[1]->getQ(j)->element(e+1, j2);
						temp = parfunc[0]->getQm(i)->element(d+1, j1)*parfunc[1]->getQ(j)->element(e+1, j2) +
							empty*parfunc[1]->getQm(j)->element(e+1, j2);
					} else if(d < j1 && e == j2)
					{
						temp1 = parfunc[0]->getQ(i)->element(d+1, j1);
						temp = parfunc[0]->getQm(i)->element(d+1, j1);
					} else if(d == j1 && e < j2)
					{
						temp2 = parfunc[1]->getQ(j)->element(e+1, j2);
						temp = parfunc[1]->getQm(j)->element(e+1, j2);
					}

					temp4 = 1.0;

					// penalize A-U stack terminal
					if((seqUp[i1] == 0 && seqDown[i2] == 3) || (seqUp[i1] == 3 && seqDown[i2] == 0))
						temp4 *= exp(-energy->Eaupenalty()/RT);
					if((seqUp[d] == 0 && seqDown[e] == 3) || (seqUp[d] == 3 && seqDown[e] == 0))
						temp4 *= exp(-energy->Eaupenalty()/RT);
					*qihh += QIh[i|j]->element(i1, d, i2, e)*temp1*temp2*temp3*temp4;
					*qihm += QIh[i|j]->element(i1, d, i2, e)*temp*temp3*temp4;
				
				}
		}

}

double PartitionFunction::capproduct(Table<double> *qir, Table<double> *qig, Table<double> *qi, int i1, int j1, int i2, int j2, int upOrDown)
{
	double result = 0.0;

	int ii = (upOrDown == 0) ? i1 : i2;
	int jj = (upOrDown == 0) ? j1 : j2;

	Table<double> *q[3] = {parfunc[upOrDown]->getQ(0), parfunc[upOrDown]->getQ(1), parfunc[upOrDown]->getQ(2)};
	Table<double> *qm[3] = {parfunc[upOrDown]->getQm(0), parfunc[upOrDown]->getQm(1), parfunc[upOrDown]->getQm(2)};

	unsigned char *seq = parfunc[upOrDown]->getSequence()->getSeq();
	unsigned int len = parfunc[upOrDown]->getSequence()->getLen();
	
	bool tempdebug = (i1 == 2 && j1 == 8 && i2 == 3 && j2 == 5 && qi == QIx_k[0][0]) || 
			(i1 == 29 && j1 == 31 && i2 == 26 && j2 == 32 && qi == QIx_k[1][0]);


	for(int d = ii; d <= jj-2; d++)
		for(int e = d+2; e <= jj; e++)
		{
			double tig = (upOrDown == 0) ? qig->element(d+1, e-1, i2, j2) : qig->element(i1, j1, d+1, e-1);
			double tir = (upOrDown == 0) ? qir->element(d+1, e-1, i2, j2) : qir->element(i1, j1, d+1, e-1);
			double ti = (upOrDown == 0) ? qi->element(d+1, e-1, i2, j2) : qi->element(i1, j1, d+1, e-1);
       			double rtemp1=1.0, rtemp2=1.0, gtemp1=1.0, gtemp2=1.0;
	
       			if (d > ii) rtemp1=q[1]->element(ii+1, d);
			if (e < jj) rtemp2=q[1]->element(e, jj-1);
       			if (d > ii) gtemp1=q[2]->element(ii+1, d);
			if (e < jj) gtemp2=q[2]->element(e, jj-1);

			result += rtemp1*rtemp2*exp(-energy->Ekissing(0, 0, 1)/RT)*tir +
				gtemp1*gtemp2*exp(-energy->Emulti(0, 0, 1)/RT)*tig;

			if (e == jj && d == ii)
			{
				double es = (upOrDown == 0) ? energy->Es(ii, jj, seq) : energy->Es(len-jj-1, len-ii-1, seq);
				result += exp(-es/RT)*ti;
			}
			else
			{
				double ebi = (upOrDown == 0) ? energy->Ebi(ii, jj, d+1, e-1, seq) : energy->Ebi(len-jj-1, len-ii-1, len-e, len-d-2, seq);
				result += exp(-ebi/RT)*ti;
				if(d > ii)
					result += qm[2]->element(ii+1, d)*exp(-energy->Emulti(1, jj-e, 1)/RT)*ti;
				if(e < jj)
					result += qm[2]->element(e, jj-1)*exp(-energy->Emulti(1, d-ii, 1)/RT)*ti;
				if(d > ii && e < jj)
					result += qm[2]->element(ii+1, d)*qm[2]->element(e, jj-1)*exp(-energy->Emulti(1, 0, 1)/RT)*ti;
			}
		}

	return result;
}

void PartitionFunction::computeXY(int i1, int j1, int i2, int j2)
{
	//QIy
	if(energy->base_pair(seqUp[i1], seqUp[j1]) && energy->base_pair(seqDown[i2], seqDown[j2]) && j1 > i1 + MIN_HAIRPIN_SIZE && j2 > i2 + MIN_HAIRPIN_SIZE)
		*QIy->estar(i1, j1, i2, j2) = capproduct(QIx_k[0][1], QIx_m[0][2], QIy, i1, j1, i2, j2, 1); 

	//QIx_m
	for(int c = 0; c < 3; c++)
	{
		if(energy->base_pair(seqUp[i1], seqUp[j1]) && j1 > i1 + MIN_HAIRPIN_SIZE)
			*QIx_m[0][c]->estar(i1, j1, i2, j2) = capproduct(QIkm[1][c], QImm[2][c], QIx_m[0][c], i1, j1, i2, j2, 0); 

		if(energy->base_pair(seqDown[i2], seqDown[j2]) && j2 > i2 + MIN_HAIRPIN_SIZE)
			*QIx_m[1][c]->estar(i1, j1, i2, j2) = capproduct(QImk[c][1], QImm[c][2], QIx_m[1][c], i1, j1, i2, j2, 1); 
	}

	//QIx_k
	for(int c = 0; c < 2; c++)
	{
		if(energy->base_pair(seqUp[i1], seqUp[j1]) && j1 > i1 + MIN_HAIRPIN_SIZE)
			*QIx_k[0][c]->estar(i1, j1, i2, j2) = capproduct(QIkk[1][c], QImk[2][c], QIx_k[0][c], i1, j1, i2, j2, 0);

		if(energy->base_pair(seqDown[i2], seqDown[j2]) && j2 > i2 + MIN_HAIRPIN_SIZE)
			*QIx_k[1][c]->estar(i1, j1, i2, j2) = capproduct(QIkk[c][1], QIkm[c][2], QIx_k[1][c], i1, j1, i2, j2, 1);
	}
}

void PartitionFunction::computeAux(int i1, int j1, int i2, int j2)
{
	if(!energy->base_pair(seqUp[j1], seqDown[j2]))
		return;
	
	
	//QI_r, QIlr, QIa_r
	for(int c = 1; c < 4; c++)
	{
		int i = c % 2;
		int j = c / 2;

		double *qi_r = QI_r[i][j]->estar(i1, j1, i2, j2);
		double *qia_r = QIa_r[i][j]->estar(i1, j1, i2, j2);
		double *qilr = QIlr[i][j]->estar(i1, j1, i2, j2);
		
		double temp1 = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);
		double temp2 = (j == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);
		double temp3 = (i+j) ? exp(-energy->Ekissing(1, 0, 0)/RT) : exp(-energy->Eintstackpenalty()/RT);
		double temp4 = 1.0;

		// penalize A-U stack terminal
		if((seqUp[i1] == 0 && seqDown[i2] == 3) || (seqUp[i1] == 3 && seqDown[i2] == 0))
			temp4 *= exp(-energy->Eaupenalty()/RT);
		if((seqUp[j1] == 0 && seqDown[j2] == 3) || (seqUp[j1] == 3 && seqDown[j2] == 0))
			temp4 *= exp(-energy->Eaupenalty()/RT);

		*qilr = QIh[1]->element(i1, j1, i2, j2)*temp3*temp4;
		*qilr += product(QIhm[i][j], QIlr[i][j], i1, j1, i2, j2) + product(QIhh[i][j], QIa_r[i][j], i1, j1, i2, j2);
		*qia_r = product(QIx[0][j], QI_r[i][j], i1, j1, i2, j2)*temp1 +
			product(QIx[1][i], QI_r[i][j], i1, j1, i2, j2)*temp2 +
			product(QIy, QI_r[i][j], i1, j1, i2, j2)*temp1*temp2;

		*qi_r = rsemiproduct(parfunc[0]->getQ(i), parfunc[1]->getQ(j), QIlr[i][j], i1, j1, i2, j2) +
			rsemiproduct(parfunc[0]->getQ(i), parfunc[1]->getQ(j), QIa_r[i][j], i1, j1, i2, j2);
	}
}

void PartitionFunction::computeImmkk(int i1, int j1, int i2, int j2)
{
	//QImm
	double temp2 = exp(-energy->Emulti(0, 0, 1)/RT);
	double temp3 = exp(-energy->Ekissing(0, 0, 1)/RT);
	for(int i = 0; i < 2; i++)
	{
		double temp1 = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);

		*QImm[i][2]->estar(i1, j1, i2, j2) = product(QIa_nn[i][2], QIx_m[0][2], i1, j1, i2, j2)*temp1 + 
							product(QIa_nn[i][2], QIx_m[1][i], i1, j1, i2, j2)*temp2 + 
							product(QIa_nn[i][2], QIy, i1, j1, i2, j2)*temp1*temp2;

		*QImm[2][i]->estar(i1, j1, i2, j2) = product(QIa_nn[2][i], QIx_m[0][i], i1, j1, i2, j2)*temp2 + 
							product(QIa_nn[2][i], QIx_m[1][2], i1, j1, i2, j2)*temp1 + 
							product(QIa_nn[2][i], QIy, i1, j1, i2, j2)*temp1*temp2;

	}
	*QImm[2][2]->estar(i1, j1, i2, j2) = product(QIa_nn[2][2], QIx_m[0][2], i1, j1, i2, j2)*temp2 + 
						product(QIa_nn[2][2], QIx_m[1][2], i1, j1, i2, j2)*temp2 + 
						product(QIa_nn[2][2], QIy, i1, j1, i2, j2)*temp2*temp2;
	//QImk
	for(int i = 0; i < 2; i++)
	{
		double temp1 = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);

		*QImk[i][1]->estar(i1, j1, i2, j2) = product(QIa_nd[i][1], QIx[0][1], i1, j1, i2, j2)*temp1 +
							product(QIa_nd[i][1], QIx_m[1][i], i1, j1, i2, j2)*temp3 + 
							product(QIa_nn[i][1], QIx_k[0][1], i1, j1, i2, j2)*temp1 +
							product(QIa_nd[i][1], QIy, i1, j1, i2, j2)*temp1*temp3;


		*QImk[2][i]->estar(i1, j1, i2, j2) = product(QIa_nd[2][i], QIx[0][i], i1, j1, i2, j2)*temp2 +
							product(QIa_nd[2][i], QIx_m[1][2], i1, j1, i2, j2)*temp1 + 
							product(QIa_nn[2][i], QIx_k[0][i], i1, j1, i2, j2)*temp2 +
							product(QIa_nd[2][i], QIy, i1, j1, i2, j2)*temp2*temp3;
	}

	//QIkm
	for(int i = 0; i < 2; i++)
	{
		double temp1 = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);

		*QIkm[1][i]->estar(i1, j1, i2, j2) = product(QIa_dn[1][i], QIx_m[0][i], i1, j1, i2, j2)*temp3 +
							product(QIa_dn[1][i], QIx[1][1], i1, j1, i2, j2)*temp1 + 
							product(QIa_nn[1][i], QIx_k[1][1], i1, j1, i2, j2)*temp1 +
							product(QIa_dn[1][i], QIy, i1, j1, i2, j2)*temp1*temp3;

		*QIkm[i][2]->estar(i1, j1, i2, j2) = product(QIa_dn[i][2], QIx_m[0][2], i1, j1, i2, j2)*temp1 +
							product(QIa_dn[i][2], QIx[1][i], i1, j1, i2, j2)*temp2 + 
							product(QIa_nn[i][2], QIx_k[1][i], i1, j1, i2, j2)*temp2 +
							product(QIa_dn[i][2], QIy, i1, j1, i2, j2)*temp1*temp2;
	}

	//QIkk
	for(int c = 1; c < 4; c++)
	{
		int i = c % 2;
		int j = c / 2;

		double temp1 = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);
		double temp4 = (j == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);

		double *qikk = QIkk[i][j]->estar(i1, j1, i2, j2);
		*qikk = QIa_r[i][j]->element(i1, j1, i2, j2) + QIlr[i][j]->element(i1, j1, i2, j2); 
		*qikk += product(auxQIa_dd[i][j], QIx[0][j], i1, j1, i2, j2)*temp1;
		*qikk += product(auxQIa_dd[i][j], QIx[1][i], i1, j1, i2, j2)*temp4;
		*qikk += product(auxQIa_dd[i][j], QIy, i1, j1, i2, j2)*temp1*temp4;
		*qikk += product(QIa_dn[i][j], QIx_k[0][j], i1, j1, i2, j2)*temp1;
		*qikk += product(QIa_nd[i][j], QIx_k[1][i], i1, j1, i2, j2)*temp4;
	}
}

void PartitionFunction::computeIa_nndd(int i1, int j1, int i2, int j2)
{
	//QIa_nn
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			if(i+j)
			{
				double *qia_nn = QIa_nn[i][j]->estar(i1, j1, i2, j2);
				double temp3 = 1.0, temp4 = 1.0;
				if(i == 1) temp3 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(i == 2) temp3 = exp(-energy->Emulti(0, 0, 1)/RT);
				if(j == 1) temp4 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(j == 2) temp4 = exp(-energy->Emulti(0, 0, 1)/RT);

				*qia_nn = product(QIx_m[0][j], QI_nn[i][j], i1, j1, i2, j2)*temp3 + 
						product(QIx_m[1][i], QI_nn[i][j], i1, j1, i2, j2)*temp4 +
						product(QIy, QI_nn[i][j], i1, j1, i2, j2)*temp3*temp4;

				*qia_nn += lsemiproduct(QIx_m[0][j], parfunc[0]->getQ(i), parfunc[1]->getQ(j),i1, j1, i2, j2)*temp3 + 
						lsemiproduct(QIx_m[1][i], parfunc[0]->getQ(i), parfunc[1]->getQ(j),i1, j1, i2, j2)*temp4 +
						lsemiproduct(QIy, parfunc[0]->getQ(i), parfunc[1]->getQ(j),i1, j1, i2, j2)*temp3*temp4;
			}
	//QIa_nd
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 2; j++)
			if(j == 1 || (i != 0 && j == 0))
			{
				double *qia_nd = QIa_nd[i][j]->estar(i1, j1, i2, j2);
				double temp3 = 1.0, temp4 = 1.0;
				if(i == 1) temp3 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(i == 2) temp3 = exp(-energy->Emulti(0, 0, 1)/RT);
				if(j == 1) temp4 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(j == 2) temp4 = exp(-energy->Emulti(0, 0, 1)/RT);

				*qia_nd = product(QIx[0][j], QI_nd[i][j], i1, j1, i2, j2)*temp3 + 
						product(QIx_k[0][j], QI_nn[i][j], i1, j1, i2, j2)*temp3 +
						product(QIx_m[1][i], QI_nd[i][j], i1, j1, i2, j2)*temp4 +
						product(QIy, QI_nd[i][j], i1, j1, i2, j2)*temp3*temp4;

				*qia_nd += lsemiproduct(QIx_k[0][j], parfunc[0]->getQ(i), parfunc[1]->getQ(j), i1, j1, i2, j2)*temp3;
			}
	//QIa_dn
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 3; j++)
			if(i == 1 || (j != 0 && i == 0))
			{
				double *qia_dn = QIa_dn[i][j]->estar(i1, j1, i2, j2);
				double temp3 = 1.0, temp4 = 1.0;
				if(i == 1) temp3 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(i == 2) temp3 = exp(-energy->Emulti(0, 0, 1)/RT);
				if(j == 1) temp4 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(j == 2) temp4 = exp(-energy->Emulti(0, 0, 1)/RT);

				*qia_dn = product(QIx_m[0][j], QI_dn[i][j], i1, j1, i2, j2)*temp3 + 
						product(QIx[1][i], QI_dn[i][j], i1, j1, i2, j2)*temp4 +
						product(QIx_k[1][i], QI_nn[i][j], i1, j1, i2, j2)*temp4 +
						product(QIy, QI_dn[i][j], i1, j1, i2, j2)*temp3*temp4;

				*qia_dn += lsemiproduct(QIx_k[1][i], parfunc[0]->getQ(i), parfunc[1]->getQ(j),i1, j1, i2, j2)*temp4;
			}

	//QIa_dd
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			if(i+j)
			{
				double *qia_dd = QIa_dd[i][j]->estar(i1, j1, i2, j2);
				double temp3 = 1.0, temp4 = 1.0;
				if(i == 1) temp3 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(j == 1) temp4 = exp(-energy->Ekissing(0, 0, 1)/RT);

				*qia_dd = product(QIx[0][j], QI_dd[i][j], i1, j1, i2, j2)*temp3 + 
						product(QIx_k[0][j], QI_dn[i][j], i1, j1, i2, j2)*temp3 +
						product(QIx_k[1][i], QI_nd[i][j], i1, j1, i2, j2)*temp4 +
						product(QIx[1][i], QI_dd[i][j], i1, j1, i2, j2)*temp4 +
						product(QIy, QI_dd[i][j], i1, j1, i2, j2)*temp3*temp4;
			}
}

void PartitionFunction::computeI_nndd(int i1, int j1, int i2, int j2)
{
	//QI_nn
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			if(i+j)
				*QI_nn[i][j]->estar(i1, j1, i2, j2) = parfunc[0]->getQ(i)->element(i1, j1)*parfunc[1]->getQ(j)->element(i2, j2)+
							rsemiproduct(parfunc[0]->getQ(i), parfunc[1]->getQ(j), QIa_nn[i][j], i1, j1, i2, j2);

	//QI_nd
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 2; j++)
			if(j == 1 || (i != 0 && j == 0))
				*QI_nd[i][j]->estar(i1, j1, i2, j2) = rsemiproduct(parfunc[0]->getQ(i), parfunc[1]->getQ(j), QIa_nd[i][j], i1, j1, i2, j2);


	//QI_dn
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 3; j++)
			if(i == 1 || (j != 0 && i == 0))
				*QI_dn[i][j]->estar(i1, j1, i2, j2) = rsemiproduct(parfunc[0]->getQ(i), parfunc[1]->getQ(j), QIa_dn[i][j], i1, j1, i2, j2);


	//QI_dd
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			if(i+j)
				*QI_dd[i][j]->estar(i1, j1, i2, j2) = rsemiproduct(parfunc[0]->getQ(i), parfunc[1]->getQ(j), QIa_dd[i][j], i1, j1, i2, j2)
									+ rsemiproduct(parfunc[0]->getQ(i), parfunc[1]->getQ(j), QIb[i][j], i1, j1, i2, j2);
}

void PartitionFunction::computeI(int i1, int j1, int i2, int j2)
{
	//QIa
	double *qia = QIa->estar(i1, j1, i2, j2);

	*qia = product(QIx[0][0], QI, i1, j1, i2, j2) + product(QIx[1][0], QI, i1, j1, i2, j2) + product(QIy, QI, i1, j1, i2, j2);
				
	*qia += lsemiproduct(QIx[0][0], parfunc[0]->getQ(0), parfunc[1]->getQ(0),i1, j1, i2, j2) + 
		lsemiproduct(QIx[1][0], parfunc[0]->getQ(0), parfunc[1]->getQ(0),i1, j1, i2, j2) +
		lsemiproduct(QIy, parfunc[0]->getQ(0), parfunc[1]->getQ(0),i1, j1, i2, j2);		

	//QIb
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
		{
			*QIb[i][j]->estar(i1, j1, i2, j2) = QIhh[i][j]->element(i1, j1, i2, j2);
			*QIb[i][j]->estar(i1, j1, i2, j2) += product(QIhm[i][j], QIb[i][j], i1, j1, i2, j2); 
			*QIb[i][j]->estar(i1, j1, i2, j2) += (i+j) ? product(QIhh[i][j], auxQIa[i][j], i1, j1, i2, j2) : product(QIhh[i][j], QIa, i1, j1, i2, j2);
		}

	//QI
	*QI->estar(i1, j1, i2, j2) = parfunc[0]->getQ(0)->element(i1, j1)*parfunc[1]->getQ(0)->element(i2, j2)+
							rsemiproduct(parfunc[0]->getQ(0), parfunc[1]->getQ(0), QIa, i1, j1, i2, j2) +
							rsemiproduct(parfunc[0]->getQ(0), parfunc[1]->getQ(0), QIb[0][0], i1, j1, i2, j2);
}

void PartitionFunction::store(FILE *outp)
{
	if(parFuncType == SINGLE)
	{
		if(!quiet) printf("Storing single partition function of %s... \n", sequence->getName());
		for(int j = 0; j < 3; j++)
		{
			Q[j]->store(outp);
			Qm[j]->store(outp);
		}
		Qb->store(outp);

		if(!quiet) printf("Tables stored. \n");
	} else
	{
		if(!quiet) printf("Storing double partition function of %s and %s... \n", parfunc[0]->getSequence()->getName(), parfunc[1]->getSequence()->getName());

		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				if(i+j) QIa_nn[i][j]->store(outp);
				if(i+j) QI_nn[i][j]->store(outp);
			}
		}

		for(int i = 0; i < 2; i++)
		{
			QIh[i]->store(outp);
			QIx_m[i][2]->store(outp);
			QImm[i][2]->store(outp);
			QImm[2][i]->store(outp);
			QImk[i][1]->store(outp);
			QImk[2][i]->store(outp);
			QIkm[1][i]->store(outp);
			QIkm[i][2]->store(outp);
			QIa_nd[i][1]->store(outp);
			QIa_nd[2][i]->store(outp);
			QI_nd[i][1]->store(outp);
			QI_nd[2][i]->store(outp);
			QIa_dn[1][i]->store(outp);
			QIa_dn[i][2]->store(outp);
			QI_dn[1][i]->store(outp);
			QI_dn[i][2]->store(outp);
			for(int j = 0; j < 2; j++)
			{
				QIb[i][j]->store(outp);
				QIhh[i][j]->store(outp);
				QIhm[i][j]->store(outp);
				QIx_m[i][j]->store(outp);
				QIx_k[i][j]->store(outp);
				if(i+j) QIkk[i][j]->store(outp);
				if(i+j) QIa_r[i][j]->store(outp);
				if(i+j) QI_r[i][j]->store(outp);
				if(i+j) QIlr[i][j]->store(outp);
				if(i+j) QIa_dd[i][j]->store(outp);
				if(i+j) QI_dd[i][j]->store(outp);
			}
		}
		QImm[2][2]->store(outp);
		QIa_nd[1][0]->store(outp);
		QI_nd[1][0]->store(outp);
		QIa_dn[0][1]->store(outp);
		QI_dn[0][1]->store(outp);
		QIy->store(outp);
		QI->store(outp);
		QIa->store(outp);
		if(!quiet) printf("Tables stored. \n");
	}
}

void PartitionFunction::retrieve(FILE *inp)
{
	if(parFuncType == SINGLE)
	{
		if(!quiet) printf("Retrieving single partition function of %s... \n", sequence->getName());
		for(int j = 0; j < 3; j++)
		{
			Q[j]->retrieve(inp);
			Qm[j]->retrieve(inp);
		}
		Qb->retrieve(inp);
		Zz = Q[0]->element(0, sequence->getLen() - 1);
		if(!quiet) printf("Tables retrieved. \n");
	} else
	{
		if(!quiet) printf("Retrieving double partition function of %s and %s... \n", parfunc[0]->getSequence()->getName(), parfunc[1]->getSequence()->getName());

		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				if(i+j) QIa_nn[i][j]->retrieve(inp);
				if(i+j) QI_nn[i][j]->retrieve(inp);
			}
		}

		for(int i = 0; i < 2; i++)
		{
			QIh[i]->retrieve(inp);
			QIx_m[i][2]->retrieve(inp);
			QImm[i][2]->retrieve(inp);
			QImm[2][i]->retrieve(inp);
			QImk[i][1]->retrieve(inp);
			QImk[2][i]->retrieve(inp);
			QIkm[1][i]->retrieve(inp);
			QIkm[i][2]->retrieve(inp);
			QIa_nd[i][1]->retrieve(inp);
			QIa_nd[2][i]->retrieve(inp);
			QI_nd[i][1]->retrieve(inp);
			QI_nd[2][i]->retrieve(inp);
			QIa_dn[1][i]->retrieve(inp);
			QIa_dn[i][2]->retrieve(inp);
			QI_dn[1][i]->retrieve(inp);
			QI_dn[i][2]->retrieve(inp);
			for(int j = 0; j < 2; j++)
			{
				QIb[i][j]->retrieve(inp);
				QIhh[i][j]->retrieve(inp);
				QIhm[i][j]->retrieve(inp);
				QIx_m[i][j]->retrieve(inp);
				QIx_k[i][j]->retrieve(inp);
				if(i+j) QIkk[i][j]->retrieve(inp);
				if(i+j) QIa_r[i][j]->retrieve(inp);
				if(i+j) QI_r[i][j]->retrieve(inp);
				if(i+j) QIlr[i][j]->retrieve(inp);
				if(i+j) QIa_dd[i][j]->retrieve(inp);
				if(i+j) QI_dd[i][j]->retrieve(inp);
			}
		}
		QImm[2][2]->retrieve(inp);
		QIa_nd[1][0]->retrieve(inp);
		QI_nd[1][0]->retrieve(inp);
		QIa_dn[0][1]->retrieve(inp);
		QI_dn[0][1]->retrieve(inp);
		QIy->retrieve(inp);
		QI->retrieve(inp);
		QIa->retrieve(inp);
		Zz = QI->element(0, parfunc[0]->getSequence()->getLen() - 1, 0, parfunc[1]->getSequence()->getLen() - 1);
		if(!quiet) printf("Tables retrieved. \n");
	}
}

// this is a test method that verifies interstrand symmetries provided that the two input sequences are identical.
#define ERR 1e20

void PartitionFunction::selftest(int i1, int j1, int i2, int j2, int l1, int l2)
{
	for(int c1 = 0; c1 < 3; c1++)
	{
		if(fabs(parfunc[0]->getQ(c1)->element(i1, j1) - parfunc[1]->getQ(c1)->element( l1 - j1 - 1, l1 - i1 - 1)) > ERR)
		{
			printf("Q[%d] fails self test: %d\t%d\n", c1, i1, j1);
			exit(1);
		}
		if(fabs(parfunc[0]->getQm(c1)->element(i1, j1) - parfunc[1]->getQm(c1)->element( l1 - j1 - 1, l1 - i1 - 1)) > ERR)
		{
			printf("Qm[%d] fails self test: %d\t%d\n", c1, i1, j1);
			exit(1);
		}
	}

	if(fabs(QI->element(i1, j1, i2, j2) - QI->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
	{
		printf("QI fails self test: %d\t%d\t%d\t%d\n", i1, j1, i2, j2);
		printf("%g\t%g\n", QI->element(i1, j1, i2, j2), QI->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1));
		exit(1);
	}
	if(fabs(QIy->element(i1, j1, i2, j2) - QIy->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
	{
		printf("QIy fails self test: %d\t%d\t%d\t%d\n", i1, j1, i2, j2);
		exit(1);
	}

	for(int c1 = 0; c1 < 2; c1++)
		for(int c2 = 0; c2 < 2; c2++)
			if(c1+c2)
			{
				if(fabs(QIb[c1][c2]->element(i1, j1, i2, j2) -  QI_r[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("QIb[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
				if(fabs(QIkk[c1][c2]->element(i1, j1, i2, j2) - QIkk[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("QIkk[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
				if(fabs(QI_dd[c1][c2]->element(i1, j1, i2, j2) -  QI_dd[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("QI_dd[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					printf("%g\t%g\n", QI_dd[c1][c2]->element(i1, j1, i2, j2), QI_dd[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1));
					exit(1);
				}
				if(fabs(QIlr[c1][c2]->element(i1, j1, i2, j2) -  QIlr[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("QIlr[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
				if(fabs(QI_r[c1][c2]->element(i1, j1, i2, j2) - QIb[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("QI_r[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
			}

	for(int c1 = 0; c1 < 3; c1++)
		for(int c2 = 0; c2 < 3; c2++)
		{
			if(c1 == 2 || c2 == 2)			
				if(fabs(QImm[c1][c2]->element(i1, j1, i2, j2) - QImm[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("QImm[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
			if(c1+c2)
				if(fabs(QI_nn[c1][c2]->element(i1, j1, i2, j2) - QI_nn[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("QI_nn[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
		}	
	for(int c1 = 0; c1 < 2; c1++)
		if(fabs(QIh[c1]->element(i1, j1, i2, j2) - QIh[c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
		{
			printf("QIh[%d] fails self test: %d\t%d\t%d\t%d\n", c1, i1, j1, i2, j2);
			exit(1);
		}
	if((i1 == 6 && j1 == 23 && i2 == 10 && j2 == 24)
		|| (i1 == 10 && j1 == 24 && i2 == 11 && j2 == 28))
	{
		printf("%d: %g\t%g\t%g\n", i1, parfunc[0]->getQ(0)->element(i1, j1)*parfunc[1]->getQ(0)->element(i2, j2),
			rsemiproduct(parfunc[0]->getQ(0), parfunc[1]->getQ(0), QIa, i1, j1, i2, j2),
			rsemiproduct(parfunc[0]->getQ(0), parfunc[1]->getQ(0), QIb[0][0], i1, j1, i2, j2));
		printf("Up: ");
		for(int tmp1 = i1; tmp1 <= j1; tmp1++)		
			printf("%d ", seqUp[tmp1]);
		printf("\nDown: ");
		for(int tmp1 = i2; tmp1 <= j2; tmp1++)		
			printf("%d ", seqDown[tmp1]);
		printf("\n");
		printf("QIa: %g\n", QIa->element(i1, j1, i2, j2));
		printf("%g\t%g\t%g\t%g\t%g\t%g\n", product(QIx[0][0], QI, i1, j1, i2, j2), 
			product(QIx[1][0], QI, i1, j1, i2, j2), 
			product(QIy, QI, i1, j1, i2, j2),
			lsemiproduct(QIx[0][0], parfunc[0]->getQ(0), parfunc[1]->getQ(0),i1, j1, i2, j2), 
			lsemiproduct(QIx[1][0], parfunc[0]->getQ(0), parfunc[1]->getQ(0),i1, j1, i2, j2),
			lsemiproduct(QIy, parfunc[0]->getQ(0), parfunc[1]->getQ(0),i1, j1, i2, j2));
		printf("QIx_k: %g\t%g\n", QIx_k[0][0]->element(i1, j1, i2, j2), QIx_k[1][0]->element(i1, j1, i2, j2));
	}
}
