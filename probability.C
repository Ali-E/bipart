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
	Desc: Probability class computes the base pair probabilities and also motif probabilities. 
It uses the Probability class to compute Q tables. 

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow, 
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: June 5, 2009
	Ver: 1.8
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
#include "probability.h"
#include "partitionfunction.h"
#include "getopt.h"
#include "sequence.h"


int inline max(int a, int b)
{
	int ret = (a > b) ? a : b;
	return ret;
}

Probability::Probability(Energy *en, PartitionFunction *pfunc, bool _debug, int _procNum, bool _quiet, bool _no_isolate)
{
	double temp;
	energy = en;
	parfunc = pfunc;
	energy->getParams(&NA, &polymer, &naConc, &mgConc, &suffix, &temp, &zip, &nodangle, &energyType);
	RT = energy->getRT();

	DEBUG = _debug;
	procNum = _procNum;
	quiet = _quiet;
	no_isolate = _no_isolate;
	finished = false;
	aProb[0] = aProb[1] = bProb = NULL;

	long int inittotal_tables_size = total_tables_size;

	if(parfunc->getType() == DUPLEX)
	{
		buildTables();
		if(!quiet) printf("Allocation of %ld MB of memory and its initialization finished. \n", (total_tables_size - inittotal_tables_size) / (1024*1024));
		computeProb(parfunc->getParFunc()[0]->getSequence()->getLen(), (time_t)0);
		finished = true;
	}
}

Probability::Probability(FILE *inp, FILE *outp, Energy *en, PartitionFunction *pfunc, time_t e_time, bool _debug, int _procNum, bool _quiet, bool _no_isolate)
{
	double temp;
	energy = en;
	parfunc = pfunc;
	energy->getParams(&NA, &polymer, &naConc, &mgConc, &suffix, &temp, &zip, &nodangle, &energyType);
	RT = energy->getRT();

	DEBUG = _debug;
	procNum = _procNum;
	quiet = _quiet;
	no_isolate = _no_isolate;
	finished = false;
	aProb[0] = aProb[1] = bProb = NULL;

	long int inittotal_tables_size = total_tables_size;
	unsigned int itr;

	if(parfunc->getType() == DUPLEX)
	{
		if(inp)
		{
			buildEmptyTables();
			if(fread(&itr, sizeof(unsigned int), 1, inp) != 1)
				perror("Reading iteration\n");
			retrieve(inp);
		}
		else
		{
			buildTables();
			itr = parfunc->getParFunc()[0]->getSequence()->getLen();
		}

		if(!quiet) printf("Allocation of %ld MB of memory and its initialization finished. \n", (total_tables_size - inittotal_tables_size) / (1024*1024));
		itr = computeProb(itr, e_time);
		if(outp)
		{
			if(fwrite(&itr, sizeof(unsigned int), 1, outp) != 1)
				perror("Writing iteration\n");
			store(outp);
		}
		finished = (itr == 0);
	}
}

void Probability::buildEmptyTables()
{
	if(parfunc->getType() == DUPLEX)
	{
		unsigned int l1 = parfunc->getParFunc()[0]->getSequence()->getLen();
		unsigned int l2 = parfunc->getParFunc()[1]->getSequence()->getLen();
		unsigned int lens[2] = {l1, l2};

		char d2 = 2, d4 = 4;
		for(int i = 0; i < 2; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				Q[i][j] = parfunc->getParFunc()[i]->getQ(j);
				Qm[i][j] = parfunc->getParFunc()[i]->getQm(j);
				P[i][j] = new Table<double>();
				Pm[i][j] = new Table<double>();
			}
			Qb[i] = parfunc->getParFunc()[i]->getQb();
			Pb[i] = new Table<double>();

			aProb[i] = new Table<double>(lens[i], &d2, 0.0);
		}

		bProb = new Table<double>(l1, l2, 1, 1, 0.0);

		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				if(i+j) QIa_nn[i][j] = parfunc->getQIa_nn(i, j);
				if(i+j) PIa_nn[i][j] = new Table<double>();
				if(i+j) QI_nn[i][j] = parfunc->getQI_nn(i, j);
				if(i+j) PI_nn[i][j] = new Table<double>();
			}
		}

		for(int i = 0; i < 2; i++)
		{
			QIh[i] = parfunc->getQIh(i);
			PIh[i] = new Table<double>();
			QIx_m[i][2] = parfunc->getQIx_m(i, 2);
			PIx_m[i][2] = new Table<double>();
			QImm[i][2] = parfunc->getQImm(i, 2);
			PImm[i][2] = new Table<double>();
			QImm[2][i] = parfunc->getQImm(2, i);
			PImm[2][i] = new Table<double>();
			QImk[i][1] = parfunc->getQImk(i, 1);
			PImk[i][1] = new Table<double>();
			QImk[2][i] = parfunc->getQImk(2, i);
			PImk[2][i] = new Table<double>();
			QIkm[1][i] = parfunc->getQIkm(1, i);
			PIkm[1][i] = new Table<double>();
			QIkm[i][2] = parfunc->getQIkm(i, 2);
			PIkm[i][2] = new Table<double>();
			QIa_nd[i][1] = parfunc->getQIa_nd(i, 1);
			PIa_nd[i][1] = new Table<double>();
			QIa_nd[2][i] = parfunc->getQIa_nd(2, i);
			PIa_nd[2][i] = new Table<double>();
			QI_nd[i][1] = parfunc->getQI_nd(i, 1);
			PI_nd[i][1] = new Table<double>();
			QI_nd[2][i] = parfunc->getQI_nd(2, i);
			PI_nd[2][i] = new Table<double>();
			QIa_dn[1][i] = parfunc->getQIa_dn(1, i);
			PIa_dn[1][i] = new Table<double>();
			QIa_dn[i][2] = parfunc->getQIa_dn(i, 2);
			PIa_dn[i][2] = new Table<double>();
			QI_dn[1][i] = parfunc->getQI_dn(1, i);
			PI_dn[1][i] = new Table<double>();
			QI_dn[i][2] = parfunc->getQI_dn(i, 2);
			PI_dn[i][2] = new Table<double>();
			for(int j = 0; j < 2; j++)
			{
				QIb[i][j] = parfunc->getQIb(i, j);
				PIb[i][j] = new Table<double>();
				QIhh[i][j] = parfunc->getQIhh(i, j);
				PIhh[i][j] = new Table<double>();
				QIhm[i][j] = parfunc->getQIhm(i, j);
				PIhm[i][j] = new Table<double>();
				QIx_m[i][j] = parfunc->getQIx_m(i, j);
				PIx_m[i][j] = new Table<double>();
				QIx_k[i][j] = parfunc->getQIx_k(i, j);
				PIx_k[i][j] = new Table<double>();
				if(i+j) QIkk[i][j] = parfunc->getQIkk(i, j);
				if(i+j) PIkk[i][j] = new Table<double>();
				if(i+j) QIa_r[i][j] = parfunc->getQIa_r(i, j);
				if(i+j) PIa_r[i][j] = new Table<double>();
				if(i+j) QI_r[i][j] = parfunc->getQI_r(i, j);
				if(i+j) PI_r[i][j] = new Table<double>();
				if(i+j) QIlr[i][j] = parfunc->getQIlr(i, j);
				if(i+j) PIlr[i][j] = new Table<double>();
				if(i+j) QIa_dd[i][j] = parfunc->getQIa_dd(i, j);
				if(i+j) PIa_dd[i][j] = new Table<double>();
				if(i+j) QI_dd[i][j] = parfunc->getQI_dd(i, j);
				if(i+j) PI_dd[i][j] = new Table<double>();
			}
		}
		QImm[2][2] = parfunc->getQImm(2, 2);
		PImm[2][2] = new Table<double>();
		QIa_nd[1][0] = parfunc->getQIa_nd(1, 0);
		PIa_nd[1][0] = new Table<double>();
		QI_nd[1][0] = parfunc->getQI_nd(1, 0);
		PI_nd[1][0] = new Table<double>();
		QIa_dn[0][1] = parfunc->getQIa_dn(0, 1);
		PIa_dn[0][1] = new Table<double>();
		QI_dn[0][1] = parfunc->getQI_dn(0, 1);
		PI_dn[0][1] = new Table<double>();
		QIy = parfunc->getQIy();
		PIy = new Table<double>();
		QI = parfunc->getQI();
		PI = new Table<double>();
		QIa = parfunc->getQIa();
		PIa = new Table<double>();

		for(int i = 0; i < 2; i++)
			for(int j = 0; j < 2; j++)
			{
				auxPIa[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*5);
				*auxPIa[i][j] = PIa_nn[i][j];
				*(auxPIa[i][j] + 1) = PIa_nd[i][j];
				*(auxPIa[i][j] + 2) = PIa_dn[i][j];
				*(auxPIa[i][j] + 3) = PIa_dd[i][j];
				*(auxPIa[i][j] + 4) = NULL;
				auxQIa[i][j] = parfunc->getauxQIa(i, j);
				auxPIa_dd[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*3);
				*auxPIa_dd[i][j] = PIa_dd[i][j];
				*(auxPIa_dd[i][j] + 1) = PIb[i][j];
				*(auxPIa_dd[i][j] + 2) = NULL;
				auxQIa_dd[i][j] = parfunc->getauxQIa_dd(i, j);
				PIx[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*3);
				*PIx[i][j] = PIx_m[i][j];
				*(PIx[i][j] + 1) = PIx_k[i][j];
				*(PIx[i][j] + 2) = NULL;
				QIx[i][j] = parfunc->getQIx(i, j);
			}
	}
}

void Probability::buildTables()
{
	if(parfunc->getType() == DUPLEX)
	{
		unsigned int l1 = parfunc->getParFunc()[0]->getSequence()->getLen();
		unsigned int l2 = parfunc->getParFunc()[1]->getSequence()->getLen();
		unsigned int lens[2] = {l1, l2};

		char d2 = 2, d4 = 4;
		for(int i = 0; i < 2; i++)
		{
			if(!quiet) printf("Allocating memory for single base pairing probability of %s... \n", parfunc->getParFunc()[i]->getSequence()->getName());

			for(int j = 0; j < 3; j++)
			{
				Q[i][j] = parfunc->getParFunc()[i]->getQ(j);
				Qm[i][j] = parfunc->getParFunc()[i]->getQm(j);
				P[i][j] = new Table<double>(lens[i], &d2, 0.0);
				Pm[i][j] = new Table<double>(lens[i], &d2, 0.0);
			}
			Qb[i] = parfunc->getParFunc()[i]->getQb();
			Pb[i] = new Table<double>(lens[i], &d2, 0.0);

			aProb[i] = new Table<double>(lens[i], &d2, 0.0);
	
			if(!quiet) printf("Allocation and initialization finished. \n");
		}

		if(!quiet) printf("Allocating memory for double base pairing probability of %s and %s... \n", parfunc->getParFunc()[0]->getSequence()->getName(), parfunc->getParFunc()[1]->getSequence()->getName());

		bProb = new Table<double>(l1, l2, 1, 1, 0.0);

		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				if(i+j) QIa_nn[i][j] = parfunc->getQIa_nn(i, j);
				if(i+j) PIa_nn[i][j] = new Table<double>(l1, l2, 0.0);
				if(i+j) QI_nn[i][j] = parfunc->getQI_nn(i, j);
				if(i+j) PI_nn[i][j] = new Table<double>(l1, l2, 0.0);
			}
		}

		for(int i = 0; i < 2; i++)
		{
			QIh[i] = parfunc->getQIh(i);
			PIh[i] = new Table<double>(l1, l2, 0.0);
			QIx_m[i][2] = parfunc->getQIx_m(i, 2);
			PIx_m[i][2] = new Table<double>(l1, l2, 0.0);
			QImm[i][2] = parfunc->getQImm(i, 2);
			PImm[i][2] = new Table<double>(l1, l2, 0.0);
			QImm[2][i] = parfunc->getQImm(2, i);
			PImm[2][i] = new Table<double>(l1, l2, 0.0);
			QImk[i][1] = parfunc->getQImk(i, 1);
			PImk[i][1] = new Table<double>(l1, l2, 0.0);
			QImk[2][i] = parfunc->getQImk(2, i);
			PImk[2][i] = new Table<double>(l1, l2, 0.0);
			QIkm[1][i] = parfunc->getQIkm(1, i);
			PIkm[1][i] = new Table<double>(l1, l2, 0.0);
			QIkm[i][2] = parfunc->getQIkm(i, 2);
			PIkm[i][2] = new Table<double>(l1, l2, 0.0);
			QIa_nd[i][1] = parfunc->getQIa_nd(i, 1);
			PIa_nd[i][1] = new Table<double>(l1, l2, 0.0);
			QIa_nd[2][i] = parfunc->getQIa_nd(2, i);
			PIa_nd[2][i] = new Table<double>(l1, l2, 0.0);
			QI_nd[i][1] = parfunc->getQI_nd(i, 1);
			PI_nd[i][1] = new Table<double>(l1, l2, 0.0);
			QI_nd[2][i] = parfunc->getQI_nd(2, i);
			PI_nd[2][i] = new Table<double>(l1, l2, 0.0);
			QIa_dn[1][i] = parfunc->getQIa_dn(1, i);
			PIa_dn[1][i] = new Table<double>(l1, l2, 0.0);
			QIa_dn[i][2] = parfunc->getQIa_dn(i, 2);
			PIa_dn[i][2] = new Table<double>(l1, l2, 0.0);
			QI_dn[1][i] = parfunc->getQI_dn(1, i);
			PI_dn[1][i] = new Table<double>(l1, l2, 0.0);
			QI_dn[i][2] = parfunc->getQI_dn(i, 2);
			PI_dn[i][2] = new Table<double>(l1, l2, 0.0);
			for(int j = 0; j < 2; j++)
			{
				QIb[i][j] = parfunc->getQIb(i, j);
				PIb[i][j] = new Table<double>(l1, l2, 0.0);
				QIhh[i][j] = parfunc->getQIhh(i, j);
				PIhh[i][j] = new Table<double>(l1, l2, 0.0);
				QIhm[i][j] = parfunc->getQIhm(i, j);
				PIhm[i][j] = new Table<double>(l1, l2, 0.0);
				QIx_m[i][j] = parfunc->getQIx_m(i, j);
				PIx_m[i][j] = new Table<double>(l1, l2, 0.0);
				QIx_k[i][j] = parfunc->getQIx_k(i, j);
				PIx_k[i][j] = new Table<double>(l1, l2, 0.0);
				if(i+j) QIkk[i][j] = parfunc->getQIkk(i, j);
				if(i+j) PIkk[i][j] = new Table<double>(l1, l2, 0.0);
				if(i+j) QIa_r[i][j] = parfunc->getQIa_r(i, j);
				if(i+j) PIa_r[i][j] = new Table<double>(l1, l2, 0.0);
				if(i+j) QI_r[i][j] = parfunc->getQI_r(i, j);
				if(i+j) PI_r[i][j] = new Table<double>(l1, l2, 0.0);
				if(i+j) QIlr[i][j] = parfunc->getQIlr(i, j);
				if(i+j) PIlr[i][j] = new Table<double>(l1, l2, 0.0);
				if(i+j) QIa_dd[i][j] = parfunc->getQIa_dd(i, j);
				if(i+j) PIa_dd[i][j] = new Table<double>(l1, l2, 0.0);
				if(i+j) QI_dd[i][j] = parfunc->getQI_dd(i, j);
				if(i+j) PI_dd[i][j] = new Table<double>(l1, l2, 0.0);
			}
		}
		QImm[2][2] = parfunc->getQImm(2, 2);
		PImm[2][2] = new Table<double>(l1, l2, 0.0);
		QIa_nd[1][0] = parfunc->getQIa_nd(1, 0);
		PIa_nd[1][0] = new Table<double>(l1, l2, 0.0);
		QI_nd[1][0] = parfunc->getQI_nd(1, 0);
		PI_nd[1][0] = new Table<double>(l1, l2, 0.0);
		QIa_dn[0][1] = parfunc->getQIa_dn(0, 1);
		PIa_dn[0][1] = new Table<double>(l1, l2, 0.0);
		QI_dn[0][1] = parfunc->getQI_dn(0, 1);
		PI_dn[0][1] = new Table<double>(l1, l2, 0.0);
		QIy = parfunc->getQIy();
		PIy = new Table<double>(l1, l2, 0.0);
		QI = parfunc->getQI();
		PI = new Table<double>(l1, l2, 0.0);
		QIa = parfunc->getQIa();
		PIa = new Table<double>(l1, l2, 0.0);

		for(int i = 0; i < 2; i++)
			for(int j = 0; j < 2; j++)
			{
				auxPIa[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*5);
				*auxPIa[i][j] = PIa_nn[i][j];
				*(auxPIa[i][j] + 1) = PIa_nd[i][j];
				*(auxPIa[i][j] + 2) = PIa_dn[i][j];
				*(auxPIa[i][j] + 3) = PIa_dd[i][j];
				*(auxPIa[i][j] + 4) = NULL;
				auxQIa[i][j] = parfunc->getauxQIa(i, j);
				auxPIa_dd[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*3);
				*auxPIa_dd[i][j] = PIa_dd[i][j];
				*(auxPIa_dd[i][j] + 1) = PIb[i][j];
				*(auxPIa_dd[i][j] + 2) = NULL;
				auxQIa_dd[i][j] = parfunc->getauxQIa_dd(i, j);
				PIx[i][j] = (Table<double> **)alloc.xmalloc(sizeof(Table<double> *)*3);
				*PIx[i][j] = PIx_m[i][j];
				*(PIx[i][j] + 1) = PIx_k[i][j];
				*(PIx[i][j] + 2) = NULL;
				QIx[i][j] = parfunc->getQIx(i, j);
			}
	}
}

Probability::~Probability()
{
	if(parfunc->getType() == DUPLEX)
	{
		delete bProb;

		for(int i = 0; i < 2; i++)
		{		
			for(int j = 0; j < 3; j++)
			{
				delete P[i][j];
				delete Pm[i][j];
			}
			delete Pb[i];

			delete aProb[i];
		}

		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				if(i+j) delete PIa_nn[i][j];
				if(i+j) delete PI_nn[i][j];
			}
		}
		delete PImm[2][2];
		delete PIa_nd[1][0];
		delete PI_nd[1][0];
		delete PIa_dn[0][1];
		delete PI_dn[0][1];
		delete PIy;
		delete PI;
		delete PIa;
		for(int i = 0; i < 2; i++)
		{
			delete PIkm[1][i];
			delete PIkm[i][2];
			delete PImk[i][1];
			delete PImk[2][i];
			delete PIa_nd[i][1];
			delete PIa_nd[2][i];
			delete PI_nd[i][1];
			delete PI_nd[2][i];
			delete PIa_dn[1][i];
			delete PIa_dn[i][2];
			delete PI_dn[1][i];
			delete PI_dn[i][2];
			delete PIh[i];
			delete PIx_m[i][2];
			delete PImm[i][2];
			delete PImm[2][i];
			for(int j = 0; j < 2; j++)
			{
				free(auxPIa[i][j]);
				free(auxPIa_dd[i][j]);
				free(PIx[i][j]);
				delete PIb[i][j];
				delete PIhh[i][j];
				delete PIhm[i][j];
				if(i+j) delete PIkk[i][j];
				if(i+j) delete PIa_r[i][j];
				if(i+j) delete PI_r[i][j];
				if(i+j) delete PIlr[i][j];
				if(i+j) delete PIa_dd[i][j];
				if(i+j) delete PI_dd[i][j];
				delete PIx_m[i][j];
				delete PIx_k[i][j];
			}
		}
	}
}

void Probability::processPbPm(int i, int j, int upOrDown)
{
	unsigned char *seq = (upOrDown == 0) ? seqUp : seqDown;
			
	double qmw = Qm[upOrDown][0]->element(i, j);
	double qmr = Qm[upOrDown][1]->element(i, j);
	double qmg = Qm[upOrDown][2]->element(i, j);

	double pmw = Pm[upOrDown][0]->element(i, j);
	double pmr = Pm[upOrDown][1]->element(i, j);
	double pmg = Pm[upOrDown][2]->element(i, j);

	if(pmw == 0.0 || qmw == 0.0)
	{
		pmw = 0.0;
		qmw = 1.0;
	}

	if(pmr == 0.0 || qmr == 0.0)
	{
		pmr = 0.0;
		qmr = 1.0;
	}

	if(pmg == 0.0 || qmg == 0.0)
	{
		pmg = 0.0;
		qmg = 1.0;
	}

	for(int d = i; d <= j - MIN_HAIRPIN_SIZE - 1; d++)
		for(int e = d + MIN_HAIRPIN_SIZE + 1; e <= j; e++)
		{
			if (d > i)
			{ 
				double deltaPw = pmw * Qm[upOrDown][0]->element(i, d-1)*Qb[upOrDown]->element(d, e) / qmw; 
			  	double deltaPr = pmr * Qm[upOrDown][1]->element(i, d-1)*Qb[upOrDown]->element(d, e)*exp(-energy->Ekissing(0, j-e, 1)/RT) / qmr; 
   				double deltaPg = pmg * Qm[upOrDown][2]->element(i, d-1)*Qb[upOrDown]->element(d, e)*exp(-energy->Emulti(0, j-e, 1)/RT) / qmg;

#pragma omp atomic
				*Pb[upOrDown]->estar(d, e) += deltaPw + deltaPr + deltaPg;
#pragma omp atomic
				*Pm[upOrDown][0]->estar(i, d-1) += deltaPw;
#pragma omp atomic
				*Pm[upOrDown][1]->estar(i, d-1) += deltaPr;
#pragma omp atomic
				*Pm[upOrDown][2]->estar(i, d-1) += deltaPg;
			}					
#pragma omp atomic
			*Pb[upOrDown]->estar(d, e) += pmw * Qb[upOrDown]->element(d, e) / qmw; 
#pragma omp atomic
			*Pb[upOrDown]->estar(d, e) += pmr * Qb[upOrDown]->element(d, e)*exp(-energy->Ekissing(0, j-i-(e-d), 1)/RT) / qmr; 
#pragma omp atomic
			*Pb[upOrDown]->estar(d, e) += pmg * Qb[upOrDown]->element(d, e)*exp(-energy->Emulti(0, j-i-(e-d), 1)/RT) / qmg; 
		}


	if(j > i + MIN_HAIRPIN_SIZE && energy->base_pair(seq[i], seq[j]) && Pb[upOrDown]->element(i, j) != 0.0 && Qb[upOrDown]->element(i, j) != 0.0)
	{
		for(int d = i+1; d <= j-MIN_HAIRPIN_SIZE-2; d++)
			for(int e = d+MIN_HAIRPIN_SIZE+1; e <= j-1; e++)
			{
				if (DEBUG) printf("Pb - b [%d,%d]\n", d, e);
				if (e == j-1 && d == i+1)
#pragma omp atomic
					*Pb[upOrDown]->estar(d, e) += Pb[upOrDown]->element(i, j) * Qb[upOrDown]->element(d, e)*exp(-energy->Es(i, j, seq)/RT) / Qb[upOrDown]->element(i, j);
				else 
				{
#pragma omp atomic
					*Pb[upOrDown]->estar(d, e) += Pb[upOrDown]->element(i, j) * Qb[upOrDown]->element(d, e)*exp(-energy->Ebi(i, j, d, e, seq)/RT) / Qb[upOrDown]->element(i, j);
					if(d > i+1)
					{
						double deltaP = Pb[upOrDown]->element(i, j) * Qb[upOrDown]->element(d, e)*Qm[upOrDown][2]->element(i + 1, d - 1)*exp(-energy->Emulti(1, j-e-1, 1)/RT) / Qb[upOrDown]->element(i, j);
#pragma omp atomic
						*Pb[upOrDown]->estar(d, e) += deltaP; 
#pragma omp atomic
						*Pm[upOrDown][2]->estar(i + 1, d - 1) += deltaP; 
					}
				}
			}
	}
}

void Probability::processP(int i, int j, int upOrDown)
{
	for(int d = i; d <= j-MIN_HAIRPIN_SIZE-1; d++)
		for(int e = d+MIN_HAIRPIN_SIZE+1; e <= j; e++)
		{
			double temp1 = 1.0, temp2 = 1.0, temp3 = 1.0;
			if (d > i)
			{
				temp1 = Q[upOrDown][0]->element(i, d-1);
				temp2 = Q[upOrDown][1]->element(i, d-1);
				temp3 = Q[upOrDown][2]->element(i, d-1);
			}

			double deltaP0 = 0.0, deltaP1 = 0.0, deltaP2 = 0.0;

			if(Q[upOrDown][0]->element(i, j) != 0.0)
				deltaP0 = P[upOrDown][0]->element(i, j) * temp1*Qb[upOrDown]->element(d, e) / Q[upOrDown][0]->element(i, j);

			if(Q[upOrDown][1]->element(i, j) != 0.0)
				deltaP1 = P[upOrDown][1]->element(i, j) * temp2*Qb[upOrDown]->element(d, e)*exp(-energy->Ekissing(0, j-e, 1)/RT) / Q[upOrDown][1]->element(i, j);

			if(Q[upOrDown][2]->element(i, j) != 0.0)
				deltaP2 = P[upOrDown][2]->element(i, j) * temp3*Qb[upOrDown]->element(d, e)*exp(-energy->Emulti(0, j-e, 1)/RT) / Q[upOrDown][2]->element(i, j);

#pragma omp atomic
			*Pb[upOrDown]->estar(d, e) += deltaP0 + deltaP1 + deltaP2;
			
			if (d > i)
			{
#pragma omp atomic
				*P[upOrDown][0]->estar(i, d-1) += deltaP0;
#pragma omp atomic
				*P[upOrDown][1]->estar(i, d-1) += deltaP1;
#pragma omp atomic
				*P[upOrDown][2]->estar(i, d-1) += deltaP2;
			}
		}
}

/* Hamid: The two RNAs interact 3'->5' and 5'->3'. But I assume that one RNA is reversed beforehand so the indices now interact in the same direction. */

void Probability::Reverse(int upOrDown)
{
	for(int j=0; j < 3; j++)
	{
		P[upOrDown][j]->Reverse();
		Pm[upOrDown][j]->Reverse(); 
	}
	Pb[upOrDown]->Reverse();
	aProb[upOrDown]->Reverse();
}


unsigned int Probability::computeProb(unsigned int itr, time_t end_time)
{
	Sequence *sequenceUp = parfunc->getParFunc()[0]->getSequence();
	Sequence *sequenceDown = parfunc->getParFunc()[1]->getSequence();

	seqUp = sequenceUp->getSeq();	
	seqDown = sequenceDown->getReverse();

	unsigned int lenUp = sequenceUp->getLen();
	unsigned int lenDown = sequenceDown->getLen();
	parfunc->getParFunc()[1]->Reverse();
	Reverse(1);

	if(!quiet) printf("Computing probabilities for %s and %s whose lengths are %d and %d... \n", sequenceUp->getName(), sequenceDown->getName(), sequenceUp->getLen(), sequenceDown->getLen());

	if(DEBUG) printf("tables: %d mem: %ld\n", big_tables_num, total_tables_size);
	
	*PI->estar(0, lenUp-1, 0, lenDown-1) = 1.0;

	for(unsigned int l1 = itr; l1 >= 1; l1--)
	{
		for(unsigned int l2 = lenDown; l2 >= 1; l2--)
		{
			if(((l1 - 1) * lenDown + l2) % 10 == 0) if(!quiet) printf("Iteration %d out of %d done.\n", (l1 - 1) * lenDown + l2, lenUp*lenDown);

/*			for(int i1 = 0; (unsigned int)i1 <= lenUp - l1; i1++)
				for(unsigned int i2 = 0; i2 <= lenDown - l2; i2++)*/
#pragma omp parallel for num_threads(procNum)
			for(int idx = 0; (unsigned int)idx < (lenUp - l1 + 1)*(lenDown - l2 + 1); idx++)
			{
				int i1 = idx / (lenDown - l2 + 1);
				int i2 = idx % (lenDown - l2 + 1);
				int j1 = i1 + l1 - 1;
				int j2 = i2 + l2 - 1;

				if(DEBUG) printf("Ib, Ia, I\t");
				processI(i1, j1, i2, j2);
				processIa(i1, j1, i2, j2);
				processIb(i1, j1, i2, j2);
				if(DEBUG) printf("xy\t");
				processXY(i1, j1, i2, j2);
				if(DEBUG) printf("mmkk\t");
				processImmkk(i1, j1, i2, j2);
				if(DEBUG) printf("Ia_nndd\t");
				processIa_nndd(i1, j1, i2, j2);
				if(DEBUG) printf("aux\t");
				processAux(i1, j1, i2, j2);
				if(DEBUG) printf("I_nndd\n");
				processI_nndd(i1, j1, i2, j2);
				if(DEBUG) printf("bands\t");
				processBands(i1, j1, i2, j2);
			}
		}
		if(end_time != 0 && end_time < time(NULL)+600 && l1 != 1)
		{
			if(!quiet) printf("Computing probabilities ran out of time... \n");

			Reverse(1);
			parfunc->getParFunc()[1]->Reverse();	

			return l1-1; 
		}
/*		for(unsigned int l2 = lenDown; l2 >= l1; l2--)
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
	for(int i1 = 0; i1 < lenUp; i1++)
		for(int i2 = 0; i2 < lenDown; i2++)
		{
			double sum = 0.0;
			for(int j1 = i1; j1 < lenUp; j1++)
				for(int j2 = i2; j2 < lenDown; j2++)
					sum += PIh[0]->element(i1, j1, i2, j2) + PIh[1]->element(i1, j1, i2, j2);
			*bProb->estar(i1, i1,  lenDown - i2 - 1, lenDown - i2 - 1) = sum;
			if(DEBUG)				
				printf("duplex: %d %d: %lf\n", i1, lenDown - i2 - 1, bProb->element(i1, i1,  lenDown - i2 - 1, lenDown - i2 - 1));
		}


	unsigned int lens[2] = {lenUp, lenDown};

	for(int upOrDown = 0; upOrDown < 2; upOrDown++)
	{
		for(unsigned int l = 1; l <= lens[upOrDown]; l++)
		{
#pragma omp parallel for num_threads(procNum)
			for(int i = 0; i <= lens[upOrDown]-l; i++)
			{
				int j = i + l - 1;
				processP(i, j, upOrDown);
				processPbPm(i, j, upOrDown);
			}
		}

		for(int i = 0; i < lens[upOrDown]; i++)
			for(int j = i; j < lens[upOrDown]; j++)
			{
				*aProb[upOrDown]->estar(i, j) += Pb[upOrDown]->element(i, j);
				for(int d = 0; d < lens[(upOrDown + 1) % 2]; d++)
					for(int e = d; e < lens[(upOrDown + 1) % 2]; e++)
					{
						for(int c = 0; c < 3; c++)
							*aProb[upOrDown]->estar(i, j) += (upOrDown == 0) ? PIx_m[upOrDown][c]->element(i, j, d, e) : PIx_m[upOrDown][c]->element(d, e, i, j);

						for(int c = 0; c < 2; c++)
							*aProb[upOrDown]->estar(i, j) += (upOrDown == 0) ? PIx_k[upOrDown][c]->element(i, j, d, e) : PIx_k[upOrDown][c]->element(d, e, i, j);

						if(upOrDown)
							*aProb[upOrDown]->estar(i, j) += PIy->element(d, e, i, j);
					}

				if(DEBUG)				
					if(upOrDown == 0)
						printf("up: %d %d: %lf\n", i, j, aProb[upOrDown]->element(i, j));
					else
						printf("down: %d %d: %lf\n", i, j, aProb[upOrDown]->element(i, j));
			}
	}
	if(!quiet) printf("Computing probabilities finished. \n");

	Reverse(1);
	parfunc->getParFunc()[1]->Reverse();
	return 0;	
}

void Probability::product(Table<double> *T1, Table<double> *T2, Table<double> *P1, Table<double> *P2, double p, double q, int i1, int j1, int i2, int j2)
{
	if(p == 0.0 || q == 0.0)
		return;

	for(int d = i1; d <= j1-1; d++)
		for(int e = i2; e <= j2-1; e++)
		{
			double deltaP = p * T1->element(i1, d, i2, e)*T2->element(d+1, j1, e+1, j2) / q;
#pragma omp atomic
			*P1->estar(i1, d, i2, e) += deltaP;
#pragma omp atomic
			*P2->estar(d+1, j1, e+1, j2) += deltaP;
		}
}

void Probability::product(Table<double> **T1, Table<double> **T2, Table<double> **P1, Table<double> **P2, double p, double q, int i1, int j1, int i2, int j2)
{
	if(p == 0.0 || q == 0.0)
		return;

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
#pragma omp atomic
				*P2[k]->estar(d+1, j1, e+1, j2) += p * T2[k]->element(d+1, j1, e+1, j2) * sum1 / q;
				sum2 += T2[k]->element(d+1, j1, e+1, j2);
				k++;
			}

			k = 0;			
			while(T1[k]) 
			{
#pragma omp atomic
				*P1[k]->estar(i1, d, i2, e) += p * T1[k]->element(i1, d, i2, e) * sum2 / q;
				k++;
			}
		}
}

void Probability::product(Table<double> **T1, Table<double> *T2, Table<double> **P1, Table<double> *P2, double p, double q, int i1, int j1, int i2, int j2)
{
	if(p == 0.0 || q == 0.0)
		return;

	for(int d = i1; d <= j1-1; d++)
		for(int e = i2; e <= j2-1; e++)
		{
			int  k = 0;
			while(T1[k]) 
			{
				double deltaP = p * T1[k]->element(i1, d, i2, e) * T2->element(d+1, j1, e+1, j2) / q;
#pragma omp atomic
				*P1[k]->estar(i1, d, i2, e) += deltaP;
#pragma omp atomic
				*P2->estar(d+1, j1, e+1, j2) += deltaP;
				k++;
			}
		}
}

void Probability::product(Table<double> *T1, Table<double> **T2, Table<double> *P1, Table<double> **P2, double p, double q, int i1, int j1, int i2, int j2)
{
	if(p == 0.0 || q == 0.0)
		return;
	
	for(int d = i1; d <= j1-1; d++)
		for(int e = i2; e <= j2-1; e++)
		{
			int k = 0;
			while(T2[k]) 
			{
				double deltaP = p * T1->element(i1, d, i2, e) * T2[k]->element(d+1, j1, e+1, j2) / q;
#pragma omp atomic
				*P1->estar(i1, d, i2, e) += deltaP;
#pragma omp atomic
				*P2[k]->estar(d+1, j1, e+1, j2) += deltaP;
				k++;
			}
		}
}

void Probability::rsemiproduct(Table<double> *Tu, Table<double> *Td, Table<double> *T2, Table<double> *Pu, Table<double> *Pd, Table<double> *P2, double p, double q, int i1, int j1, int i2, int j2)
{
	if(p == 0.0 || q == 0.0)
		return;

	for(int d = i1; d <= j1; d++)
		for(int e = i2; e <= j2; e++)
		{
			double temp1 = 1.0, temp2 = 1.0;
			if(d > i1) temp1 = Tu->element(i1, d-1);
			if(e > i2) temp2 = Td->element(i2, e-1);
			double deltaP = p * temp1 * temp2 * T2->element(d, j1, e, j2) / q;
			if(d > i1) 
#pragma omp atomic
				*Pu->estar(i1, d-1) += deltaP;
			if(e > i2)
#pragma omp atomic
				*Pd->estar(i2, e-1) += deltaP;
#pragma omp atomic
			*P2->estar(d, j1, e, j2) += deltaP;
		}
}

void Probability::lsemiproduct(Table<double> *T1, Table<double> *Tu, Table<double> *Td, Table<double> *P1, Table<double> *Pu, Table<double> *Pd, double p, double q, int i1, int j1, int i2, int j2)
{
	if(p == 0.0 || q == 0.0)
		return;

	for(int d = i1; d <= j1; d++){
		double temp1 = 1.0;		
		if(d < j1) temp1 = Tu->element(d+1, j1);
		double deltaP = p * T1->element(i1, d, i2, j2) * temp1 / q;
#pragma omp atomic
		*P1->estar(i1, d, i2, j2) += deltaP;
		if(d < j1) 
#pragma omp atomic
			*Pu->estar(d+1, j1) += deltaP;
	}

	for(int e = i2; e < j2; e++){
		double deltaP = p * T1->element(i1, j1, i2, e) * Td->element(e+1, j2) / q;
#pragma omp atomic
		*P1->estar(i1, j1, i2, e) += deltaP;
#pragma omp atomic
		*Pd->estar(e+1, j2) += deltaP;
	}
}

void Probability::lsemiproduct(Table<double> **T1, Table<double> *Tu, Table<double> *Td, Table<double> **P1, Table<double> *Pu, Table<double> *Pd, double p, double q, int i1, int j1, int i2, int j2)
{
	if(p == 0.0 || q == 0.0)
		return;

	for(int d = i1; d <= j1; d++){
		double temp1 = 1.0;
		if(d < j1) temp1 = Tu->element(d+1, j1);
		int  k = 0;
		while(T1[k]) 
		{
			double deltaP = p * T1[k]->element(i1, d, i2, j2) * temp1 / q;
#pragma omp atomic
			*P1[k]->estar(i1, d, i2, j2) += deltaP;
			if(d < j1) 
#pragma omp atomic
				*Pu->estar(d+1, j1) += deltaP;
			k++;
		}
	}

	for(int e = i2; e < j2; e++){
		int  k = 0;
		while(T1[k]) 
		{
			double deltaP = p * T1[k]->element(i1, j1, i2, e) * Td->element(e+1, j2) / q;				
#pragma omp atomic
			*P1[k]->estar(i1, j1, i2, e) += deltaP;
#pragma omp atomic
			*Pd->estar(e+1, j2) += deltaP;
			k++;
		}
	}
}

void Probability::processBands(int i1, int j1, int i2, int j2)
{
	if(!energy->base_pair(seqUp[i1], seqDown[i2]))
		return;

	//QIhh and QIhm
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
		{
			double pihh = PIhh[i][j]->element(i1, j1, i2, j2);
			double pihm = PIhm[i][j]->element(i1, j1, i2, j2);
			double qihh = QIhh[i][j]->element(i1, j1, i2, j2);
			double qihm = QIhm[i][j]->element(i1, j1, i2, j2);

			if((pihh == 0.0 && pihm == 0.0) || (qihh == 0.0 && qihm == 0.0))
				continue;

			for(int d = i1; d <= j1; d++)
				for(int e = i2; e <= j2; e++)
				{
					double temp1 = 1.0, temp2 = 1.0, temp5 = 0.0, temp6 = 0.0, 
					temp3 = (i+j) ? exp(-energy->Ekissing(1, 0, 0)/RT) : exp(-energy->Eintstackpenalty()/RT), temp4;
					double empty = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, j1 - d, 0)/RT);
					if(d < j1 && e < j2)
					{
						temp1 = Q[0][i]->element(d+1, j1);
						temp2 = Q[1][j]->element(e+1, j2);
						temp5 = Qm[0][i]->element(d+1, j1)*Q[1][j]->element(e+1, j2);
						temp6 = empty*Qm[1][j]->element(e+1, j2);
					} else if(d < j1 && e == j2)
					{
						temp1 = Q[0][i]->element(d+1, j1);
						temp5 = Qm[0][i]->element(d+1, j1);
					} else if(d == j1 && e < j2)
					{
						temp2 = Q[1][j]->element(e+1, j2);
						temp6 = Qm[1][j]->element(e+1, j2);
					}

					temp4 = 1.0;

					// penalize A-U stack terminal
					if((seqUp[i1] == 0 && seqDown[i2] == 3) || (seqUp[i1] == 3 && seqDown[i2] == 0))
						temp4 *= exp(-energy->Eaupenalty()/RT);
					if((seqUp[d] == 0 && seqDown[e] == 3) || (seqUp[d] == 3 && seqDown[e] == 0))
						temp4 *= exp(-energy->Eaupenalty()/RT);

					double deltaPh = 0.0, deltaPm1 = 0.0, deltaPm2 = 0.0;

					if (qihh != 0.0) 
						deltaPh = pihh * QIh[i|j]->element(i1, d, i2, e)*temp1*temp2*temp3*temp4 / qihh;

					
					if (qihm != 0.0) 
					{
						deltaPm1 = pihm * QIh[i|j]->element(i1, d, i2, e)*temp5*temp3*temp4 / qihm;
						deltaPm2 = pihm * QIh[i|j]->element(i1, d, i2, e)*temp6*temp3*temp4 / qihm;
					}
					
#pragma omp atomic
					*PIh[i|j]->estar(i1, d, i2, e) += deltaPh + deltaPm1 + deltaPm2;
					/* here we take care of the bonds that appear on the right of Ih in bProb */
					if(d != i1 || e != i2) 
#pragma omp atomic
						*PIh[i|j]->estar(d, d, e, e) += deltaPh + deltaPm1 + deltaPm2; 

					if(d < j1 && e < j2)
					{
#pragma omp atomic
						*P[0][i]->estar(d+1, j1) += deltaPh;
#pragma omp atomic
						*P[1][j]->estar(e+1, j2) += deltaPh;
#pragma omp atomic
						*Pm[0][i]->estar(d+1, j1) += deltaPm1;
#pragma omp atomic
						*P[1][j]->estar(e+1, j2) += deltaPm1;
#pragma omp atomic
						*Pm[1][j]->estar(e+1, j2) += deltaPm2;
					} else if(d < j1 && e == j2)
					{
#pragma omp atomic
						*P[0][i]->estar(d+1, j1) += deltaPh;
#pragma omp atomic
						*Pm[0][i]->estar(d+1, j1) += deltaPm1;
					} else if(d == j1 && e < j2)
					{
#pragma omp atomic
						*P[1][j]->estar(e+1, j2) += deltaPh;
#pragma omp atomic
						*Pm[1][j]->estar(e+1, j2) += deltaPm2;
					}
				}
		}

	//QIh
	if(energy->base_pair(seqUp[j1], seqDown[j2]))
	{
		for(int i = 0; i < 2; i++)
		{
			double pih = PIh[i]->element(i1, j1, i2, j2);
			double qih = QIh[i]->element(i1, j1, i2, j2);

			if(pih == 0.0 || qih == 0.0)
				continue;

			for(int d = i1+1; d < j1 && d < i1 + MAX_HYBRID_LEN; d++)
				for(int e = i2+1; e < j2 && e < i2 + MAX_HYBRID_LEN; e++)
				{
					double temp = pih * QIh[i]->element(d, j1, e, j2) / qih;
					double *pihs = PIh[i]->estar(d, j1, e, j2);

					if(d == i1+1 && e == i2+1)
						if(i == 0)
#pragma omp atomic
							*pihs += temp*exp(-energy->Eintstack(0, i1, i2, seqUp, seqDown)/RT);
						else
#pragma omp atomic
							*pihs += temp*exp(-energy->Ekissingstack(i1, i2, seqUp, seqDown)/RT);
					else
						if(i == 0)
#pragma omp atomic
							*pihs += temp*exp(-energy->Eintbi(i1, i2, d, e, seqUp, seqDown)/RT);
						else
#pragma omp atomic
							*pihs += temp*exp(-energy->Ekissingbi(i1, i2, d, e, seqUp, seqDown)/RT);

				}


/*			for(int d = j1-1; d > i1 && d > j1 - MAX_HYBRID_LEN; d--)
				for(int e = j2-1; e > i2 && e > j2 - MAX_HYBRID_LEN; e--)
				{
					double temp = pih * QIh[i]->element(i1, d, i2, e) / qih;
					double *pihs = PIh[i]->estar(i1, d, i2, e);

					if(d == j1-1 && e == j2-1)
						if(i == 0)
#pragma omp atomic
							*pihs += temp*exp(-energy->Eintstack(0, d, e, seqUp, seqDown)/RT);
						else
#pragma omp atomic
							*pihs += temp*exp(-energy->Ekissingstack(d, e, seqUp, seqDown)/RT);
					else
						if(i == 0)
#pragma omp atomic
							*pihs += temp*exp(-energy->Eintbi(d, e, j1, j2, seqUp, seqDown)/RT);
						else
#pragma omp atomic
							*pihs += temp*exp(-energy->Ekissingbi(d, e, j1, j2, seqUp, seqDown)/RT);

				} */
		}
	}
}


void Probability::capproduct(Table<double> *qir, Table<double> *qig, Table<double> *qi, 
Table<double> *pir, Table<double> *pig, Table<double> *pi, int i1, int j1, int i2, int j2, int upOrDown)
{
	double totalp = pi->element(i1, j1, i2, j2);
	double totalq = qi->element(i1, j1, i2, j2);

	if(totalp == 0.0 || totalq == 0.0)
		return;

	int ii = (upOrDown == 0) ? i1 : i2;
	int jj = (upOrDown == 0) ? j1 : j2;

	Table<double> *q[3] = {Q[upOrDown][0], Q[upOrDown][1], Q[upOrDown][2]};
	Table<double> *qm[3] = {Qm[upOrDown][0], Qm[upOrDown][1], Qm[upOrDown][2]};

	Table<double> *p[3] = {P[upOrDown][0], P[upOrDown][1], P[upOrDown][2]};
	Table<double> *pm[3] = {Pm[upOrDown][0], Pm[upOrDown][1], Pm[upOrDown][2]};

	unsigned char *seq = parfunc->getParFunc()[upOrDown]->getSequence()->getSeq();
	unsigned int len = parfunc->getParFunc()[upOrDown]->getSequence()->getLen();

	for(int d = ii; d <= jj-2; d++)
		for(int e = d+2; e <= jj; e++)
		{
			double tig = (upOrDown == 0) ? qig->element(d+1, e-1, i2, j2) : qig->element(i1, j1, d+1, e-1);
			double tir = (upOrDown == 0) ? qir->element(d+1, e-1, i2, j2) : qir->element(i1, j1, d+1, e-1);
			double ti = (upOrDown == 0) ? qi->element(d+1, e-1, i2, j2) : qi->element(i1, j1, d+1, e-1);
			double *ptig = (upOrDown == 0) ? pig->estar(d+1, e-1, i2, j2) : pig->estar(i1, j1, d+1, e-1);
			double *ptir = (upOrDown == 0) ? pir->estar(d+1, e-1, i2, j2) : pir->estar(i1, j1, d+1, e-1);
			double *pti = (upOrDown == 0) ? pi->estar(d+1, e-1, i2, j2) : pi->estar(i1, j1, d+1, e-1);
       			double rtemp1=1.0, rtemp2=1.0, gtemp1=1.0, gtemp2=1.0;
	
       			if (d > ii) rtemp1=q[1]->element(ii+1, d);
			if (e < jj) rtemp2=q[1]->element(e, jj-1);
       			if (d > ii) gtemp1=q[2]->element(ii+1, d);
			if (e < jj) gtemp2=q[2]->element(e, jj-1);


			double deltaPk = totalp*rtemp1*rtemp2*exp(-energy->Ekissing(0, 0, 1)/RT)*tir / totalq;
			double deltaPm = totalp*gtemp1*gtemp2*exp(-energy->Emulti(0, 0, 1)/RT)*tig / totalq;

#pragma omp atomic
			*ptir += deltaPk;
#pragma omp atomic
			*ptig += deltaPm;
       			if (d > ii) 
#pragma omp atomic
				*p[1]->estar(ii+1, d) += deltaPk;
			if (e < jj) 
#pragma omp atomic
				*p[1]->estar(e, jj-1) += deltaPk;
       			if (d > ii) 
#pragma omp atomic
				*p[2]->estar(ii+1, d) += deltaPm;
			if (e < jj) 
#pragma omp atomic
				*p[2]->estar(e, jj-1) += deltaPm;


			double deltaPi = 0.0;

			if (e == jj && d == ii)
			{
				double es = (upOrDown == 0) ? energy->Es(ii, jj, seq) : energy->Es(len-jj-1, len-ii-1, seq);
				deltaPi += totalp*exp(-es/RT)*ti / totalq;
			}
			else
			{
				double ebi = (upOrDown == 0) ? energy->Ebi(ii, jj, d+1, e-1, seq) : energy->Ebi(len-jj-1, len-ii-1, len-e, len-d-2, seq);
   				deltaPi += totalp*exp(-ebi/RT)*ti / totalq;
				
				if(d > ii)
				{
					double tmp = totalp*qm[2]->element(ii+1, d)*exp(-energy->Emulti(1, jj-e, 1)/RT)*ti / totalq;
#pragma omp atomic
					*pm[2]->estar(ii+1, d) += tmp;
					deltaPi += tmp;
				}
				if(e < jj)
				{
					double tmp = totalp*qm[2]->element(e, jj-1)*exp(-energy->Emulti(1, d-ii, 1)/RT)*ti / totalq;
#pragma omp atomic
					*pm[2]->estar(e, jj-1) += tmp;
					deltaPi += tmp;
				}
				if(d > ii && e < jj)
				{
					double tmp = totalp*qm[2]->element(ii+1, d)*qm[2]->element(e, jj-1)*exp(-energy->Emulti(1, 0, 1)/RT)*ti / totalq;
#pragma omp atomic
					*pm[2]->estar(ii+1, d) += tmp;
#pragma omp atomic
					*pm[2]->estar(e, jj-1) += tmp;
					deltaPi += tmp;
				}
			}

#pragma omp atomic
			*pti += deltaPi;
		}
}

void Probability::processXY(int i1, int j1, int i2, int j2)
{
	//QIy
	if(energy->base_pair(seqUp[i1], seqUp[j1]) && energy->base_pair(seqDown[i2], seqDown[j2]) && j1 > i1 + MIN_HAIRPIN_SIZE && j2 > i2 + MIN_HAIRPIN_SIZE)
		capproduct(QIx_k[0][1], QIx_m[0][2], QIy, PIx_k[0][1], PIx_m[0][2], PIy, i1, j1, i2, j2, 1); 

	//QIx_m
	for(int c = 0; c < 3; c++)
	{
		if(energy->base_pair(seqUp[i1], seqUp[j1]) && j1 > i1 + MIN_HAIRPIN_SIZE)
			capproduct(QIkm[1][c], QImm[2][c], QIx_m[0][c], PIkm[1][c], PImm[2][c], PIx_m[0][c], i1, j1, i2, j2, 0); 

		if(energy->base_pair(seqDown[i2], seqDown[j2]) && j2 > i2 + MIN_HAIRPIN_SIZE)
			capproduct(QImk[c][1], QImm[c][2], QIx_m[1][c], PImk[c][1], PImm[c][2], PIx_m[1][c], i1, j1, i2, j2, 1); 
	}

	//QIx_k
	for(int c = 0; c < 2; c++)
	{
		if(energy->base_pair(seqUp[i1], seqUp[j1]) && j1 > i1 + MIN_HAIRPIN_SIZE)
			capproduct(QIkk[1][c], QImk[2][c], QIx_k[0][c], PIkk[1][c], PImk[2][c], PIx_k[0][c], i1, j1, i2, j2, 0);

		if(energy->base_pair(seqDown[i2], seqDown[j2]) && j2 > i2 + MIN_HAIRPIN_SIZE)
			capproduct(QIkk[c][1], QIkm[c][2], QIx_k[1][c], PIkk[c][1], PIkm[c][2], PIx_k[1][c], i1, j1, i2, j2, 1);
	}
}

void Probability::processAux(int i1, int j1, int i2, int j2)
{
	if(!energy->base_pair(seqUp[j1], seqDown[j2]))
		return;
	
	
	//QI_r, QIlr, QIa_r
	for(int c = 1; c < 4; c++)
	{
		int i = c % 2;
		int j = c / 2;

		double temp1 = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);
		double temp2 = (j == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);
		double temp3 = (i+j) ? exp(-energy->Ekissing(1, 0, 0)/RT) : exp(-energy->Eintstackpenalty()/RT);
		double temp4 = 1.0;

		// penalize A-U stack terminal
		if((seqUp[i1] == 0 && seqDown[i2] == 3) || (seqUp[i1] == 3 && seqDown[i2] == 0))
			temp4 *= exp(-energy->Eaupenalty()/RT);
		if((seqUp[j1] == 0 && seqDown[j2] == 3) || (seqUp[j1] == 3 && seqDown[j2] == 0))
			temp4 *= exp(-energy->Eaupenalty()/RT);

		double deltaP = 0.0;

		if(QIlr[i][j]->element(i1, j1, i2, j2) != 0.0)
			deltaP = PIlr[i][j]->element(i1, j1, i2, j2) * QIh[1]->element(i1, j1, i2, j2)*temp3*temp4 / QIlr[i][j]->element(i1, j1, i2, j2); 
#pragma omp atomic
		*PIh[1]->estar(i1, j1, i2, j2) += deltaP;
		// for bProb
		if(i1 != j1 || i2 != j2) 
#pragma omp atomic
			*PIh[1]->estar(j1, j1, j2, j2) += deltaP;

		product(QIhm[i][j], QIlr[i][j], PIhm[i][j], PIlr[i][j], PIlr[i][j]->element(i1, j1, i2, j2), QIlr[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
		product(QIhh[i][j], QIa_r[i][j], PIhh[i][j], PIa_r[i][j], PIlr[i][j]->element(i1, j1, i2, j2), QIlr[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);

		product(QIx[0][j], QI_r[i][j], PIx[0][j], PI_r[i][j], PIa_r[i][j]->element(i1, j1, i2, j2)*temp1, QIa_r[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIx[1][i], QI_r[i][j], PIx[1][i], PI_r[i][j], PIa_r[i][j]->element(i1, j1, i2, j2)*temp2, QIa_r[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIy, QI_r[i][j], PIy, PI_r[i][j], PIa_r[i][j]->element(i1, j1, i2, j2)*temp1*temp2, QIa_r[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);


		rsemiproduct(Q[0][i], Q[1][j], QIlr[i][j], P[0][i], P[1][j], PIlr[i][j], PI_r[i][j]->element(i1, j1, i2, j2), QI_r[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		rsemiproduct(Q[0][i], Q[1][j], QIa_r[i][j], P[0][i], P[1][j], PIa_r[i][j], PI_r[i][j]->element(i1, j1, i2, j2), QI_r[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
	}
}

void Probability::processImmkk(int i1, int j1, int i2, int j2)
{
	//QImm
	double temp2 = exp(-energy->Emulti(0, 0, 1)/RT);
	double temp3 = exp(-energy->Ekissing(0, 0, 1)/RT);
	for(int i = 0; i < 2; i++)
	{
		double temp1 = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);

		product(QIa_nn[i][2], QIx_m[0][2], PIa_nn[i][2], PIx_m[0][2], PImm[i][2]->element(i1, j1, i2, j2)*temp1, QImm[i][2]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
		product(QIa_nn[i][2], QIx_m[1][i], PIa_nn[i][2], PIx_m[1][i], PImm[i][2]->element(i1, j1, i2, j2)*temp2, QImm[i][2]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
		product(QIa_nn[i][2], QIy, PIa_nn[i][2], PIy, PImm[i][2]->element(i1, j1, i2, j2)*temp1*temp2, QImm[i][2]->element(i1, j1, i2, j2), i1, j1, i2, j2);


		product(QIa_nn[2][i], QIx_m[0][i], PIa_nn[2][i], PIx_m[0][i], PImm[2][i]->element(i1, j1, i2, j2)*temp2, QImm[2][i]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
		product(QIa_nn[2][i], QIx_m[1][2], PIa_nn[2][i], PIx_m[1][2], PImm[2][i]->element(i1, j1, i2, j2)*temp1, QImm[2][i]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
		product(QIa_nn[2][i], QIy, PIa_nn[2][i], PIy, PImm[2][i]->element(i1, j1, i2, j2)*temp1*temp2, QImm[2][i]->element(i1, j1, i2, j2), i1, j1, i2, j2);
	}

	product(QIa_nn[2][2], QIx_m[0][2], PIa_nn[2][2], PIx_m[0][2], PImm[2][2]->element(i1, j1, i2, j2)*temp2, QImm[2][2]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
	product(QIa_nn[2][2], QIx_m[1][2], PIa_nn[2][2], PIx_m[1][2], PImm[2][2]->element(i1, j1, i2, j2)*temp2, QImm[2][2]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
	product(QIa_nn[2][2], QIy, PIa_nn[2][2], PIy, PImm[2][2]->element(i1, j1, i2, j2)*temp2*temp2, QImm[2][2]->element(i1, j1, i2, j2), i1, j1, i2, j2);
	
	//QImk
	for(int i = 0; i < 2; i++)
	{
		double temp1 = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);

		product(QIa_nd[i][1], QIx[0][1], PIa_nd[i][1], PIx[0][1], PImk[i][1]->element(i1, j1, i2, j2)*temp1, QImk[i][1]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIa_nd[i][1], QIx_m[1][i], PIa_nd[i][1], PIx_m[1][i], PImk[i][1]->element(i1, j1, i2, j2)*temp3, QImk[i][1]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
		product(QIa_nn[i][1], QIx_k[0][1], PIa_nn[i][1], PIx_k[0][1], PImk[i][1]->element(i1, j1, i2, j2)*temp1, QImk[i][1]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIa_nd[i][1], QIy, PIa_nd[i][1], PIy, PImk[i][1]->element(i1, j1, i2, j2)*temp1*temp3, QImk[i][1]->element(i1, j1, i2, j2), i1, j1, i2, j2);


		product(QIa_nd[2][i], QIx[0][i], PIa_nd[2][i], PIx[0][i], PImk[2][i]->element(i1, j1, i2, j2)*temp2, QImk[2][i]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIa_nd[2][i], QIx_m[1][2], PIa_nd[2][i], PIx_m[1][2], PImk[2][i]->element(i1, j1, i2, j2)*temp1, QImk[2][i]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
		product(QIa_nn[2][i], QIx_k[0][i], PIa_nn[2][i], PIx_k[0][i], PImk[2][i]->element(i1, j1, i2, j2)*temp2, QImk[2][i]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIa_nd[2][i], QIy, PIa_nd[2][i], PIy, PImk[2][i]->element(i1, j1, i2, j2)*temp2*temp3, QImk[2][i]->element(i1, j1, i2, j2), i1, j1, i2, j2);
	}

	//QIkm
	for(int i = 0; i < 2; i++)
	{
		double temp1 = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);

		product(QIa_dn[1][i], QIx_m[0][i], PIa_dn[1][i], PIx_m[0][i], PIkm[1][i]->element(i1, j1, i2, j2)*temp3, QIkm[1][i]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIa_dn[1][i], QIx[1][1], PIa_dn[1][i], PIx[1][1], PIkm[1][i]->element(i1, j1, i2, j2)*temp1, QIkm[1][i]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
		product(QIa_nn[1][i], QIx_k[1][1], PIa_nn[1][i], PIx_k[1][1], PIkm[1][i]->element(i1, j1, i2, j2)*temp1, QIkm[1][i]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIa_dn[1][i], QIy, PIa_dn[1][i], PIy, PIkm[1][i]->element(i1, j1, i2, j2)*temp1*temp3, QIkm[1][i]->element(i1, j1, i2, j2), i1, j1, i2, j2);

		product(QIa_dn[i][2], QIx_m[0][2], PIa_dn[i][2], PIx_m[0][2], PIkm[i][2]->element(i1, j1, i2, j2)*temp1, QIkm[i][2]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIa_dn[i][2], QIx[1][i], PIa_dn[i][2], PIx[1][i], PIkm[i][2]->element(i1, j1, i2, j2)*temp2, QIkm[i][2]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
		product(QIa_nn[i][2], QIx_k[1][i], PIa_nn[i][2], PIx_k[1][i], PIkm[i][2]->element(i1, j1, i2, j2)*temp2, QIkm[i][2]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIa_dn[i][2], QIy, PIa_dn[i][2], PIy, PIkm[i][2]->element(i1, j1, i2, j2)*temp1*temp2, QIkm[i][2]->element(i1, j1, i2, j2), i1, j1, i2, j2);
	}

	//QIkk
	for(int c = 1; c < 4; c++)
	{
		int i = c % 2;
		int j = c / 2;

		double temp1 = (i == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);
		double temp4 = (j == 0) ? 1.0 : exp(-energy->Ekissing(0, 0, 1)/RT);

		if(PIkk[i][j]->element(i1, j1, i2, j2) == 0.0 || QIkk[i][j]->element(i1, j1, i2, j2) == 0.0)
			continue;

#pragma omp atomic
		*PIa_r[i][j]->estar(i1, j1, i2, j2) += PIkk[i][j]->element(i1, j1, i2, j2) * QIa_r[i][j]->element(i1, j1, i2, j2) / QIkk[i][j]->element(i1, j1, i2, j2); 

#pragma omp atomic
		*PIlr[i][j]->estar(i1, j1, i2, j2) += PIkk[i][j]->element(i1, j1, i2, j2) * QIlr[i][j]->element(i1, j1, i2, j2) / QIkk[i][j]->element(i1, j1, i2, j2); 

		product(auxQIa_dd[i][j], QIx[0][j], auxPIa_dd[i][j], PIx[0][j], PIkk[i][j]->element(i1, j1, i2, j2)*temp1, QIkk[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(auxQIa_dd[i][j], QIx[1][i], auxPIa_dd[i][j], PIx[1][i], PIkk[i][j]->element(i1, j1, i2, j2)*temp4, QIkk[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(auxQIa_dd[i][j], QIy, auxPIa_dd[i][j], PIy, PIkk[i][j]->element(i1, j1, i2, j2)*temp1*temp4, QIkk[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIa_dn[i][j], QIx_k[0][j], PIa_dn[i][j], PIx_k[0][j], PIkk[i][j]->element(i1, j1, i2, j2)*temp1, QIkk[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		product(QIa_nd[i][j], QIx_k[1][i], PIa_nd[i][j], PIx_k[1][i], PIkk[i][j]->element(i1, j1, i2, j2)*temp4, QIkk[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
	}
}

void Probability::processIa_nndd(int i1, int j1, int i2, int j2)
{
	//QIa_nn
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			if(i+j)
			{
				double temp3 = 1.0, temp4 = 1.0;
				if(i == 1) temp3 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(i == 2) temp3 = exp(-energy->Emulti(0, 0, 1)/RT);
				if(j == 1) temp4 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(j == 2) temp4 = exp(-energy->Emulti(0, 0, 1)/RT);

				product(QIx_m[0][j], QI_nn[i][j], PIx_m[0][j], PI_nn[i][j], PIa_nn[i][j]->element(i1, j1, i2, j2)*temp3, QIa_nn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
				product(QIx_m[1][i], QI_nn[i][j], PIx_m[1][i], PI_nn[i][j], PIa_nn[i][j]->element(i1, j1, i2, j2)*temp4, QIa_nn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				product(QIy, QI_nn[i][j], PIy, PI_nn[i][j], PIa_nn[i][j]->element(i1, j1, i2, j2)*temp3*temp4, QIa_nn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);

				lsemiproduct(QIx_m[0][j], Q[0][i], Q[1][j], PIx_m[0][j], P[0][i], P[1][j], PIa_nn[i][j]->element(i1, j1, i2, j2)*temp3, QIa_nn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				lsemiproduct(QIx_m[1][i], Q[0][i], Q[1][j], PIx_m[1][i], P[0][i], P[1][j], PIa_nn[i][j]->element(i1, j1, i2, j2)*temp4, QIa_nn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				lsemiproduct(QIy, Q[0][i], Q[1][j], PIy, P[0][i], P[1][j], PIa_nn[i][j]->element(i1, j1, i2, j2)*temp3*temp4, QIa_nn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
			}

	//QIa_nd
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 2; j++)
			if(j == 1 || (i != 0 && j == 0))
			{
				double temp3 = 1.0, temp4 = 1.0;
				if(i == 1) temp3 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(i == 2) temp3 = exp(-energy->Emulti(0, 0, 1)/RT);
				if(j == 1) temp4 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(j == 2) temp4 = exp(-energy->Emulti(0, 0, 1)/RT);

				product(QIx[0][j], QI_nd[i][j], PIx[0][j], PI_nd[i][j], PIa_nd[i][j]->element(i1, j1, i2, j2)*temp3, QIa_nd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
				product(QIx_k[0][j], QI_nn[i][j], PIx_k[0][j], PI_nn[i][j], PIa_nd[i][j]->element(i1, j1, i2, j2)*temp3, QIa_nd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				product(QIx_m[1][i], QI_nd[i][j], PIx_m[1][i], PI_nd[i][j], PIa_nd[i][j]->element(i1, j1, i2, j2)*temp4, QIa_nd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				product(QIy, QI_nd[i][j], PIy, PI_nd[i][j], PIa_nd[i][j]->element(i1, j1, i2, j2)*temp3*temp4, QIa_nd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);

				lsemiproduct(QIx_k[0][j], Q[0][i], Q[1][j], PIx_k[0][j], P[0][i], P[1][j], PIa_nd[i][j]->element(i1, j1, i2, j2)*temp3, QIa_nd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
			}

	//QIa_dn
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 3; j++)
			if(i == 1 || (j != 0 && i == 0))
			{
				double temp3 = 1.0, temp4 = 1.0;
				if(i == 1) temp3 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(i == 2) temp3 = exp(-energy->Emulti(0, 0, 1)/RT);
				if(j == 1) temp4 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(j == 2) temp4 = exp(-energy->Emulti(0, 0, 1)/RT);

				product(QIx_m[0][j], QI_dn[i][j], PIx_m[0][j], PI_dn[i][j], PIa_dn[i][j]->element(i1, j1, i2, j2)*temp3, QIa_dn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
				product(QIx[1][i], QI_dn[i][j], PIx[1][i], PI_dn[i][j], PIa_dn[i][j]->element(i1, j1, i2, j2)*temp4, QIa_dn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				product(QIx_k[1][i], QI_nn[i][j], PIx_k[1][i], PI_nn[i][j], PIa_dn[i][j]->element(i1, j1, i2, j2)*temp4, QIa_dn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				product(QIy, QI_dn[i][j], PIy, PI_dn[i][j], PIa_dn[i][j]->element(i1, j1, i2, j2)*temp3*temp4, QIa_dn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);

				lsemiproduct(QIx_k[1][i], Q[0][i], Q[1][j], PIx_k[1][i], P[0][i], P[1][j], PIa_dn[i][j]->element(i1, j1, i2, j2)*temp4, QIa_dn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
			}

	//QIa_dd
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			if(i+j)
			{
				double temp1, temp2, temp3, temp4;
				if(i == 0) temp1 = 1.0, temp3 = 1.0;
				if(i == 1) temp1 = exp(-energy->Ekissing(0, 1, 0)/RT), temp3 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(i == 2) temp1 = exp(-energy->Emulti(0, 1, 0)/RT), temp3 = exp(-energy->Emulti(0, 0, 1)/RT);
				if(j == 0) temp2 = 1.0, temp4 = 1.0;
				if(j == 1) temp2 = exp(-energy->Ekissing(0, 1, 0)/RT), temp4 = exp(-energy->Ekissing(0, 0, 1)/RT);
				if(j == 2) temp2 = exp(-energy->Emulti(0, 1, 0)/RT), temp4 = exp(-energy->Emulti(0, 0, 1)/RT);

				product(QIx[0][j], QI_dd[i][j], PIx[0][j], PI_dd[i][j], PIa_dd[i][j]->element(i1, j1, i2, j2)*temp3, QIa_dd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2); 
				product(QIx_k[0][j], QI_dn[i][j], PIx_k[0][j], PI_dn[i][j], PIa_dd[i][j]->element(i1, j1, i2, j2)*temp3, QIa_dd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				product(QIx_k[1][i], QI_nd[i][j], PIx_k[1][i], PI_nd[i][j], PIa_dd[i][j]->element(i1, j1, i2, j2)*temp4, QIa_dd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				product(QIx[1][i], QI_dd[i][j], PIx[1][i], PI_dd[i][j], PIa_dd[i][j]->element(i1, j1, i2, j2)*temp4, QIa_dd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				product(QIy, QI_dd[i][j], PIy, PI_dd[i][j], PIa_dd[i][j]->element(i1, j1, i2, j2)*temp3*temp4, QIa_dd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
			}
}

void Probability::processI_nndd(int i1, int j1, int i2, int j2)
{
	//QI_nn
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			if(i+j)
			{
				if(QI_nn[i][j]->element(i1, j1, i2, j2) == 0.0)
					continue;
 
				double deltaP = PI_nn[i][j]->element(i1, j1, i2, j2) * Q[0][i]->element(i1, j1)*Q[1][j]->element(i2, j2) / QI_nn[i][j]->element(i1, j1, i2, j2);
#pragma omp atomic
				*P[0][i]->estar(i1, j1) += deltaP;
#pragma omp atomic
				*P[1][j]->estar(i2, j2) += deltaP;

				rsemiproduct(Q[0][i], Q[1][j], QIa_nn[i][j], P[0][i], P[1][j], PIa_nn[i][j], PI_nn[i][j]->element(i1, j1, i2, j2), QI_nn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
			}

	//QI_nd
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 2; j++)
			if(j == 1 || (i != 0 && j == 0))
				rsemiproduct(Q[0][i], Q[1][j], QIa_nd[i][j], P[0][i], P[1][j], PIa_nd[i][j], PI_nd[i][j]->element(i1, j1, i2, j2), QI_nd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);


	//QI_dn
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 3; j++)
			if(i == 1 || (j != 0 && i == 0))
				rsemiproduct(Q[0][i], Q[1][j], QIa_dn[i][j], P[0][i], P[1][j], PIa_dn[i][j], PI_dn[i][j]->element(i1, j1, i2, j2), QI_dn[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);


	//QI_dd
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			if(i+j)
			{
				rsemiproduct(Q[0][i], Q[1][j], QIa_dd[i][j], P[0][i], P[1][j], PIa_dd[i][j], PI_dd[i][j]->element(i1, j1, i2, j2), QI_dd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
				rsemiproduct(Q[0][i], Q[1][j], QIb[i][j], P[0][i], P[1][j], PIb[i][j], PI_dd[i][j]->element(i1, j1, i2, j2), QI_dd[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
			}
}

void Probability::processI(int i1, int j1, int i2, int j2)
{
	//QI
	double deltaP; 

	if (QI->element(i1, j1, i2, j2) != 0.0)
		deltaP = PI->element(i1, j1, i2, j2) * Q[0][0]->element(i1, j1)*Q[1][0]->element(i2, j2) / QI->element(i1, j1, i2, j2);

#pragma omp atomic
	*P[0][0]->estar(i1, j1) += deltaP; 
#pragma omp atomic
	*P[1][0]->estar(i2, j2) += deltaP; 

	rsemiproduct(Q[0][0], Q[1][0], QIa, P[0][0], P[1][0], PIa, PI->element(i1, j1, i2, j2), QI->element(i1, j1, i2, j2), i1, j1, i2, j2);
	rsemiproduct(Q[0][0], Q[1][0], QIb[0][0], P[0][0], P[1][0], PIb[0][0], PI->element(i1, j1, i2, j2), QI->element(i1, j1, i2, j2), i1, j1, i2, j2);
}

void Probability::processIb(int i1, int j1, int i2, int j2)
{
	//QIb
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
		{
#pragma omp atomic
			*PIhh[i][j]->estar(i1, j1, i2, j2) += PIb[i][j]->element(i1, j1, i2, j2) * QIhh[i][j]->element(i1, j1, i2, j2) / QIb[i][j]->element(i1, j1, i2, j2);


			if(i+j == 0)
				product(QIhh[i][j], QIa, PIhh[i][j], PIa, PIb[i][j]->element(i1, j1, i2, j2), QIb[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
			else
				product(QIhh[i][j], auxQIa[i][j], PIhh[i][j], auxPIa[i][j], PIb[i][j]->element(i1, j1, i2, j2), QIb[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);


			product(QIhm[i][j], QIb[i][j], PIhm[i][j], PIb[i][j], PIb[i][j]->element(i1, j1, i2, j2), QIb[i][j]->element(i1, j1, i2, j2), i1, j1, i2, j2);
		}
}
void Probability::processIa(int i1, int j1, int i2, int j2)
{
	//QIa
	double qia = QIa->element(i1, j1, i2, j2);
	double pia = PIa->element(i1, j1, i2, j2);

	product(QIx[0][0], QI, PIx[0][0], PI, pia, qia, i1, j1, i2, j2); 
	product(QIx[1][0], QI, PIx[1][0], PI, pia, qia, i1, j1, i2, j2);  
	product(QIy, QI, PIy, PI, pia, qia, i1, j1, i2, j2);
				
	lsemiproduct(QIx[0][0], Q[0][0], Q[1][0], PIx[0][0], P[0][0], P[1][0], pia, qia, i1, j1, i2, j2); 
	lsemiproduct(QIx[1][0], Q[0][0], Q[1][0], PIx[1][0], P[0][0], P[1][0], pia, qia, i1, j1, i2, j2); 
	lsemiproduct(QIy, Q[0][0], Q[1][0], PIy, P[0][0], P[1][0], pia, qia, i1, j1, i2, j2); 		
}

void Probability::print(FILE *outp)
{
	for(int i = 0; i < 2; i++)
	{
		fprintf(outp, "# probabilities of arcs in %s \n", parfunc->getParFunc()[i]->getSequence()->getName());
		for(int i1 = 0; i1 < parfunc->getParFunc()[i]->getSequence()->getLen(); i1++)
			for(int j1 = i1 + MIN_HAIRPIN_SIZE + 1; j1 < parfunc->getParFunc()[i]->getSequence()->getLen(); j1++)
				fprintf(outp, "%d\t%d\t%g\n", i1+1, j1+1, aProb[i]->element(i1, j1));
	}

	fprintf(outp, "# probabilities of bonds between %s and %s \n", parfunc->getParFunc()[0]->getSequence()->getName(), 
									parfunc->getParFunc()[1]->getSequence()->getName());
	for(int i1 = 0; i1 < parfunc->getParFunc()[0]->getSequence()->getLen(); i1++)
		for(int i2 = 0; i2 < parfunc->getParFunc()[1]->getSequence()->getLen(); i2++)
			fprintf(outp, "%d\t%d\t%g\n", i1+1, i2+1, bProb->element(i1, i1, i2, i2));
}

void Probability::store(FILE *outp)
{
	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			P[i][j]->store(outp);
			Pm[i][j]->store(outp);
		}
		Pb[i]->store(outp);
	}

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(i+j) PIa_nn[i][j]->store(outp);
			if(i+j) PI_nn[i][j]->store(outp);
		}
	}

	for(int i = 0; i < 2; i++)
	{
		PIh[i]->store(outp);
		PIx_m[i][2]->store(outp);
		PImm[i][2]->store(outp);
		PImm[2][i]->store(outp);
		PImk[i][1]->store(outp);
		PImk[2][i]->store(outp);
		PIkm[1][i]->store(outp);
		PIkm[i][2]->store(outp);
		PIa_nd[i][1]->store(outp);
		PIa_nd[2][i]->store(outp);
		PI_nd[i][1]->store(outp);
		PI_nd[2][i]->store(outp);
		PIa_dn[1][i]->store(outp);
		PIa_dn[i][2]->store(outp);
		PI_dn[1][i]->store(outp);
		PI_dn[i][2]->store(outp);
		for(int j = 0; j < 2; j++)
		{
			PIb[i][j]->store(outp);
			PIhh[i][j]->store(outp);
			PIhm[i][j]->store(outp);
			PIx_m[i][j]->store(outp);
			PIx_k[i][j]->store(outp);
			if(i+j) PIkk[i][j]->store(outp);
			if(i+j) PIa_r[i][j]->store(outp);
			if(i+j) PI_r[i][j]->store(outp);
			if(i+j) PIlr[i][j]->store(outp);
			if(i+j) PIa_dd[i][j]->store(outp);
			if(i+j) PI_dd[i][j]->store(outp);
		}
	}
	PImm[2][2]->store(outp);
	PIa_nd[1][0]->store(outp);
	PI_nd[1][0]->store(outp);
	PIa_dn[0][1]->store(outp);
	PI_dn[0][1]->store(outp);
	PIy->store(outp);
	PI->store(outp);
	PIa->store(outp);
}

void Probability::retrieve(FILE *inp)
{
	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			P[i][j]->retrieve(inp);
			Pm[i][j]->retrieve(inp);
		}
		Pb[i]->retrieve(inp);
	}

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(i+j) PIa_nn[i][j]->retrieve(inp);
			if(i+j) PI_nn[i][j]->retrieve(inp);
		}
	}

	for(int i = 0; i < 2; i++)
	{
		PIh[i]->retrieve(inp);
		PIx_m[i][2]->retrieve(inp);
		PImm[i][2]->retrieve(inp);
		PImm[2][i]->retrieve(inp);
		PImk[i][1]->retrieve(inp);
		PImk[2][i]->retrieve(inp);
		PIkm[1][i]->retrieve(inp);
		PIkm[i][2]->retrieve(inp);
		PIa_nd[i][1]->retrieve(inp);
		PIa_nd[2][i]->retrieve(inp);
		PI_nd[i][1]->retrieve(inp);
		PI_nd[2][i]->retrieve(inp);
		PIa_dn[1][i]->retrieve(inp);
		PIa_dn[i][2]->retrieve(inp);
		PI_dn[1][i]->retrieve(inp);
		PI_dn[i][2]->retrieve(inp);
		for(int j = 0; j < 2; j++)
		{
			PIb[i][j]->retrieve(inp);
			PIhh[i][j]->retrieve(inp);
			PIhm[i][j]->retrieve(inp);
			PIx_m[i][j]->retrieve(inp);
			PIx_k[i][j]->retrieve(inp);
			if(i+j) PIkk[i][j]->retrieve(inp);
			if(i+j) PIa_r[i][j]->retrieve(inp);
			if(i+j) PI_r[i][j]->retrieve(inp);
			if(i+j) PIlr[i][j]->retrieve(inp);
			if(i+j) PIa_dd[i][j]->retrieve(inp);
			if(i+j) PI_dd[i][j]->retrieve(inp);
		}
	}
	PImm[2][2]->retrieve(inp);
	PIa_nd[1][0]->retrieve(inp);
	PI_nd[1][0]->retrieve(inp);
	PIa_dn[0][1]->retrieve(inp);
	PI_dn[0][1]->retrieve(inp);
	PIy->retrieve(inp);
	PI->retrieve(inp);
	PIa->retrieve(inp);
}

// this is a test method that verifies interstrand symmetries provided that the two input sequences are identical.
#define ERR 1e-2
void Probability::selftest(int i1, int j1, int i2, int j2, int l1, int l2)
{
	if(fabs(PI->element(i1, j1, i2, j2) - PI->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
	{
		printf("PI fails self test: %d\t%d\t%d\t%d\n", i1, j1, i2, j2);
		exit(1);
	}
	if(fabs(PIy->element(i1, j1, i2, j2) - PIy->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
	{
		printf("PIy fails self test: %d\t%d\t%d\t%d\n", i1, j1, i2, j2);
		printf("%g\t%g\n", PIy->element(i1, j1, i2, j2), PIy->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1));
		exit(1);
	}

	for(int c1 = 0; c1 < 2; c1++)
		for(int c2 = 0; c2 < 2; c2++)
			if(c1+c2)
			{
				if(fabs(PIb[c1][c2]->element(i1, j1, i2, j2) -  PI_r[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("PIb[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
				if(fabs(PIkk[c1][c2]->element(i1, j1, i2, j2) - PIkk[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("PIkk[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
				if(fabs(PI_dd[c1][c2]->element(i1, j1, i2, j2) - PI_dd[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("PI_dd[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					printf("%g\t%g\n", QI_dd[c1][c2]->element(i1, j1, i2, j2), QI_dd[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1));
					exit(1);
				}
				if(fabs(PIlr[c1][c2]->element(i1, j1, i2, j2) - PIlr[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("PIlr[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					printf("%g\t%g\n", PIlr[c1][c2]->element(i1, j1, i2, j2), PIlr[c1][c2]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1));
					exit(1);
				}
				if(fabs(PI_r[c1][c2]->element(i1, j1, i2, j2) - PIb[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("PI_r[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
			}

	for(int c1 = 0; c1 < 3; c1++)
		for(int c2 = 0; c2 < 3; c2++)
		{
			if(c1 == 2 || c2 == 2)			
				if(fabs(PImm[c1][c2]->element(i1, j1, i2, j2) - PImm[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("PImm[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
			if(c1+c2)
				if(fabs(PI_nn[c1][c2]->element(i1, j1, i2, j2) - PI_nn[c2][c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
				{
					printf("PI_nn[%d][%d] fails self test: %d\t%d\t%d\t%d\n", c1, c2, i1, j1, i2, j2);
					exit(1);
				}
		}	
	for(int c1 = 0; c1 < 2; c1++)
		if(fabs(PIh[c1]->element(i1, j1, i2, j2) - PIh[c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1)) > ERR)
		{
			printf("PIh[%d] fails self test: %d\t%d\t%d\t%d\n", c1, i1, j1, i2, j2);
			printf("%g\t%g\n", PIh[c1]->element(i1, j1, i2, j2), PIh[c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1));
			printf("%g\t%g\n", QIh[c1]->element(i1, j1, i2, j2), QIh[c1]->element(l2 - j2 - 1, l2 - i2 - 1, l1 - j1 - 1, l1 - i1 - 1));
			exit(1);
		}
}
