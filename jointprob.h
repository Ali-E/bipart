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
	Desc: JointProb class computes the conditional and joint probabilities of base pairs and unpaired regions in a single sequence.

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow,
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: April 16, 2009
	Ver: 1.3
*/

#ifndef JOINTPROB_H
#define JOINTPROB_H

#include "energy.h"
#include "sequence.h"
#include "table.h"
#include "alloc.h"
#include "collection.h"
#include "mtree.h"


class JointProb {
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


	Energy *energy;

	Sequence *sequence;

	int sublen;

	Table<double> *Q, *Qm, *Qb, *Qm2;
	Table<double> *P, *Pm, *Pb, *Pm2;
	Table<double> *Pfree;
	Table<double> *Phairpin;

	double RT;

	RegionCollection *bpCollection; // collection of base pairs
	RegionCollection *frCollection; // collection of free regions
	RegionCollection *hpCollection; // collection of hairpins
	RegionCollection *sitesCollection; // collection of site regions for which joint probability is computed
	Table<double> *bp;
	Table<double> *fr;
	Table<double> *hp;
	Table<double> *jointfree;
	Table<double> *mutualinfo;

	RegionCollection *forcedbp; // collection of forced base pairs
	RegionCollection *forcednotbp; // collection of forced not base pairs
	RegionCollection *forcedup; // collection of forced unpaired regions

	MTree *markovTree;

	void computeQ(int i, int j);
	void computeQbQm(int i, int j);

	void processP(int i, int j);
	void processPbPm(int i, int j);

	void filter(RegionCollection *, Table<double> *, Table<double> *, double threshold);
	void computeParFunc(bool);
	void computeProb(bool);
	void buildBayesianNet();
public:
	JointProb(Energy *energy, Sequence *, int, bool _debug = false, int _procNum = 16, bool _quiet = false, bool _no_isolate = true);
	JointProb(JointProb *, Region<double> *);
	~JointProb();
	void setSites(RegionCollection *);
	void reportProbs(FILE *);
	void dumpTree(FILE *f) {markovTree->dump(f);};
	double basepairProb(int, int);
	double freeProb(int, int);
	double hairpinProb(int, int);
	double basepairEnergy(int, int);
	double freeEnergy(int, int);
	double hairpinEnergy(int, int);
	double freeProb(int);
	double jointFreeProb(int *list) {return markovTree->inference(list);};
	double jointFreeEnergy(int *list) {return -RT*log(markovTree->inference(list));};
	double jointFreeProb(int, int);
	double jointFreeEnergy(int, int);
	double mutualInfo(int, int);
	Energy *getEnergy() {return energy;}
	bool getDebug() {return DEBUG;}
	int getProcNum() {return procNum;}
	bool getNo_isolate() {return no_isolate;}
	Sequence *getSequence() {return sequence;}
	int getSubLen() {return sublen;}
};

#endif
