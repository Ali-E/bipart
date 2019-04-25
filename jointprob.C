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

	Last Update by Hamid Reza Chitsaz: April 17, 2009
	Ver: 1.5
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
#include "jointprob.h"
#include "getopt.h"

#define ZERO -1e-7

int inline max(int a, int b)
{
	int ret = (a > b) ? a : b;
	return ret;
}

JointProb::JointProb(Energy *en, Sequence *seq, int slen, bool _debug, int _procNum, bool _quiet, bool _no_isolate)
{
	double temp;
	energy = en;
	energy->getParams(&NA, &polymer, &naConc, &mgConc, &suffix, &temp, &zip, &nodangle, &energyType);
	RT = energy->getRT();

	DEBUG = _debug;
	procNum = _procNum;
	quiet = _quiet;
	no_isolate = _no_isolate;

	sequence = seq;
	sublen = slen;
	char d2 = 2;

	if(!quiet) printf("Allocating memory for single partition function of %s... \n", sequence->getName());
	Q = new Table<double>(sequence->getLen(), &d2, 0.0);
	Qm = new Table<double>(sequence->getLen(), &d2, 0.0);
	Qm2 = new Table<double>(sequence->getLen(), &d2, 0.0);
	Qb = new Table<double>(sequence->getLen(), &d2, 0.0);
	P = new Table<double>(sequence->getLen(), &d2, 0.0);
	Pm = new Table<double>(sequence->getLen(), &d2, 0.0);
	Pm2 = new Table<double>(sequence->getLen(), &d2, 0.0);
	Pb = new Table<double>(sequence->getLen(), &d2, 0.0);
	Pfree = new Table<double>(sequence->getLen(), sublen, &d2, 0.0);
	Phairpin = new Table<double>(sequence->getLen(), &d2, 0.0);

	forcedbp = new RegionCollection();
	forcednotbp = new RegionCollection();
	forcedup = new RegionCollection();

	frCollection = new RegionCollection();
	bpCollection = new RegionCollection();
	hpCollection = new RegionCollection();
	fr = new Table<double>(sequence->getLen(), sublen, &d2, 0.0);
	bp = new Table<double>(sequence->getLen(), &d2, 0.0);
	hp = new Table<double>(sequence->getLen(), &d2, 0.0);

	jointfree = NULL;
	mutualinfo = NULL;
	markovTree = NULL;
	if(!quiet) printf("Allocation and initialization finished. \n");

	computeParFunc(true);
	computeProb(true);
}

JointProb::JointProb(JointProb *jp, Region<double> *r)
{
	double temp;
	energy = jp->getEnergy();
	energy->getParams(&NA, &polymer, &naConc, &mgConc, &suffix, &temp, &zip, &nodangle, &energyType);
	RT = energy->getRT();

	DEBUG = jp->getDebug();
	procNum = jp->getProcNum();
	quiet = true;
	no_isolate = jp->getNo_isolate();

	sequence = jp->getSequence();
	sublen = jp->getSubLen();
	char d2 = 2;

	if(!quiet) printf("Allocating memory for single partition function of %s... \n", sequence->getName());
	Q = new Table<double>(sequence->getLen(), &d2, 0.0);
	Qm = new Table<double>(sequence->getLen(), &d2, 0.0);
	Qm2 = new Table<double>(sequence->getLen(), &d2, 0.0);
	Qb = new Table<double>(sequence->getLen(), &d2, 0.0);
	P = new Table<double>(sequence->getLen(), &d2, 0.0);
	Pm = new Table<double>(sequence->getLen(), &d2, 0.0);
	Pm2 = new Table<double>(sequence->getLen(), &d2, 0.0);
	Pb = new Table<double>(sequence->getLen(), &d2, 0.0);
	Pfree = new Table<double>(sequence->getLen(), sublen, &d2, 0.0);
	Phairpin = new Table<double>(sequence->getLen(), &d2, 0.0);

	forcedbp = new RegionCollection();
	forcednotbp = new RegionCollection();
	forcedup = new RegionCollection();

	frCollection = new RegionCollection();
	bpCollection = new RegionCollection();
	hpCollection = new RegionCollection();
	fr = new Table<double>(sequence->getLen(), sublen, &d2, 0.0);
	bp = new Table<double>(sequence->getLen(), &d2, 0.0);
	hp = new Table<double>(sequence->getLen(), &d2, 0.0);

	jointfree = NULL;
	mutualinfo = NULL;
	markovTree = NULL;
	if(!quiet) printf("Allocation and initialization finished. \n");

	forcedup->add(r);
	computeParFunc(false);
	computeProb(false);
}


JointProb::~JointProb()
{
	delete Q;
	delete Qm;
	delete Qm2;
	delete Qb;
	delete P;
	delete Pm;
	delete Pm2;
	delete Pb;
	delete Pfree;
	delete Phairpin;

	delete forcedbp;
	delete forcednotbp;
	delete forcedup;

	delete bpCollection;
	delete frCollection;
	delete hpCollection;
	delete bp;
	delete fr;
	delete hp;

	if(jointfree) delete jointfree;
	if(mutualinfo) delete mutualinfo;
	if(markovTree) delete markovTree;
}

void JointProb::computeParFunc(bool show)
{
	if(!quiet && show) 
		printf("Computing partition function for %s whose length is %d... \n", sequence->getName(), sequence->getLen());

	unsigned int len = sequence->getLen();
	for(unsigned int l = 1; l <= len; l++)
	{
#pragma omp parallel for num_threads(procNum)
		for(int i = 0; (unsigned int)i <= len - l; i++)
		{
			
			int j = i + l - 1;
			if(DEBUG) printf("i:%d\tj:%d\n", i, j);
			computeQbQm(i, j);
			computeQ(i, j);
		}
		
		if(!quiet && show && !(l % 10)) printf("Length %d done.\n", l);		
	}

	if(!quiet && show) 
		printf("Computing partition function finished. \n");
}

void JointProb::computeQbQm(int i, int j)
{
	double temp = 0.0;
	unsigned char *seq = sequence->getSeq();

	if(j > i + MIN_HAIRPIN_SIZE && energy->base_pair(seq[i], seq[j]) && !forcednotbp->contains(Region<double>(i, j)) && !forcedup->includes(Region<double>(i, i)) && !forcedup->includes(Region<double>(j, j)))
	{
		temp = exp(-energy->Eh(i, j, seq)/RT) + Qm2->element(i+1, j-1)*exp(-energy->Emulti(1, 0, 0)/RT);

		for(int d = i + 1; d <= j - MIN_HAIRPIN_SIZE - 2; d++)
			for(int e = d + MIN_HAIRPIN_SIZE + 1; e <= j-1; e++)
			{
				if (e == j-1 && d == i+1)
					temp += Qb->element(d, e)*exp(-energy->Es(i, j, seq)/RT);
				else 
					temp += Qb->element(d, e)*exp(-energy->Ebi(i, j, d, e, seq)/RT);
			}

		*Qb->estar(i, j) = temp;
	}


	double temp1 = 0.0, temp2 = 0.0;

	for(int d = i; d <= j - MIN_HAIRPIN_SIZE - 1; d++)
		for(int e = d + MIN_HAIRPIN_SIZE + 1; e <= j; e++)
		{
			if (d > i)
				temp1 += Qm->element(i, d-1)*Qb->element(d, e)*exp(-energy->Emulti(0, j-e, 1)/RT);

			temp2 += Qb->element(d, e)*exp(-energy->Emulti(0, j-i-(e-d), 1)/RT); 
		}

	*(Qm->estar(i, j)) = temp1 + temp2;
	*(Qm2->estar(i, j)) = temp1;
}

void JointProb::computeQ(int i, int j)
{
	*(Q->estar(i, j)) = 1.0;

	for(int d = i; d <= j - MIN_HAIRPIN_SIZE - 1; d++)
		for(int e = d + MIN_HAIRPIN_SIZE + 1; e <= j; e++)
		{
			double temp1 = 1.0;
			if (d > i)
				temp1 = Q->element(i, d-1);

			*(Q->estar(i, j)) += temp1*(Qb->element(d, e));
		}
}

void JointProb::computeProb(bool show)
{
	if(!quiet && show) 
		printf("Computing probabilities for %s whose length is %d... \n", sequence->getName(), sequence->getLen());

	unsigned int len = sequence->getLen();

	*P->estar(0, len - 1) = 1.0;

	for(unsigned int l = len; l >= 1; l--)
	{
		for(int i = 0; (unsigned int)i <= len - l; i++)
		{
			int j = i + l - 1;
			if(DEBUG) printf("i:%d\tj:%d\n", i, j);
			processP(i, j);
			processPbPm(i, j);

			if(l <= sublen)
			{			
				double temp1 = 1.0, temp2 = 1.0, temp;

				if(i > 0)
					temp1 = Q->element(0, i-1);	

				if(j < len-1)
					temp2 = Q->element(j+1, len-1);
	
				temp = temp1*temp2 / Q->element(0, len-1);

				*Pfree->estar(i, j) += temp;
			}
		}
		
		if(!quiet && show && !(l % 10)) printf("Length %d done.\n", l);		
	}

	filter(bpCollection, bp, Pb, ZERO);
	filter(hpCollection, hp, Phairpin, ZERO);
	filter(frCollection, fr, Pfree, ZERO);
	
	if(!quiet && show) 
		printf("Computing probabilities finished. \n");
}


void JointProb::processPbPm(int i, int j)
{
	unsigned char *seq = sequence->getSeq();

	double qm = Qm->element(i, j);
	double pm = Pm->element(i, j);
	double qm2 = Qm2->element(i, j);
	double pm2 = Pm2->element(i, j);

	if(pm == 0.0 || qm == 0.0)
	{
		pm = 0.0;
		qm = 1.0;
	}

	if(pm2 == 0.0 || qm2 == 0.0)
	{
		pm2 = 0.0;
		qm2 = 1.0;
	}

	for(int d = i; d <= j - MIN_HAIRPIN_SIZE - 1; d++)
		for(int e = d + MIN_HAIRPIN_SIZE + 1; e <= j; e++)
		{
			double deltaP;
			if (d > i)
			{
   				deltaP = ((pm / qm) + (pm2 / qm2)) * Qm->element(i, d-1)*Qb->element(d, e)*exp(-energy->Emulti(0, j-e, 1)/RT);
				*Pb->estar(d, e) += deltaP;
				*Pm->estar(i, d-1) += deltaP;
			}

			deltaP = pm * Qb->element(d, e)*exp(-energy->Emulti(0, j-i-(e-d), 1)/RT) / qm;
			*Pb->estar(d, e) += deltaP;
		}


	if(j > i + MIN_HAIRPIN_SIZE && energy->base_pair(seq[i], seq[j]) && Pb->element(i, j) != 0.0 && Qb->element(i, j) != 0.0)
	{
		double dp = Pb->element(i, j) * exp(-energy->Eh(i, j, seq)/RT) / Qb->element(i, j);

		*Phairpin->estar(i, j) += dp;

		for(int k = i+1; k <= j-1; k++)
			for(int l = k; l <= j-1 && l <= k + sublen-1; l++)
			{
				double deltaP = dp;
				if(k > i+1) 
					deltaP += Pb->element(i, j) * Qm2->element(i+1, k-1) * exp(-energy->Emulti(1, j-k, 0)/RT) / Qb->element(i, j); 
				if(l < j-1) 
					deltaP += Pb->element(i, j) * Qm2->element(l+1, j-1) * exp(-energy->Emulti(1, l-i, 0)/RT) / Qb->element(i, j); 
				if(k > i+1 && l < j-1) 
					deltaP += Pb->element(i, j) * Qm->element(i+1, k-1) * Qm->element(l+1, j-1) * exp(-energy->Emulti(1, l-k+1, 0)/RT) / Qb->element(i, j); 
				*Pfree->estar(k, l) += deltaP;
			}

		double dp2 = Pb->element(i, j) * Qm2->element(i+1, j-1)*exp(-energy->Emulti(1, 0, 0)/RT) / Qb->element(i, j);
		*Pm2->estar(i+1, j-1) += dp2;

		for(int d = i+1; d <= j-MIN_HAIRPIN_SIZE-2; d++)
			for(int e = d+MIN_HAIRPIN_SIZE+1; e <= j-1; e++)
			{
				if (e == j-1 && d == i+1)
				{
					double dp3 = Pb->element(i, j) * Qb->element(d, e)*exp(-energy->Es(i, j, seq)/RT) / Qb->element(i, j);
					*Pb->estar(d, e) += dp3;
				}
				else
				{
					double deltaP = Pb->element(i, j) * Qb->element(d, e)*exp(-energy->Ebi(i, j, d, e, seq)/RT) / Qb->element(i, j);
					*Pb->estar(d, e) += deltaP;

					if(deltaP != 0.0)
					{
						for(int k = i+1; k <= d-1; k++)
							for(int l = k; l <= d-1 && l <= k + sublen-1; l++)
								*Pfree->estar(k, l) += deltaP;  
					
						for(int k = e+1; k <= j-1; k++)
							for(int l = k; l <= j-1 && l <= k + sublen-1; l++)
								*Pfree->estar(k, l) += deltaP;  
					}
				}
			}
	}
}

void JointProb::processP(int i, int j)
{
	for(int d = i; d <= j-MIN_HAIRPIN_SIZE-1; d++)
		for(int e = d+MIN_HAIRPIN_SIZE+1; e <= j; e++)
		{
			double temp1 = 1.0;
			if (d > i)
				temp1 = Q->element(i, d-1);

			double deltaP0 = P->element(i, j) * temp1*Qb->element(d, e) / Q->element(i, j);

			*Pb->estar(d, e) += deltaP0;
			if (d > i)
				*P->estar(i, d-1) += deltaP0;
		}
}

void JointProb::buildBayesianNet()
{
	if(sitesCollection == NULL)
		return;

	if(!quiet) printf("Computing conditional probabilities for pairs of regions and building the graphical model for %s whose length is %d... \n", sequence->getName(), sequence->getLen());
	if(!quiet) printf("Computing for %d regions... \n", sitesCollection->size());
	char d2 = 2;
	jointfree = new Table<double>(sitesCollection->size(), &d2, 0.0);
	mutualinfo = new Table<double>(sitesCollection->size(), &d2, 0.0);

	int count = 0;

#pragma omp parallel for num_threads(procNum)
	for(int k = 0; k < sitesCollection->size(); k++)
	{
		Region<double> *r = new Region<double>(sitesCollection->item(k));
		JointProb *jp = new JointProb(this, r);
/*		forcedup->add(r);
		Q->reset();
		Qm->reset();
		Qm2->reset();
		Qb->reset();
		P->reset();
		Pm->reset();
		Pm2->reset();
		Pb->reset();
		Pfree->reset();
		Phairpin->reset();

		computeParFunc(false);
		computeProb(false); */
		for(int l = k+1; l < sitesCollection->size(); l++)
		{
			*jointfree->estar(k, l) = freeProb(r->i(), r->j()) * jp->freeProb(sitesCollection->item(l)->i(), sitesCollection->item(l)->j());
			*mutualinfo->estar(k, l) = mutualInfo(k, l);
		}

#pragma omp atomic
		count++;
		if(!quiet) printf("Region %d out of %d regions done.\n", count, sitesCollection->size());
//		forcedup->flush();
		delete jp;
	}

	if(!quiet) printf("Building Markov tree...\n");
	markovTree = new MTree(sitesCollection, mutualinfo, jointfree);
	if(!quiet) printf("Markov tree was built.\n");

	int list[3];
	list[0] = 2;

	double sum = 0.0;
	double max = 0.0;
	double maxpercent = 0.0;
	count = 0;
	for(int k = 0; k < sitesCollection->size(); k++)
		for(int l = k+1; l < sitesCollection->size(); l++)
		{
			count++;
			list[1] = k;
			list[2] = l;
			double temp = RT*fabs(log(jointFreeProb(k, l))-log(markovTree->inference((int *)&list)));
			if(DEBUG)
				if(markovTree->isEdge(k, l))
					printf("> %d\t%d\n", k, l);

			sum += temp;
			if(temp > max) max = temp;
			if(-temp/(RT*log(jointFreeProb(k, l))) > maxpercent) maxpercent = -temp/(RT*log(jointFreeProb(k, l)));
		}

	printf("error: avg %g, max %g, max%% %g\n", sum/count, max, maxpercent*100);

	if(!quiet) printf("Computing conditional probabilities and building the graphical model finished. \n");
}

void JointProb::reportProbs(FILE *outp)
{
	fprintf(outp, "# Base pairs\n# i\tj\tP(i,j)\n");
	for(int i = 0; i < sequence->getLen(); i++)
		for(int j = i; j < sequence->getLen(); j++)
			if(bp->element(i, j) > ZERO) fprintf(outp, "%d\t%d\t%g\n", i+1, j+1, bp->element(i, j));

	fprintf(outp, "# Hairpins\n# i\tj\tP(i,j)\n");
	for(int i = 0; i < sequence->getLen(); i++)
		for(int j = i; j < sequence->getLen(); j++)
			if(hp->element(i, j) > ZERO) fprintf(outp, "%d\t%d\t%g\n", i+1, j+1, hp->element(i, j));

	fprintf(outp, "# Free regions\n# i\tj\tP(i,j)\t%d\n", frCollection->size());
	for(int k = 0; k < frCollection->size(); k++)
	{
		int i = frCollection->item(k)->i();
		int j = frCollection->item(k)->j();
		fprintf(outp, "%d\t%d\t%g\n", i+1, j+1, fr->element(i, j));
	}

	if(DEBUG)	
		for(int k = 0; k < frCollection->size(); k++)
		{
			int i = frCollection->item(k)->i();
			int j = frCollection->item(k)->j();
			if(i == j)
			{
				double sum = 0.0;
				for(int d = 0; d < sequence->getLen(); d++)
					if(d < i)
						sum += bp->element(d, i);
					else
						sum += bp->element(i, d);

				if(fabs(fr->element(i, j) - 1 + sum) > 1e-6)		
					printf("%d\t%g\t%g\t%g\n", i+1, fr->element(i, j), 1 - sum, -fr->element(i, j) + 1 - sum);
				else
					printf("%d\t%g\n", i+1, fr->element(i, j));
			}
		} 
}


void JointProb::filter(RegionCollection *regions, Table<double> *dest, Table<double> *p, double threshold)
{
	regions->flush();

	for(int i = 0; i < sequence->getLen(); i++)
		for(int j = i; j < sequence->getLen() && j <= i + sublen - 1; j++)
		{
			*dest->estar(i, j) = p->element(i, j);
			if(p->element(i, j) > threshold)
			{
				Region<double> *newr = new Region<double>(i, j);
				newr->addAttrib(p->element(i, j));
				regions->add(newr);
			}
		}

	regions->sort(0, false);
}

double JointProb::basepairProb(int i, int j)
{
	if(i < 0 || i >= sequence->getLen() || j < 0 || j >= sequence->getLen() || i >= j) 
		return 0.0;

	return bp->element(i, j);
}

double JointProb::hairpinProb(int i, int j)
{
	if(i < 0 || i >= sequence->getLen() || j < 0 || j >= sequence->getLen() || i >= j) 
		return 0.0;

	return hp->element(i, j);
}

double JointProb::freeProb(int i, int j)
{
	if(i < 0 || i >= sequence->getLen() || j < 0 || j >= sequence->getLen() || i > j) 
		return 0.0;

	return fr->element(i, j);
}

double JointProb::basepairEnergy(int i, int j)
{
	if(i < 0 || i >= sequence->getLen() || j < 0 || j >= sequence->getLen() || i >= j) 
		return 1.0/0.0;

	return -RT*log(bp->element(i, j)*Q->element(0, sequence->getLen()-1));
}

double JointProb::hairpinEnergy(int i, int j)
{
	if(i < 0 || i >= sequence->getLen() || j < 0 || j >= sequence->getLen() || i >= j) 
		return 0.0;

	return -RT*log(hp->element(i, j)*Q->element(0, sequence->getLen()-1));
}

double JointProb::freeEnergy(int i, int j)
{
	if(i < 0 || i >= sequence->getLen() || j < 0 || j >= sequence->getLen() || i > j) 
		return 0.0;

	return -RT*log(fr->element(i, j));
}

double JointProb::freeProb(int i)
{
	if(i < 0 || i >= sitesCollection->size()) 
		return 0.0;

	return fr->element(sitesCollection->item(i)->i(), sitesCollection->item(i)->j());
}

double JointProb::jointFreeProb(int i, int j)
{
	if(i < 0 || i >= sitesCollection->size() || j < 0 || j >= sitesCollection->size()) 
		return 0.0;

	return (i <= j) ? jointfree->element(i, j) : jointfree->element(j, i);
}

double JointProb::jointFreeEnergy(int i, int j)
{
	if(i < 0 || i >= sitesCollection->size() || j < 0 || j >= sitesCollection->size() || i > j) 
		return 0.0;

	return -RT*log(jointfree->element(i, j));
}

double JointProb::mutualInfo(int i, int j)
{
	if(i < 0 || i >= sitesCollection->size() || j < 0 || j >= sitesCollection->size()) 
		return 0.0;

	double a = freeProb(i);
	double b = freeProb(j);
	double na = 1 - a;
	double nb = 1 - b;
	double ab = jointFreeProb(i, j);
	double anb = a - ab;
	double nab = b - ab;
	double nanb = nb - anb;

	if(anb < 0.0)
		anb = 0.0;
	if(nab < 0.0)
		nab = 0.0;
	if(nanb < 0.0)
		nanb = 0.0;

	double ablogab = finite(log(ab / (a * b))) ? ab*log(ab / (a * b)) : 0.0;
	double anbloganb = finite(log(anb / (a * nb))) ? anb*log(anb / (a * nb)) : 0.0;
	double nablognab = finite(log(nab / (na * b))) ? nab*log(nab / (na * b)) : 0.0;
	double nanblognanb = finite(log(nanb / (na * nb))) ? nanb*log(nanb / (na * nb)) : 0.0;
	double ret = ablogab + anbloganb + nablognab + nanblognanb;

	return (ret > 0.0) ? ret : 0.0;
}

void JointProb::setSites(RegionCollection *s)
{
	sitesCollection = s;
	for(int i=0; i < sitesCollection->size(); i++)
		sitesCollection->item(i)->addAttrib(freeProb(i));

	buildBayesianNet();
}

