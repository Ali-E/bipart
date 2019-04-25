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
	Desc: Interaction class drives the whole program. main() is here. This is a program to predict the interaction site of an sRNA with a target RNA.

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow, 
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: April 25, 2009
	Ver: 1.6
*/

#include <math.h>
#include <time.h>

#include "birna.h"

#define INT_THRESH 0.1

long int total_tables_size;
int big_tables_num;

int main(int argc, char** argv)
{
	Interaction inter(argc, argv);
	inter.computeBindingSites();
}

Interaction::Interaction(int argc, char** argv)
{
	opts = new GetOpt(argc, argv, OPTIONS);
	total_tables_size = 0;
	big_tables_num = 0;

	NA = 0;
	naConc = 1.0;
	mgConc = 0.0;
	polymer = 0;
	nodangle = 0;
	tMin = 37;
	tInc = 1;
	tMax = 37;
	suffix = NULL;
	zip = 0;
	debug = false;
	quiet = false;
	procNum = 16;
	files = 0;
	no_isolate = true;
	energyType = 1;
	seq_num = 0;
	max_sublen = 25;
	maxCandidates = 1000;
	dump_sites = true;


	while (opts->hasNext())
	{
		Option *current = opts->next();
		char count = current->getShortForm();

		if (count == FREE_ARG)
    			filenames[files++] = current->getArg();
		else if (count == 2)
    			++nodangle;
		else if (count == 3)
    			dump_sites = false;
		else if (count == 10)
    			no_isolate = false;
      		else if (count == 'V')
			version("birna");
      		else if (count == 'h')
		{
			printf("Usage: birna [options] filename1 filename2\n");
			printf("%s\n", opts->help());
			exit(0);
		}
      		else if (count == 'n')
		{
	  		if (!strcmp(current->getArg(), "RNA"))
	    			NA = 0;
	  		else if (!strcmp(current->getArg(), "DNA"))
	    			NA = 1;
		}
      		else if (count == 'e')
		{
	  		if (!strcmp(current->getArg(), "|sig|"))
	    			energyType = 1;
	  		else if (!strcmp(current->getArg(), "sig"))
	    			energyType = 0;
		}
      		else if (count == 't')
			tMin = atof(current->getArg());
      		else if (count == 'i')
			tInc = atof(current->getArg());
      		else if (count == 'T')
			tMax = atof(current->getArg());
      		else if (count == 's')
			suffix = current->getArg();
      		else if (count == 'd')
			debug = true;
      		else if (count == 'q')
			quiet = true;
      		else if (count == 'N')
			naConc = atof(current->getArg());
      		else if (count == 'M')
			mgConc = atof(current->getArg());
      		else if (count == 'm')
			maxCandidates = atoi(current->getArg());
      		else if (count == 'w')
			max_sublen = atoi(current->getArg());
      		else if (count == 'P')
			++polymer;
      		else if (count == 'p')
			procNum = atoi(current->getArg());
      		else if (count == 'z')
			++zip;
	}

	if (NA == 0 && (naConc != 1.0 || mgConc != 0.0 || polymer))
		fputs("Warning: salt concentrations ignored for RNA\n", stderr);

	if (suffix && (naConc != 1 || mgConc != 0 || polymer))
		fputs("Warning: salt concentrations ignored with suffix\n", stderr);

	if (!suffix && tMin > tMax)
	{
		fputs("Error: tMax must be greater than or equal to tMin.\n", stderr);
		exit(1);
	}

	if (tMin + tInc == tMin)
	{
		fputs("Error: tInc is too small compared to tMin.\n", stderr);
		exit(1);
	}

	if(strrchr(argv[0],'/')) 
	{
		*strrchr(argv[0],'/') = 0;
		strcpy(dataDir, argv[0]);
		strcat(dataDir, "/data/");
	}
	else 
		strcpy(dataDir, "../../data/");

	energy = new Energy(dataDir, NA, polymer, naConc, mgConc, suffix, tMin, zip, nodangle, energyType);

	if (files < 1)
	{
		fputs("Error: data not specified, try 'birna -h' for help.\n", stderr);
		exit(1);
	}
	
	if (files == 2)
	{
		seq[0] = new Sequence(filenames[0]);
		seq[1] = new Sequence(filenames[1]);
		seq_num = 2;
		if(!quiet) printf("Sequences %s and %s are read.\n", seq[0]->getName(), seq[1]->getName());
	} else
	{
		seq_num = 0;
		bool more_seq = true;
		FILE *fp = fopen(filenames[0], "rt");
		if(!fp)
		{
			perror("Input file error:");
			exit(1);
		}
		while(more_seq)
		{
			Sequence *motherseq = seq[seq_num] = new Sequence(fp, filenames[0]);
			more_seq = motherseq->isLoaded();
			if(more_seq)
			{
				seq_num++;
				while((seq[seq_num] = motherseq->split())) seq_num++;
			}
		}
		fclose(fp);
		if(seq_num % 2)
		{
			fputs("Error: odd number of sequences. Please provide pairs.\n", stderr);
			exit(1);
		}
	}
}

Interaction::~Interaction()
{
	delete opts;
	delete energy;
	for(int i = 0; i < seq_num; i++)
		delete seq[i];
}

void Interaction::computeBindingSites()
{
	char *buffer = (char *)alloc.xmalloc(10000);
	FILE *logfile = NULL, *sitesFileAB = NULL;

	for(int i = 0; i < seq_num; i += 2)
	{
		time_t now = time(NULL);
		if(dump_sites)
		{		
			fname2(&buffer, seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1);
			strcat(buffer, ".z");
			if (!(sitesFileAB = fopen(buffer, "wt")))
			{
				perror(buffer);
				exit(EXIT_FAILURE);
			}
			fputs("#T\tdG\tdGI\tdGu\tdGd\tdGfu\tdGfd\ti1\tj1\ti2\tj2\n", sitesFileAB);
		}

		fname2(&buffer, seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1);
		strcat(buffer, ".z.run");
		if (!(logfile = fopen(buffer, "wt")))
		{
			perror(buffer);
			exit(EXIT_FAILURE);
		}
		fprintf(logfile, "birna %s ran on %s and %s at %s\n", PACKAGE_VERSION, seq[i]->getFileName(), seq[i+1]->getFileName(), ctime(&now));
		if (suffix)
			fprintf(logfile, "suffix = %s\n", suffix);
		else
		{
			fprintf(logfile, "NA = %s\n", NA ? "DNA" : "RNA");
			fprintf(logfile, "tMin = %g\n", tMin);
			fprintf(logfile, "tInc = %g\n", tInc);
			fprintf(logfile, "tMax = %g\n", tMax);
			fprintf(logfile, "[Na+] = %g\n", naConc);
			fprintf(logfile, "[Mg++] = %g\n", mgConc);
		}
		if (nodangle)
			fputs("no dangle\n", logfile);
		if (polymer)
			fputs("polymer mode\n", logfile);
		if (!no_isolate)
			fputs("isolates allowed\n", logfile);
		fprintf(logfile, "energy type = %d\n", energyType);
		fprintf(logfile, "no. cpus = %d\n", procNum);


		PartitionFunction *ab[2];
		JointProb *abProb[2];

		unsigned int sublen[2] = {seq[i]->getLen(), seq[i+1]->getLen()};

		if(sublen[0] > max_sublen) sublen[0] = max_sublen;
		if(sublen[1] > max_sublen) sublen[1] = max_sublen;

		for (double t = tMin; t <= tMax; t += tInc)
		{
			RegionCollection *hpIntRegions = new RegionCollection(maxCandidates);
			hpIntRegions->sort(0, false);
			energy->setTemperature(t);
			if(!quiet) printf("Calculating for %s and %s, temperature = %lf\n", seq[i]->getName(), seq[i+1]->getName(), t);
	
			for(int j = 0; j < 2; j++)
			{
				ab[j] = new PartitionFunction(energy, seq[i+j], sublen[j], debug, procNum, quiet, no_isolate);
				abProb[j] = new JointProb(energy, seq[i+j], sublen[j], debug, procNum, quiet, no_isolate);
			}

			unsigned int l1 = seq[i]->getLen();
			unsigned int l2 = seq[i+1]->getLen();
			PartitionFunction *duplex = new PartitionFunction(energy, ab, debug, procNum, quiet, no_isolate);

			filterHighProbIntRegions(hpIntRegions, ab, duplex, abProb, l1, l2, ab[0]->getSubLen(), ab[1]->getSubLen());

			if(dump_sites) 
			{
				hpIntRegions->sort(1, true);
				outputSites(hpIntRegions, sitesFileAB, t);
			}
			fflush(NULL);

			RegionCollection *sites[2];
			for(int s = 0; s < 2; s++)
			{
				sites[s] = hpIntRegions->component(s);
				abProb[s]->setSites(sites[s]);
				if(dump_sites) 
				{
					fprintf(sitesFileAB, "#Free sites for %s\n", seq[i+s]->getName());
					for(int k = 0; k < sites[s]->size(); k++)
						fprintf(sitesFileAB, "%d\t%d\t%d\t%g\n", k+1, sites[s]->item(k)->i()+1, sites[s]->item(k)->j()+1, sites[s]->item(k)->getAttrib(0));

					fputs("#Markov tree\n", sitesFileAB);
					abProb[s]->dumpTree(sitesFileAB);

					fputs("#Joint probs and mutual info\n", sitesFileAB);
					for(int k = 0; k < sites[s]->size(); k++)
						for(int l = k+1; l < sites[s]->size(); l++)
							fprintf(sitesFileAB, "%d\t%d\t%g\t%g\n", k+1, l+1, abProb[s]->jointFreeProb(k, l), abProb[s]->mutualInfo(k, l));

				}
				
			}

			if(dump_sites)
			{
				computeBestMatches(sitesFileAB, sites, ab, duplex, abProb, l2);
				computeBestMatches(sitesFileAB, 3, sites, ab, duplex, abProb, l2);
			} else
			{
				computeBestMatches(stdout, sites, ab, duplex, abProb, l2);
				computeBestMatches(stdout, 3, sites, ab, duplex, abProb, l2);
			}


			delete hpIntRegions;
			delete duplex;
			for(int j = 0; j < 2; j++)
			{
				delete ab[j];
				delete abProb[j];
				delete sites[j];
			}
		}
 
		if(!quiet) printf("Running time: %ld seconds.\n", time(NULL) - now);
		fprintf(logfile, "Running time: %ld seconds.\n", time(NULL) - now);
		if(dump_sites) fclose(sitesFileAB);
		fclose(logfile);
	}
	free(buffer);
}

#define SCALE 1.3

void Interaction::filterHighProbIntRegions(RegionCollection *hpIntRegions, PartitionFunction *ab[2], PartitionFunction *duplex, JointProb *abProb[2], int l1, int l2, int sub1, int sub2)
{
	for(int i1 = 0; i1 <= l1 - sub1; i1++)
		for(int j1 = i1; j1 <= i1 + sub1 - 1; j1++)
			for(int i2 = 0; i2 <= l2 - sub2; i2++)
				for(int j2 = i2; j2 <= i2 + sub2 - 1; j2++)
				{
					double RT = energy->getRT();
					double qi = duplex->getQI()->element(i1, j1, l2 - j2 -1, l2 - i2 -1);
					double qu = ab[0]->getQ(0)->element(i1, j1);
					double qd = ab[1]->getQ(0)->element(i2, j2);
					double ei = -RT*log(qi)*SCALE;
					double eu = -RT*log(qu)*SCALE;
					double ed = -RT*log(qd)*SCALE;
					double efu = abProb[0]->freeEnergy(i1, j1);
					double efd = abProb[1]->freeEnergy(i2, j2);
					double pfu = abProb[0]->freeProb(i1, j1);
					double pfd = abProb[1]->freeProb(i2, j2);
/*					r->addAttrib(-ei + (eu+ed));
					r->addAttrib(-RT*log(pfu) + -RT*log(pfd) + ei - (eu + ed));*/
					if((-efu -efd - ei + (eu + ed)) / (j1-i1+1 + j2-i2+1) > INT_THRESH)
					{
						Region<double> *r = new Region<double>(i1, j1, i2, j2);
						r->addAttrib(-efu -efd - ei + (eu + ed));
						r->addAttrib(ei);
						r->addAttrib(eu);
						r->addAttrib(ed);
						r->addAttrib(efu);
						r->addAttrib(efd);
						if(!hpIntRegions->sortadd(r)) delete r;
					}
				}
}

void Interaction::computeBestMatches(FILE *out, RegionCollection *sites[2], PartitionFunction *ab[2], PartitionFunction *duplex, JointProb *abProb[2], int l2)
{
	double min = 1e10;
	double RT = energy->getRT();

	int maxr[2] = {-1, -1}, maxs[2] = {-1, -1};
	int r[2], s[2];
	for(r[0] = 0; r[0] < sites[0]->size(); r[0]++)
		for(s[0] = r[0]+1; s[0] < sites[0]->size(); s[0]++)
			if(!sites[0]->item(r[0])->overlaps(sites[0]->item(s[0])))
			{
				for(s[1] = 0; s[1] < sites[1]->size(); s[1]++)
					for(r[1] = s[1]+1; r[1] < sites[1]->size(); r[1]++)
						if(!sites[1]->item(r[1])->overlaps(sites[1]->item(s[1])))
						{
							double intdG1 = 0.0, intdG2 = 0.0;
							intdG1 += 
-RT*log(duplex->getQI()->element(sites[0]->item(r[0])->i(), sites[0]->item(r[0])->j(), l2 - sites[1]->item(r[1])->j() -1, 
l2 - sites[1]->item(r[1])->i() -1))*SCALE;
							intdG2 += 
-RT*log(duplex->getQI()->element(sites[0]->item(s[0])->i(), sites[0]->item(s[0])->j(), l2 - sites[1]->item(s[1])->j() -1, 
l2 - sites[1]->item(s[1])->i() -1))*SCALE;
							intdG1 -= 
-RT*log(ab[0]->getQ(0)->element(sites[0]->item(r[0])->i(), sites[0]->item(r[0])->j()) * 
ab[1]->getQ(0)->element(sites[1]->item(r[1])->i(), sites[1]->item(r[1])->j()))*SCALE;

							intdG2 -= 
-RT*log(ab[0]->getQ(0)->element(sites[0]->item(s[0])->i(), sites[0]->item(s[0])->j()) * 
ab[1]->getQ(0)->element(sites[1]->item(s[1])->i(), sites[1]->item(s[1])->j()))*SCALE; 

							double efu = abProb[0]->jointFreeEnergy(r[0], s[0]);
							double efd = abProb[1]->jointFreeEnergy(s[1], r[1]);
							if(efu + efd + intdG1 + intdG2 < min)
							{
								min = efu + efd + intdG1 + intdG2;
								maxr[0] = r[0];
								maxr[1] = r[1];
								maxs[0] = s[0];
								maxs[1] = s[1];
							}

							intdG1 = intdG2 = 0.0;
							intdG1 += 
-RT*log(duplex->getQI()->element(sites[0]->item(r[0])->i(), sites[0]->item(r[0])->j(), l2 - sites[1]->item(s[1])->j() -1, 
l2 - sites[1]->item(s[1])->i() -1))*SCALE;
							intdG2 += 
-RT*log(duplex->getQI()->element(sites[0]->item(s[0])->i(), sites[0]->item(s[0])->j(), l2 - sites[1]->item(r[1])->j() -1, 
l2 - sites[1]->item(r[1])->i() -1))*SCALE;
							intdG1 -= 
-RT*log(ab[0]->getQ(0)->element(sites[0]->item(r[0])->i(), sites[0]->item(r[0])->j()) * 
ab[1]->getQ(0)->element(sites[1]->item(s[1])->i(), sites[1]->item(s[1])->j()))*SCALE;
							intdG2 -= 
-RT*log(ab[0]->getQ(0)->element(sites[0]->item(s[0])->i(), sites[0]->item(s[0])->j()) * 
ab[1]->getQ(0)->element(sites[1]->item(r[1])->i(), sites[1]->item(r[1])->j()))*SCALE; 

							if(efu + efd + intdG1 + intdG2 < min)
							{
								min = efu + efd + intdG1 + intdG2;
								maxr[0] = r[0];
								maxr[1] = s[1];
								maxs[0] = s[0];
								maxs[1] = r[1];
							}
						}				
			}
	if(min < 1e10)
	{
		fprintf(out, "#Double binding sites\ndG = %g\t", min);
		for(int idx = 0; idx < 2; idx++)
		{
			if(maxr[idx] >= 0 && maxs[idx] >= 0)
			{
				fprintf(out, "%d\t%d\t%d\t%d", sites[idx]->item(maxr[idx])->i()+1, sites[idx]->item(maxr[idx])->j()+1, sites[idx]->item(maxs[idx])->i()+1, sites[idx]->item(maxs[idx])->j()+1);
				if(idx)
					fprintf(out, "\n");
				else		
					fprintf(out, ",\t");
			}
		}
	}
}

void Interaction::computeBestMatches(FILE* out, int sites_num, RegionCollection *sites[2], PartitionFunction *ab[2], PartitionFunction *duplex, JointProb *abProb[2], int l2)
{
	double min = 1e10;
	double RT = energy->getRT();

	int *bestMatch[2];

	for(int upOrDown = 0; upOrDown < 2; upOrDown++)
	{
		bestMatch[upOrDown] = (int *)alloc.xmalloc(sizeof(int)*sites_num);
		for(int i = 0; i < sites_num; i++)
			bestMatch[upOrDown][i] = -1;
	}

	int *s[2];

	for(int upOrDown = 0; upOrDown < 2; upOrDown++)
	{
		s[upOrDown] = (int *)alloc.xmalloc(sizeof(int)*sites_num);
		for(int i = 0; i < sites_num; i++)
			s[upOrDown][i] = 0;
	}
	
	int *assignment = (int *)alloc.xmalloc(sizeof(int)*sites_num);

	double *weights = (double *)alloc.xmalloc(sizeof(double)*sites_num*sites_num);

	for(s[0][0] = 0; s[0][0] < sites[0]->size(); s[0][0]++)
	for(s[0][1] = s[0][0]+1; s[0][1] < sites[0]->size(); s[0][1]++)
	for(s[0][2] = s[0][1]+1; s[0][2] < sites[0]->size(); s[0][2]++)
		if(!sites[0]->item(s[0][0])->overlaps(sites[0]->item(s[0][1])) && !sites[0]->item(s[0][0])->overlaps(sites[0]->item(s[0][2])) && !sites[0]->item(s[0][1])->overlaps(sites[0]->item(s[0][2])))
			for(s[1][0] = 0; s[1][0] < sites[1]->size(); s[1][0]++)
			for(s[1][1] = s[1][0]+1; s[1][1] < sites[1]->size(); s[1][1]++)
			for(s[1][2] = s[1][1]+1; s[1][2] < sites[1]->size(); s[1][2]++)
				if(!sites[1]->item(s[1][0])->overlaps(sites[1]->item(s[1][1])) && !sites[1]->item(s[1][1])->overlaps(sites[1]->item(s[1][2])) && !sites[1]->item(s[1][0])->overlaps(sites[1]->item(s[1][2])))
				{
					int list[2][4];
					list[0][0] = list[1][0] = 3;
					for(int k = 0; k < sites_num; k++)
					{
						list[0][k+1] = s[0][k];
						list[1][k+1] = s[1][k];
					}
					double efu = abProb[0]->jointFreeEnergy(list[0]);
					double efd = abProb[1]->jointFreeEnergy(list[1]);

					for(int k = 0; k < sites_num; k++)
						for(int l = 0; l < sites_num; l++)
						{
							Region<double> *upr = sites[0]->item(s[0][k]), *downr = sites[1]->item(s[1][l]);
							weights[k*sites_num + l] = -RT*log(duplex->getQI()->element(upr->i(), upr->j(), l2 - downr->j() -1, 
l2 - downr->i() -1))*SCALE + RT*log(ab[0]->getQ(0)->element(upr->i(), upr->j()) * 
ab[1]->getQ(0)->element(downr->i(), downr->j()))*SCALE;
						}					

					double intdG = match(sites_num, weights, assignment);

					if(efu + efd + intdG < min)
					{
						min = efu + efd + intdG;
						for(int k = 0; k < sites_num; k++)
						{
							bestMatch[0][k] = s[0][k];
							bestMatch[1][k] = s[1][assignment[k]];
						}
					}
				}				
	if(min < 1e10)
	{
		fprintf(out, "#Triple binding sites\ndG = %g\t", min);
		for(int idx = 0; idx < 2; idx++)
		{
			if(bestMatch[idx][0] >= 0 && bestMatch[idx][1] >= 0 && bestMatch[idx][2] >= 0)
			{
				fprintf(out, "%d\t%d\t%d\t%d\t%d\t%d", sites[idx]->item(bestMatch[idx][0])->i()+1, sites[idx]->item(bestMatch[idx][0])->j()+1, 
sites[idx]->item(bestMatch[idx][1])->i()+1, sites[idx]->item(bestMatch[idx][1])->j()+1,
sites[idx]->item(bestMatch[idx][2])->i()+1, sites[idx]->item(bestMatch[idx][2])->j()+1);
				if(idx)
					fprintf(out, "\n");
				else		
					fprintf(out, ",\t");
			}
		}
	}
	for(int upOrDown = 0; upOrDown < 2; upOrDown++)
	{
		free(bestMatch[upOrDown]);
		free(s[upOrDown]);
	}
	free(assignment);
	free(weights);
}

double Interaction::match(int sites_num, double *weights, int *assign)
{
	double min = 1e10;
	if(weights[0*sites_num + 2] + weights[1*sites_num + 0] + weights[2*sites_num + 1] < min)
	{
		min = weights[0*sites_num + 2] + weights[1*sites_num + 0] + weights[2*sites_num + 1];
		assign[0] = 2;
		assign[1] = 0;
		assign[2] = 1;
	}
	if(weights[0*sites_num + 2] + weights[1*sites_num + 1] + weights[2*sites_num + 0] < min)
	{
		min = weights[0*sites_num + 2] + weights[1*sites_num + 1] + weights[2*sites_num + 0];
		assign[0] = 2;
		assign[1] = 1;
		assign[2] = 0;
	}
	if(weights[0*sites_num + 0] + weights[1*sites_num + 1] + weights[2*sites_num + 2] < min)
	{
		min = weights[0*sites_num + 0] + weights[1*sites_num + 1] + weights[2*sites_num + 2];
		assign[0] = 0;
		assign[1] = 1;
		assign[2] = 2;
	}
	if(weights[0*sites_num + 0] + weights[1*sites_num + 2] + weights[2*sites_num + 1] < min)
	{
		min = weights[0*sites_num + 0] + weights[1*sites_num + 2] + weights[2*sites_num + 1];
		assign[0] = 0;
		assign[1] = 2;
		assign[2] = 1;
	}
	if(weights[0*sites_num + 1] + weights[1*sites_num + 0] + weights[2*sites_num + 2] < min)
	{
		min = weights[0*sites_num + 1] + weights[1*sites_num + 0] + weights[2*sites_num + 2];
		assign[0] = 1;
		assign[1] = 0;
		assign[2] = 2;
	}
	if(weights[0*sites_num + 1] + weights[1*sites_num + 2] + weights[2*sites_num + 0] < min)
	{
		min = weights[0*sites_num + 1] + weights[1*sites_num + 2] + weights[2*sites_num + 0];
		assign[0] = 1;
		assign[1] = 2;
		assign[2] = 0;
	}
	return min;
}

void Interaction::outputSites(RegionCollection *regions, FILE *outp, double Temp)
{
	for(int i = 0; i < regions->size(); i++)
	{
		Region<double> *r = regions->item(i);
		fprintf(outp, "%g\t", Temp);
		for(int j = 0; j < 6; j++)
			fprintf(outp, "%0.5lf\t", r->getAttrib(j));
		fprintf(outp, "%d\t%d\t%d\t%d\n", r->i1()+1, r->j1()+1, r->i2()+1, r->j2()+1);
	}
}

void Interaction::fname1(char **res, char *fn, int i)
{
	if(seq_num > 2)
		sprintf(*res, "%s%d", fn, i);
	else
		sprintf(*res, "%s", fn);
}

void Interaction::fname2(char **res, char *fn1, char *fn2, int i, int j)
{
	if(seq_num > 2)
		sprintf(*res, "%s%d-%s%d", fn1, i, fn2, j);
	else
		sprintf(*res, "%s-%s", fn1, fn2);
}


