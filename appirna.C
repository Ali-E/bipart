/*
Copyright 2012 Hamid Reza Chitsaz (chitsaz@wayne.edu)

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
	Desc: APPartition class drives the whole program. main() is here.

	Author: Hamid Reza Chitsaz, Elmirasadat Forouzmand
		Wayne State University
		Algorithmic Biology Lab

	Last Update by Hamid Reza Chitsaz: Sep 21, 2012
	Ver: 1.0
*/

#include <math.h>
#include <time.h>

#include "appirna.h"
#include "partitionfunction.h"

#define Euler_Mascheroni (0.577215664901532)
#define MAX_CONST_CNT 5
#define CONST_THRESH 0.5


long int total_tables_size;
int big_tables_num;

int main(int argc, char** argv)
{
	srand(time(NULL));

	APPartition p(argc, argv);
	p.computeUpperbound();
}

APPartition::APPartition(int argc, char** argv)
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
	store = false;
	retrieve = false;

	while (opts->hasNext())
	{
		Option *current = opts->next();
		char count = current->getShortForm();

		if (count == FREE_ARG)
    			filenames[files++] = current->getArg();
		else if (count == 2)
    			++nodangle;
		else if (count == 10)
    			no_isolate = false;
      		else if (count == 17)
		{
			store = true;
		}
      		else if (count == 18)
		{
			retrieve = true;
		}
      		else if (count == 'V')
			version("appirna");
      		else if (count == 'h')
		{
			printf("Usage: appirna [options] filename\n");
			printf("    or appirna [options] filename1 filename2\n");
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

	if (files < 1 || files > 2)
	{
		fputs("Error: data not specified, try 'appirna -h' for help.\n", stderr);
		exit(1);
	}
	
	if (files == 2)
	{
		seq.push_back(new Sequence(filenames[0]));
		seq.push_back(new Sequence(filenames[1]));
		if(!quiet) printf("Sequences %s and %s are read.\n", seq[0]->getName(), seq[1]->getName());
	} else if (files == 1)
	{
		bool more_seq = true;
		FILE *fp = fopen(filenames[0], "rt");
		if(!fp)
		{
			perror("Input file error:");
			exit(1);
		}
		while(more_seq)
		{
			Sequence *motherseq = new Sequence(fp, filenames[0]);
			more_seq = motherseq->isLoaded();
			if(more_seq) 
			{
				if(!quiet) printf("Sequence %s is read.\n", motherseq->getName());
				seq.push_back(motherseq); 
			}
			else 
				delete motherseq;
		}
		fclose(fp);
	}
}

APPartition::~APPartition()
{
	for(unsigned long int i = 0; i < seq.size(); i++)
		delete seq[i];

	delete opts;
	delete energy;
}

FILE * APPartition::openfile(char *fn, int i, char *ext, char *type)
{
	FILE *out;
	char res[10000];

	sprintf(res, "%s%s", fn, ext);
	if (!(out = fopen(res, type)))
	{
		perror(res);
		exit(EXIT_FAILURE);
	}
	return out;
}

FILE * APPartition::openfile(char *fn1, char *fn2, int i, int j, char *ext, char *type)
{
	FILE *out;
	char res[10000];

	sprintf(res, "%s-%s%s", fn1, fn2, ext);
	if (!(out = fopen(res, type)))
	{
		perror(res);
		exit(EXIT_FAILURE);
	}
	return out;
}

void APPartition::computeUpperbound()
{
	FILE *logfile = NULL, *dGFile = NULL;

	time_t now = time(NULL);

	dGFile = openfile(seq[0]->getFileName(), 0, (char *)".dG", (char *)"wt");
	logfile = openfile(seq[0]->getFileName(), 0, (char *)".run", (char *)"wt");
	fprintf(logfile, "appirna %s ran on %s at %s\n", PACKAGE_VERSION, seq[0]->getFileName(), ctime(&now));

//	fputs("#T\t-RT ln Z\tZ\tIterations\n", dGFile);

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


	char dim = 2;

	int count = 0;
	double avgenerr = 0.0, avgperr = 0.0, avgenper = 0.0, avgpper = 0.0, avgiter = 0.0;

	for(unsigned long int sqn = 0; sqn < seq.size(); sqn++)
	{
		unsigned int len = seq[sqn]->getLen();

		Table<double> *perturbation[2] = {new Table<double>(len, &dim, 0.0), new Table<double>(len, &dim, 0.0)};

		for (double t = tMin; t <= tMax; t += tInc)
		{
			energy->setTemperature(t);

			if(!quiet) printf("Calculating for %s, temperature = %lf\n", seq[sqn]->getName(), t);
	
			int iter, maxiter = len;

			double esum = 0.0, expectation = 0.0, parf;
			int constcnt = 0;
	
			for(iter = 0; iter < maxiter && constcnt < MAX_CONST_CNT; iter++)
			{
				sample(perturbation, len);
				UBPartitionFunction *app = new UBPartitionFunction(energy, seq[sqn], perturbation, debug, procNum, quiet, no_isolate);
				esum += app->freeEnsembleEnergy();
			
				if(fabs((esum / (iter+1)) - expectation) < CONST_THRESH)
					constcnt++;
				else
					constcnt = 0;

				expectation = esum / (iter+1);
				if(!quiet) printf("iteration: %i\tenergy: %g\texpectation: %g\n", iter+1,  app->freeEnsembleEnergy(), expectation);
				parf = app->parFunc(expectation);
				delete app;
			}

			PartitionFunction p(energy, seq[sqn], debug, procNum, quiet, no_isolate);

			avgiter += (iter+1.0)/len;
			avgenerr += fabs(expectation - p.freeEnsembleEnergy(1.0));
			avgperr += fabs(parf - p.Z(1.0));
			avgenper += fabs((expectation/p.freeEnsembleEnergy(1.0)) - 1.0);
			avgpper += fabs((parf/p.Z(1.0)) - 1.0);
			count++;

			fprintf(dGFile, "%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%s\n", t, expectation, parf, iter+1, p.freeEnsembleEnergy(1.0), p.Z(1.0), avgenerr, avgperr, avgenper, avgpper, avgiter, count, seq[sqn]->getName());
			fflush(NULL);
		}

		delete perturbation[0], perturbation[1];
	}

	printf("avg_iter: %g\tavg_en_err: %g\tavg_p_err: %g\tavg_en_percentage: %g\tavg_p_percentage: %g\n", (avgiter*1.0)/count, avgenerr/count, avgperr/count, avgenper*100.0/count, avgpper*100.0/count);

	if(!quiet) printf("Running time: %ld seconds.\n", time(NULL) - now);
	fprintf(logfile, "Running time: %ld seconds.\n", time(NULL) - now);
	fclose(dGFile);
	fclose(logfile);
}

void APPartition::sample(Table<double> *samp[2], unsigned int len)
{
	for(int sign = 0; sign < 2; sign++)
	{
		for(int i = 0; i < len; i++)
			for(int j = i; j < len; j++)
			{
				double u = (rand()*1.0) / RAND_MAX;
				(*samp[sign])(i, j) = -log(-log(u))-Euler_Mascheroni;
			}
	}
}

