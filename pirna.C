/*
Copyright 2008  Hamidreza Chitsaz (chitsaz@wayne.edu)

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
	Desc: Partition class drives the whole program. main() is here.

	Author: Hamidreza Chitsaz
		Wayne State University
		Algorithmic Biology Lab

	Last Update by Hamidreza Chitsaz: Oct 25, 2018
	Ver: 3.9
*/

#include <math.h>
#include <time.h>

#include "pirna.h"
#include "probability.h"

long int total_tables_size;
int big_tables_num;

int main(int argc, char** argv)
{
	Partition p(argc, argv);
	if(p.doPartition()) p.computeThermoValues();
	if(p.doConc()) p.computeConcentrations();
	if(p.doEns()) p.computedG(); 
	if(p.doTm()) p.computeTm();
}

Partition::Partition(int argc, char** argv)
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
	AA = true;
	BB = true;
	AB = true;
	no_isolate = true;
	energyType = 1;
	Na0 = Nb0 = 0.0;
	seq_num = 0;
	ensemble = false;
	tm = false;
	prob = false;
	storeP = false;
	retrieveP = false;
	storeQ = false;
	retrieveQ = false;
	end_time = time(NULL) + 84000;

	while (opts->hasNext())
	{
		Option *current = opts->next();
		char count = current->getShortForm();

		if (count == FREE_ARG)
    			filenames.push_back(current->getArg());
		else if (count == 2)
    			++nodangle;
		else if (count == 3)
    			AA = false;
		else if (count == 4)
    			BB = false;
		else if (count == 5)
    			AB = false;
		else if (count == 10)
    			no_isolate = false;
      		else if (count == 12)
			Na0 = atof(current->getArg());
      		else if (count == 13)
			Nb0 = atof(current->getArg());
      		else if (count == 14)
			ensemble = true;
      		else if (count == 15)
			tm = true;
      		else if (count == 16)
			prob = true;
      		else if (count == 17)
		{
			storeP = !strcmp(current->getArg(), "P") || !strcmp(current->getArg(), "PQ");
			storeQ = !strcmp(current->getArg(), "Q") || !strcmp(current->getArg(), "PQ");
		}
      		else if (count == 18)
		{
			retrieveP = !strcmp(current->getArg(), "P") || !strcmp(current->getArg(), "PQ");
			retrieveQ = !strcmp(current->getArg(), "Q") || !strcmp(current->getArg(), "PQ");
		}
      		else if (count == 'V')
			version("pirna");
      		else if (count == 'h')
		{
			printf("Usage: pirna [options] filename1 filename2\n");
			printf("    or pirna [options] filename\n");
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

	files = filenames.size();


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
		fputs("Error: data not specified, try 'pirna -h' for help.\n", stderr);
		exit(1);
	}
	
	if (files == 2 && Na0 == 0.0 && Nb0 == 0.0 && !ensemble && !tm)
	{
		seq.push_back(new Sequence(filenames[0]));
		seq.push_back(new Sequence(filenames[1]));
		seq_num = seq.size();
		if(!quiet) printf("Sequences %s and %s are read.\n", seq[0]->getName(), seq[1]->getName());
	} else if (Na0 == 0.0 && Nb0 == 0.0  && !ensemble && !tm)
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
			Sequence *motherseq = new Sequence(fp, filenames[0]);
			more_seq = motherseq->isLoaded();
			if(more_seq)
			{
				seq.push_back(motherseq);
				Sequence *s_temp;
				while((s_temp = motherseq->split())) seq.push_back(s_temp);
			}
		}
		seq_num = seq.size();

		fclose(fp);
		if(seq_num % 2)
		{
			fputs("Error: odd number of sequences. Please provide pairs.\n", stderr);
			exit(1);
		}
	}
}

FILE * Partition::openfile(char *fn, int i, char *ext, char *type)
{
	FILE *out;
	char res[10000];

	if(seq_num > 2)
		sprintf(res, "%s%d", fn, i);
	else
		sprintf(res, "%s", fn);

	strcat(res, ext);
	if (!(out = fopen(res, type)))
	{
		perror(res);
		exit(EXIT_FAILURE);
	}
	return out;
}

FILE * Partition::openfile(char *fn1, char *fn2, int i, int j, char *ext, char *type)
{
	FILE *out;
	char res[10000];

	if(seq_num > 2)
		sprintf(res, "%s%d-%s%d", fn1, i, fn2, j);
	else
		sprintf(res, "%s-%s", fn1, fn2);

	strcat(res, ext);
	if (!(out = fopen(res, type)))
	{
		perror(res);
		exit(EXIT_FAILURE);
	}
	return out;
}

void Partition::computeThermoValues()
{
	char *buffer = (char *)alloc.xmalloc(10000);
	char *buffer2 = (char *)alloc.xmalloc(10000);
	FILE *logfile = NULL, *dGFileA = NULL, *dGFileB = NULL, *dGFileAA = NULL, *dGFileBB = NULL, *dGFileAB = NULL;
	FILE *probFileAA = NULL, *probFileBB = NULL, *probFileAB = NULL;

	for(int i = 0; i < seq_num; i += 2)
	{
		time_t now = time(NULL);
      		dGFileA = openfile(seq[i]->getFileName(), i, (char *)".dG", (char *)"wt");
      		dGFileB = openfile(seq[i+1]->getFileName(), i+1, (char *)".dG", (char *)"wt");
		fputs("#T\t-RT ln Z\tZ\n", dGFileA);
		fputs("#T\t-RT ln Z\tZ\n", dGFileB);

		if(AA)
		{
			dGFileAA = openfile(seq[i]->getFileName(), seq[i]->getFileName(), i, i, (char *)".dG", (char *)"wt");
			fputs("#T\t-RT ln Z\tZ\n", dGFileAA);
			if(prob)
				probFileAA = openfile(seq[i]->getFileName(), seq[i]->getFileName(), i, i, (char *)".prob", (char *)"wt");
		}
		if(BB)
		{
			dGFileBB = openfile(seq[i+1]->getFileName(), seq[i+1]->getFileName(), i+1, i+1, (char *)".dG", (char *)"wt");
			fputs("#T\t-RT ln Z\tZ\n", dGFileBB);
			if(prob)
				probFileBB = openfile(seq[i+1]->getFileName(), seq[i+1]->getFileName(), i+1, i+1, (char *)".prob", (char *)"wt");
		}
		if(AB)
		{	
			dGFileAB = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, (char *)".dG", (char *)"wt");
			fputs("#T\t-RT ln Z\tZ\n", dGFileAB);
			if(prob)
				probFileAB = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, (char *)".prob", (char *)"wt");
		}

		logfile = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, (char *)".run", (char *)"wt");
		fprintf(logfile, "pirna %s ran on %s and %s at %s\n", PACKAGE_VERSION, seq[i]->getFileName(), seq[i+1]->getFileName(), ctime(&now));
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
		PartitionFunction *aa[2];
		PartitionFunction *bb[2];


		for (double t = tMin; t <= tMax; t += tInc)
		{
			energy->setTemperature(t);
			if(retrieveQ)
			{
				FILE *inpA, *inpB;
				sprintf(buffer, "a.tables.%g", t);
		      		inpA = openfile(seq[i]->getFileName(), i, buffer, (char *)"rb");
				sprintf(buffer, "b.tables.%g", t);
		      		inpB = openfile(seq[i+1]->getFileName(), i+1, buffer, (char *)"rb");
				ab[0] = new PartitionFunction(inpA, energy, seq[i], debug, procNum, quiet, no_isolate);
				ab[1] = new PartitionFunction(inpB, energy, seq[i+1], debug, procNum, quiet, no_isolate);
				fclose(inpA);
				fclose(inpB);
			} else
			{
				if(!quiet) printf("Calculating for %s and %s, temperature = %lf\n", seq[i]->getName(), seq[i+1]->getName(), t);
	
				for(int j = 0; j < 2; j++)
					ab[j] = new PartitionFunction(energy, seq[i+j], debug, procNum, quiet, no_isolate);

				if(storeQ)
				{
					FILE *outpA, *outpB;
					sprintf(buffer, "a.tables.%g", t);
			      		outpA = openfile(seq[i]->getFileName(), i, buffer, (char *)"wb");
					sprintf(buffer, "b.tables.%g", t);
			      		outpB = openfile(seq[i+1]->getFileName(), i+1, buffer, (char *)"wb");
					ab[0]->store(outpA);
					ab[1]->store(outpB);
					fclose(outpA);
					fclose(outpB);
				}
			}

			aa[0] = ab[0];
			aa[1] = new PartitionFunction(ab[0]);
			bb[0] = ab[1];
			bb[1] = new PartitionFunction(ab[1]);
			fprintf(dGFileA, "%g\t%g\t%g\n", t, ab[0]->freeEnsembleEnergy(1.0), ab[0]->Z(1.0));
			fflush(NULL);
			fprintf(dGFileB, "%g\t%g\t%g\n", t, ab[1]->freeEnsembleEnergy(1.0), ab[1]->Z(1.0));
			fflush(NULL);

			sprintf(buffer, ".tables.%g", t);
			if(AA)
			{
				PartitionFunction *duplex;
				if(retrieveQ)
				{
					FILE *inp;
					inp = openfile(seq[i]->getFileName(), seq[i]->getFileName(), i, i, buffer, (char *)"rb");
					duplex = new PartitionFunction(inp, energy, aa, debug, procNum, quiet, no_isolate);
					fclose(inp);
				}
				else
				{
					duplex = new PartitionFunction(energy, aa, debug, procNum, quiet, no_isolate);
					if(storeQ)
					{
						FILE *outp;
						outp = openfile(seq[i]->getFileName(), seq[i]->getFileName(), i, i, buffer, (char *)"wb");
						duplex->store(outp);
						fclose(outp);
					}
				}
				fprintf(dGFileAA, "%g\t%g\t%g\n", t, duplex->freeEnsembleEnergy(0.5), duplex->Z(0.5));
				fflush(NULL);
				if(prob)
				{
					FILE *inp = NULL, *outp = NULL;
					Probability *dupprob;
					if(retrieveP)
					{
						sprintf(buffer2, ".pin%s", buffer);
						inp = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, buffer2, (char *)"rb");
					}
					if(storeP)
					{
						sprintf(buffer2, ".pout%s", buffer);
						outp = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, buffer2, (char *)"wb");
					}
					if(storeP || retrieveP)
						dupprob = new Probability(inp, outp, energy, duplex, end_time, debug, procNum, quiet, no_isolate);
					else
						dupprob = new Probability(energy, duplex, debug, procNum, quiet, no_isolate);
					if(dupprob->isFinished()) dupprob->print(probFileAA);
					if(inp) fclose(inp);
					if(outp) fclose(outp);
					delete dupprob;
				}
				delete duplex;
			}
			if(BB)
			{
				PartitionFunction *duplex;
				if(retrieveQ)
				{
					FILE *inp;
					inp = openfile(seq[i+1]->getFileName(), seq[i+1]->getFileName(), i+1, i+1, buffer, (char *)"rb");
					duplex = new PartitionFunction(inp, energy, bb, debug, procNum, quiet, no_isolate);
					fclose(inp);
				}
				else
				{
					duplex = new PartitionFunction(energy, bb, debug, procNum, quiet, no_isolate);
					if(storeQ)
					{
						FILE *outp;
						outp = openfile(seq[i+1]->getFileName(), seq[i+1]->getFileName(), i+1, i+1, buffer, (char *)"wb");
						duplex->store(outp);
						fclose(outp);
					}
				}
				fprintf(dGFileBB, "%g\t%g\t%g\n", t, duplex->freeEnsembleEnergy(0.5), duplex->Z(0.5));
				fflush(NULL);
				if(prob)
				{
					FILE *inp = NULL, *outp = NULL;
					Probability *dupprob;
					if(retrieveP)
					{
						sprintf(buffer2, ".pin%s", buffer);
						inp = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, buffer2, (char *)"rb");
					}
					if(storeP)
					{
						sprintf(buffer2, ".pout%s", buffer);
						outp = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, buffer2, (char *)"wb");
					}
					if(storeP || retrieveP)
						dupprob = new Probability(inp, outp, energy, duplex, end_time, debug, procNum, quiet, no_isolate);
					else
						dupprob = new Probability(energy, duplex, debug, procNum, quiet, no_isolate);
					if(dupprob->isFinished()) dupprob->print(probFileBB);
					if(inp) fclose(inp);
					if(outp) fclose(outp);
					delete dupprob;
				}
				delete duplex;
			}
			if(AB)
			{
				PartitionFunction *duplex;
				if(retrieveQ)
				{
					FILE *inp;
					inp = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, buffer, (char *)"rb");
					duplex = new PartitionFunction(inp, energy, ab, debug, procNum, quiet, no_isolate);
					fclose(inp);
				}
				else
				{
					duplex = new PartitionFunction(energy, ab, debug, procNum, quiet, no_isolate);
					if(storeQ)
					{
						FILE *outp;
						outp = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, buffer, (char *)"wb");
						duplex->store(outp);
						fclose(outp);
					}
				}
				fprintf(dGFileAB, "%g\t%g\t%g\n", t, duplex->freeEnsembleEnergy(1.0), duplex->Z(1.0));
				fflush(NULL);
				if(prob)
				{
					FILE *inp = NULL, *outp = NULL;
					Probability *dupprob;
					if(retrieveP)
					{
						sprintf(buffer2, ".pin%s", buffer);
						inp = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, buffer2, (char *)"rb");
					}
					if(storeP)
					{
						sprintf(buffer2, ".pout%s", buffer);
						outp = openfile(seq[i]->getFileName(), seq[i+1]->getFileName(), i, i+1, buffer2, (char *)"wb");
					}
					if(storeP || retrieveP)
						dupprob = new Probability(inp, outp, energy, duplex, end_time, debug, procNum, quiet, no_isolate);
					else
						dupprob = new Probability(energy, duplex, debug, procNum, quiet, no_isolate);
					if(dupprob->isFinished()) dupprob->print(probFileAB);
					if(inp) fclose(inp);
					if(outp) fclose(outp);
					delete dupprob;
				}
				delete duplex;
			}
			delete aa[0];
			delete aa[1];
			delete bb[0];
			delete bb[1];
		}
 
		if(!quiet) printf("Running time: %ld seconds.\n", time(NULL) - now);
		fprintf(logfile, "Running time: %ld seconds.\n", time(NULL) - now);
		fclose(dGFileA);
		fclose(dGFileB);
		if(AA) fclose(dGFileAA);
		if(BB) fclose(dGFileBB);
		if(AB) fclose(dGFileAB);
		if(prob && AA) fclose(probFileAA);
		if(prob && BB) fclose(probFileBB);
		if(prob && AB) fclose(probFileAB);
		fclose(logfile);
	}
	free(buffer);
	free(buffer2);
}

Partition::~Partition()
{
	delete opts;
	delete energy;
}

double inline ln1pex(double x)
{
  /* calculates ln(1 + exp(x)) in numerically stable manner */

  if (x >= 35.0)
    return x;

  return log(1.0 + exp(x));
}

bool inline readZ(FILE *f, double *t, double *z)
{
	char buffer[10000];
	strcpy(buffer, "#");

	while(buffer[0] == '#') if(!fgets(buffer, 10000, f)) return false;
	return sscanf(buffer, "%lg%*g%lg", t, z) == 2;
}

FILE *openF1(char *n, char *ext, const char *rw)
{
	FILE *ret;
	char buffer[1000];

	sprintf(buffer, "%s.%s", n, ext);
	if (!(ret = fopen(buffer, rw)))
	{
		perror(buffer);
		exit(EXIT_FAILURE);
	}
	return ret;
}

FILE *openF2(char *n1, char *n2, char *ext, const char *rw)
{
	FILE *ret;
	char buffer[1000];

	sprintf(buffer, "%s-%s.%s", n1, n2, ext);
	if (!(ret = fopen(buffer, rw)))
	{
		perror(buffer);
		exit(EXIT_FAILURE);
	}
	return ret;
}

void Partition::computeConcentrations()
{
	if(files < 2)
	{
		fputs("Error: you need to provide the two prefixes, e.g. pirna --A0=1.0 --B0=1.0 rna1 rna2.\n", stderr);
		exit(EXIT_FAILURE);
	}

	if(!quiet) printf("Initial concentrations A0: %g B0: %g\n", Na0, Nb0);

	FILE *dGFileA, *dGFileB, *dGFileAA, *dGFileBB, *dGFileAB, *concFile;

	dGFileA = openF1(filenames[0], (char *)"dG", "rt");
	dGFileB = openF1(filenames[1], (char *)"dG", "rt");
	dGFileAA = openF2(filenames[0], filenames[0], (char *)"dG", "rt");
	dGFileBB = openF2(filenames[1], filenames[1], (char *)"dG", "rt");
	dGFileAB = openF2(filenames[0], filenames[1], (char *)"dG", "rt");
	concFile = openF2(filenames[0], filenames[1], (char *)"conc", "wt");

	fputs("T\t[Af]\t[Bf]\t[A]\t[B]\t[AA]\t[BB]\t[AB]\n", concFile);
	
	while(!feof(dGFileAB))
	{
		double Ta, Tb, Taa, Tbb, Tab;
		double Za, Zb, Zaa, Zbb, Zab;

		if (!readZ(dGFileA, &Ta, &Za)) break;
		if (!readZ(dGFileB, &Tb, &Zb)) break;
		if (!readZ(dGFileAA, &Taa, &Zaa)) break;
		if (!readZ(dGFileBB, &Tbb, &Zbb)) break;
		if (!readZ(dGFileAB, &Tab, &Zab)) break;

		if(Ta != Tb || Ta != Taa || Ta != Tbb || Ta != Tab || Tb != Taa || Tb != Tbb || Tb != Tab || Taa != Tbb || Taa != Tab || Tbb != Tab)
		{
			fprintf(stderr, "Error: temperatures do not match across different files %g\t%g\t%g\t%g\t%g\n", Ta, Tb, Taa, Tbb, Tab);
			exit(EXIT_FAILURE);
		}

		if(!quiet) printf("Calculating concentrations for temperature %g\n", Tab);

		double lKa = log(Zaa) - 2.0*log(Za);
		double lKb = log(Zbb) - 2.0*log(Zb);
		double lKab = log(Zab) - log(Za) - log(Zb);
		double lA0 = log(Na0);
		double lB0 = log(Nb0);


		double lAleft = log(2.0) + lA0 - ln1pex(lKab + lB0) - ln1pex(ln1pex(log(8.0) + lKa + lA0 - 2.0 * ln1pex(lKab + lB0)) / 2.0);
		double lAright = lA0;
		double lA = (lAleft + lAright) / 2.0, lAold;
		do {
			lAold = lA;
			double lB1 = log((Na0 - 2.0 * exp(lKa + 2.0 * lA) - exp(lA)) / exp(lKab + lA));
			double lB2 = log(2.0) + lB0 - ln1pex(lKab + lA) - ln1pex(ln1pex(log(8.0) + lKb + lB0 - 2.0 * ln1pex(lKab + lA)) / 2.0);
			if (lB2 < lB1)
			{
				lAleft = lA;
				lA = (lA + lAright) / 2.0;
			}
			else
			{
				lAright = lA;
				lA = (lA + lAleft) / 2.0;
			}
		} while (fabs(1.0 - exp(lA - lAold)) > 1e-7);
		
		double lB = log(2.0) + lB0 - ln1pex(lKab + lA) - ln1pex(ln1pex(log(8.0) + lKb + lB0 - 2.0 * ln1pex(lKab + lA)) / 2.0);


		double Na = exp(lA);		
		double Nb = exp(lB);		
		double Naa = exp(lKa + 2*lA);
		double Nbb = exp(lKb + 2*lB);
		double Nab = exp(lKab + lA + lB);

		fprintf(concFile, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", Tab, Na*(1.0 - 1.0 / Za), Nb*(1.0 - 1.0 / Zb), Na, Nb, Naa, Nbb, Nab);
	}
	fclose(dGFileA);
	fclose(dGFileB);
	fclose(dGFileAA);
	fclose(dGFileBB);
	fclose(dGFileAB);
	fclose(concFile);
}


void Partition::computedG()
{
	FILE *dGFileA, *dGFileB, *dGFileAA, *dGFileBB, *dGFileAB, *concFile;
	FILE *outFileA, *outFileB, *outFileAA, *outFileBB, *outFileAB, *outFile;

	if(files < 2)
	{
		fputs("Error: you need to provide the two prefixes, e.g. pirna --ensemble rna1 rna2.\n", stderr);
		exit(EXIT_FAILURE);
	}

	dGFileA = openF1(filenames[0], (char *)"dG", "rt");
	dGFileB = openF1(filenames[1], (char *)"dG", "rt");
	dGFileAA = openF2(filenames[0], filenames[0], (char *)"dG", "rt");
	dGFileBB = openF2(filenames[1], filenames[1], (char *)"dG", "rt");
	dGFileAB = openF2(filenames[0], filenames[1], (char *)"dG", "rt");
	concFile = openF2(filenames[0], filenames[1], (char *)"conc", "rt");

	outFileA = openF1(filenames[0], (char *)"A.dG", "wt");
	outFileB = openF1(filenames[1], (char *)"B.dG", "wt");
	outFileAA = openF2(filenames[0], filenames[0], (char *)"AA.dG", "wt");
	outFileBB = openF2(filenames[1], filenames[1], (char *)"BB.dG", "wt");
	outFileAB = openF2(filenames[0], filenames[1], (char *)"AB.dG", "wt");
	outFile = openF2(filenames[0], filenames[1], (char *)"ens.dG", "wt");

	/* ignore first line */
	int tmp;
	for (char c = 0; c != '\n'; tmp = fscanf(concFile, "%c", &c));

	/* read each line, compute ensemble free energy, print */
	fputs("#T\tFree energy\n", outFileA);
	fputs("#T\tFree energy\n", outFileB);
	fputs("#T\tFree energy\n", outFileAA);
	fputs("#T\tFree energy\n", outFileBB);
	fputs("#T\tFree energy\n", outFileAB);
	fputs("#T\tFree energy\n", outFile);

	double Za, Zb, Zaa, Zbb, Zab;
	double Ta, Tb, Taa, Tbb, Tab, Tconc;
	double Na, Nb, Naa, Nbb, Nab;
	double A0 = 0.0, B0 = 0.0, maxA0B0 = 0.0;
	double muA, muB;
	const double TOLERANCE = 2e-4;
	const double R = 0.0019872;

	while (1)
	{
		if (!readZ(dGFileA, &Ta, &Za)) break;
		if (!readZ(dGFileB, &Tb, &Zb)) break;
		if (!readZ(dGFileAA, &Taa, &Zaa)) break;
		if (!readZ(dGFileBB, &Tbb, &Zbb)) break;
		if (!readZ(dGFileAB, &Tab, &Zab)) break;

		if(fscanf(concFile, "%lg%*g%*g%lg%lg%lg%lg%lg", &Tconc, &Na, &Nb, &Naa, &Nbb, &Nab) < 6)
			break;

		if (Ta != Tconc)
			fprintf(stderr, "Warning: temperature %g in %s.dG doesn't match temperature %g in %s-%s.conc\n", Ta, filenames[0], Tconc, filenames[0], filenames[1]);
		if (Tb != Tconc)
			fprintf(stderr, "Warning: temperature %g in %s.dG doesn't match temperature %g in %s-%s.conc\n", Tb, filenames[0], Tconc, filenames[0], filenames[1]);
		if (Taa != Tconc)
			fprintf(stderr, "Warning: temperature %g in %s-%s.dG doesn't match temperature %g in %s-%s.conc\n", Taa, filenames[0], filenames[1], Tconc, filenames[0], filenames[1]);
		if (Tbb != Tconc)
			fprintf(stderr, "Warning: temperature %g in %s-%s.dG doesn't match temperature %g in %s-%s.conc\n", Tbb, filenames[0], filenames[1], Tconc, filenames[0], filenames[1]);
		if (Tab != Tconc)
			fprintf(stderr, "Warning: temperature %g in %s-%s.dG doesn't match temperature %g in %s-%s.conc\n", Tab, filenames[0], filenames[1], Tconc, filenames[0], filenames[1]);


		if(A0 == 0.0 && B0 == 0.0)
		{
			A0 = Na + Nab + 2 * Naa;
			B0 = Nb + Nab + 2 * Nbb;
			maxA0B0 = A0 > B0 ? A0 : B0;
		}
		else
		{
			if (fabs(Na + Nab + 2.0 * Naa - A0) / A0 > TOLERANCE && !quiet)
				fprintf(stderr, "Warning: at %g degrees the relative error of [A]+2[AA]+[AB] is %g\n", Tconc, fabs(Na + Nab + 2.0 * Naa - A0) / A0);
			if (fabs(Nb + Nab + 2.0 * Nbb - B0) / B0 > TOLERANCE && !quiet)
				fprintf(stderr, "Warning: at %g degrees the relative error of [B]+2[BB]+[AB] is %g\n", Tconc, fabs(Nb + Nab + 2.0 * Nbb - B0) / B0);
		}

		muA = (Na == 0.0) ? 0.0 : R * (Tconc + 273.15) * (log(Na / A0) - log(Za));
		muB = (Nb == 0.0) ? 0.0 : R * (Tconc + 273.15) * (log(Nb / B0) - log(Zb));

		fprintf(outFileA, "%g\t%g\n", Tconc, muA * Na / maxA0B0);
		fprintf(outFileB, "%g\t%g\n", Tconc, muB * Nb / maxA0B0);
		fprintf(outFileAA, "%g\t%g\n", Tconc, 2.0 * muA * Naa / maxA0B0);
		fprintf(outFileBB, "%g\t%g\n", Tconc, 2.0 * muB * Nbb / maxA0B0);
		fprintf(outFileAB, "%g\t%g\n", Tconc, (muA + muB) * Nab / maxA0B0);
		fprintf(outFile, "%g\t%g\n", Tconc, (muA * A0 + muB * B0) / maxA0B0);
	}

	fclose(dGFileA);
	fclose(dGFileB);
	fclose(dGFileAA);
	fclose(dGFileBB);
	fclose(dGFileAB);
	fclose(concFile);
	fclose(outFileA);
	fclose(outFileB);
	fclose(outFileAA);
	fclose(outFileBB);
	fclose(outFileAB);
	fclose(outFile);
}

void Partition::computeTm()
{
	double dx, fpp, z, cp[3], cpmax, Tm, global_max, global_tm;
	double *T, *G;
	int capacity;
	FILE *dGfile, *CpFile, *tmCp;

	int m = 1, n;

	/* local maxima of less than CUTOFF times the global maximum are not output */
//	const double CUTOFF = 0.1;

	global_max = 0;
	global_tm = 0;

	dGfile = openF2(filenames[0], filenames[1], (char *)"ens.dG", "rt");

	capacity = 1024;
	T = (double *)alloc.xmalloc(capacity * sizeof(double));
	G = (double *)alloc.xmalloc(capacity * sizeof(double));

	int tmp;
	for (char c = 0; c != '\n'; tmp = fscanf(dGfile, "%c", &c))
		;
	for (n = 0; fscanf(dGfile, "%lg %lg", &T[n], &G[n]) == 2; ++n)
	{
		for (char c = 0; c != '\n' && !feof(dGfile); tmp = fscanf(dGfile, "%c", &c))
			;
		if (n == capacity - 1)
		{
			capacity += 1024;
			T = (double *)alloc.xrealloc(T, capacity * sizeof(double));
			G = (double *)alloc.xrealloc(G, capacity * sizeof(double));
		}
	}
	fclose(dGfile);

	if (2 * m + 1 > n)
	{
		fputs("Too few points for computation.\n", stderr);
		return;
	}

	CpFile = openF2(filenames[0], filenames[1], (char *)"ens.Cp", "wt");
	fputs("#T\tHeat Capacity\n", CpFile);
	tmCp = openF2(filenames[0], filenames[1], (char *)"Tm", "wt");

	dx = (T[n - 1] - T[0]) / (n - 1);
	z = m * (m + 1);
	cp[0] = cp[1] = cp[2] = cpmax = 0.0;
	Tm = T[0];
	for (int i = m; i < n - m; ++i)
	{
		fpp = 0.0;
		for (int j = -m; j <= m; ++j)
			fpp += (3.0 * j * j - z) * G[i + j];
		fpp = 30.0 * fpp / (dx * dx * m * (m + 1) * (4 * m * m) * (2 * m + 3));
		fprintf(CpFile, "%g\t%g\n", T[i], -(T[i] + 273.15) * fpp);

		cp[0] = cp[1];
		cp[1] = cp[2];
		cp[2] = -fpp * (273.15 + T[i]);
		if (cp[1] - cp[0] > 1e-6 && cp[1] - cp[2] > 1e-6)
		{
			Tm = T[i - 1] + 0.5 * dx * (cp[0] - cp[2]) / (cp[0] - 2.0 * cp[1] + cp[2]);
			cpmax = cp[1] - 0.125 * (cp[2] - cp[0]) * (cp[2] - cp[0]) / (cp[0] - 2.0 * cp[1] + cp[2]);
			if (cpmax > global_max)
			{
				global_max = cpmax;
				global_tm = Tm;
			}
/*			if (cpmax >= global_max * CUTOFF)
				fprintf(tmCp, "%g\t%g\n", Tm, cpmax);*/
		}
	}

	fprintf(tmCp, "%g\t%g\n", global_tm, global_max);

	free(T);
	free(G);
	fclose(CpFile);
	fclose(tmCp);
}


void tableTest(unsigned int len1, unsigned int len2, unsigned int sublen1, unsigned int sublen2)
{
//	printf("2d %d %d\n", sublen1, sublen2);
	char d2 = 2, d4 = 4;
	Table<double> *Q2 = new Table<double>(len1, sublen1, &d2, 0.0);

	for(unsigned int i1 = 0; i1 < len1; i1++)
		for(unsigned int j1 = i1; j1 < len1 && j1 < i1 + sublen1; j1++)
			*Q2->estar(i1, j1) += 1.0;

	for(unsigned int i1 = 0; i1 < len1; i1++)
		for(unsigned int j1 = i1; j1 < len1 && j1 < i1 + sublen1; j1++)
			if(Q2->element(i1, j1) != 1.0)
				printf("%d\t%d\t%d\n", i1, j1, sublen1);

	if(Q2->hasAny(0.0))
		printf("2d redundant\n");

	delete Q2;

//	printf("4d\n");
	Table<double> *Q4 = new Table<double>(len1, sublen1, &d4, 0.0);

	for(unsigned int i1 = 0; i1 < len1; i1++)
		for(unsigned int j1 = i1; j1 < len1 && j1 < i1 + sublen1; j1++)
			for(unsigned int d = i1; d <= j1; d++)
				for(unsigned int e = d; e <= j1 ; e++)
					*Q4->estar(i1, j1, d, e) += 1.0;

	for(unsigned int i1 = 0; i1 < len1; i1++)
		for(unsigned int j1 = i1; j1 < len1 && j1 < i1 + sublen1; j1++)
			for(unsigned int d = i1; d <= j1; d++)
				for(unsigned int e = d; e <= j1 ; e++)
					if(Q4->element(i1, j1, d, e) != 1.0)
						printf("%d\t%d\t%d\t%d\t%d\n", i1, j1, d, e, sublen1);

	if(Q2->hasAny(0.0))
		printf("4d redundant\n");

	delete Q4;
//	printf("duplex\n");
	Table<double> *Qd = new Table<double>(len1, len2, sublen1, sublen2, 0.0);

	for(unsigned int i1 = 0; i1 < len1; i1++)
		for(unsigned int j1 = i1; j1 < len1 && j1 < i1 + sublen1; j1++)
			for(unsigned int i2 = 0; i2 < len2; i2++)
				for(unsigned int j2 = i2; j2 < len2 && j2 < i2 + sublen2; j2++)
					*Qd->estar(i1, j1, i2, j2) += 1.0;

	for(unsigned int i1 = 0; i1 < len1; i1++)
		for(unsigned int j1 = i1; j1 < len1 && j1 < i1 + sublen1; j1++)
			for(unsigned int i2 = 0; i2 < len2; i2++)
				for(unsigned int j2 = i2; j2 < len2 && j2 < i2 + sublen2; j2++)
					if(Qd->element(i1, j1, i2, j2) != 1.0)
						printf("%d\t%d\t%d\t%d\n", i1, j1, i2, j2);

	if(Q2->hasAny(0.0))
		printf("duplex redundant\n");

	delete Qd;
}

/*int main(int argc, char** argv)
{
	for(int i=1; i <= 20; i++)
		for(int j=1; j <= 20; j++)
			tableTest(20, 20, i, j);
}*/

