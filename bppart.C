/*
Copyright 2018  Hamidreza Chitsaz (chitsaz@chitsazlab.org)

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
	Desc: BPPart class drives the whole program. main() is here.

	Authors: Hamidreza Chitsaz and Ali Ebrahimpour Boroojeny
		Colorado State University
		Algorithmic Biology Lab

*/

#include <math.h>
#include <time.h>
#include <string.h>
#include <sstream>
#include "bppart.h"
#include <limits>

int main(int argc, char** argv)
{
	BPPart bp(argc, argv);
	while(bp.more_pairs())
	{
		bp.advance();
		bp.allocate();
		bp.compute();
		bp.release();
	}
}

BPPart::BPPart(int argc, char** argv)
{
	opts = new GetOpt(argc, argv, OPTIONS);
	total_tables_size = 0;
	big_tables_num = 0;

	suffix = NULL;
	debug = false;
	quiet = false;
	procNum = 32;
	files = 0;
	seq_num = 0;
	current_seq = 0;
	var2 = 2;
	var3 = 1;

	while (opts->hasNext())
	{
		Option *current = opts->next();
		char count = current->getShortForm();

		if (count == FREE_ARG)
    			filenames.push_back(current->getArg());
      		else if (count == 'V')
			version("bppart");
      		else if (count == 'h')
		{
			printf("Usage: bppart [options] filename1 filename2\n");
			printf("    or bppart [options] filename\n");
			printf("%s\n", opts->help());
			exit(0);
		}
      		else if (count == 's')
			suffix = current->getArg();
      		else if (count == 'd')
			debug = true;
      		else if (count == 'q')
			quiet = true;
      		else if (count == 'p')
			procNum = atoi(current->getArg());
		else if (count == 'M')
			var2 = atof(current->getArg());
		else if (count == 'N')
			var3 = atof(current->getArg());
	}

	files = filenames.size();

	if (files < 1)
	{
		fputs("Error: data not specified, try 'bppart -h' for help.\n", stderr);
		exit(1);
	}

	if (files > 2)
		fputs("Warning: extra filenames are being ignored.\n", stderr);

	if (files >= 2)
	{
		seqs.push_back(new Sequence(filenames[0]));
		seqs.push_back(new Sequence(filenames[1]));
		seq[0] = seqs[0];
		seq[1] = seqs[1];
		if(!quiet) printf("Sequences %s and %s are read.\n", seq[0]->getName(), seq[1]->getName());
	} else
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
				seqs.push_back(motherseq);
				Sequence *s_temp;
				while((s_temp = motherseq->split())) seqs.push_back(s_temp);
			}
		}
		fclose(fp);
	}

	seq_num = seqs.size();

	if(seq_num % 2 || seq_num < 2)
	{
		fputs("Error: odd number of sequences. Please provide pairs.\n", stderr);
		exit(1);
	}

	time_t now = time(NULL);

	std::ostringstream var2_str;
        var2_str << var2;
        std::string str_1 = var2_str.str();
        char *var2_s = new char[str_1.length() + 1];
        strcpy(var2_s, str_1.c_str());
        //const char *var2_s = var2_str.str().c_str();

        std::ostringstream var3_str;
        var3_str << var3;
        std::string str_2 = var3_str.str();
        char *var3_s = new char[str_2.length() + 1];
        strcpy(var3_s, str_2.c_str());

	if(files >= 2) 
		logfile = openfile(seq[0]->getFileName(), seq[1]->getFileName(), (char *)".bppart.run.", (char *) var2_s, (char *) var3_s, (char *)"wt");
	else
		logfile = openfile(filenames[0], (char *)".bppart.run.", (char *) var2_s, (char *) var3_s, (char *)"wt");

	if(files >= 2) 
		outfile = openfile(seq[0]->getFileName(), seq[1]->getFileName(), (char *)".parts.", (char *) var2_s, (char *) var3_s, (char *)"wt");
	else
		outfile = openfile(filenames[0], (char *)".parts.", (char *) var2_s, (char *) var3_s, (char *)"wt");

	if(files >= 2) 
		fprintf(logfile, "bppart %s ran on %s and %s at %s\n", PACKAGE_VERSION, seq[0]->getFileName(), seq[1]->getFileName(), ctime(&now));
	else
		fprintf(logfile, "bppart %s ran on %s at %s\n", PACKAGE_VERSION, filenames[0], ctime(&now));

	if (suffix)
		fprintf(logfile, "suffix = %s\n", suffix);
	fprintf(logfile, "no. cpus = %d\n", procNum);

	fprintf(outfile, "#seq1\tseq2\tQI\tQ0\tQ1\t-log(QI)\t-log(QI-Q0*Q1)\tlog_interaction_normalized\tlog_full_normalized\n");
}

FILE * BPPart::openfile(char *fn, char *ext, char *var2_s, char *var3_s, char *type)
{
	FILE *out;
	char res[10000];

	sprintf(res, "%s", fn);

	strcat(res, ext);
	strcat(res, var2_s);
        strcat(res, "_");
        strcat(res, var3_s);
	if (!(out = fopen(res, type)))
	{
		perror(res);
		exit(EXIT_FAILURE);
	}
	return out;
}

FILE * BPPart::openfile(char *fn1, char *fn2, char *ext, char *var2_s, char *var3_s, char *type)
{
	FILE *out;
	char res[10000];

	sprintf(res, "%s-%s", fn1, fn2);
	strcat(res, ext);
	strcat(res, var2_s);
        strcat(res, "_");
        strcat(res, var3_s);
	if (!(out = fopen(res, type)))
	{
		perror(res);
		exit(EXIT_FAILURE);
	}
	return out;
}

double BPPart::score(int a, int b, double var2, double var3)
{
	//int sc = scorer.intra_score(a, b);

	double sc = scorer.intra_score(a, b, var2, var3);

	if(sc != -std::numeric_limits<double>::infinity())
		return exp(sc);
	else
		return 0;
}

double BPPart::iscore(int a, int b, double var2, double var3)
{
	//int sc = scorer.inter_score(a, b);

	double sc = scorer.inter_score(a, b, var2, var3);

	if(sc != -std::numeric_limits<double>::infinity())
		return exp(sc);
	else
		return 0;
}


void BPPart::allocate()
{
	for(int i = 0; i < 2; i++)
	{
		if(!quiet) printf("Allocating memory for single basepair partition function of %s... \n", seq[i]->getName());
		int sublen = seq[i]->getLen();
		Q[i] = new Table<double>(seq[i]->getLen()+1, sublen+1, &d2, 1);
		Qz[i] = new Table<double>(seq[i]->getLen()+1, sublen+1, &d2, 0);
		if(!quiet) printf("Allocation and initialization finished. \n");
	}

	if(!quiet) printf("Allocating memory for joint basepair partition function of %s and %s... \n", seq[0]->getName(), seq[1]->getName());
	len1 = seq[0]->getLen();
	len2 = seq[1]->getLen();
	sublen1 = len1;
	sublen2 = len2;

	

	QI = new Table<double>(len1+1, len2+1, sublen1+1, sublen2+1, 0);
	QIa = new Table<double>(len1+1, len2+1, sublen1+1, sublen2+1, 0);
	QIac = new Table<double>(len1+1, len2+1, sublen1+1, sublen2+1, 0);
	for(int i = 0; i < 2; i++)
	{
		QIs[i] = new Table<double>(len1+1, len2+1, sublen1+1, sublen2+1, 0);
		QIaux[i] = new Table<double>(len1+1, len2+1, sublen1+1, sublen2+1, 0);
	}
	QIe = new Table<double>(len1+1, len2+1, sublen1+1, sublen2+1, 0);
	QIm = new Table<double>(len1+1, len2+1, sublen1+1, sublen2+1, 0);
	if(!quiet) printf("Allocation and initialization finished. \n");
}

void BPPart::compute()
{
	time_t now = time(NULL);

	fprintf(logfile, "Calculating for %s and %s\n", seq[0]->getName(), seq[1]->getName());
	if(!quiet) printf("Calculating for %s and %s\n", seq[0]->getName(), seq[1]->getName());

	fprintf(outfile, "%s\t%s", seq[0]->getName(), seq[1]->getName());	

	for(int s = 0; s < 2; s++)
	{
		int len = seq[s]->getLen();
		int sublen = len;
		unsigned char *sq = seq[s]->getSeq();

		for(int l = 5; l <= sublen; l++)
#pragma omp parallel for num_threads(procNum)
			for(int i = 0; i <= len-l; i++)
			{
				int j = i + l - 1;
				double *res = Qz[s]->estar(i, j+1);

				for(int d = i+4; d <= j; d++)
				{
					double sc = score(sq[(s == 0) ? i : (len2 - i - 1)], sq[(s == 0) ? d : (len2 - d - 1)], var2, var3);

					double n = Q[s]->element(i+1,d) * sc * Q[s]->element(d+1, j+1);

					*res += n;
				}

				res = Q[s]->estar(i, j+1);

				for(int d = i; d <= j-4; d++)
					*res += Qz[s]->element(d,j+1);
			}
	}


	sq1 = seq[0]->getSeq();
	sq2 = seq[1]->getSeq();

	int it = 0;

	for(int l1 = 0; l1 <= sublen1; l1++)
		for(int l2 = 0; l2 <= sublen2; l2++)
		{
			it++;

			if(!(it % 1000) && !quiet)
				printf("%d \n", it);	
	
#pragma omp parallel for num_threads(procNum)
			for(int i1 = 0; i1 <= len1 - l1; i1++)
				for(int i2 = 0; i2 <= len2 - l2; i2++)
				{
					register int j1 = i1 + l1 - 1;
					register int j2 = i2 + l2 - 1;

					register double *res;


					res = QIs[0]->estar(i1, j1+1, i2, j2+1);
					if(i1<=j1-4 && i2<=j2)
					{
						for(int k1 = i1+1; k1 <= j1-1; k1++)
							*res += Q[0]->element(i1+1,k1) * QIaux[0]->element(k1, j1, i2, j2+1);

						*res *= score(sq1[i1], sq1[j1], var2, var3);
					}

					res = QIs[1]->estar(i1, j1+1, i2, j2+1);
					if(i1<=j1 && i2<=j2-4)
					{
						for(int k2 = i2+1; k2 <= j2-1; k2++)
							*res += Q[1]->element(i2+1,k2) * QIaux[1]->element(i1, j1+1, k2, j2);

						*res *= score(sq2[len2-i2-1], sq2[len2-j2-1], var2, var3);
					}

					res = QIm->estar(i1, j1+1, i2, j2+1);
					if(i1 == j1 && i2 == j2)
						*res = iscore(sq1[i1], sq2[len2-i2-1], var2, var3);
					else
						if(i1 < j1 && i2 < j2)
						{
							for(int k1 = i1+1; k1 <= j1; k1++)
								for(int k2 = i2+1; k2 <= j2; k2++)
									*res += (QIa->element(i1, k1, i2, k2) + QI->element(i1+1, k1, i2+1, k2) * iscore(sq1[i1], sq2[len2-i2-1], var2, var3)) * QIac->element(k1, j1+1, k2, j2+1);

							*res += QI->element(i1+1, j1, i2+1, j2) * iscore(sq1[i1], sq2[len2-i2-1], var2, var3) * iscore(sq1[j1], sq2[len2-j2-1], var2, var3) + 
								iscore(sq1[j1], sq2[len2-j2-1], var2, var3) * QIa->element(i1, j1, i2, j2);					
						}

					res = QIaux[0]->estar(i1, j1+1, i2, j2+1);
					for(int k1 = i1; k1 <= j1; k1++)
						*res += (QIs[0]->element(i1, k1+1, i2, j2+1) + QIm->element(i1, k1+1, i2, j2+1)) * Q[0]->element(k1+1, j1+1); 

					res = QIaux[1]->estar(i1, j1+1, i2, j2+1);
					for(int k2 = i2; k2 <= j2; k2++)
						*res += (QIs[1]->element(i1, j1+1, i2, k2+1) + QIm->element(i1, j1+1, i2, k2+1)) * Q[1]->element(k2+1, j2+1); 
					

					res = QIe->estar(i1, j1+1, i2, j2+1);
					if(i1 <= j1-4 && i2 <= j2-4)
						*res = (QI->element(i1+1, j1, i2+1, j2) - Q[0]->element(i1+1, j1)*Q[1]->element(i2+1, j2)) * score(sq1[i1], sq1[j1], var2, var3) * score(sq2[len2-i2-1], sq2[len2-j2-1], var2, var3);

					res = QIac->estar(i1, j1+1, i2, j2+1);
					*res = QIs[0]->element(i1, j1+1, i2, j2+1) + QIs[1]->element(i1, j1+1, i2, j2+1) + QIe->element(i1, j1+1, i2, j2+1);
				
					res = QIa->estar(i1, j1+1, i2, j2+1);
					// (i1 <= j1 && i2 <= j2) 
					for(int k1 = i1; k1 <= j1; k1++)
						for(int k2 = i2; k2 <= j2; k2++)
							*res += QIac->element(i1, k1+1, i2, k2+1) * QI->element(k1+1,j1+1,k2+1,j2+1);


					res = QI->estar(i1, j1+1, i2, j2+1);
					*res = Q[0]->element(i1,j1+1) * Q[1]->element(i2,j2+1);
					// (i1 <= j1 && i2 <= j2) 
					for(int k1 = i1; k1 <= j1; k1++)
						for(int k2 = i2; k2 <= j2; k2++)
							*res += Q[0]->element(i1, k1) * Q[1]->element(i2, k2) * 
								(iscore(sq1[k1], sq2[len2-k2-1], var2, var3) * QI->element(k1+1,j1+1,k2+1,j2+1) + QIa->element(k1,j1+1,k2,j2+1));

				}

		}		

	if(!quiet)
		printf("\t%e\t%e\t%e\t%g\t%g\t%g\t%g\n", QI->element(0, len1, 0, len2), Q[0]->element(0, len1), Q[1]->element(0, len2), -log(QI->element(0, len1, 0, len2)), -log(QI->element(0, len1, 0, len2) - Q[0]->element(0, len1)*Q[1]->element(0, len2)), log(QI->element(0, len1, 0, len2) - Q[0]->element(0, len1)*Q[1]->element(0, len2))/MIN(len1, len2), log(QI->element(0, len1, 0, len2))/(len1+len2));

	fprintf(logfile, "\t%e\t%e\t%e\t%g\t%g\t%g\t%g\n", QI->element(0, len1, 0, len2), Q[0]->element(0, len1), Q[1]->element(0, len2), -log(QI->element(0, len1, 0, len2)), -log(QI->element(0, len1, 0, len2) - Q[0]->element(0, len1)*Q[1]->element(0, len2)), log(QI->element(0, len1, 0, len2) - Q[0]->element(0, len1)*Q[1]->element(0, len2))/MIN(len1, len2), log(QI->element(0, len1, 0, len2))/(len1+len2));

	fprintf(outfile, "\t%e\t%e\t%e\t%g\t%g\t%g\t%g\n", QI->element(0, len1, 0, len2), Q[0]->element(0, len1), Q[1]->element(0, len2), -log(QI->element(0, len1, 0, len2)), -log(QI->element(0, len1, 0, len2) - Q[0]->element(0, len1)*Q[1]->element(0, len2)), log(QI->element(0, len1, 0, len2) - Q[0]->element(0, len1)*Q[1]->element(0, len2))/MIN(len1, len2), log(QI->element(0, len1, 0, len2))/(len1+len2));


	if(!quiet) printf("Running time: %ld seconds.\n", time(NULL) - now);
	fprintf(logfile, "Running time: %ld seconds.\n", time(NULL) - now);
}

bool BPPart::more_pairs()
{
	return(current_seq < seq_num);
}

void BPPart::advance()
{
	seq[0] = seqs[current_seq++];
	seq[1] = seqs[current_seq++];
}

void BPPart::release()
{
	for(int i = 0; i < 2; i++)
	{
		delete Q[i];
		delete Qz[i];
	}

	delete QI;
	delete QIa;
	delete QIac;
	for(int i = 0; i < 2; i++)
	{
		delete QIs[i];
		delete QIaux[i];
	}
	delete QIe;
	delete QIm;
}


BPPart::~BPPart()
{
	delete opts;
	delete seq[0], seq[1];

	fclose(logfile);
}


