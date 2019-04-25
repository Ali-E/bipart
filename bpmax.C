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
	Desc: BPMax class drives the whole program. main() is here.

	Author: Hamidreza Chitsaz
		Colorado State University
		Algorithmic Biology Lab

	Last Update by Hamidreza Chitsaz: Oct 26, 2018
	Ver: 1.4
*/

#include <math.h>
#include <time.h>

#include "bpmax.h"

int main(int argc, char** argv)
{
	BPMax bp(argc, argv);
	while(bp.more_pairs())
	{
		bp.advance();
		bp.forward();
		bp.backtrace();
		bp.release();
	}
}

BPMax::BPMax(int argc, char** argv)
{
	opts = new GetOpt(argc, argv, OPTIONS);
	total_tables_size = 0;
	big_tables_num = 0;

	suffix = NULL;
	debug = false;
	quiet = false;
	backtr = true;
	procNum = 32;
	seq_num = 0;
	current_seq = 0;

	while (opts->hasNext())
	{
		Option *current = opts->next();
		char count = current->getShortForm();

		if (count == FREE_ARG)
    			filenames.push_back(current->getArg());
      		else if (count == 'V')
			version("bpmax");
      		else if (count == 'h')
		{
			printf("Usage: bpmax [options] filename1 filename2\n");
			printf("    or bpmax [options] filename\n");
			printf("%s\n", opts->help());
			exit(0);
		}
      		else if (count == 's')
			suffix = current->getArg();
      		else if (count == 'd')
			debug = true;
      		else if (count == 'n')
			backtr = false;
      		else if (count == 'q')
			quiet = true;
      		else if (count == 'p')
			procNum = atoi(current->getArg());
	}

	files = filenames.size();

	if (files < 1)
	{
		fputs("Error: data not specified, try 'bpmax -h' for help.\n", stderr);
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

	if(files >= 2) 
		logfile = openfile(seq[0]->getFileName(), seq[1]->getFileName(), (char *)".bpmax.run", (char *)"wt");
	else
		logfile = openfile(filenames[0], (char *)".bpmax.run", (char *)"wt");

	if(files >= 2) 
		outfile = openfile(seq[0]->getFileName(), seq[1]->getFileName(), (char *)".scores", (char *)"wt");
	else
		outfile = openfile(filenames[0], (char *)".scores", (char *)"wt");

	if(files >= 2) 
		fprintf(logfile, "bpmax %s ran on %s and %s at %s\n", PACKAGE_VERSION, seq[0]->getFileName(), seq[1]->getFileName(), ctime(&now));
	else
		fprintf(logfile, "bpmax %s ran on %s at %s\n", PACKAGE_VERSION, filenames[0], ctime(&now));

	if (suffix)
		fprintf(logfile, "suffix = %s\n", suffix);
	fprintf(logfile, "no. cpus = %d\n", procNum);
	if(!backtr) fprintf(logfile, "no backtracing\n");

	fprintf(outfile, "#seq1\tseq2\tfull\tnormalized_full\tinteraction\tnormalized_interaction\n");
}

FILE * BPMax::openfile(char *fn, char *ext, char *type)
{
	FILE *out;
	char res[10000];

	sprintf(res, "%s", fn);

	strcat(res, ext);
	if (!(out = fopen(res, type)))
	{
		perror(res);
		exit(EXIT_FAILURE);
	}
	return out;
}

FILE * BPMax::openfile(char *fn1, char *fn2, char *ext, char *type)
{
	FILE *out;
	char res[10000];

	sprintf(res, "%s-%s", fn1, fn2);

	strcat(res, ext);
	if (!(out = fopen(res, type)))
	{
		perror(res);
		exit(EXIT_FAILURE);
	}
	return out;
}

void BPMax::forward()
{
	time_t now = time(NULL);

	fprintf(logfile, "Calculating for %s and %s\n", seq[0]->getName(), seq[1]->getName());
	if(!quiet) printf("Calculating for %s and %s\n", seq[0]->getName(), seq[1]->getName());

	fprintf(outfile, "%s\t%s", seq[0]->getName(), seq[1]->getName());	

	for(int i = 0; i < 2; i++)
	{
		if(!quiet) printf("Allocating memory for single basepair maximization of %s... \n", seq[i]->getName());
		int sublen = seq[i]->getLen();
		S[i] = new Table<uint16_t>(seq[i]->getLen(), sublen, &d2, 0);
		if(!quiet) printf("Allocation and initialization finished. \n");
	}

	if(!quiet) printf("Allocating memory for joint basepair maximization of %s and %s... \n", seq[0]->getName(), seq[1]->getName());
	len1 = seq[0]->getLen();
	len2 = seq[1]->getLen();
	sublen1 = len1;
	sublen2 = len2;
	F = new Table<uint16_t>(len1, len2, sublen1, sublen2, 0);
	if(!backtr) C = new Table<uint16_t>(len1, len2, sublen1, sublen2, 0);
	if(!quiet) printf("Allocation and initialization finished. \n");

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
				uint16_t *res = S[s]->estar(i, j);

				for(int d = i; d < j; d++)
				{
					uint16_t n = S[s]->element(i, d) + S[s]->element(d+1, j);
					if (n > *res)
						*res = n;
				}

				uint16_t n = scorer.intra_score(sq[(s == 0) ? i : (len2 - i - 1)], sq[(s == 0) ? j : (len2 - j - 1)]) + S[s]->element(i+1, j-1);
				if(n > *res)
					*res = n;
			}
	}


	sq1 = seq[0]->getSeq();
	sq2 = seq[1]->getSeq();

	int it = 0;

	// (l1 == 1 && l2 == 1)

	for(int i1 = 0; i1 <= len1 - 1; i1++)
		for(int i2 = 0; i2 <= len2 - 1; i2++)
		{
			int j1 = i1;
			int j2 = i2;

			uint16_t *res = F->estar(i1, j1, i2, j2);
			uint16_t *c = (!backtr) ? C->estar(i1, j1, i2, j2) : NULL;

			*res = scorer.inter_score(sq1[i1], sq2[len2 - i2 - 1]);
			if(!backtr) *c = *res;
		}


	for(int l1 = 1; l1 <= sublen1; l1++)
		for(int l2 = 1; l2 <= sublen2; l2++)
		{
			it++;

			if(!(it % 1000) && !quiet)
				printf("%d \n", it);	
	
#pragma omp parallel for num_threads(procNum)
			for(int i1 = 0; i1 <= len1 - l1; i1++)
				for(int i2 = 0; i2 <= len2 - l2; i2++)
				{
					int j1 = i1 + l1 - 1;
					int j2 = i2 + l2 - 1;

					uint16_t register *res = F->estar(i1, j1, i2, j2);
					uint16_t register *c = (!backtr) ? C->estar(i1, j1, i2, j2) : NULL;

					for(int k1 = i1; k1 < j1; k1++)
						for(int k2 = i2; k2 < j2; k2++)
						{
							uint16_t n = F->element(i1, k1, i2, k2) + F->element(k1+1, j1, k2+1, j2);

							if(n > *res)
							{
								*res = n;
								if(!backtr) *c = C->element(i1, k1, i2, k2) + C->element(k1+1, j1, k2+1, j2);
							}
						}

					for(int k1 = i1; k1 < j1; k1++)
					{
						uint16_t n = F->element(i1, k1, i2, j2) + S[0]->element(k1+1, j1);

						if(n > *res)
						{
							*res = n;
							if(!backtr) *c = C->element(i1, k1, i2, j2);
						}

						n = S[0]->element(i1, k1) + F->element(k1+1, j1, i2, j2);

						if(n > *res)
						{
							*res = n;
							if(!backtr) *c = C->element(k1+1, j1, i2, j2);
						}
					}
					
					for(int k2 = i2; k2 < j2; k2++)
					{
						uint16_t n = F->element(i1, j1, i2, k2) + S[1]->element(k2+1, j2);

						if(n > *res)
						{
							*res = n;
							if(!backtr) *c = C->element(i1, j1, i2, k2);
						}

						n = S[1]->element(i2, k2) + F->element(i1, j1, k2+1, j2);

						if(n > *res)
						{
							*res = n;
							if(!backtr) *c = C->element(i1, j1, k2+1, j2);
						}
					}

					if(j1 >= i1 + 4)
					{
						uint16_t n = scorer.intra_score(sq1[i1], sq1[j1]) + F->element(i1+1, j1-1, i2, j2);						

						if(n > *res)
						{
							*res = n;
							if(!backtr) *c = C->element(i1+1, j1-1, i2, j2);
						}
					}
				
					if(j2 >= i2 + 4)
					{
						uint16_t n = scorer.intra_score(sq2[len2 - i2 - 1], sq2[len2 - j2 - 1]) + F->element(i1, j1, i2+1, j2-1);						

						if(n > *res)
						{
							*res = n;
							if(!backtr) *c = C->element(i1, j1, i2+1, j2-1);
						}
					}
				
				}
		}		

	if(!quiet)
	{ 
		printf("Total score: %d\t", S[0]->element(0, len1-1));
		printf("%d\t", S[1]->element(0, len2-1));
 		printf("%d\n", F->element(0, len1-1, 0, len2-1));
 		if(!backtr) printf("Interaction score: %d\n", C->element(0, len1-1, 0, len2-1));
	}	
	fprintf(logfile, "Total score: %d\t", S[0]->element(0, len1-1));
	fprintf(logfile, "%d\t", S[1]->element(0, len2-1));
 	fprintf(logfile, "%d\n", F->element(0, len1-1, 0, len2-1));
 	if(!backtr) fprintf(logfile, "Interaction score: %d\n", C->element(0, len1-1, 0, len2-1));

	fprintf(outfile, "\t%d\t%g", F->element(0, len1-1, 0, len2-1), F->element(0, len1-1, 0, len2-1)*1.0/(len1+len2));


	if(!quiet) printf("Running time: %ld seconds.\n", time(NULL) - now);
	fprintf(logfile, "Running time: %ld seconds.\n", time(NULL) - now);
}


void BPMax::backtrace(int i, int j, int s, vector<Pair> *res)
{
	if(j < i + 4)
		return;

	if(debug)
		printf("%d\t%d\tstrand %d\n", i, j, s);

	Table<uint16_t> *R = S[s];

	int len = seq[s]->getLen();
	int sublen = len;
	unsigned char *sq = seq[s]->getSeq();


	uint16_t r = R->element(i, j);

	for(int d = i; d < j; d++)
	{
		uint16_t n = R->element(i, d) + R->element(d+1, j);
		if (n == r)
		{
			backtrace(i, d, s, res);
			backtrace(d+1, j, s, res);

			return;
		}
	}


	int idx1 = (s == 0) ? i : (len2 - j - 1);
	int idx2 = (s == 0) ? j : (len2 - i - 1);

	if(debug)
		printf("idx: %d\t%d\t%d\n", idx1, idx2, len2);

	uint16_t n = scorer.intra_score(sq[idx1], sq[idx2]) + R->element(i+1, j-1);
	if(debug)
		printf("val: %d\t%d\n", n, r);

	if(n == r)
	{
		if(scorer.intra_score(sq[idx1], sq[idx2]) > 0)
		{
			Pair p(idx1, idx2, (s == 0) ? INTRA1 : INTRA2);	
			res->push_back(p);
		}

		backtrace(i+1, j-1, s, res);

		return;
	}

	printf("Panic: did not expect an s value without a matching case: i: %d j: %d strand: %d s: %d \n", i, j, s, r);

	if(debug)
		printf("end: %d\t%d\tstrand %d\n", i, j, s);
}

uint16_t BPMax::backtrace(int i1, int j1, int i2, int j2, vector<Pair> *res)
{
	if(debug)
		printf("%d\t%d\t%d\t%d\n", i1, j1, i2, j2);

	uint16_t f = F->element(i1, j1, i2, j2);

	if(i1 == j1 && i2 == j2)
	{
		Pair p(i1, len2 - i2 - 1, INTER);

		if (f > 0)
			res->push_back(p);

		return f;
	}

	for(int k1 = i1; k1 < j1; k1++)
		for(int k2 = i2; k2 < j2; k2++)
		{
			uint16_t n = F->element(i1, k1, i2, k2) + F->element(k1+1, j1, k2+1, j2);

			if(n == f)
				return backtrace(i1, k1, i2, k2, res) + backtrace(k1+1, j1, k2+1, j2, res);
		}

	for(int k1 = i1; k1 < j1; k1++)
	{
		uint16_t n = F->element(i1, k1, i2, j2) + S[0]->element(k1+1, j1);

		if(n == f)
		{
			backtrace(k1+1, j1, 0, res);
			return backtrace(i1, k1, i2, j2, res);
		}

		n = S[0]->element(i1, k1) + F->element(k1+1, j1, i2, j2);

		if(n == f)
		{
			backtrace(i1, k1, 0, res);
			return backtrace(k1+1, j1, i2, j2, res);
		}
	}
					
	for(int k2 = i2; k2 < j2; k2++)
	{
		uint16_t n = F->element(i1, j1, i2, k2) + S[1]->element(k2+1, j2);

		if(n == f)
		{
			backtrace(k2+1, j2, 1, res);
			return backtrace(i1, j1, i2, k2, res);
		}

		n = S[1]->element(i2, k2) + F->element(i1, j1, k2+1, j2);

		if(n == f)
		{
			backtrace(i2, k2, 1, res);
			return backtrace(i1, j1, k2+1, j2, res);
		}
	}

	if(j1 >= i1 + 4)
	{
		uint16_t n = scorer.intra_score(sq1[i1], sq1[j1]) + F->element(i1+1, j1-1, i2, j2);						

		if(n == f)
		{
			if(scorer.intra_score(sq1[i1], sq1[j1]) > 0)
			{
				Pair p(i1, j1, INTRA1);
				res->push_back(p);
			}

			return backtrace(i1+1, j1-1, i2, j2, res);
		}
	}
				
	if(j2 >= i2 + 4)
	{
		uint16_t n = scorer.intra_score(sq2[len2 - i2 - 1], sq2[len2 - j2 - 1]) + F->element(i1, j1, i2+1, j2-1);						

		if(n == f)
		{
			if(scorer.intra_score(sq2[len2 - i2 - 1], sq2[len2 - j2 - 1]) > 0)
			{
				Pair p(len2 - i2 - 1,  len2 - j2 - 1, INTRA2);
				res->push_back(p);
			}
			return backtrace(i1, j1, i2+1, j2-1, res);
		}
	}

	printf("Panic: did not expect an f value without a matching case: i1: %d j1: %d i2: %d j2: %d  f: %d \n", i1, j1, i2, j2, f);

	return 0;
}


void BPMax::backtrace()
{
	if(!backtr)
		return;

	fprintf(logfile, "Back tracing for %s and %s\n", seq[0]->getName(), seq[1]->getName());
	if(!quiet) printf("Back tracing for %s and %s\n", seq[0]->getName(), seq[1]->getName());
	

	time_t now = time(NULL);


	FILE *struc = openfile(seq[0]->getFileName(), seq[1]->getFileName(), (char *)".struct", (char *)"wt");


#pragma omp parallel for num_threads(procNum)
	for(int i1 = 0; i1 <= len1 - sublen1; i1++)
		for(int i2 = 0; i2 <= len2 - sublen2; i2++)
		{
			int j1 = i1 + sublen1 - 1;
			int j2 = i2 + sublen2 - 1;

			vector<Pair> res;

			int sco = backtrace(i1, j1, i2, j2, &res);

			if(!quiet)
			{
				printf("C(%d, %d, %d, %d) = %d\n", i1, j1, i2, j2, sco);
				for(int i = 0; i < res.size(); i++)
					res[i].print((char *)"\n");
			}

			fprintf(logfile, "C(%d, %d, %d, %d) = %d\n", i1, j1, i2, j2, sco);
			for(int i = 0; i < res.size(); i++)
				res[i].print(logfile, (char *)"\n");


			fprintf(struc, "C(%d, %d, %d, %d) = %d\n", i1, j1, i2, j2, sco);
			for(int i = 0; i < res.size(); i++)
				res[i].print(struc, (char *)"\n");
			
			fprintf(outfile, "\t%d\t%g\n", sco, sco*1.0/MIN(sublen1, sublen2));
		}		


	fclose(struc);

	if(!quiet) printf("Running time: %ld seconds.\n", time(NULL) - now);
	fprintf(logfile, "Running time: %ld seconds.\n", time(NULL) - now);
}

bool BPMax::more_pairs()
{
	return(current_seq < seq_num);
}

void BPMax::advance()
{
	seq[0] = seqs[current_seq++];
	seq[1] = seqs[current_seq++];
}

void BPMax::release()
{
	for(int i = 0; i < 2; i++)
		delete S[i];
	delete F;
	if(!backtr) delete C;
}

BPMax::~BPMax()
{
	delete opts;
	for(int s = 0; s < seq_num; s++)
		delete seqs[s];

	fclose(logfile);
	fclose(outfile);
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

