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
	Desc: Energy class computes energies.

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow, 
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: June 2, 2009
	Ver: 1.6
*/

#include "config.h"

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#if HAVE_IEEEFP_H
# include <ieeefp.h>
#endif

#ifdef __SGI
#include <sigfpe.h>
#include <sys/fpu.h>
#endif

#include <ctype.h>
#include <string.h>

#include "energy.h"
#include "alloc.h"


#ifndef PRECISION
# define PRECISION 1
#endif

#ifdef INTEGER
# define isFinite(x) (x < INFINITY / 2)
#else
# define isFinite(x) finite(x)
#endif

#ifdef INFINITY
# undef INFINITY
#endif

#ifndef isinf
# define isinf(x) (!finite(x) && x == x)
#endif

#ifdef NO_GU_BASEPAIRS
const int BPI[6][6] = {{6, 6, 6, 0, 6, 6},
		       {6, 6, 1, 6, 6, 6},
		       {6, 2, 6, 6, 6, 6},
		       {3, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6}};
#else
const int BPI[6][6] = {{6, 6, 6, 0, 6, 6},
		       {6, 6, 1, 6, 6, 6},
		       {6, 2, 6, 4, 6, 6},
		       {3, 6, 5, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6}};
#endif
#define basePairIndex(a, b) BPI[a][b]

#define scale(d) ((d) * PRECISION)
const double INFINITY = 1.0 / 0.0;


#ifndef __SGI
#ifndef __SUNOS
#ifndef __LINUX
bool finite(double x)
{
  return (x <= DBL_MAX && x >= -DBL_MAX);
}
#endif
#endif
#endif

int min(int a, int b)
{
  return (a < b) ? a : b;
}

int min3(int a, int b, int c)
{
  if (a <= b && a <= c)
    return a;
  if (b <= c)
    return b;
  return c;
}

int triloopcmp(const void* loop1, const void* loop2)
{
  int i;
  const unsigned char* h1 = (unsigned char*) loop1;
  const struct triloop *h2 = (triloop *) loop2;

  for (i = 0; i < 5; ++i)
    if (h1[i] < h2->loop[i])
      return -1;
    else if (h1[i] > h2->loop[i])
      return 1;

  return 0;
}

int tloopcmp(const void* loop1, const void* loop2)
{
  int i;
  const unsigned char* h1 = (unsigned char*)loop1;
  const struct tloop *h2 = (tloop *)loop2;

  for (i = 0; i < 6; ++i)
    if (h1[i] < h2->loop[i])
      return -1;
    else if (h1[i] > h2->loop[i])
      return 1;

  return 0;
}


int hexaloopcmp(const void* loop1, const void* loop2)
{
  int i;
  const unsigned char* h1 = (unsigned char*) loop1;
  const struct hexaloop *h2 = (hexaloop *) loop2;

  for (i = 0; i < 8; ++i)
    if (h1[i] < h2->loop[i])
      return -1;
    else if (h1[i] > h2->loop[i])
      return 1;

  return 0;
}


Energy::Energy(char *dataDir, int _NA, int _polymer, double _naConc, double _mgConc, char *suf, double temperature, int z, int nodangle, int itype) 
{
	R = .0019872; 
	BASES[0] = 'A'; 
	BASES[1] = 'C'; 
	BASES[2] = 'G'; 
	BASES[3] = 'U'; 
	BASES[4] = 'N'; 
	strcpy(BASE_PAIRS[0], "A-U"); 
	strcpy(BASE_PAIRS[1], "C-G"); 
	strcpy(BASE_PAIRS[2], "G-C"); 
	strcpy(BASE_PAIRS[3], "U-A"); 
	strcpy(BASE_PAIRS[4], "G-U"); 
	strcpy(BASE_PAIRS[5], "U-G"); 
	strcpy(PKGDATADIR, dataDir);
	NA = _NA;
	polymer = _polymer;
	naConc = _naConc;
	mgConc = _mgConc;
	suffix = suf;
	tRatio = (temperature + 273.15) / 310.15;
	scaleFactor = 1.0;
	inttype = itype;
	zip = z;
	g_nodangle = nodangle;
	saltCorrection = ion();
	loadEntropiesEnthalpies();
	computeEnergies();
}

Energy::~Energy()
{
	if(!suffix)
	{
		free(g_triloop);
		free(g_tloop);
		free(g_hexaloop);
	}
}

bool Energy::base_pair(unsigned char a, unsigned char b) { return basePairIndex(a, b) < 6; }

unsigned char Energy::toNum(char c)
{
  c = toupper(c);
  switch (c)
    {
    case 'A': case '0':
      return 0;
    case 'C': case '1':
      return 1;
    case 'G': case '2':
      return 2;
    case 'T': case 'U': case '3':
      return 3;
    }
  return 4;
}


void Energy::getParams(int *_NA, int *_polymer, double *_naConc, double *_mgConc, char **suf, double *temperature, int *z, int *nodangle, int *itype)
{
	*_NA = NA;
	*_polymer = polymer;
	*_naConc = naConc;
	*_mgConc = mgConc;
	*suf = suffix;
	*itype = inttype;
	*z = zip;
	*nodangle = g_nodangle;
	*temperature = tRatio*310.15 - 273.15;
}

void Energy::setTemperature(double temperature)
{
	tRatio = (temperature + 273.15) / 310.15;
	computeEnergies();
}


double Energy::ion()
{
  if (NA == 0)
    return 0;
  else
    if (polymer)
      {
	if (mgConc != 0.0)
	  fputs("Warning: [Mg++] correction ignored for polymer mode\n", stderr);
	return -(0.2 + 0.175 * log(naConc));
      }
    else
      return -0.114 * log(naConc + 3.3 * sqrt(mgConc));
}

FILE* Energy::openFile(char* name)
{
  FILE* file;
  char* buffer;

  file = fopen(name, "rt");

  if (!file)
    {
      buffer = (char *)alloc.xmalloc(strlen(PKGDATADIR) + strlen(name) + 1);
      strcpy(buffer, PKGDATADIR);
      strcat(buffer, name);
      if (!(file = fopen(buffer, "rt")))
	{
	  perror(name);
	  exit(EXIT_FAILURE);
	}
      free(buffer);
    }

  return file;
}

double readFloat(FILE *f)
{
	char buffer[10000];

	int tmp = fscanf(f, "%s", buffer);

	for(int i = 0; i < strlen(buffer); i++)
		buffer[i] = toupper(buffer[i]);

	double ret;
	if(!strcmp(buffer, "INF"))
		ret = INFINITY;
	else
		ret = atof(buffer);

	return ret;
} 

void Energy::loadStack()
{
  int i, j, ii, jj;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"stack.DGD");
      hFile = openFile((char *)"stack.DHD");
    }
  else
    {
      gFile = openFile((char *)"stack.DG");
      hFile = openFile((char *)"stack.DH");
    }

  for (i = 0; i < 5; ++i)
    for (ii = 0; ii < 5; ++ii)
      for (j = 0; j < 5; ++j)
	for (jj = 0; jj < 5; ++jj)
	  if (i == 4 || j == 4 || ii == 4 || jj == 4)
	    stackEnthalpies[i][j][ii][jj] = INFINITY;
	  else
	    {
	      stackEnergies[i][j][ii][jj] = readFloat(gFile);
	      stackEnergies[i][j][ii][jj] += saltCorrection;
	      stackEnthalpies[i][j][ii][jj] = readFloat(hFile);
	    }

  fclose(gFile);
  fclose(hFile);
}

void Energy::combineStack()
{
  int i, j, ii, jj;

  for (i = 0; i < 5; ++i)
    for (ii = 0; ii < 5; ++ii)
      for (j = 0; j < 5; ++j)
	for (jj = 0; jj < 5; ++jj)
	  if (i == 4 || j == 4 || ii == 4 || jj == 4)
	    g_stack[i][j][ii][jj] = INFINITY;
	  else if (!isFinite(stackEnergies[i][j][ii][jj]) || !isFinite(stackEnthalpies[i][j][ii][jj]))
	    g_stack[i][j][ii][jj] = INFINITY;
	  else
	    g_stack[i][j][ii][jj] = scale(tRatio * stackEnergies[i][j][ii][jj] + (1.0 - tRatio) * stackEnthalpies[i][j][ii][jj]);
}

void Energy::calculateStack(double stack[5][5][5][5], double estack[5][5][5][5])
{
  int i, j, ii, jj;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (ii = 0; ii < 5; ++ii)
	for (jj = 0; jj < 5; ++jj)
	  stack[i][j][ii][jj] = exp(-estack[i][j][ii][jj] / RT) / scaleFactor / scaleFactor;
}

void Energy::calculateStack2(double stack[5][5][6][6], double estack[5][5][6][6])
{
  int i, j, ii, jj;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (ii = 0; ii < 6; ++ii)
	for (jj = 0; jj < 6; ++jj)
	  stack[i][j][ii][jj] = exp(-estack[i][j][ii][jj] / RT) / scaleFactor / scaleFactor;
}

void Energy::calculateZipStack2(double stack[5][5][6][6])
{
  int i, j, ii, jj;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (ii = 0; ii < 6; ++ii)
	for (jj = 0; jj < 6; ++jj)
	  if (ii == 5 || jj == 5)
	    stack[i][j][ii][jj] = exp(-stack[i][j][ii][jj] / RT) / scaleFactor / scaleFactor;
	  else if (!isFinite(stack[i][j][ii][jj]))
	    stack[i][j][ii][jj] = (-g_dangle3[i][j][ii] - g_dangle5[i][j][jj]) / scaleFactor / scaleFactor;
	  else
	    stack[i][j][ii][jj] = (exp(-stack[i][j][ii][jj] / RT) - g_dangle3[i][j][ii] - g_dangle5[i][j][jj] - 1) / scaleFactor / scaleFactor;
}

void Energy::calculateZeroStack2(double stack[5][5][6][6])
{
  int i, j, ii, jj;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (ii = 0; ii < 6; ++ii)
	for (jj = 0; jj < 6; ++jj)
	  stack[i][j][ii][jj] = 0;
}

void Energy::calculateInfStack2(double stack[5][5][6][6])
{
  int i, j, ii, jj;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (ii = 0; ii < 6; ++ii)
	for (jj = 0; jj < 6; ++jj)
	  stack[i][j][ii][jj] = INFINITY;
}

void Energy::loadStackSuffix()
{
  int i, j, ii, jj;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(6 + strlen(suffix) + 1);
  strcpy(buffer, "stack.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i = 0; i < 5; ++i)
    for (ii = 0; ii < 5; ++ii)
      for (j = 0; j < 5; ++j)
	for (jj = 0; jj < 5; ++jj)
	  if (i == 4 || j == 4 || ii == 4 || jj == 4)
	    g_stack[i][j][ii][jj] = INFINITY;
	  else
	    {
	      d = readFloat(file);
	      g_stack[i][j][ii][jj] = scale(d);
	    }

  fclose(file);
}

void Energy::symmetryCheckStack(double stack[4][4][4][4], char* which)
{
  int i, j, ii, jj;

  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      for (ii = 0; ii < 4; ++ii)
	for (jj = 0; jj < 4; ++jj)
	  if (stack[i][j][ii][jj] != stack[jj][ii][j][i])
	    fprintf(stderr, "Warning: %c-%c/%c-%c stack %s is %g; %c-%c/%c-%c stack %s is %g\n", BASES[i], BASES[j], BASES[ii], BASES[jj], which, stack[i][j][ii][jj], BASES[jj], BASES[ii], BASES[j], BASES[i], which, stack[jj][ii][j][i]);
}

double Energy::estimateScale(double stack[5][5][5][5])
{
  int i1, i2;
  double avg = 0.0;

  for (i1 = 0; i1 < 4; ++i1)
    for (i2 = 0; i2 < 4; ++i2)
      avg += stack[i1][3 - i1][i2][3 - i2];

  return avg / 80;
}

void Energy::loadDangle()
{
  int i, j, k;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"dangle.DGD");
      hFile = openFile((char *)"dangle.DHD");
    }
  else
    {
      gFile = openFile((char *)"dangle.DG");
      hFile = openFile((char *)"dangle.DH");
    }

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	if (i == 4 || j == 4)
	  dangleEnthalpies3[i][j][k] = INFINITY;
	else if (k == 4)
	  dangleEnthalpies3[i][j][k] = INFINITY;
	else if (k == 5)
	  dangleEnthalpies3[i][j][k] = INFINITY;
	else
	  {
	    dangleEnergies3[i][j][k] = readFloat(gFile);
	    dangleEnergies3[i][j][k] += 0.5 * saltCorrection;
	    dangleEnthalpies3[i][j][k] = readFloat(hFile);
	  }

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	if (i == 4 || j == 4)
	  dangleEnthalpies5[i][j][k] = INFINITY;
	else if (k == 4)
	  dangleEnthalpies5[i][j][k] = INFINITY;
	else if (k == 5)
	  dangleEnthalpies5[i][j][k] = INFINITY;
	else
	{
	  dangleEnergies5[i][j][k] = readFloat(gFile);
	  dangleEnergies5[i][j][k] += 0.5 * saltCorrection;
	  dangleEnthalpies5[i][j][k] = readFloat(hFile);
	}

  fclose(gFile);
  fclose(hFile);
}

void Energy::combineDangle()
{
  int i, j, k;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  if (i == 4 || j == 4)
	    {
	      g_dangle3[i][j][k] = INFINITY;
	      g_dangle5[i][j][k] = INFINITY;
	    }
	  else if (k == 4)
	    {
	      g_dangle3[i][j][k] = INFINITY;
	      g_dangle5[i][j][k] = INFINITY;
	    }
	  else if (k == 5)
	    {
	      g_dangle3[i][j][k] = INFINITY;
	      g_dangle5[i][j][k] = INFINITY;
	    }
	  else
	    {
	      if (!isFinite(dangleEnergies3[i][j][k]) || !isFinite(dangleEnthalpies3[i][j][k]))
		g_dangle3[i][j][k] = INFINITY;
	      else
		g_dangle3[i][j][k] = scale(tRatio * dangleEnergies3[i][j][k] + (1.0 - tRatio) * dangleEnthalpies3[i][j][k]);
	      if (!isFinite(dangleEnergies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k]))
		g_dangle5[i][j][k] = INFINITY;
	      else
		g_dangle5[i][j][k] = scale(tRatio * dangleEnergies5[i][j][k] + (1.0 - tRatio) * dangleEnthalpies5[i][j][k]);
	    }
	}
}

void Energy::combineDangleNew()
{
  int i, j, k;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  if (i == 4 || j == 4)
	    {
	      g_dangle3[i][j][k] = 0;
	      g_dangle5[i][j][k] = 0;
	    }
	  else if (k == 4)
	    {
	      g_dangle3[i][j][k] = 0;
	      g_dangle5[i][j][k] = 0;
	    }
	  else if (k == 5)
	    {
	      g_dangle3[i][j][k] = 0;
	      g_dangle5[i][j][k] = 0;
	    }
	  else
	    {
	      if (!isFinite(dangleEnergies3[i][j][k]) || !isFinite(dangleEnthalpies3[i][j][k]))
		g_dangle3[i][j][k] = 0;
	      else
		{
		  double ea, b;
		  b = -dangleEnthalpies3[i][j][k] / R * exp(-dangleEnergies3[i][j][k] / R / 310.15) / (exp(-dangleEnergies3[i][j][k] / R / 310.15) - 1);
		  ea = (exp(-dangleEnergies3[i][j][k] / R / 310.15) - 1) / exp(b / 310.15);
		  g_dangle3[i][j][k] = ea * exp(b / tRatio / 310.15);
		}
	      if (!isFinite(dangleEnergies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k]))
		g_dangle5[i][j][k] = 0;
	      else
		{
		  double ea, b;
		  b = -dangleEnthalpies5[i][j][k] / R * exp(-dangleEnergies5[i][j][k] / R / 310.15) / (exp(-dangleEnergies5[i][j][k] / R / 310.15) - 1);
		  ea = (exp(-dangleEnergies5[i][j][k] / R / 310.15) - 1) / exp(b / 310.15);
		  g_dangle5[i][j][k] = ea * exp(b / tRatio / 310.15);
		}
	    }
	}
}

void Energy::calculateDangle()
{
  int i, j, k;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  z_dangle3[i][j][k] = exp(-g_dangle3[i][j][k] / RT) / scaleFactor;
	  z_dangle5[i][j][k] = exp(-g_dangle5[i][j][k] / RT) / scaleFactor;
	}
}

void Energy::calculateZipDangle()
{
  int i, j, k;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	if (k == 5)
	  {
	    z_dangle3[i][j][k] = exp(-g_dangle3[i][j][k] / RT) / scaleFactor;
	    z_dangle5[i][j][k] = exp(-g_dangle5[i][j][k] / RT) / scaleFactor;
	  }
	else
	  {
	    z_dangle3[i][j][k] = (exp(-g_dangle3[i][j][k] / RT) - 1.0) / scaleFactor;
	    z_dangle5[i][j][k] = (exp(-g_dangle5[i][j][k] / RT) - 1.0) / scaleFactor;
	  }
}

void Energy::calculateZeroDangle()
{
  int i, j, k;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  z_dangle3[i][j][k] = 0;
	  z_dangle5[i][j][k] = 0;
	}
}

void Energy::calculateInfDangle()
{
  int i, j, k;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  z_dangle3[i][j][k] = INFINITY;
	  z_dangle5[i][j][k] = INFINITY;
	}
}

void Energy::loadDangleSuffix()
{
  int i, j, k;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(7 + strlen(suffix) + 1);
  strcpy(buffer, "dangle.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	if (i == 4 || j == 4)
	  g_dangle3[i][j][k] = INFINITY;
	else if (k == 4)
	  g_dangle3[i][j][k] = INFINITY;
	else if (k == 5)
	  g_dangle3[i][j][k] = INFINITY;
	else
	  {
	    d = readFloat(file);
	    g_dangle3[i][j][k] = scale(d);
	  }
    
  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	if (i == 4 || j == 4)
	  g_dangle5[i][j][k] = INFINITY;
	else if (k == 4)
	  g_dangle5[i][j][k] = INFINITY;
	else if (k == 5)
	  g_dangle5[i][j][k] = INFINITY;
	else
	  {
	    d = readFloat(file);
	    g_dangle5[i][j][k] = scale(d);
	  }
    
  fclose(file);
}

void Energy::zipDangle()
{
  int i, j;
 
  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      g_dangle5[i][j][5] = g_dangle3[i][j][5] = 0.0;
}

void Energy::addZeroDangle()
{
  int i, j, k;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  g_dangle3[i][j][k] = -RT * log(1 + exp(-g_dangle3[i][j][k] / RT));
	  g_dangle5[i][j][k] = -RT * log(1 + exp(-g_dangle5[i][j][k] / RT));
	}
}

void Energy::minZeroDangle()
{
  int i, j, k;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  if (g_dangle3[i][j][k] > 0.0)
	    g_dangle3[i][j][k] = 0.0;
	  if (g_dangle5[i][j][k] > 0.0)
	    g_dangle5[i][j][k] = 0.0;
	}
}

void Energy::loadLoop()
{
  int k;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"loop.DGD");
      hFile = openFile((char *)"loop.DHD");
    }
  else
    {
      gFile = openFile((char *)"loop.DG");
      hFile = openFile((char *)"loop.DH");
    }

  for (k = 0; k < 30; ++k)
    {
      int tmp = fscanf(gFile, "%*f"); 
      interiorLoopEnergies[k] = readFloat(gFile);
      bulgeLoopEnergies[k] = readFloat(gFile);
      hairpinLoopEnergies[k] = readFloat(gFile);
      bulgeLoopEnergies[k] += saltCorrection * (1.0 + 0.5 * min(k, 10));
      interiorLoopEnergies[k] += saltCorrection * (1.0 + 0.5 * min(k, 10));
      tmp = fscanf(hFile, "%*f"); 
      interiorLoopEnthalpies[k] = readFloat(hFile);
      bulgeLoopEnthalpies[k] = readFloat(hFile);
      hairpinLoopEnthalpies[k] = readFloat(hFile);
    }

  fclose(gFile);
  fclose(hFile);
}

void Energy::combineLoop()
{
  int k;

  for (k = 0; k < 30; ++k)
    {
      g_hairpinLoop[k] = scale(tRatio * hairpinLoopEnergies[k] + (1.0 - tRatio) * hairpinLoopEnthalpies[k]);
      g_interiorLoop[k] = scale(tRatio * interiorLoopEnergies[k] + (1.0 - tRatio) * interiorLoopEnthalpies[k]);
      g_bulgeLoop[k] = scale(tRatio * bulgeLoopEnergies[k] + (1.0 - tRatio) * bulgeLoopEnthalpies[k]);
    }
}

void Energy::calculateLoop()
{
  int k;
  const double RT = tRatio * 310.15 * R;

  for (k = 0; k < 30; ++k)
    {
      z_hairpinLoop[k] = exp(-g_hairpinLoop[k] / RT) / pow(scaleFactor, k + 3);
      z_interiorLoop[k] = exp(-g_interiorLoop[k] / RT) / pow(scaleFactor, k + 3);
      z_bulgeLoop[k] = exp(-g_bulgeLoop[k] / RT) / pow(scaleFactor, k + 3);
    }
}

void Energy::loadLoopSuffix()
{
  int k;
  double d1, d2, d3;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(5 + strlen(suffix) + 1);
  strcpy(buffer, "loop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (k = 0; k < 30; ++k)
    {
      int tmp = fscanf(file, "%*f"); 
      d1 = readFloat(file); 
      d2 = readFloat(file); 
      d3 = readFloat(file);
      g_interiorLoop[k] = scale(d1);
      g_bulgeLoop[k] = scale(d2);
      g_hairpinLoop[k] = scale(d3);
    }

  fclose(file);
}

void Energy::loadSint2()
{
  int b, c, i, j;
  FILE* gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"sint2.DGD");
      hFile = openFile((char *)"sint2.DHD");
    }
  else
    {
      gFile = openFile((char *)"sint2.DG");
      hFile = openFile((char *)"sint2.DH");
    }

  for (b = 0; b < 7; ++b)
    for (i = 0; i < 5; ++i)
      for (c = 0; c < 7; ++c)
	for (j = 0; j < 5; ++j)
	  if (b == 6 || c == 6)
	    sint2Enthalpies[b][c][i][j] = INFINITY;
	  else if (i == 4 || j == 4)
	    sint2Enthalpies[b][c][i][j] = 0;
	  else
	    {
	      sint2Energies[b][c][i][j] = readFloat(gFile);
	      sint2Energies[b][c][i][j] += 1.0 * saltCorrection;
	      sint2Enthalpies[b][c][i][j] = readFloat(hFile);
	    }

  fclose(gFile);
  fclose(hFile);
}

void Energy::combineSint2()
{
  int b, c, i, j;

  for (b = 0; b < 6; ++b)
    for (c = 0; c < 6; ++c)
      for (i = 0; i < 5; ++i)
	for (j = 0; j < 5; ++j)
	  {
	    if (b == 6 || c == 6)
	      g_sint2[b][c][i][j] = INFINITY;
	    else if (i == 4 || j == 4)
	      g_sint2[b][c][i][j] = 0;
	    else if (!isFinite(sint2Energies[b][c][i][j]) || !isFinite(sint2Enthalpies[b][c][i][j]))
	      g_sint2[b][c][i][j] = INFINITY;
	    else
	      g_sint2[b][c][i][j] = scale(tRatio * sint2Energies[b][c][i][j] + (1.0 - tRatio) * sint2Enthalpies[b][c][i][j]);
	  }
}

void Energy::calculateSint2()
{
  int b, c, i, j;
  const double RT = tRatio * 310.15 * R;

  for (b = 0; b < 6; ++b)
    for (c = 0; c < 6; ++c)
      for (i = 0; i < 5; ++i)
	for (j = 0; j < 5; ++j)
	  z_sint2[b][c][i][j] = exp(-g_sint2[b][c][i][j] / RT) / pow(scaleFactor, 4);
}

void Energy::loadSint2Suffix()
{
  int b, c, i, j;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(6 + strlen(suffix) + 1);
  strcpy(buffer, "sint2.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (b = 0; b < 6; ++b)
    for (i = 0; i < 5; ++i)
      for (c = 0; c < 6; ++c)
	for (j = 0; j < 5; ++j)
	  if (b == 6 || c == 6)
	    g_sint2[b][c][i][j] = INFINITY;
	  else if (i == 4 || j == 4)
	    g_sint2[b][c][i][j] = 0;
	  else
	    {
	      d = readFloat(file);
	      g_sint2[b][c][i][j] = scale(d);
	    }

  fclose(file);
}

void Energy::symmetryCheckSint2(double sint2[6][6][4][4], char* which)
{
  int b, c, i, j;

  for (b = 0; b < 6; ++b)
    for (c = 0; c < 6; ++c)
      for (i = 0; i < 4; ++i)
	for (j = 0; j < 4; ++j)
	  {
	    int bb, cc;
	    bb = (b > 3) ? 9 - b : 3 - b;
	    cc = (c > 3) ? 9 - c : 3 - c;

	    if (sint2[b][c][i][j] != sint2[cc][bb][j][i])
	      fprintf(stderr, "Warning: %s/%s/%c/%c sint2 %s is %g; %s/%s/%c/%c sint2 %s is %g\n", BASE_PAIRS[b], BASE_PAIRS[c], BASES[i], BASES[j], which, sint2[b][c][i][j], BASE_PAIRS[cc], BASE_PAIRS[bb], BASES[j], BASES[i], which, sint2[cc][bb][j][i]);
	  }
}

void Energy::loadAsint1x2()
{
  int b, c, i, j, k;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"asint1x2.DGD");
      hFile = openFile((char *)"asint1x2.DHD");
    }
  else
    {
      gFile = openFile((char *)"asint1x2.DG");
      hFile = openFile((char *)"asint1x2.DH");
    }

  for (b = 0; b < 7; ++b)
    for (k = 0; k < 5; ++k)
      for (i = 0; i < 5; ++i)
	for (c = 0; c < 7; ++c)
	  for (j = 0; j < 5; ++j)
	    if (b == 6 || c == 6)
	      asint1x2Enthalpies[b][c][i][j][k] = INFINITY;
	    else if (i == 4 || j == 4 || k == 4)
	      asint1x2Enthalpies[b][c][i][j][k] = 0;
	    else
	      {
		asint1x2Energies[b][c][i][j][k] = readFloat(gFile);
		asint1x2Energies[b][c][i][j][k] += 2.5 * saltCorrection;
		asint1x2Enthalpies[b][c][i][j][k] = readFloat(hFile);
	      }

  fclose(gFile);
  fclose(hFile);
}

void Energy::combineAsint1x2()
{
  int b, c, i, j, k;

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i = 0; i < 5; ++i)
	for (j = 0; j < 5; ++j)
	  for (k = 0; k < 5; ++k)
	    {
	      if (b == 6 || c == 6)
		g_asint1x2[b][c][i][j][k] = INFINITY;
	      else if (i == 4 || j == 4 || k == 4)
		g_asint1x2[b][c][i][j][k] = 0;
	      else if (!isFinite(asint1x2Energies[b][c][i][j][k]) || !isFinite(asint1x2Enthalpies[b][c][i][j][k]))
		g_asint1x2[b][c][i][j][k] = INFINITY;
	      else
		g_asint1x2[b][c][i][j][k] = scale(tRatio * asint1x2Energies[b][c][i][j][k] + (1.0 - tRatio) * asint1x2Enthalpies[b][c][i][j][k]);
	    }
}

void Energy::calculateAsint1x2()
{
  int b, c, i, j, k;
  const double RT = tRatio * 310.15 * R;

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i = 0; i < 5; ++i)
	for (j = 0; j < 5; ++j)
	  for (k = 0; k < 5; ++k)
	    z_asint1x2[b][c][i][j][k] = exp(-g_asint1x2[b][c][i][j][k] / RT) / pow(scaleFactor, 5);
}

void Energy::loadAsint1x2Suffix()
{
  int b, c, i, j, k;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "asint1x2.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (b = 0; b < 7; ++b)
    for (k = 0; k < 5; ++k)
      for (i = 0; i < 5; ++i)
	for (c = 0; c < 7; ++c)
	  for (j = 0; j < 5; ++j)
	    if (b == 6 || c == 6)
	      g_asint1x2[b][c][i][j][k] = INFINITY;
	    else if (i == 4 || j == 4 || k == 4)
	      g_asint1x2[b][c][i][j][k] = 0;
	    else
	      {
		d = readFloat(file);
		g_asint1x2[b][c][i][j][k] = scale(d);
	      }

  fclose(file);
}

void Energy::loadSint4()
{
  int b, c, i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"sint4.DGD");
      hFile = openFile((char *)"sint4.DHD");
    }
  else
    {
      gFile = openFile((char *)"sint4.DG");
      hFile = openFile((char *)"sint4.DH");
    }

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i1 = 0; i1 < 5; ++i1)
	for (j1 = 0; j1 < 5; ++j1)
	  for (i2 = 0; i2 < 5; ++i2)
	    for (j2 = 0; j2 < 5; ++j2)
	      if (b == 6 || c == 6)
		sint4Enthalpies[b][c][i1][j1][i2][j2] = INFINITY;
	      else if (i1 == 4 || j1 == 4 || i2 == 4 || j2 == 4)
		sint4Enthalpies[b][c][i1][j1][i2][j2] = 0;
	      else
		{
		  sint4Energies[b][c][i1][j1][i2][j2] = readFloat(gFile);
		  sint4Energies[b][c][i1][j1][i2][j2] += 3.0 * saltCorrection;
		  sint4Enthalpies[b][c][i1][j1][i2][j2] = readFloat(hFile);
		}

  fclose(gFile);
  fclose(hFile);
}

void Energy::combineSint4()
{
  int b, c, i1, j1, i2, j2;

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i1 = 0; i1 < 5; ++i1)
	for (j1 = 0; j1 < 5; ++j1)
	  for (i2 = 0; i2 < 5; ++i2)
	    for (j2 = 0; j2 < 5; ++j2)
	      {
		if (b == 6 || c == 6)
		  g_sint4[b][c][i1][j1][i2][j2] = INFINITY;
		else if (i1 == 4 || j1 == 4 || i2 == 4 || j2 == 4)
		  g_sint4[b][c][i1][j1][i2][j2] = 0;
		else if (!isFinite(sint4Energies[b][c][i1][j1][i2][j2]) || !isFinite(sint4Enthalpies[b][c][i1][j1][i2][j2]))
		  g_sint4[b][c][i1][j1][i2][j2] = INFINITY;
		else
		  g_sint4[b][c][i1][j1][i2][j2] = scale(tRatio * sint4Energies[b][c][i1][j1][i2][j2] + (1.0 - tRatio) * sint4Enthalpies[b][c][i1][j1][i2][j2]);
	      }
}

void Energy::calculateSint4()
{
  int b, c, i1, j1, i2, j2;
  const double RT = tRatio * 310.15 * R;

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i1 = 0; i1 < 5; ++i1)
	for (j1 = 0; j1 < 5; ++j1)
	  for (i2 = 0; i2 < 5; ++i2)
	    for (j2 = 0; j2 < 5; ++j2)
	      z_sint4[b][c][i1][j1][i2][j2] = exp(-g_sint4[b][c][i1][j1][i2][j2] / RT) / pow(scaleFactor, 6);
}

void Energy::loadSint4Suffix()
{
  int b, c, i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(6 + strlen(suffix) + 1);
  strcpy(buffer, "sint4.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i1 = 0; i1 < 5; ++i1)
	for (j1 = 0; j1 < 5; ++j1)
	  for (i2 = 0; i2 < 5; ++i2)
	    for (j2 = 0; j2 < 5; ++j2)
	      if (b == 6 || c == 6)
		g_sint4[b][c][i1][j1][i2][j2] = INFINITY;
	      else if (i1 == 4 || j1 == 4 || i2 == 4 || j2 == 4)
		g_sint4[b][c][i1][j1][i2][j2] = 0;
	      else
		{
		  int tmp = fscanf(file, "%lf", &d);
		  g_sint4[b][c][i1][j1][i2][j2] = scale(d);
		}

  fclose(file);
}

void Energy::symmetryCheckSint4(double sint4[6][6][4][4][4][4], char* which)
{
  int b, c, i1, j1,i2, j2;

  for (b = 0; b < 6; ++b)
    for (c = 0; c < 6; ++c)
      for (i1 = 0; i1 < 4; ++i1)
	for (j1 = 0; j1 < 4; ++j1)
	  for (i2 = 0; i2 < 4; ++i2)
	    for (j2 = 0; j2 < 4; ++j2)
	      {
		int bb, cc;
		bb = (b > 3) ? 9 - b : 3 - b;
		cc = (c > 3) ? 9 - c : 3 - c;

		if (sint4[b][c][i1][j1][i2][j2] != sint4[cc][bb][j2][i2][j1][i1])
		  fprintf(stderr, "Warning: %s/%s/%c/%c/%c/%c sint4 %s is %g; %s/%s/%c/%c/%c/%c sint4 %s is %g\n", BASE_PAIRS[b], BASE_PAIRS[c], BASES[i1], BASES[j1], BASES[i2], BASES[j2], which, sint4[b][c][i1][j1][i2][j2], BASE_PAIRS[cc], BASE_PAIRS[bb], BASES[j2], BASES[i2], BASES[j1], BASES[i1], which, sint4[cc][bb][j2][i2][j1][i1]);
	      }
}

void Energy::loadTstacki()
{
  int i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"tstacki.DGD");
      hFile = openFile((char *)"tstacki.DHD");
    }
  else
    {
      gFile = openFile((char *)"tstacki.DG");
      hFile = openFile((char *)"tstacki.DH");
    }

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackiEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackiEnthalpies[i1][j1][i2][j2] = 0;
	  else
	    {
	      tstackiEnergies[i1][j1][i2][j2] = readFloat(gFile);
	      tstackiEnthalpies[i1][j1][i2][j2] = readFloat(hFile);;
	    }

  fclose(gFile);
  fclose(hFile);
}

void Energy::loadTstackh()
{
  int i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"tstackh.DGD");
      hFile = openFile((char *)"tstackh.DHD");
    }
  else
    {
      gFile = openFile((char *)"tstackh.DG");
      hFile = openFile((char *)"tstackh.DH");
    }

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackhEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackhEnthalpies[i1][j1][i2][j2] = 0;
	  else
	    {
	      tstackhEnergies[i1][j1][i2][j2] = readFloat(gFile);
	      tstackhEnthalpies[i1][j1][i2][j2] = readFloat(hFile);
	    }

  fclose(gFile);
  fclose(hFile);
}

void Energy::loadTstackm()
{
  int i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"tstackm.DGD");
      hFile = openFile((char *)"tstackm.DHD");
    }
  else
    {
      gFile = openFile((char *)"tstackm.DG");
      hFile = openFile((char *)"tstackm.DH");
    }

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 6; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 6; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackmEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 5 || j2 == 5)
	    tstackmEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackmEnthalpies[i1][j1][i2][j2] = 0;
	  else
	    {
	      tstackmEnergies[i1][j1][i2][j2] = readFloat(gFile);
	      tstackmEnergies[i1][j1][i2][j2] += saltCorrection;
	      tstackmEnthalpies[i1][j1][i2][j2] = readFloat(hFile);
	    }

  fclose(gFile);
  fclose(hFile);
}

void Energy::loadTstacke()
{
  int i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"tstacke.DGD");
      hFile = openFile((char *)"tstacke.DHD");
    }
  else
    {
      gFile = openFile((char *)"tstacke.DG");
      hFile = openFile((char *)"tstacke.DH");
    }

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 6; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 6; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackeEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 5 || j2 == 5)
	    tstackeEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackeEnthalpies[i1][j1][i2][j2] = 0;
	  else
	    {
	      tstackeEnergies[i1][j1][i2][j2] = readFloat(gFile);
	      tstackeEnergies[i1][j1][i2][j2] += saltCorrection;
	      tstackeEnthalpies[i1][j1][i2][j2] = readFloat(hFile);
	    }

  fclose(gFile);
  fclose(hFile);
}

void Energy::combineTstack(double tstackEnergies[4][4][4][4], double tstackEnthalpies[5][5][5][5], double tstack[5][5][5][5])
{
  int i1, j1, i2, j2;

  for (i1 = 0; i1 < 5; ++i1)
    for (j1 = 0; j1 < 5; ++j1)
      for (i2 = 0; i2 < 5; ++i2)
	for (j2 = 0; j2 < 5; ++j2)
	  {
	    if (i1 == 4 || j1 == 4)
	      tstack[i1][j1][i2][j2] = INFINITY;
	    else if (i2 == 4 || j2 == 4)
	      tstack[i1][j1][i2][j2] = 0;
	    else if (!isFinite(tstackEnergies[i1][j1][i2][j2]) || !isFinite(tstackEnthalpies[i1][j1][i2][j2]))
	      tstack[i1][j1][i2][j2] = INFINITY;
	    else
	      tstack[i1][j1][i2][j2] = scale(tRatio * tstackEnergies[i1][j1][i2][j2] + (1.0 - tRatio) * tstackEnthalpies[i1][j1][i2][j2]);
	  }
}

void Energy::combineTstack2(double tstackEnergies[4][4][4][4], double tstackEnthalpies[5][5][6][6], double tstack[5][5][6][6])
{
  int i1, j1, i2, j2;

  for (i1 = 0; i1 < 5; ++i1)
    for (j1 = 0; j1 < 5; ++j1)
      for (i2 = 0; i2 < 6; ++i2)
	for (j2 = 0; j2 < 6; ++j2)
	  {
	    if (i1 == 4 || j1 == 4)
	      tstack[i1][j1][i2][j2] = INFINITY;
	    else if (i2 == 5 || j2 == 5)
	      tstack[i1][j1][i2][j2] = INFINITY;
	    else if (i2 == 4 || j2 == 4)
	      tstack[i1][j1][i2][j2] = 0;
	    else if (!isFinite(tstackEnergies[i1][j1][i2][j2]) || !isFinite(tstackEnthalpies[i1][j1][i2][j2]))
	      tstack[i1][j1][i2][j2] = INFINITY;
	    else
	      tstack[i1][j1][i2][j2] = scale(tRatio * tstackEnergies[i1][j1][i2][j2] + (1.0 - tRatio) * tstackEnthalpies[i1][j1][i2][j2]);
	  }
}

void Energy::loadTstackiSuffix()
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "tstacki.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    g_tstacki[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    g_tstacki[i1][j1][i2][j2] = 0;
	  else
	    {
	      d = readFloat(file);
	      g_tstacki[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void Energy::loadTstackhSuffix()
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "tstackh.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    g_tstackh[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    g_tstackh[i1][j1][i2][j2] = 0;
	  else
	    {
	      d = readFloat(file);
	      g_tstackh[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void Energy::loadTstackmSuffix()
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "tstackm.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 6; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 6; ++j2)
	  if (i1 == 4 || j1 == 4)
	    g_tstackm[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 5 || j2 == 5)
	    g_tstackm[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    g_tstackm[i1][j1][i2][j2] = 0;
	  else
	    {
	      d = readFloat(file);
	      g_tstackm[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void Energy::loadTstackeSuffix()
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "tstacke.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 6; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 6; ++j2)
	  if (i1 == 4 || j1 == 4)
	    g_tstacke[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 5 || j2 == 5)
	    g_tstacke[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    g_tstacke[i1][j1][i2][j2] = 0;
	  else
	    {
	      d = readFloat(file);
	      g_tstacke[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void Energy::loadCoaxialSuffix()
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "coaxial.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    g_coaxial[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    g_coaxial[i1][j1][i2][j2] = 0;
	  else
	    {
	      d = readFloat(file);
	      g_coaxial[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void Energy::loadTstackcoaxSuffix()
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(11 + strlen(suffix) + 1);
  strcpy(buffer, "tstackcoax.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    g_tstackcoax[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    g_tstackcoax[i1][j1][i2][j2] = 0;
	  else
	    {
	      d = readFloat(file);
	      g_tstackcoax[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void Energy::loadCoaxstackSuffix()
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(10 + strlen(suffix) + 1);
  strcpy(buffer, "coaxstack.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    g_coaxstack[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    g_coaxstack[i1][j1][i2][j2] = 0;
	  else
	    {
	      d = readFloat(file);
	      g_coaxstack[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void Energy::loadMulti()
{
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"miscloop.DGD");
      hFile = openFile((char *)"miscloop.DHD");
    }
  else
    {
      gFile = openFile((char *)"miscloop.DG");
      hFile = openFile((char *)"miscloop.DH");
    }

  int tmp = fscanf(gFile, "%*g%*g%*g%*g%*g%*g%lg%lg%lg", &multiEnergies[0], &multiEnergies[1], &multiEnergies[2]);
  tmp = fscanf(hFile, "%*g%*g%*g%*g%*g%*g%lg%lg%lg", &multiEnthalpies[0], &multiEnthalpies[1], &multiEnthalpies[2]);

  fclose(gFile);
  fclose(hFile);
}

void Energy::loadMulti2()
{
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"miscloop.DGD");
      hFile = openFile((char *)"miscloop.DHD");
    }
  else
    {
      gFile = openFile((char *)"miscloop.DG");
      hFile = openFile((char *)"miscloop.DH");
    }

  int tmp = fscanf(gFile, "%*g%*g%*g%*g%*g%*g%*g%*g%*g%lg%lg%lg", &multiEnergies2[0], &multiEnergies2[1], &multiEnergies2[2]);
  tmp = fscanf(hFile, "%*g%*g%*g%*g%*g%*g%*g%*g%*g%lg%lg%lg", &multiEnthalpies2[0], &multiEnthalpies2[1], &multiEnthalpies2[2]);

  fclose(gFile);
  fclose(hFile);
}

void Energy::combineMulti()
{
  int i;

  for (i = 0; i < 3; ++i)
  {
    g_multi[i] = scale(tRatio * multiEnergies[i] + (1.0 - tRatio) * multiEnthalpies[i]);
    g_multi2[i] = scale(tRatio * multiEnergies2[i] + (1.0 - tRatio) * multiEnthalpies2[i]);
  }
}

void Energy::calculateMulti()
{
  int i;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 3; ++i)
  {
    z_multi[i] = exp(-g_multi[i] / RT);
    z_multi2[i] = exp(-g_multi2[i] / RT);
  }
  z_multi[0] /= scaleFactor * scaleFactor;
  z_multi2[0] /= scaleFactor * scaleFactor;
}

void Energy::loadMultiSuffix()
{
  double d1, d2, d3;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "miscloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  int tmp = fscanf(file, "%*g%*g%*g%*g%*g%*g%lg%lg%lg", &d1, &d2, &d3);
  g_multi[0] = scale(d1);
  g_multi[1] = scale(d2);
  g_multi[2] = scale(d3);

  fclose(file);
}

void Energy::loadMulti2Suffix()
{
  double d1, d2, d3;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "miscloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  int tmp = fscanf(file, "%*g%*g%*g%*g%*g%*g%*g%*g%*g%lg%lg%lg", &d1, &d2, &d3);
  g_multi2[0] = scale(d1);
  g_multi2[1] = scale(d2);
  g_multi2[2] = scale(d3);

  fclose(file);
}

void Energy::loadMisc()
{
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"miscloop.DGD");
      hFile = openFile((char *)"miscloop.DHD");
    }
  else
    {
      gFile = openFile((char *)"miscloop.DG");
      hFile = openFile((char *)"miscloop.DH");
    }

  int tmp = fscanf(gFile, "%lg%lg%lg%lg%lg%lg%*g%*g%*g%*g%*g%*g%lg%lg%lg%lg%lg%lg%lg", &miscEnergies[12], &miscEnergies[4], &miscEnergies[0], &miscEnergies[1], &miscEnergies[2], &miscEnergies[3], &miscEnergies[6], &miscEnergies[8], &miscEnergies[9], &miscEnergies[10], &miscEnergies[11], &miscEnergies[5], &miscEnergies[7]);
  tmp = fscanf(hFile, "%lg%lg%lg%lg%lg%lg%*g%*g%*g%*g%*g%*g%lg%lg%lg%lg%lg%lg%lg", &miscEnthalpies[12], &miscEnthalpies[4], &miscEnthalpies[0], &miscEnthalpies[1], &miscEnthalpies[2], &miscEnthalpies[3], &miscEnthalpies[6], &miscEnthalpies[8], &miscEnthalpies[9], &miscEnthalpies[10], &miscEnthalpies[11], &miscEnthalpies[5], &miscEnthalpies[7]);

  fclose(gFile);
  fclose(hFile);
}

void Energy::combineMisc()
{
  int i;

  for (i = 0; i < 7; ++i)
    g_misc[i] = scale(tRatio * miscEnergies[i] + (1.0 - tRatio) * miscEnthalpies[i]);

  g_misc[7] = miscEnergies[7] == 1 || miscEnthalpies[7] == 1;
  for (i = 8; i < 13; ++i)
    g_misc[i] = scale(tRatio * miscEnergies[i] + (1.0 - tRatio) * miscEnthalpies[i]);
}

void Energy::calculateMisc()
{
  int i;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 7; ++i)
    z_misc[i] = exp(-g_misc[i] / RT);
  for (i = 8; i < 13; ++i)
    z_misc[i] = exp(-g_misc[i] / RT);
}

void Energy::loadMiscSuffix()
{
  int i;
  double d[13];
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "miscloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  int tmp = fscanf(file, "%lg%lg%lg%lg%lg%lg%*g%*g%*g%*g%*g%*g%lg%lg%lg%lg%lg%lg%lg", &d[12], &d[4], &d[0], &d[1], &d[2], &d[3], &d[6], &d[8], &d[9], &d[10], &d[11], &d[5], &d[7]);
  for (i = 0; i < 13; ++i)
    g_misc[i] = scale(d[i]);

  fclose(file);
}

void Energy::makeAUPenalty()
{
  int i, j;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      g_aup[i][j] = 0.0;

  g_aup[0][3] = g_aup[3][0] = g_aup[2][3] = g_aup[3][2] = g_misc[6];
}

void Energy::makeAUPenaltyH()
{
  int i, j;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      g_aup[i][j] = 0.0;

  g_aup[0][3] = g_aup[3][0] = g_aup[2][3] = g_aup[3][2] = g_misc[6];
}

void Energy::loadTriloop()
{
  FILE *gFile, *hFile;
  int i, size;
  double energy;

  if (NA)
    gFile = openFile((char *)"triloop.DGD");
  else
    gFile = openFile((char *)"triloop.DG");

  numTriloops = 0;
  size = 16;
  triloopEnergies = (triloopE *)calloc(16, sizeof(struct triloopE));

  while (fscanf(gFile, "%5s %lg", triloopEnergies[numTriloops].loop, &energy) == 2)
    {
      for (i = 0; i < 5; ++i)
	triloopEnergies[numTriloops].loop[i] = toNum(triloopEnergies[numTriloops].loop[i]);
      triloopEnergies[numTriloops].energy = energy;
      ++numTriloops;
      if (numTriloops == size)
	{
	  size *= 2;
	  triloopEnergies = (triloopE *)realloc(triloopEnergies, size * sizeof(struct triloopE));
	}
    }

  triloopEnergies = (triloopE *)realloc(triloopEnergies, numTriloops * sizeof(struct triloopE));

  fclose(gFile);

  if (NA)
    hFile = openFile((char *)"triloop.DHD");
  else
    hFile = openFile((char *)"triloop.DH");

  numTriloops = 0;
  size = 16;
  triloopEnthalpies = (triloopE *)calloc(16, sizeof(struct triloopE));

  while (fscanf(hFile, "%5s %lg", triloopEnthalpies[numTriloops].loop, &energy) == 2)
    {
      for (i = 0; i < 5; ++i)
	triloopEnthalpies[numTriloops].loop[i] = toNum(triloopEnthalpies[numTriloops].loop[i]);
      triloopEnthalpies[numTriloops].energy = energy;
      ++numTriloops;
      if (numTriloops == size)
	{
	  size *= 2;
	  triloopEnthalpies = (triloopE *)realloc(triloopEnthalpies, size * sizeof(struct triloopE));
	}
    }

  triloopEnthalpies = (triloopE *)realloc(triloopEnthalpies, numTriloops * sizeof(struct triloopE));

  fclose(hFile);
}

void Energy::combineTriloop()
{
  int i;

  for (i = 0; i < numTriloops; ++i)
    {
      memcpy(g_triloop[i].loop, triloopEnergies[i].loop, 5);
      g_triloop[i].energy = scale(tRatio * triloopEnergies[i].energy + (1.0 - tRatio) * triloopEnthalpies[i].energy);
    }
}

void Energy::calculateTriloop()
{
  int i;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < numTriloops; ++i)
    z_triloop[i].energy = exp(-g_triloop[i].energy / RT);
}

void Energy::loadTriloopSuffix(triloop** tri, int* num)
{
  int i, size;
  double energy;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "triloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  *num = 0;
  size = 16;
  *tri = (triloop *)calloc(16, sizeof(struct triloop));

  while (fscanf(file, "%5s %lg", (*tri)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 5; ++i)
	(*tri)[*num].loop[i] = toNum((*tri)[*num].loop[i]);
      (*tri)[*num].energy = scale(energy);
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *tri = (triloop *)realloc(*tri, size * sizeof(struct triloop));
	}
    }

  *tri = (triloop *)realloc(*tri, *num * sizeof(struct triloop));

  fclose(file);
}

void Energy::loadTloop()
{
  FILE *gFile, *hFile;
  int i, size;
  double energy;

  if (NA)
    gFile = openFile((char *)"tloop.DGD");
  else
    gFile = openFile((char *)"tloop.DG");

  numTloops = 0;
  size = 16;
  tloopEnergies = (tloopE *)calloc(16, sizeof(struct tloopE));

  while (fscanf(gFile, "%6s %lg", tloopEnergies[numTloops].loop, &energy) == 2)
    {
      for (i = 0; i < 6; ++i)
	tloopEnergies[numTloops].loop[i] = toNum(tloopEnergies[numTloops].loop[i]);
      tloopEnergies[numTloops].energy = energy;
      ++numTloops;
      if (numTloops == size)
	{
	  size *= 2;
	  tloopEnergies = (tloopE *)realloc(tloopEnergies, size * sizeof(struct tloopE));
	}
    }

  tloopEnergies = (tloopE *)realloc(tloopEnergies, numTloops * sizeof(struct tloopE));

  fclose(gFile);

  if (NA)
    hFile = openFile((char *)"tloop.DHD");
  else
    hFile = openFile((char *)"tloop.DH");

  numTloops = 0;
  size = 16;
  tloopEnthalpies = (tloopE *)calloc(16, sizeof(struct tloopE));

  while (fscanf(hFile, "%6s %lg", tloopEnthalpies[numTloops].loop, &energy) == 2)
    {
      for (i = 0; i < 6; ++i)
	tloopEnthalpies[numTloops].loop[i] = toNum(tloopEnthalpies[numTloops].loop[i]);
      tloopEnthalpies[numTloops].energy = energy;
      ++numTloops;
      if (numTloops == size)
	{
	  size *= 2;
	  tloopEnthalpies = (tloopE *)realloc(tloopEnthalpies, size * sizeof(struct tloopE));
	}
    }

  tloopEnthalpies = (tloopE *)realloc(tloopEnthalpies, numTloops * sizeof(struct tloopE));

  fclose(hFile);
}

void Energy::combineTloop()
{
  int i;

  for (i = 0; i < numTloops; ++i)
    {
      memcpy(g_tloop[i].loop, tloopEnergies[i].loop, 6);
      g_tloop[i].energy = scale(tRatio * tloopEnergies[i].energy + (1.0 - tRatio) * tloopEnthalpies[i].energy);
    }
}

void Energy::calculateTloop()
{
  int i;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < numTloops; ++i)
    z_tloop[i].energy = exp(-g_tloop[i].energy / RT);
}

void Energy::loadTloopSuffix(struct tloop** tl, int* num)
{
  int i, size;
  double energy;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(6 + strlen(suffix) + 1);
  strcpy(buffer, "tloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  *num = 0;
  size = 16;
  *tl = (tloop *)calloc(16, sizeof(struct tloop));

  while (fscanf(file, "%6s %lg", (*tl)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 6; ++i)
	(*tl)[*num].loop[i] = toNum((*tl)[*num].loop[i]);
      (*tl)[*num].energy = scale(energy);
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *tl = (tloop *)realloc(*tl, size * sizeof(struct tloop));
	}
    }

  *tl = (tloop *)realloc(*tl, *num * sizeof(struct tloop));

  fclose(file);
}

void Energy::loadHexaloop()
{
  FILE *gFile, *hFile;
  int i, size;
  double energy;

  if (NA)
    gFile = openFile((char *)"hexaloop.DGD");
  else
    gFile = openFile((char *)"hexaloop.DG");

  numHexaloops = 0;
  size = 16;
  hexaloopEnergies = (hexaloopE *)calloc(16, sizeof(struct hexaloopE));

  while (fscanf(gFile, "%8s %lg", hexaloopEnergies[numHexaloops].loop, &energy) == 2)
    {
      for (i = 0; i < 8; ++i)
	hexaloopEnergies[numHexaloops].loop[i] = toNum(hexaloopEnergies[numHexaloops].loop[i]);
      hexaloopEnergies[numHexaloops].energy = energy;
      ++numHexaloops;
      if (numHexaloops == size)
	{
	  size *= 2;
	  hexaloopEnergies = (hexaloopE *)realloc(hexaloopEnergies, size * sizeof(struct hexaloopE));
	}
    }

  hexaloopEnergies = (hexaloopE *)realloc(hexaloopEnergies, numHexaloops * sizeof(struct hexaloopE));

  fclose(gFile);

  if (NA)
    hFile = openFile((char *)"hexaloop.DHD");
  else
    hFile = openFile((char *)"hexaloop.DH");

  numHexaloops = 0;
  size = 16;
  hexaloopEnthalpies = (hexaloopE *)calloc(16, sizeof(struct hexaloopE));

  while (fscanf(hFile, "%8s %lg", hexaloopEnthalpies[numHexaloops].loop, &energy) == 2)
    {
      for (i = 0; i < 8; ++i)
	hexaloopEnthalpies[numHexaloops].loop[i] = toNum(hexaloopEnthalpies[numHexaloops].loop[i]);
      hexaloopEnthalpies[numHexaloops].energy = energy;
      ++numHexaloops;
      if (numHexaloops == size)
	{
	  size *= 2;
	  hexaloopEnthalpies = (hexaloopE *)realloc(hexaloopEnthalpies, size * sizeof(struct hexaloopE));
	}
    }

  hexaloopEnthalpies = (hexaloopE *)realloc(hexaloopEnthalpies, numHexaloops * sizeof(struct hexaloopE));

  fclose(hFile);
}

void Energy::combineHexaloop()
{
  int i;

  for (i = 0; i < numHexaloops; ++i)
    {
      memcpy(g_hexaloop[i].loop, hexaloopEnergies[i].loop, 8);
      g_hexaloop[i].energy = scale(tRatio * hexaloopEnergies[i].energy + (1.0 - tRatio) * hexaloopEnthalpies[i].energy);
    }
}

void Energy::calculateHexaloop()
{
  int i;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < numHexaloops; ++i)
    z_hexaloop[i].energy = exp(-g_hexaloop[i].energy / RT);
}

void Energy::loadHexaloopSuffix(struct hexaloop** hexal, int* num)
{
  int i, size;
  double energy;
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "hexaloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  *num = 0;
  size = 16;
  *hexal = (hexaloop *)calloc(16, sizeof(struct hexaloop));

  while (fscanf(file, "%8s %lg", (*hexal)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 8; ++i)
	(*hexal)[*num].loop[i] = toNum((*hexal)[*num].loop[i]);
      (*hexal)[*num].energy = scale(energy);
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *hexal = (hexaloop *)realloc(*hexal, size * sizeof(struct hexaloop));
	}
    }

  *hexal = (hexaloop *)realloc(*hexal, *num * sizeof(struct hexaloop));

  fclose(file);
}

void Energy::loadRTSuffix(double* RT)
{
  char* buffer;
  FILE* file;

  buffer = (char *)alloc.xmalloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "miscloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  if (fscanf(file, "%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%lf", RT) < 20)
    *RT = 37;

  *RT = R * (*RT + 273.15);

  fclose(file);
}

void Energy::loadInteraction()
{
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile((char *)"interaction.DGD");
      hFile = openFile((char *)"interaction.DHD");
    }
  else
    {
      gFile = openFile((char *)"interaction.DG");
      hFile = openFile((char *)"interaction.DH");
    }

  int tmp = fscanf(gFile, "%lg%lg%lg%lg", &kissingEnergies[0], &kissingEnergies[1], &kissingEnergies[2], &kissingEnergies[3]);
  tmp = fscanf(hFile, "%lg%lg%lg%lg", &kissingEnthalpies[0], &kissingEnthalpies[1], &kissingEnthalpies[2], &kissingEnthalpies[3]);

  tmp = fscanf(gFile, "%lg%lg%lg", &intstackEnergies[0], &intstackEnergies[1], &aupenEnergy);
  tmp = fscanf(hFile, "%lg%lg%lg", &intstackEnthalpies[0], &intstackEnthalpies[1], &aupenEnthalpy);

  fclose(gFile);
  fclose(hFile);
}

void Energy::combineInteraction()
{
  for(int i = 0; i < 4; i++)
     g_kissing[i] = tRatio*kissingEnergies[i] + (1-tRatio)*kissingEnthalpies[i];

  for(int i = 0; i < 2; i++)
     g_intstack[i] = tRatio*intstackEnergies[i] + (1-tRatio)*intstackEnthalpies[i];

  g_aupen = tRatio*aupenEnergy + (1-tRatio)*aupenEnthalpy;
}

void Energy::loadEntropiesEnthalpies()
{
	if (suffix)
	{
		double RT;		
		loadRTSuffix(&RT);
      		double t = RT / R - 273.15;
      		setTemperature(t);

      		loadStackSuffix();
		loadDangleSuffix();
		loadLoopSuffix();
		loadSint2Suffix();
		loadAsint1x2Suffix();
		loadSint4Suffix();
		loadTstackhSuffix();
		loadTstackiSuffix();
  		loadTstackmSuffix();
		loadTstackeSuffix();
		loadMiscSuffix();
		loadTriloopSuffix(&g_triloop, &numTriloops);
		loadTloopSuffix(&g_tloop, &numTloops);
      		loadHexaloopSuffix(&g_hexaloop, &numHexaloops);
		loadMultiSuffix();
		loadMulti2Suffix();
		loadCoaxialSuffix();
		loadTstackcoaxSuffix();
		loadCoaxstackSuffix();
	}
	else
	{
		loadStack();
		symmetryCheckStack(stackEnergies, (char *)"energy");
		/* symmetryCheckStack(stackEnthalpies, "enthalpy"); */
		loadDangle();
		loadLoop();
		loadSint2();
		symmetryCheckSint2(sint2Energies, (char *)"energy");
		/* symmetryCheckSint2(sint2Enthalpies, "enthalpy"); */
		loadAsint1x2();
		loadSint4();
		symmetryCheckSint4(sint4Energies, (char *)"energy");
		/* symmetryCheckSint4(sint4Enthalpies, "enthalpy"); */
		loadTstackh();
		loadTstacki();
		loadTstackm();
		loadTstacke();
		loadTriloop();
		loadTloop();
		loadHexaloop();
		loadMulti();
		loadMulti2();
		loadMisc();
		loadInteraction();
	}
	g_triloop = (struct triloop*) alloc.xcalloc(numTriloops, sizeof(struct triloop));
	g_tloop = (struct tloop*) alloc.xcalloc(numTloops, sizeof(struct tloop));
	g_hexaloop = (struct hexaloop*) alloc.xcalloc(numHexaloops, sizeof(struct hexaloop));

	z_triloop = (struct triloop*) alloc.xcalloc(numTriloops, sizeof(struct triloop));
	z_tloop = (struct tloop*) alloc.xcalloc(numTloops, sizeof(struct tloop));
	z_hexaloop = (struct hexaloop*) alloc.xcalloc(numHexaloops, sizeof(struct hexaloop));
}

void Energy::computeEnergies()
{
	if (!suffix)
	{
		combineStack();
		combineDangle();
		combineLoop();
		combineSint2();
		combineAsint1x2();
		combineSint4();
		combineTstack(tstackhEnergies, tstackhEnthalpies, g_tstackh);
		combineTstack(tstackiEnergies, tstackiEnthalpies, g_tstacki);
		combineTstack2(tstackmEnergies, tstackmEnthalpies, g_tstackm);
		combineTstack2(tstackeEnergies, tstackeEnthalpies, g_tstacke);
		combineMulti();
		combineMisc();
      		combineTriloop();
		combineTloop();
		combineHexaloop();
		combineInteraction();
		makeAUPenalty();
	}
/*	memcpy(z_triloop, g_triloop, numTriloops*sizeof(struct triloop));
	memcpy(z_tloop, g_tloop, numTloops*sizeof(struct tloop));
	memcpy(z_hexaloop, g_hexaloop, numHexaloops*sizeof(struct hexaloop));
	calculateZOfEnergies();*/
}

void Energy::calculateZOfEnergies()
{
	
	calculateStack(z_stack, g_stack);
	if (g_nodangle)
		calculateZeroDangle();
	else if (zip)
		calculateZipDangle();
	else
		calculateDangle();
	calculateLoop();
	calculateSint2();
	calculateAsint1x2();
	calculateSint4();
	calculateStack(z_tstackh, g_tstackh);
	calculateStack(z_tstacki, g_tstacki);
	if (g_nodangle)
		calculateZeroStack2(z_tstackm);
	else if (zip)
		calculateZipStack2(z_tstackm);
	else
		calculateStack2(z_tstackm, g_tstackm);
	if (g_nodangle)
		calculateZeroStack2(z_tstacke);
	else if (zip)
		calculateZipStack2(z_tstacke);
	else
		calculateStack2(z_tstacke, g_tstacke);
	calculateMulti();
	calculateMisc();
	calculateTriloop();
	calculateTloop();
	calculateHexaloop();
//	makeAUPenalty(1);
}

double Energy::Ed5(int i, int j, int k, unsigned char *seq)
{
  if (g_nodangle)
    return INFINITY;

  if (k == j - 1)
    return g_dangle5[seq[i]][seq[j]][seq[k]];
  else if (k == i - 1)
    return g_dangle5[seq[j]][seq[i]][seq[k]];
  else
    {
      fprintf(stderr, "Error: Ed5(%d, %d, %d)\n", i, j, k);
      return 0.0;
    }
}

double Energy::Ed3(int i, int j, int k, unsigned char *seq)
{
  if (g_nodangle)
    return INFINITY;

  if (k == i + 1)
    return g_dangle3[seq[i]][seq[j]][seq[k]];
  else if (k == j + 1)
    return g_dangle3[seq[j]][seq[i]][seq[k]];
  else
    {
      fprintf(stderr, "Error: Ed3(%d, %d, %d)\n", i, j, k);
      return 0.0;
    }
}

double Energy::Etstackm(int i, int j, unsigned char *seq)
{
  if (g_nodangle)
    return INFINITY;

  return g_tstackm[seq[i]][seq[j]][seq[i+1]][seq[j]];
}

double Energy::Etstacke(int i, int j, unsigned char *seq)
{
  if (g_nodangle)
    return INFINITY;

  return g_tstacke[seq[i]][seq[j]][seq[i+1]][seq[j]];
}


double Energy::Eh(int i, int j, unsigned char *seq)
{
  int loopSize = j - i - 1;
  double energy = 0.0;
  int k;

  if (loopSize < MIN_HAIRPIN_SIZE)
    return INFINITY;

  if (loopSize <= 30)
    energy = g_hairpinLoop[loopSize - 1];
  else
    energy = g_hairpinLoop[29]  + g_misc[12]*log((double) loopSize / 30);

  if (loopSize > 3)
    energy += g_tstackh[seq[i]][seq[j]][seq[i + 1]][seq[j-1]];
  else
    energy += auPenalty(i, j, seq);

  if (loopSize == 3)
    {
      triloop* loop;
      if (numTriloops)
	if ((loop = (triloop *)bsearch(seq + i, g_triloop, numTriloops, sizeof(struct triloop), triloopcmp)))
	  energy += loop->energy;
    }
  else if (loopSize == 4)
    {
      struct tloop* loop;
      if (numTloops)
	if ((loop = (tloop *)bsearch(seq + i, g_tloop, numTloops, sizeof(struct tloop), tloopcmp)))
	  energy += loop->energy;
    }
  else if (loopSize == 6)
    {
      struct hexaloop* loop;
      if (numHexaloops)
	if ((loop = (hexaloop *)bsearch(seq + i, g_hexaloop, numHexaloops, sizeof(struct hexaloop), hexaloopcmp)))
	  energy += loop->energy;
    }

  /* GGG */
  if (i >= 3 && seq[i - 2] == 2 && seq[i - 1] == 2 && seq[i] == 2 && seq[j] == 3)
    energy += g_misc[8];

  /* poly-C */
  if (loopSize == 3 && seq[i + 1] == 1 && seq[i + 2] == 1 && seq[i + 3] == 1)
    energy += g_misc[11];
  else
    {
      for (k = 1; k <= loopSize; ++k)
	if (seq[i + k] != 1)
	  return energy;
      energy += g_misc[9] * loopSize + g_misc[10];
    }

  return energy;
}

double Energy::Es(int i, int j, unsigned char *seq)
{
  if (i >= j)
    return INFINITY;

  return g_stack[seq[i]][seq[j]][seq[i + 1]][seq[j - 1]];
}


double Energy::Ebi(int i, int j, int ii, int jj, unsigned char *seq)
{
  int loopSize1, loopSize2;
  double loopEnergy, asPenalty;

    if(basePairIndex(seq[i], seq[j]) == 6 || basePairIndex(seq[ii], seq[jj]) == 6)
      return INFINITY;  

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;


  if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	return g_bulgeLoop[0] + g_stack[seq[i]][seq[j]][seq[ii]][seq[jj]];
      else if (loopSize2 <= 30)
	return g_bulgeLoop[loopSize2 - 1] + auPenalty(i, j, seq) + auPenalty(ii, jj, seq);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize2 / 30) + auPenalty(i, j, seq) + auPenalty(ii, jj, seq);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return g_bulgeLoop[0] + g_stack[seq[i]][seq[j]][seq[ii]][seq[jj]];
      else if (loopSize1 <= 30)
	return g_bulgeLoop[loopSize1 - 1] + auPenalty(i, j, seq) + auPenalty(ii, jj, seq);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize1 / 30) + auPenalty(i, j, seq) + auPenalty(ii, jj, seq);
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    return g_sint2[basePairIndex(seq[i], seq[j])][basePairIndex(seq[ii], seq[jj])][seq[i + 1]][seq[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return g_asint1x2[basePairIndex(seq[i], seq[j])][basePairIndex(seq[ii], seq[jj])][seq[i + 1]][seq[j - 1]][seq[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return g_asint1x2[basePairIndex(seq[jj], seq[ii])][basePairIndex(seq[j], seq[i])][seq[jj + 1]][seq[ii - 1]][seq[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return g_sint4[basePairIndex(seq[i], seq[j])][basePairIndex(seq[ii], seq[jj])][seq[i + 1]][seq[j - 1]][seq[i + 2]][seq[j - 2]];
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = g_interiorLoop[29] + g_misc[12] * log((double) (loopSize1 + loopSize2) / 30);
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy += g_tstacki[seq[i]][seq[j]][0][0];
	  loopEnergy += g_tstacki[seq[jj]][seq[ii]][0][0];
	}
      else
	{
	  loopEnergy += g_tstacki[seq[i]][seq[j]][seq[i + 1]][seq[j - 1]];
	  loopEnergy += g_tstacki[seq[jj]][seq[ii]][seq[jj + 1]][seq[ii - 1]];
	}
      asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4])
	asPenalty = g_misc[4];
      loopEnergy += asPenalty;

      return loopEnergy;
    }
}

double Energy::Eibi(int loopSize1, int loopSize2, unsigned char seqi, unsigned char seqj, unsigned char seqii, unsigned char seqjj, 
unsigned char seqip1, unsigned char seqjp1, unsigned char seqip2, unsigned char seqjp2, 
unsigned char seqiim1, unsigned char seqjjm1, unsigned char seqiim2, unsigned char seqjjm2)
{
  double loopEnergy, asPenalty;

    if(basePairIndex(seqi, seqj) == 6 || basePairIndex(seqii, seqjj) == 6)
      return INFINITY;  

  if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	return g_bulgeLoop[0] + g_stack[seqi][seqj][seqii][seqjj];
      else if (loopSize2 <= 30)
	return g_bulgeLoop[loopSize2 - 1] + auPenalty(seqi, seqj) + auPenalty(seqii, seqjj);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize2 / 30) + auPenalty(seqi, seqj) + auPenalty(seqii, seqjj);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return g_bulgeLoop[0] + g_stack[seqi][seqj][seqii][seqjj];
      else if (loopSize1 <= 30)
	return g_bulgeLoop[loopSize1 - 1] + auPenalty(seqi, seqj) + auPenalty(seqii, seqjj);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize1 / 30) + auPenalty(seqi, seqj) + auPenalty(seqii, seqjj);
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    return g_sint2[basePairIndex(seqi, seqj)][basePairIndex(seqii, seqjj)][seqip1][seqjp1];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return g_asint1x2[basePairIndex(seqi, seqj)][basePairIndex(seqii, seqjj)][seqip1][seqjp1][seqjp2];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return g_asint1x2[basePairIndex(seqjj, seqii)][basePairIndex(seqj, seqi)][seqjjm1][seqiim1][seqiim2];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return g_sint4[basePairIndex(seqi, seqj)][basePairIndex(seqii, seqjj)][seqip1][seqjp1][seqip2][seqjp2];
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = g_interiorLoop[29] + g_misc[12] * log((double) (loopSize1 + loopSize2) / 30);
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy += g_tstacki[seqi][seqj][0][0];
	  loopEnergy += g_tstacki[seqjj][seqii][0][0];
	}
      else
	{
	  loopEnergy += g_tstacki[seqi][seqj][seqip1][seqjp1];
	  loopEnergy += g_tstacki[seqjj][seqii][seqjjm1][seqiim1];
	}
      asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4])
	asPenalty = g_misc[4];
      loopEnergy += asPenalty;

      return loopEnergy;
    }
}

double Energy::Emulti(int penalty, int unpaired, int loops)
{
	return g_multi[0]*penalty + g_multi[1] * unpaired + g_multi[2] * loops;
}

/* Warning: Emulti2 cannot be used with Dynamic Programming */
double Energy::Emulti2(int penalty, int unpaired, int loops)
{
	if (unpaired > 6)
		return g_multi2[0]*penalty + 6 * g_multi2[1] + g_misc[12] * log((double) (unpaired) / 6) + g_multi2[2] * loops;
	else
		return g_multi[0]*penalty + g_multi[1] * unpaired + g_multi[2] * loops;
}

double Energy::Ekissing(int penalty, int unpaired, int loops)
{
	return g_kissing[0]*penalty + g_kissing[1]*unpaired + g_kissing[2]*loops;
}

double Energy::Ekissingbi(int i, int j, int ii, int jj, unsigned char *seq1, unsigned char *seq2)
{
	int loopSize1 = ii-i-1;
	int loopSize2 = jj-j-1;

	if(loopSize1 < 0 || loopSize2 < 0)
		return INFINITY;

	unsigned char seqi = seq1[i]; 
	unsigned char seqj = seq2[j];  
	unsigned char seqii = seq1[ii];  
	unsigned char seqjj = seq2[jj];  
	unsigned char seqip1 = (i + 1 <= ii) ? seq1[i + 1] : 6;  
	unsigned char seqjp1 = (j + 1 <= jj) ? seq2[j + 1] : 6; 
	unsigned char seqip2 = (i + 2 <= ii) ? seq1[i + 2] : 6;   
	unsigned char seqjp2 = (j + 2 <= jj) ? seq2[j + 2] : 6;  
	unsigned char seqiim1 = (i + 1 <= ii) ? seq1[ii - 1] : 6;   
	unsigned char seqjjm1 = (j + 1 <= jj) ? seq2[jj - 1] : 6;  
	unsigned char seqiim2 = (i + 2 <= ii) ? seq1[ii - 2] : 6;    
	unsigned char seqjjm2 = (j + 2 <= jj) ? seq2[jj - 2] : 6;   

	if(inttype == 0)
		return g_kissing[3]*Eibi(loopSize1, loopSize2, seqi, seqj, seqii, seqjj, seqip1, seqjp1, seqip2, seqjp2, seqiim1, seqjjm1, seqiim2, seqjjm2);
	else
	{
		double t = Eibi(loopSize1, loopSize2, seqi, seqj, seqii, seqjj, seqip1, seqjp1, seqip2, seqjp2, seqiim1, seqjjm1, seqiim2, seqjjm2);
		return t + g_kissing[3]*fabs(t);
	}

}

double Energy::Erawstack(int i1, int i2, unsigned char *seq1, unsigned char *seq2)
{
	return g_stack[seq1[i1]][seq2[i2]][seq1[i1 + 1]][seq2[i2 + 1]];
}

double Energy::Erawstack(unsigned char s1, unsigned char s2, unsigned char n1, unsigned char n2)
{
	return g_stack[s1][s2][n1][n2];
}

double Energy::Erawbi(int i, int j, int ii, int jj, unsigned char *seq1, unsigned char *seq2)
{
	int loopSize1 = ii-i-1;
	int loopSize2 = jj-j-1;

	if(loopSize1 < 0 || loopSize2 < 0)
		return INFINITY;

	unsigned char seqi = seq1[i]; 
	unsigned char seqj = seq2[j];  
	unsigned char seqii = seq1[ii];  
	unsigned char seqjj = seq2[jj];  
	unsigned char seqip1 = (i + 1 <= ii) ? seq1[i + 1] : 6;  
	unsigned char seqjp1 = (j + 1 <= jj) ? seq2[j + 1] : 6; 
	unsigned char seqip2 = (i + 2 <= ii) ? seq1[i + 2] : 6;   
	unsigned char seqjp2 = (j + 2 <= jj) ? seq2[j + 2] : 6;  
	unsigned char seqiim1 = (i + 1 <= ii) ? seq1[ii - 1] : 6;   
	unsigned char seqjjm1 = (j + 1 <= jj) ? seq2[jj - 1] : 6;  
	unsigned char seqiim2 = (i + 2 <= ii) ? seq1[ii - 2] : 6;    
	unsigned char seqjjm2 = (j + 2 <= jj) ? seq2[jj - 2] : 6;   

	return Eibi(loopSize1, loopSize2, seqi, seqj, seqii, seqjj, seqip1, seqjp1, seqip2, seqjp2, seqiim1, seqjjm1, seqiim2, seqjjm2);
}

double Energy::Ekissingstack(int i1, int i2, unsigned char *seq1, unsigned char *seq2)
{
	if(inttype == 0)
		return g_kissing[3]*(g_stack[seq1[i1]][seq2[i2]][seq1[i1 + 1]][seq2[i2 + 1]]);
	else
	{
		double t = g_stack[seq1[i1]][seq2[i2]][seq1[i1 + 1]][seq2[i2 + 1]];
		return  t + g_kissing[3]*fabs(t);
	}
}

double Energy::Eintstack(int penalty, int i1, int i2, unsigned char *seq1, unsigned char *seq2)
{

	if(inttype == 0)
		return g_intstack[0]*penalty + g_intstack[1]*(g_stack[seq1[i1]][seq2[i2]][seq1[i1 + 1]][seq2[i2 + 1]]);
	else
	{
		double t = g_stack[seq1[i1]][seq2[i2]][seq1[i1 + 1]][seq2[i2 + 1]];
		return g_intstack[0]*penalty + t + g_intstack[1]*fabs(t);
	}
	
}

double Energy::Eintstackpenalty()
{
	return g_intstack[0];
}

double Energy::Eintbi(int i, int j, int ii, int jj, unsigned char *seq1, unsigned char *seq2)
{
	int loopSize1 = ii-i-1;
	int loopSize2 = jj-j-1;

	if(loopSize1 < 0 || loopSize2 < 0)
		return INFINITY;

	unsigned char seqi = seq1[i]; 
	unsigned char seqj = seq2[j];  
	unsigned char seqii = seq1[ii];  
	unsigned char seqjj = seq2[jj];  
	unsigned char seqip1 = (i + 1 <= ii) ? seq1[i + 1] : 6;  
	unsigned char seqjp1 = (j + 1 <= jj) ? seq2[j + 1] : 6; 
	unsigned char seqip2 = (i + 2 <= ii) ? seq1[i + 2] : 6;   
	unsigned char seqjp2 = (j + 2 <= jj) ? seq2[j + 2] : 6;  
	unsigned char seqiim1 = (i + 1 <= ii) ? seq1[ii - 1] : 6;   
	unsigned char seqjjm1 = (j + 1 <= jj) ? seq2[jj - 1] : 6;  
	unsigned char seqiim2 = (i + 2 <= ii) ? seq1[ii - 2] : 6;    
	unsigned char seqjjm2 = (j + 2 <= jj) ? seq2[jj - 2] : 6;   

	if(inttype == 0)
		return g_intstack[1]*Eibi(loopSize1, loopSize2, seqi, seqj, seqii, seqjj, seqip1, seqjp1, seqip2, seqjp2, seqiim1, seqjjm1, seqiim2, seqjjm2);
	else
	{
		double t = Eibi(loopSize1, loopSize2, seqi, seqj, seqii, seqjj, seqip1, seqjp1, seqip2, seqjp2, seqiim1, seqjjm1, seqiim2, seqjjm2);
		return t + g_intstack[1]*fabs(t);
	}
}

double Energy::Eaupenalty()
{
  return g_aupen;
}

double Energy::auPenalty(int i, int j, unsigned char *seq)
{
  return auPenalty(seq[i], seq[j]);
}

double Energy::auPenalty(unsigned char seqi, unsigned char seqj)
{
  return g_aup[seqi][seqj];
}

void Energy::getIntParams(double *p)
{
	p[0] = g_kissing[0];
	p[1] = g_kissing[1];
	p[2] = g_kissing[2];
	p[3] = g_kissing[3];
	p[4] = g_intstack[0];
	p[5] = g_intstack[1];
}

void Energy::setIntParams(double *p)
{
	g_kissing[0] = p[0];
	g_kissing[1] = p[1];
	g_kissing[2] = p[2];
	g_kissing[3] = p[3];
	g_intstack[0] = p[4];
	g_intstack[1] = p[5];
}
