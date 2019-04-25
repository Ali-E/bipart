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
	Desc: Energy class computes energies based on UNAFOLD 3.6 code.

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow, 
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: June 2, 2009
	Ver: 1.6
*/


#ifndef ENERGY_H
#define ENERGY_H

#define MIN_HAIRPIN_SIZE 3 /* minimum size of hairpin loop */

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


#include "config.h"
#include "alloc.h"

struct triloopE { char loop[5]; double energy; };
struct triloop { char loop[5]; double energy; };

struct tloopE { char loop[6]; double energy; };
struct tloop { char loop[6]; double energy; };

struct hexaloopE { char loop[8]; double energy; };
struct hexaloop { char loop[8]; double energy; };

#ifndef __SGI
#ifndef __SUNOS
#ifndef __LINUX
bool finite(double x);
#endif
#endif
#endif


class Energy {
public:
	void setTemperature(double);
	Energy(char *dataDir, int _NA, int _polymer, double _naConc, double _mgConc, char *suf, double temperature, int z, int nodangle, int itype);
	~Energy();
	void computeEnergies();
	void calculateZOfEnergies();
	double Ed3(int i, int j, int k, unsigned char *seq);
	double Ed5(int i, int j, int k, unsigned char *seq);
	double Eh(int i, int j, unsigned char *seq);
	double Es(int i, int j, unsigned char *seq);
	double Ebi(int i, int j, int ii, int jj, unsigned char *seq);
	double Eibi(int loopSize1, int loopSize2, unsigned char seqi, unsigned char seqj, unsigned char seqii, unsigned char seqjj, 
unsigned char seqip1, unsigned char seqjp1, unsigned char seqip2, unsigned char seqjp2, 
unsigned char seqiim1, unsigned char seqjjm1, unsigned char seqiim2, unsigned char seqjjm2);
	double Emulti(int, int ,int);
	/* Warning: Emulti2 cannot be used with Dynamic Programming */
	double Emulti2(int, int ,int);
	double Ekissing(int, int, int);
	double Ekissingstack(int, int, unsigned char *, unsigned char *);
	double Ekissingbi(int i, int j, int ii, int jj, unsigned char *seq1, unsigned char *seq2);	
	double Erawstack(int, int, unsigned char *, unsigned char *);
	double Erawstack(unsigned char, unsigned char, unsigned char, unsigned char);
	double Erawbi(int i, int j, int ii, int jj, unsigned char *seq1, unsigned char *seq2);	
	double Eintstack(int, int, int, unsigned char *, unsigned char *);
	double Eintstackpenalty();
	double Eintbi(int i, int j, int ii, int jj, unsigned char *seq1, unsigned char *seq2);	
	double Eaupenalty();
	double Etstackm(int i, int j, unsigned char *seq);
	double Etstacke(int i, int j, unsigned char *seq);

	double auPenalty(int i, int j, unsigned char *seq);
	double auPenalty(unsigned char, unsigned char);

	void getIntParams(double *);
	void setIntParams(double *);

	double getRT() {return R * (tRatio * 310.15);};
	void getParams(int *_NA, int *_polymer, double *_naConc, double *_mgConc, char **suf, double *temperature, int *z, int *nodangle, int *itype);
	bool base_pair(unsigned char a, unsigned char b);

private:
	double R;
	char BASES[5];
	char BASE_PAIRS[6][4];
	char PKGDATADIR[100];
	Alloc alloc;

	int NA, polymer;
	double tRatio;
	double naConc, mgConc;
	double saltCorrection;
	char *suffix;
	double scaleFactor;
	int zip, g_nodangle;
	int inttype;

	double g_dangle3[5][5][6];
	double g_dangle5[5][5][6];
	double g_stack[5][5][5][5];
	double g_interiorLoop[30];
	double g_bulgeLoop[30];
	double g_hairpinLoop[30];
	double g_sint2[7][7][5][5];
	double g_asint1x2[7][7][5][5][5];
	double g_sint4[7][7][5][5][5][5];
	double g_tstackh[5][5][5][5];
	double g_tstacki[5][5][5][5];
	double g_tstackm[5][5][6][6];
	double g_tstacke[5][5][6][6];
	double g_misc[13];
	double g_aup[5][5];
	struct triloop* g_triloop; int numTriloops;
	struct tloop* g_tloop; int numTloops;
	struct hexaloop* g_hexaloop; int numHexaloops;
	double g_multi[3];
	double g_multi2[3];
	double g_coaxial[5][5][5][5];
	double g_tstackcoax[5][5][5][5];
	double g_coaxstack[5][5][5][5];
	double g_kissing[4]; // {beta1, beta2, beta3, sigma} please refer to the paper for more info
	double g_intstack[2]; // {beta1, sigma} please refer to the paper for more info
	double g_aupen; // A-U stack terminal penalty, please refer to Biochem. Vol. 37, No. 42, 1998 by Xia et al.

	double z_dangle3[5][5][6];
	double z_dangle5[5][5][6];
	double z_stack[5][5][5][5];
	double z_interiorLoop[30];
	double z_bulgeLoop[30];
	double z_hairpinLoop[30];
	double z_sint2[7][7][5][5];
	double z_asint1x2[7][7][5][5][5];
	double z_sint4[7][7][5][5][5][5];
	double z_tstackh[5][5][5][5];
	double z_tstacki[5][5][5][5];
	double z_tstackm[5][5][6][6];
	double z_tstacke[5][5][6][6];
	double z_misc[13];
	double z_aup[5][5];
	struct triloop* z_triloop;
	struct tloop* z_tloop;
	struct hexaloop* z_hexaloop;
	double z_multi[3];
	double z_multi2[3];
	double z_coaxial[5][5][5][5];
	double z_tstackcoax[5][5][5][5];
	double z_coaxstack[5][5][5][5];
	double z_kissing[4]; // {beta1, beta2, beta3, sigma} please refer to the paper for more info
	double z_intstack[2]; // {beta1, sigma} please refer to the paper for more info

	double dangleEnergies3[4][4][4];
	double dangleEnthalpies3[5][5][6];
	double dangleEnergies5[4][4][4];
	double dangleEnthalpies5[5][5][6];
	double stackEnergies[4][4][4][4];
	double stackEnthalpies[5][5][5][5];
	double interiorLoopEnergies[30];
	double bulgeLoopEnergies[30];
	double hairpinLoopEnergies[30];
	double interiorLoopEnthalpies[30];
	double bulgeLoopEnthalpies[30];
	double hairpinLoopEnthalpies[30];
	double sint2Energies[6][6][4][4];
	double sint2Enthalpies[7][7][5][5];
	double asint1x2Energies[6][6][4][4][4];
	double asint1x2Enthalpies[7][7][5][5][5];
	double sint4Energies[6][6][4][4][4][4];
	double sint4Enthalpies[7][7][5][5][5][5];
	double tstackiEnergies[4][4][4][4];
	double tstackiEnthalpies[5][5][5][5];
	double tstackhEnergies[4][4][4][4];
	double tstackhEnthalpies[5][5][5][5];
	double tstackmEnergies[4][4][4][4];
	double tstackmEnthalpies[5][5][6][6];
	double tstackeEnergies[4][4][4][4];
	double tstackeEnthalpies[5][5][6][6];
	double miscEnergies[13];
	double miscEnthalpies[13];
	struct triloopE* triloopEnergies;
	struct triloopE* triloopEnthalpies;
	struct tloopE* tloopEnergies;
	struct tloopE* tloopEnthalpies;
	struct hexaloopE* hexaloopEnergies;
	struct hexaloopE* hexaloopEnthalpies;
	double multiEnergies[3];
	double multiEnthalpies[3];
	double multiEnergies2[3];
	double multiEnthalpies2[3];
	double kissingEnergies[4];
	double kissingEnthalpies[4];
	double intstackEnergies[2];
	double intstackEnthalpies[2];
	double aupenEnergy, aupenEnthalpy; // A-U stack terminal penalty, please refer to Biochem. Vol. 37, No. 42, 1998 by Xia et al.

	void loadEntropiesEnthalpies();

	double ion();

	FILE* openFile(char*);

	void loadInteraction();

	void combineInteraction();

	void loadStack();

	void combineStack();

	void calculateStack(double stack[5][5][5][5], double estack[5][5][5][5]);

	void calculateStack2(double stack[5][5][6][6], double estack[5][5][6][6]);

	void calculateZipStack2(double stack[5][5][6][6]);

	void calculateZeroStack2(double stack[5][5][6][6]);

	void calculateInfStack2(double stack[5][5][6][6]);

	void loadStackSuffix();

	void symmetryCheckStack(double stack[4][4][4][4], char* which);

	double estimateScale(double stack[5][5][5][5]);

	void loadDangle();

	void combineDangle();

	void combineDangleNew();

	void calculateDangle();

	void calculateZipDangle();

	void calculateZeroDangle();

	void calculateInfDangle();

	void loadDangleSuffix();

	void zipDangle();

	void addZeroDangle();

	void minZeroDangle();

	void loadLoop();

	void combineLoop();

	void calculateLoop();

	void loadLoopSuffix();

	void loadSint2();

	void combineSint2();

	void calculateSint2();

	void loadSint2Suffix();

	void symmetryCheckSint2(double sint2[6][6][4][4], char* which);

	void loadAsint1x2();

	void combineAsint1x2();

	void calculateAsint1x2();

	void loadAsint1x2Suffix();

	void loadSint4();

	void combineSint4();

	void calculateSint4();

	void loadSint4Suffix();

	void symmetryCheckSint4(double sint4[6][6][4][4][4][4], char* which);

	void loadTstacki();

	void loadTstackh();

	void loadTstackm();

	void loadTstacke();

	void combineTstack(double tstackEnergies[4][4][4][4], double tstackEnthalpies[5][5][5][5], double tstack[5][5][5][5]);

	void combineTstack2(double tstackEnergies[4][4][4][4], double tstackEnthalpies[5][5][6][6], double tstack[5][5][6][6]);

	void loadTstackiSuffix();

	void loadTstackhSuffix();

	void loadTstackmSuffix();

	void loadTstackeSuffix();

	void loadCoaxialSuffix();

	void loadTstackcoaxSuffix();

	void loadCoaxstackSuffix();

	void loadMulti();

	void loadMulti2();

	void combineMulti();

	void calculateMulti();

	void loadMultiSuffix();

	void loadMulti2Suffix();

	void loadMisc();

	void combineMisc();

	void calculateMisc();

	void loadMiscSuffix();

	void makeAUPenalty();

	void makeAUPenaltyH();

	void loadTriloop();

	void combineTriloop();

	void calculateTriloop();

	void loadTriloopSuffix(struct triloop**, int*);

	void loadTloop();

	void combineTloop();

	void calculateTloop();

	void loadTloopSuffix(struct tloop**, int*);

	void loadHexaloop();

	void combineHexaloop();

	void calculateHexaloop();

	void loadHexaloopSuffix(struct hexaloop**, int*);

	void loadRTSuffix(double* RT);

	unsigned char toNum(char c);
};
#endif
