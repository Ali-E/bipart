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
	Desc: BPScore class calculates base pairing scores.

	Author: Hamidreza Chitsaz and Ali Ebrahimpour Boroojeny
		Colorado State University
		Algorithmic Biology Lab

	Last Update by Hamidreza Chitsaz: Oct 20, 2018
	Ver: 1.0
*/

#include "bpscore.h"

uint16_t BPScore::intra_score(int a, int b)
{
	//[0-4] for [A,C,G,TU,N]

	if((a == 0 && b == 3) || (a == 3 && b == 0))
		return 2;

	if((a == 1 && b == 2) || (a == 2 && b == 1))
		return 3;

	if((a == 2 && b == 3) || (a == 3 && b == 2))
		return 1;

	return 0;
}

uint16_t BPScore::inter_score(int a, int b)
{
	//[0-4] for [A,C,G,TU,N]

	if((a == 0 && b == 3) || (a == 3 && b == 0))
		return 2;

	if((a == 1 && b == 2) || (a == 2 && b == 1))
		return 3;

	if((a == 2 && b == 3) || (a == 3 && b == 2))
		return 1;

	return 0;
}



