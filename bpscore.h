/*
Copyright 2018 Hamid Reza Chitsaz (chitsaz@chitsazlab.org)

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

	Authors: Hamidreza Chitsaz and Ali Ebrahimpour Boroojeny
		Colorado State University
		Algorithmic Biology Lab

*/

#ifndef BPSCORE_H
#define BPSCORE_H

#include <stdint.h>

#include "config.h"


class BPScore {
public:
	BPScore() {};
	~BPScore() {};

	double intra_score(int a, int b, double var2, double var3, double no_interact_value);
	double inter_score(int a, int b, double var2, double var3, double no_interact_value);
};

#endif
