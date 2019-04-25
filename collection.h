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
	Desc: Region and RegionCollection maintain a region and a collection of regions in a sequence or pair of sequences.

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow,
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: April 17, 2009
	Ver: 1.2
*/

#ifndef COLLECTION_H
#define COLLECTION_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "alloc.h"

typedef enum {D2 = 0, D4 = 1} RegionType;

template<class T>
class Region {
private:
	Alloc alloc;
	int I1, J1; 
	int I2, J2;
	RegionType regionType;
	T *attributes;
	int attrib_num;
	void createFrom(Region *r)
	{
		regionType = r->getRegionType(); 
		if(regionType == D2)
		{
			I1 = r->i();
			J1 = r->j();
		}
		if(regionType == D4)
		{
			I1 = r->i1();
			J1 = r->j1();
			I2 = r->i2();
			J2 = r->j2();
		}
		attrib_num = r->getAttribNum();
		attributes = (T *)alloc.xrealloc(attributes, attrib_num*sizeof(T));
		for(int i = 0; i < attrib_num; i++)
			attributes[i] = r->getAttrib(i);
	}
public:
	Region(int i, int j) {I1 = i; J1 = j; regionType = D2; attrib_num = 0; attributes = NULL;};
	Region(int _i1, int _j1, int _i2, int _j2) {I1 = _i1; J1 = _j1; I2 = _i2; J2 = _j2; regionType = D4; attrib_num = 0; attributes = NULL;};
	Region(Region *r) {attributes = NULL; createFrom(r);};
	~Region() {if(attributes) free(attributes);};
	int i() {return I1;}
	int j() {return J1;}
	int i1() {return I1;}
	int j1() {return J1;}
	int i2() {return I2;}
	int j2() {return J2;}
	int i(int idx) {return (idx == 0) ? I1 : I2;}
	int j(int idx) {return (idx == 0) ? J1 : J2;}
	bool equals(Region *r) 
	{
		if(regionType != r->getRegionType())
			return false;
		if(regionType == D2)
			return (I1 == r->i() && J1 == r->j());
		else
			return (I1 == r->i1() && J1 == r->j1() && I2 == r->i2() && J2 == r->j2());
	}
	bool includes(Region *r) 
	{
		if(regionType != r->getRegionType())
			return false;
		if(regionType == D2)
			return (I1 <= r->i() && J1 >= r->j());
		else
			return (I1 <= r->i1() && J1 >= r->j1() && I2 <= r->i2() && J2 >= r->j2());
	}
	bool overlaps(Region *r)
	{
		if(regionType != r->getRegionType())
			return false;

		if(regionType == D2)
			return !((I1 < r->i() && J1 < r->i()) || (I1 > r->j() && J1 > r->j()));
		else
			return false;
	}
	void operator=(Region r) {createFrom(&r);}
	RegionType getRegionType() {return regionType;};
	int getAttribNum() {return attrib_num;}
	T getAttrib(int i) {T temp; if(i >= attrib_num || i < 0) return temp; else return attributes[i];}
	void addAttrib(T a)
	{
		attrib_num++;
		attributes = (T *)alloc.xrealloc(attributes, sizeof(T)*attrib_num);
		attributes[attrib_num-1] = a;
	}
};

class RegionCollection {
private:
	Alloc alloc;
	Region<double> **list;
	int num;
	int max_num;
	int sorted_attrib;
	bool sorted_dir;
	bool sorted;
	void swap(int, int);
	void bubble(int, int, bool);
public:
	RegionCollection();
	RegionCollection(int);
	RegionCollection(Region<double> **l);
	~RegionCollection();
	void add(Region<double> *);
	void replace(int, Region<double> *);
	bool sortadd(Region<double> *);
	void sort(int, bool);
	Region<double> *item(int i) {if(i >= num || i < 0 || list == NULL) return NULL; else return list[i];};
	void flush();
	int size() { return num;};
	bool contains(Region<double> *);
	bool includes(Region<double> *);
	bool contains(Region<double>);
	bool includes(Region<double>);
	RegionCollection *component(int s);
	RegionCollection *component(int s, int, double);
};

#endif
