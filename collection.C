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

#include "alloc.h"
#include "collection.h"

#define ZERO -1e-7

RegionCollection::RegionCollection(Region<double> **l)
{
	max_num = 0;
	sorted = false;
	list = NULL;
	int count = 0;
	while(l[count]) count++;
	num = count;
	list = (Region<double> **)alloc.xmalloc(sizeof(Region<double> *)*num);
	count = 0;
	while(l[count])
	{ 
		list[count] = l[count];
		count++;
	}
}

RegionCollection::RegionCollection(int m)
{
	sorted = false;
	list = NULL;
	num = 0;
	max_num = m;
}

RegionCollection::RegionCollection()
{
	sorted = false;
	list = NULL;
	num = 0;
	max_num = 0;
}

void RegionCollection::swap(int a, int b)
{
	Region<double> *temp = list[a];
	list[a] = list[b];
	list[b] = temp;
}

void RegionCollection::bubble(int itm, int attr, bool dir)
{
	double cf = (dir) ? 1.0 : -1.0;

	int cur = itm;

	while((cur <= 0) ? false : (cf*list[cur]->getAttrib(attr) < cf*list[cur-1]->getAttrib(attr)))
	{
		swap(cur, cur-1);
		cur--;
	}
}

void RegionCollection::sort(int attr, bool dir)
{
	double cf = (dir) ? -1.0 : 1.0;

	for(int i = 0; i < num; i++)
	{
		int max_idx = i;
		double max_p = cf*list[i]->getAttrib(attr);

		for(int j = i; j < num; j++)
		{
			if(cf*list[j]->getAttrib(attr) > max_p)
			{
				max_p = cf*list[j]->getAttrib(attr);
				max_idx = j;
			}
		}
		if(max_idx != i)
			swap(i, max_idx);
	}
	sorted = true;
	sorted_attrib = attr;
	sorted_dir = dir;
}

RegionCollection::~RegionCollection()
{
	flush();
}

void RegionCollection::add(Region<double> *r)
{
	num++;
	list = (Region<double> **)alloc.xrealloc(list, sizeof(Region<double> *)*num);
	list[num-1] = r;
}

void RegionCollection::replace(int idx, Region<double> *r)
{
	if(idx < 0 || idx >= num)
		return;

	delete list[idx];
	list[idx] = r;
}


bool RegionCollection::sortadd(Region<double> *r)
{
	if(!sorted)
		return false;

	double cf = sorted_dir ? 1.0 : -1.0;

	if(num < max_num)
	{
		num++;
		list = (Region<double> **)alloc.xrealloc(list, sizeof(Region<double> *)*num);
		list[num-1] = r;
		bubble(num-1, sorted_attrib, sorted_dir);
	}
	else
	{
		if(cf*list[num-1]->getAttrib(sorted_attrib) > cf*r->getAttrib(sorted_attrib))
		{
			delete list[num-1];
			list[num-1] = r;
			bubble(num-1, sorted_attrib, sorted_dir);
		} else
			return false;
	}
	return true;
}

void RegionCollection::flush()
{
	if(list == NULL)
	{
		num = 0;
		return;
	}

	for(int i=0; i < num; i++)
		delete list[i];

	num = 0;

	free(list);
	list = NULL;
}

bool RegionCollection::contains(Region<double> *r)
{
	if(num == 0 || list == NULL)
		return false;

	for(int i=0; i < num; i++)
		if(list[i]->equals(r))
			return true;
	
	return false;
}

bool RegionCollection::contains(Region<double> r)
{
	return contains(&r);
}

bool RegionCollection::includes(Region<double> *r)
{
	if(num == 0 || list == NULL)
		return false;

	for(int i=0; i < num; i++)
		if(list[i]->includes(r))
			return true;
	
	return false;
}

bool RegionCollection::includes(Region<double> r)
{
	return includes(&r);
}

RegionCollection *RegionCollection::component(int s)
{
	RegionCollection *ret = new RegionCollection();
	for(int i = 0; i < num; i++)
	{
		Region<double> *r = item(i);
		Region<double> *newr = new Region<double>(r->i(s), r->j(s));
		if(!ret->contains(newr))
			ret->add(newr);
		else
			delete newr;
	}
	return ret;
}

RegionCollection *RegionCollection::component(int s, int attr, double threshold)
{
	RegionCollection *ret = new RegionCollection();
	for(int i = 0; i < num; i++)
	{
		Region<double> *r = item(i);
		Region<double> *newr = new Region<double>(r->i(s), r->j(s));
		double pf = r->getAttrib(attr + s);
		newr->addAttrib(pf);
		bool ad = true;
		int replace = -1;
		for(int k = 0; k < ret->size() && ad; k++)
		{
			if(ret->item(k)->includes(newr) && (ret->item(k)->getAttrib(0) / pf >= threshold))
			{
				ad = false;
				break;
			}

			if(newr->includes(ret->item(k)) && (pf / ret->item(k)->getAttrib(0) >= threshold))
			{
				ad = false;
				replace = k;
				break;
			}	
		}

		if(!ad) 
			if(replace == -1)			
				delete newr;
			else
				ret->replace(replace, newr);
		else
			ret->add(newr);
	}
	return ret;
}


