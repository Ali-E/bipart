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
	Desc: MTree class builds a Markov Tree (A Probabilistic Graphical Model) and performs inference and max a posteriori.

	Author: Hamid Reza Chitsaz
		Postdoctoral Fellow, 
		SFU, Computational Biology Lab

	Last Update by Hamid Reza Chitsaz: May 1, 2009
	Ver: 1.2
*/

#ifndef MTREE_H
#define MTREE_H

#include "alloc.h"
#include "table.h"
#include "collection.h"

class Node {
private:
	int id;
	Alloc alloc;
	Node **neighbors;
	int neighbors_num;
	double prob;
	double *jprobs[2][2];
	bool flag;
	int color;
public:
	Node(double p, int _id) {neighbors_num = 0; neighbors = NULL; 
		jprobs[0][0] = jprobs[0][1] = jprobs[1][0] = jprobs[1][1] = NULL; prob = p; id = _id;}
	~Node() 
	{
		if(neighbors) free(neighbors); 
		for(int i=0; i < 2; i++)
			for(int j=0; j < 2; j++)
				if(jprobs[i][j]) free(jprobs[i][j]);
	}
	void add_neighbor(Node *n, double jp) 
	{
		neighbors_num++;
		neighbors = (Node **)alloc.xrealloc(neighbors, sizeof(Node *)*neighbors_num);
		for(int i=0; i < 2; i++)
			for(int j=0; j < 2; j++)
				jprobs[i][j] = (double *)alloc.xrealloc(jprobs[i][j], sizeof(double)*neighbors_num);
		neighbors[neighbors_num-1] = n;
		jprobs[0][0][neighbors_num-1] = jp;
		jprobs[1][0][neighbors_num-1] = n->self_prob(0) - jp;
		if(jprobs[1][0][neighbors_num-1] < 0.0) jprobs[1][0][neighbors_num-1] = 0.0;
		jprobs[0][1][neighbors_num-1] = self_prob(0) - jp;
		if(jprobs[0][1][neighbors_num-1] < 0.0) jprobs[0][1][neighbors_num-1] = 0.0;
		jprobs[1][1][neighbors_num-1] = 1 - self_prob(0) - n->self_prob(0) + jp;
		if(jprobs[1][1][neighbors_num-1] < 0.0) jprobs[1][1][neighbors_num-1] = 0.0;
	}
	bool get_flag() {return flag;};
	int get_color() {return color;};
	void toggle() {flag = !flag;};
	void set_flag(bool f) {flag = f;};
	void set_color(int c) {color = c;};
	Node *neighbor(int i) {return neighbors[i];};
	double joint_prob(int cself, int cn, int i) {return jprobs[cself][cn][i];};
	double self_prob(int c) {return (c == 0) ? prob : 1 - prob;};
	int get_id() {return id;};
	int degree() {return neighbors_num;};
};


class MTree {
private:
	Alloc alloc;
	Node **nodes;
	int nodes_num;
	bool causeCycle(int i, int j);
	void maximumSpanningTree(Table<double> *weights, Table<double> *jprobs);
	void connectedComponent(Node *r);
	void paintConnectedComponent(Node *r, int);
	void message(Node *r, Node *father, double *m);
public:
	MTree(RegionCollection *, Table<double> *, Table<double> *);
	~MTree(); 
	double inference(int *list);
	bool isEdge(int i, int j);
	void dump(FILE *out);
};

#endif
