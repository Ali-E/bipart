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

#include "mtree.h"

MTree::MTree(RegionCollection *sites, Table<double> *weights, Table<double> *jprobs)
{
	nodes_num = sites->size(); 
	nodes = (Node **)alloc.xmalloc(sizeof(Node *)*nodes_num);

	for(int i=0; i < nodes_num; i++)
	{
		nodes[i] = new Node(sites->item(i)->getAttrib(0), i);
		nodes[i]->set_color(i);
	}
		
	maximumSpanningTree(weights, jprobs);
}

MTree::~MTree() 
{
	if(nodes) 
	{
		for(int i=0; i < nodes_num; i++)
			delete nodes[i];

		free(nodes);
	}
}

void MTree::maximumSpanningTree(Table<double> *weights, Table<double> *jprobs)
{
	int i, j;
	bool done = false;
	int components = nodes_num;
	
	while(!done)
	{
		if(weights->max(&i, &j) == 0.0)
			break;

		if(!causeCycle(i, j))
		{
			paintConnectedComponent(nodes[j], nodes[i]->get_color());
			nodes[i]->add_neighbor(nodes[j], jprobs->element(i, j));
			nodes[j]->add_neighbor(nodes[i], jprobs->element(i, j));
			components--;
		}
		*weights->estar(i, j) = 0.0;
		done = (components == 1);
	}
	for(i=0; i < nodes_num; i++)
		nodes[i]->set_flag(false);

	connectedComponent(nodes[0]);
	for(i=0; i < nodes_num; i++)
		if(!nodes[i]->get_flag())
			perror("Markov Tree Computation Error");

}

bool MTree::causeCycle(int ni, int nj)
{
	return (nodes[ni]->get_color() == nodes[nj]->get_color());
}

void MTree::connectedComponent(Node *r)
{
	if(r->get_flag())
		return;
	else
	{
		r->toggle();
		for(int i = 0; i < r->degree(); i++)
			connectedComponent(r->neighbor(i));
	}
}

void MTree::paintConnectedComponent(Node *r, int c)
{
	if(r->get_color() == c)
		return;
	else
	{
		r->set_color(c);
		for(int i = 0; i < r->degree(); i++)
			paintConnectedComponent(r->neighbor(i), c);
	}
}

double MTree::inference(int *list)
{
	int num = list[0];

	for(int i=0; i < nodes_num; i++)
		nodes[i]->set_flag(false);

	for(int i=0; i < num; i++)
		nodes[list[i+1]]->toggle();
	
	double childm[2];
	double prod = 1.0;

	Node *r = nodes[list[1]];
	for(int j=0; j < r->degree(); j++)
	{
		message(r->neighbor(j), r, childm);
		prod *= childm[0];
	}
	return prod*r->self_prob(0);
}

void MTree::message(Node *r, Node *father, double *m)
{
	double childm[2];
	double prod[2];
	double sum;
	int fatheri;

	prod[0] = prod[1] = 1.0;
	for(int j = 0; j < r->degree(); j++)
		if(r->neighbor(j) != father)
		{
			message(r->neighbor(j), r, childm);
			for(int c2 = 0; c2 < (r->get_flag() ? 1 : 2); c2++)
				prod[c2] *= childm[c2];
		} else
			fatheri = j;

	for(int c1 = 0; c1 < 2; c1++)
	{
		sum = 0.0;
		if(father->self_prob(c1) != 0.0)
			for(int c2=0; c2 < (r->get_flag() ? 1 : 2); c2++)
				sum += prod[c2]*r->joint_prob(c2, c1, fatheri) / father->self_prob(c1);
		m[c1] = sum;
	}	
}

bool MTree::isEdge(int i, int j)
{
	bool res = false;

	for(int k=0; k < nodes[i]->degree(); k++)
		if(nodes[i]->neighbor(k) == nodes[j]) res = true;

	return res;
}

void MTree::dump(FILE *out)
{
	for(int t=0; t < nodes_num; t++)
	{
		fprintf(out, "%d <-> ", nodes[t]->get_id());
		for(int s=0; s < nodes[t]->degree(); s++)
			if(nodes[t]->neighbor(s)->get_id() > nodes[t]->get_id())
				fprintf(out, "%d ", nodes[t]->neighbor(s)->get_id());
		fprintf(out, "\n");
	}
}
