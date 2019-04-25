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

#include <string.h>
#include <ctype.h>
#include "sequence.h"


Sequence::Sequence(char *fn)
{
	file_name = (char *)alloc.xmalloc(strlen(fn) + 1);
	strcpy(file_name, fn);
	len = 0;
	name = NULL;
	string = NULL;
	seq = NULL;
	reverse_seq = NULL;
	FILE *f = fopen(file_name, "rt");
	if(!f)
		perror("Sequence: error opening sequence file: ");

	loaded = (load(f) != -1);
	if(loaded) 
		fclose(f);
	computeReverse();

	char *temp = file_name;
	char *res;

	for (res = file_name; *temp; ++temp)
		if (*temp == '/' || *temp == '\\')
			res = temp + 1;

	strncpy(short_filename, res, 99);

	if (strlen(short_filename) > 4 && !strcmp(short_filename + strlen(short_filename) - 4, ".seq"))
		short_filename[strlen(short_filename) - 4] = 0;

	if (!name)
	{
		name = (char *)alloc.xmalloc(strlen(short_filename) + 1);
		strcpy(name, short_filename);
	}
}

Sequence::Sequence(FILE *f, char *fn)
{
	file_name = (char *)alloc.xmalloc(strlen(fn) + 1);
	strcpy(file_name, fn);
	len = 0;
	name = NULL;
	string = NULL;
	seq = NULL;
	reverse_seq = NULL;
	loaded = (load(f) != -1);
	computeReverse();

	char *temp = file_name;
	char *res;

	for (res = file_name; *temp; ++temp)
		if (*temp == '/' || *temp == '\\')
			res = temp + 1;

	strncpy(short_filename, res, 99);

	if (strlen(short_filename) > 4 && !strcmp(short_filename + strlen(short_filename) - 4, ".seq"))
		short_filename[strlen(short_filename) - 4] = 0;

	if (!name)
	{
		name = (char *)alloc.xmalloc(strlen(short_filename) + 1);
		strcpy(name, short_filename);
	}
}


Sequence::Sequence(char *fn, char *n, char *st, unsigned char *sq, unsigned int l)
{
	file_name = fn;
	len = l;
	name = n;
	string = st;
	seq = sq;
	reverse_seq = NULL;
	loaded = false;
	computeReverse();

	char *temp = file_name;
	char *res;

	for (res = file_name; *temp; ++temp)
		if (*temp == '/' || *temp == '\\')
			res = temp + 1;

	strncpy(short_filename, res, 99);

	if (strlen(short_filename) > 4 && !strcmp(short_filename + strlen(short_filename) - 4, ".seq"))
		short_filename[strlen(short_filename) - 4] = 0;

	if (!name)
	{
		name = (char *)alloc.xmalloc(strlen(short_filename) + 1);
		strcpy(name, short_filename);
	}
}

Sequence::~Sequence()
{
	if(loaded)
	{
		free(file_name);
		free(name);
	}
	free(string);
	free(seq);
	free(reverse_seq);
}

int Sequence::load(FILE *f)
{
	if(readSequence(f) < 0)
		return -1;

	return 1;
}

char* Sequence::getName()
{
	return name;
}

char* Sequence::getFileName()
{
	return short_filename;
}

char* Sequence::getString()
{
	return string;
}

unsigned char* Sequence::getSeq()
{
	return seq;
}

unsigned char* Sequence::getReverse()
{
	return reverse_seq;
}

unsigned int Sequence::getLen()
{
	return len;
}

int Sequence::readSequence(FILE *f)
{
	int c = input(f);

	if(c < 0)
		return -1;

	len = strlen(string);

	seq = (unsigned char *)alloc.xmalloc(len);
	for (unsigned int i = 0; i < len; ++i)
		seq[i] = toNum(string[i]);

	return 1;
}

unsigned char Sequence::toNum(char c)
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

int Sequence::input(FILE* file)
{
  int current, last, state;
  unsigned int availableSeq, usedSeq, availableName, usedName;

  string = (char *)alloc.xmalloc(1024);
  availableSeq = 1024;
  usedSeq = 0;
  name = (char *)alloc.xmalloc(80);
  availableName = 80;
  usedName = 0;
  state = 0;
  last = '\n';

  while ((current = fgetc(file)) != EOF)
    {
      if (last == '\n' && current == '>')
	{
	  if (usedSeq)
	    {
	      ungetc('>', file);
	      break;
	    }
	  else
	    state = 1;
	}
      else if (current == '\n')
	state = 0;
      else if (state == 0)
	{
	  if (('a' <= current && current <= 'z') || ('A' <= current && current <= 'Z') || ('0' <= current && current <= '9') || current == '&')
	    string[usedSeq++] = toupper(current);
	  else if (current == ';')
	    break;
	}
      else if (state == 1)
	name[usedName++] = current;

      checkArray(&name, &availableName, usedName, 80);
      checkArray(&string, &availableSeq, usedSeq, 1024);
      last = current;
    }
  name[usedName] = 0;
  string[usedSeq] = 0;
  
  if (strlen(name) == 0)
    {
      free(name);
      name = NULL;
    }
  else
    name = (char *)alloc.xrealloc(name, strlen(name) + 1);

  if (strlen(string) == 0)
    {
      free(string);
      string = NULL;
      return -1;
    }
  else
    string = (char *)alloc.xrealloc(string, strlen(string) + 1);
  return 1;
}

void Sequence::checkArray(char** array, unsigned int* available, unsigned int used, unsigned int increment)
{
  if (used == *available)
    {
      *array = (char *) alloc.xrealloc(*array, *available + increment);
      *available += increment;
    }
}

void Sequence::computeReverse()
{
	reverse_seq = (unsigned char *) alloc.xrealloc(reverse_seq, len);
	for(unsigned int i = 0; i < len; i++)
		reverse_seq[i] = seq[len - i - 1];
}

Sequence *Sequence::split()
{
	int i = len - 1;
	char *new_str;
	unsigned char *new_seq;

	while(i >= 0 && string[i] != '&') i--;

	if(string[i] != '&' || i == len - 1)
		return NULL;

	new_str = (char *)alloc.xmalloc(len - i);
	new_seq = (unsigned char *)alloc.xmalloc(len - i - 1);
	for(unsigned int j = i+1; j <= len; j++)
	{
		new_str[j-i-1] = string[j];
		if(j < len) new_seq[j-i-1] = seq[j];
	}

	Sequence *ret = new Sequence(file_name, name, new_str, new_seq, len-i-1);

	len = i;
	string[len] = 0;	
	string = (char *) alloc.xrealloc(string, len+1);
	seq = (unsigned char *) alloc.xrealloc(seq, len);
	computeReverse();
	return ret;	
}







