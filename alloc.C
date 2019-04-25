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


#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
# include <stdlib.h>
#include "alloc.h"

void* Alloc::xcalloc(unsigned long m, unsigned long n)
{
  void* ptr;

  if(m == 0 || n == 0)
    return NULL;

  if (!(ptr = calloc(m, n)))
    {
      fputs("Error in calloc()\n", stderr);
      exit(EXIT_FAILURE);
    }

  return ptr;
}

void* Alloc::xmalloc(unsigned long n)
{
  void* ptr;

  if (!(ptr = malloc(n)))
    {
      fputs("Error in malloc()\n", stderr);
      exit(EXIT_FAILURE);
    }

  return ptr;
}

void* Alloc::xrealloc(void* ptr, unsigned long n)
{
  if (!(ptr = realloc(ptr, n)) && n != 0)
    {
      fputs("Error in realloc()\n", stderr);
      exit(EXIT_FAILURE);
    }

  return ptr;
}
