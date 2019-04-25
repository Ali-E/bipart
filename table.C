/*
Copyright 2008 Hamid Reza Chitsaz (chitsaz@wayne.edu)

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
	Desc: Table class stores dynamic programming values.

	Author: Hamid Reza Chitsaz
		Wayne State University, Algorithmic Biology Lab

	Last Update by Hamid Reza Chitsaz: Sep 14, 2012
	Ver: 2.6
*/

//#include "table.h"

template <typename T> Table<T>::Table(unsigned int l1, unsigned int l2, T iv)
{
	Q = NULL;
	len1 = l1;
	len2 = l2;
	sublen1 = l1;
	sublen2 = l2;
	duplex = true;
	dim = 0;
	reverse = false;
	allocate();
	init(iv);
}

template <typename T> Table<T>::Table(unsigned int l1, unsigned int l2, unsigned int subl1, unsigned int subl2, T iv)
{
	Q = NULL;
	len1 = l1;
	len2 = l2;
	sublen1 = (subl1 <= l1) ? subl1 : l1;
	sublen2 = (subl2 <= l2) ? subl2 : l2;
	duplex = true;
	dim = 0;
	reverse = false;
	allocate();
	init(iv);
}

template <typename T> Table<T>::Table(unsigned int l, char *d, T iv)
{
	Q = NULL;
	len = l;
	sublen = l;
	duplex = false;
	reverse = false;
	dim =*d;
	allocate();
	init(iv);
}

template <typename T> Table<T>::Table(unsigned int l, unsigned int subl, char *d, T iv)
{
	Q = NULL;
	len = l;
	sublen = (subl <= l) ? subl : l;
	duplex = false;
	reverse = false;
	dim =*d;
	allocate();
	init(iv);
}

template <typename T> Table<T>::Table(Table<T> *t)
{
	Q = NULL;
	t->getParams(&duplex, &len1, &len2, &sublen1, &sublen2, &len, &sublen, &reverse, &dim);
	allocate();
	init(t->canvas());
}

template <typename T> Table<T>::Table()
{
	Q = NULL;
}

template <typename T> Table<T>::~Table()
{
	if(Q)
	{
		free(Q);
		total_tables_size -= sizeof(T)*tableSize;
		if(duplex || dim == 4)
			big_tables_num--;
	}

//	printf("table destructed, remaining: %ld\n", total_tables_size);
}

template <typename T> long int Table<T>::allocate()
{
	if(duplex)
	{
		tableSize1 = len1*(len1 + 1)/2 - (len1 - sublen1)*(len1 - sublen1 + 1)/2;
		tableSize2 = len2*(len2 + 1)/2 - (len2 - sublen2)*(len2 - sublen2 + 1)/2;
		tableSize = tableSize1*tableSize2;
		big_tables_num++;
	}
	else
	{
		if(dim == 2)
			tableSize = len*(len + 1)/2 - (len - sublen)*(len - sublen + 1)/2;
		else if(dim == 4)
		{
			tableSize = sublen*(sublen + 1)*(sublen + 2)*(4*len - 3*(sublen - 1))/24;
			big_tables_num++;
		}
	}

	if(Q) free(Q);

	Q = (T *)alloc.xmalloc(sizeof(T)*tableSize);
	if(!Q)
		exit(EXIT_FAILURE);

	total_tables_size += sizeof(T)*tableSize;

	return tableSize;
}

template <typename T> inline long unsigned int Table<T>::addr_dup(unsigned int i1, unsigned int j1, unsigned int i2, unsigned int j2) 
{
	long unsigned int idx1 = (i1 < len1 - sublen1) ? (i1 * sublen1 + (j1 - i1)) : tableSize1 - (len1 - i1)*(len1 - i1 + 1)/2 + (j1 - i1);
	long unsigned int idx2 = (i2 < len2 - sublen2) ? (i2 * sublen2 + (j2 - i2)) : tableSize2 - (len2 - i2)*(len2 - i2 + 1)/2 + (j2 - i2);

#ifdef DEBUG_TABLE
	if(idx1*tableSize2 + idx2 >= tableSize)
	{
		printf("table panic %d\t%d\t%d\t%d\tlen1:%d\tlen2:%d\n", i1, j1, i2, j2, len1, len2);
		exit(1);
	}

	if(j1 < i1 || j2 < i2)
	{
		printf("table panic %d\t%d\t%d\t%d\tlen1:%d\tlen2:%d\n", i1, j1, i2, j2, len1, len2);
		exit(1);
	}
#endif

	return idx1*tableSize2 + idx2;
}


template <typename T> inline long unsigned int Table<T>::singaddr2(unsigned int i, unsigned int j) 
{

	unsigned int x, y;
	if(reverse)
	{
		x = len - 1 - j;
		y = len - 1 - i;
	}
	else
	{
		x = i;
		y = j;
	}

#ifdef DEBUG_TABLE
	if(((x < len - sublen) ? (x * sublen + (y - x)) : tableSize - (len - x)*(len - x + 1)/2 + (y - x)) >= tableSize)
	{
		printf("table panic %d\t%d\tlen:%d\ttablesize:%ld\taddr:%ld\n", i, j, len
		, tableSize, (x < len - sublen) ? (x * sublen + (y - x)) : tableSize - (len - x)*(len - x + 1)/2 + (y - x));
		exit(1);
	}
	if(i >= len || j >= len || j < i)
	{
		printf("table panic %d\t%d\tlen:%d\n", i, j, len);
		exit(1);
	}
#endif

	return (x < len - sublen) ? (x * sublen + (y - x)) : tableSize - (len - x)*(len - x + 1)/2 + (y - x);
}

template <typename T> inline long unsigned int Table<T>::singaddr4(unsigned int i, unsigned int j, unsigned int d, unsigned int e) 
{
	unsigned int x, y, z, t;
	if(reverse)
	{
		x = len - 1 - j;
		y = len - 1 - i;
		z = len - 1 - e;
		t = len - 1 - d;
	} else
	{
		x = i;
		y = j;
		z = d;
		t = e;
	}

#ifdef DEBUG_TABLE
	if(((x < len - sublen) ? (x * sublen * (sublen + 1) * (sublen + 2)/6 + (y - x + 1) * (y - x + 2) * (y - x + 3)/6 - (y - z + 1)*(y - z + 2)/2 + (t - z)) : tableSize - (len - x)*(len - x + 1)*(len - x + 2)*(len - x + 3)/24 + (y - x + 1)*(y - x + 2)*(y - x + 3)/6 - (y - z + 1)*(y - z + 2)/2 + (t - z)) >= tableSize)
	{
		printf("table panic %d\t%d\t%d\t%d\tlen:%d\n", i, j, d, e, len);
		exit(1);
	}
	if(i >= len || j >= len || d >= len || e >= len || j < i || d < i || d > j || e < d || e > j)
	{
		printf("table panic %d\t%d\t%d\t%d\tlen:%d\n", i, j, d, e, len);
		exit(1);
	}
#endif
	
	return (x < len - sublen) ? (x * sublen * (sublen + 1) * (sublen + 2)/6 + (y - x + 1) * (y - x + 2) * (y - x + 3)/6 - (y - z + 1)*(y - z + 2)/2 + (t - z)) : tableSize - (len - x)*(len - x + 1)*(len - x + 2)*(len - x + 3)/24 + (y - x + 1)*(y - x + 2)*(y - x + 3)/6 - (y - z + 1)*(y - z + 2)/2 + (t - z);
}

template <typename T> void Table<T>::init(T *t)
{
	memcpy(Q, t, sizeof(T)*tableSize);
}

template <typename T> void Table<T>::init(T iv)
{
	IV = iv;
	reset();

/*	if (duplex)
	{
		for(unsigned int i1 = 0; i1 < len1; i1++)
			for(unsigned int j1 = i1; j1 < len1 && j1 < i1 + sublen1; j1++)
				for(unsigned int i2 = 0; i2 < len2; i2++)
					for(unsigned int j2 = i2; j2 < len2 && j2 < i2 + sublen2; j2++)
						Q[addr_dup(i1, j1, i2, j2)] = iv;		
	}
	else
	{
		if(dim == 2)
		{
			for(unsigned int i = 0; i < len; i++)
				for(unsigned int j = i; j < len && j < i + sublen; j++)
					Q[singaddr2(i, j)] = iv;
		}
		else if(dim == 4)
		{
			for(unsigned int i = 0; i < len; i++)
				for(unsigned int j = i; j < len && j < i + sublen; j++)
					for(unsigned int d = i; d <= j; d++)
						for(unsigned int e = d; e <= j; e++)
							Q[singaddr4(i, j, d, e)] = iv;		
		}
	}*/
}

template <typename T> void Table<T>::reset()
{
	for(int i = 0; i < tableSize; i++)
		Q[i] = IV;
}

template <typename T> inline T Table<T>::element(int i, int j, int d, int e)
{
	if(duplex)
		return Q[addr_dup(i, j, d, e)];
	else
		return Q[singaddr4(i, j, d, e)];
}

template <typename T> inline T Table<T>::element(int i, int j)
{
	return Q[singaddr2(i, j)];
}

template <typename T> inline T* Table<T>::estar(int i, int j, int d, int e)
{
	if(duplex)
		return &(Q[addr_dup(i, j, d, e)]);
	else
		return &(Q[singaddr4(i, j, d, e)]);
}

template <typename T> inline T* Table<T>::estar(int i, int j)
{
	return &(Q[singaddr2(i, j)]);
}

template <typename T> inline T & Table<T>::operator()(int i, int j, int d, int e)
{
	if(duplex)
		return Q[addr_dup(i, j, d, e)];
	else
		return Q[singaddr4(i, j, d, e)];
}

template <typename T> inline T & Table<T>::operator()(int i, int j)
{
	return Q[singaddr2(i, j)];
}

template <typename T> void Table<T>::Reverse()
{
	reverse = !reverse;	
}

template <typename T> bool Table<T>::hasAny(T iv)
{
	for(int i = 0; i < tableSize; i++)
		if (Q[i] == iv)
			return true;

	return false;
}

template <typename T> bool Table<T>::hasOtherThan(T iv)
{
	for(int i = 0; i < tableSize; i++)
		if (Q[i] != iv)
			return true;

	return false;
}

template <typename T> T Table<T>::max(int *mi, int *mj)
{
	T maxval;
	if(dim != 2)
		return maxval;

	maxval = Q[singaddr2(0, 0)];
	unsigned int maxi = 0, maxj = 0;

	for(unsigned int i = 0; i < len; i++)
		for(unsigned int j = i; j < len && j < i + sublen; j++)
			if(Q[singaddr2(i, j)] > maxval)
			{
				maxval = Q[singaddr2(i, j)];
				maxi = i;
				maxj = j;
			};

	*mi = maxi;
	*mj = maxj;
	return maxval;
}

template <typename T> T Table<T>::max(int *i, int *j, int *d, int *e)
{
// to be fixed later
/*	if (duplex)
	{
		for(unsigned int i1 = 0; i1 < len1; i1++)
			for(unsigned int j1 = i1; j1 < len1 && j1 < i1 + sublen1; j1++)
				for(unsigned int i2 = 0; i2 < len2; i2++)
					for(unsigned int j2 = i2; j2 < len2 && j2 < i2 + sublen2; j2++)
						Q[addr_dup(i1, j1, i2, j2)] = iv;		
	}
	else
	{
		if(dim == 2)
		{
			for(unsigned int i = 0; i < len; i++)
				for(unsigned int j = i; j < len && j < i + sublen; j++)
					Q[singaddr2(i, j)] = iv;
		}
		else if(dim == 4)
		{
			for(unsigned int i = 0; i < len; i++)
				for(unsigned int j = i; j < len && j < i + sublen; j++)
					for(unsigned int d = i; d <= j; d++)
						for(unsigned int e = d; e <= j; e++)
							Q[singaddr4(i, j, d, e)] = iv;		
		}
	}*/

	T maxval;
	return maxval;
}

template <typename T> void Table<T>::getParams(bool *_duplex, unsigned int *_len1, unsigned int *_len2, unsigned int *_sublen1, unsigned int *_sublen2, 
unsigned int *_len, unsigned int *_sublen, bool *_reverse, char *_dim)
{
	*_duplex = duplex;
	*_len1 = len1; 
	*_len2 = len2;
	*_sublen1 = sublen1; 
	*_sublen2 = sublen2, 
	*_len = len;
	*_sublen = sublen;
	*_reverse = reverse;
	*_dim = dim;
}

template <typename T> void Table<T>::store(FILE *outp)
{
	char cduplex = 0, creverse = 0;
	if(duplex) cduplex = 1;
	if(reverse) creverse = 1;

	if(fwrite(&IV, sizeof(T), 1, outp) != 1)
		perror("Storing table");
	if(fwrite(&cduplex, 1, 1, outp) != 1)
		perror("Storing table");
	if(fwrite(&creverse, 1, 1, outp) != 1)
		perror("Storing table");
	if(fwrite(&dim, 1, 1, outp) != 1)
		perror("Storing table");
	if(fwrite(&len1, sizeof(unsigned int), 1, outp) != 1)
		perror("Storing table");
	if(fwrite(&len2, sizeof(unsigned int), 1, outp) != 1)
		perror("Storing table");
	if(fwrite(&sublen1, sizeof(unsigned int), 1, outp) != 1)
		perror("Storing table");
	if(fwrite(&sublen2, sizeof(unsigned int), 1, outp) != 1)
		perror("Storing table");
	if(fwrite(&len, sizeof(unsigned int), 1, outp) != 1)
		perror("Storing table");
	if(fwrite(&sublen, sizeof(unsigned int), 1, outp) != 1)
		perror("Storing table");
	if((long unsigned int)fwrite(Q, sizeof(T), tableSize, outp) != tableSize)
		perror("Storing table");
}

template <typename T> void Table<T>::retrieve(FILE *inp)
{
	char cduplex = 0, creverse = 0;

	if(fread(&IV, sizeof(T), 1, inp) != 1)
		perror("Retrieving table");
	if(fread(&cduplex, 1, 1, inp) != 1)
		perror("Retrieving table");
	if(fread(&creverse, 1, 1, inp) != 1)
		perror("Retrieving table");
	if(fread(&dim, 1, 1, inp) != 1)
		perror("Retrieving table");
	if(fread(&len1, sizeof(unsigned int), 1, inp) != 1)
		perror("Retrieving table");
	if(fread(&len2, sizeof(unsigned int), 1, inp) != 1)
		perror("Retrieving table");
	if(fread(&sublen1, sizeof(unsigned int), 1, inp) != 1)
		perror("Retrieving table");
	if(fread(&sublen2, sizeof(unsigned int), 1, inp) != 1)
		perror("Retrieving table");
	if(fread(&len, sizeof(unsigned int), 1, inp) != 1)
		perror("Retrieving table");
	if(fread(&sublen, sizeof(unsigned int), 1, inp) != 1)
		perror("Retrieving table");
	duplex = (cduplex != 0);
	reverse = (creverse != 0);
	allocate();
	if((long unsigned int)fread(Q, sizeof(T), tableSize, inp) != tableSize)
		perror("Retrieving table");
}


