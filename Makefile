#saigon (SunOS), gcc without openmp
#C = /usr/sfw/bin/g++
#CPP = /usr/sfw/bin/g++
#CFLAGS = -m64 -Wall -pedantic -O3 

#saigon (SunOS), cc with openmp
#C = cc
#CPP = CC
#CFLAGS = -fast -xarch=amd64 -xopenmp 

#IRIX, cc with openmp
#C = cc
#CPP = CC
#CFLAGS = -mp -O3 -64 -woffall

#C = gcc
#CPP = g++
#CFLAGS = -maix64 -Wall -O3 

#AIX, xlc_r with openmp
#C = xlc_r
#CPP = xlC_r
#CFLAGS = -q64 -O5 -qsmp=omp 

#Linux, gcc with openmp
C = gcc
CPP = g++
CFLAGS = -m64 -fopenmp -O3

CSRC = 
COBJECTS = $(CSRC:.c=.o)

PCPPSRC = pirna.C alloc.C partitionfunction.C probability.C energy.C getopt.C sequence.C
PCPPOBJECTS = $(PCPPSRC:.C=.o)

ICPPSRC = birna.C partitionfunction.C jointprob.C energy.C getopt.C sequence.C alloc.C collection.C mtree.C
ICPPOBJECTS = $(ICPPSRC:.C=.o)

I2CPPSRC = birna2.C getopt.C sequence.C alloc.C bpscore.C
I2CPPOBJECTS = $(I2CPPSRC:.C=.o)

APPCPPSRC = appirna.C alloc.C ubpartitionfunction.C partitionfunction.C energy.C getopt.C sequence.C
APPCPPOBJECTS = $(APPCPPSRC:.C=.o)

BPMAXCPPSRC = bpmax.C alloc.C getopt.C sequence.C bpscore.C
BPMAXCPPOBJECTS = $(BPMAXCPPSRC:.C=.o)

BPPARTCPPSRC = bppart.C alloc.C getopt.C sequence.C bpscore.C
BPPARTCPPOBJECTS = $(BPPARTCPPSRC:.C=.o)

PART = pirna
INTERACT = birna
INTERACT2 = birna2
APPART = appirna
BPMAX = bpmax
BPPART = bppart

all : $(PART) $(INTERACT) $(INTERACT2) $(APPART) $(BPMAX) $(BPPART)

$(PART): $(PCPPOBJECTS) $(COBJECTS)
	$(CPP) $(CFLAGS) -o $(PART) $(PCPPOBJECTS) $(COBJECTS) -lm

$(INTERACT): $(ICPPOBJECTS) $(COBJECTS)
	$(CPP) $(CFLAGS) -o $(INTERACT) $(ICPPOBJECTS) $(COBJECTS) -lm

$(INTERACT2): $(I2CPPOBJECTS) $(COBJECTS)
	$(CPP) -g $(CFLAGS) -o $(INTERACT2) $(I2CPPOBJECTS) $(COBJECTS) -lm

$(APPART): $(APPCPPOBJECTS) $(COBJECTS)
	$(CPP) $(CFLAGS) -o $(APPART) $(APPCPPOBJECTS) $(COBJECTS) -lm

$(BPMAX): $(BPMAXCPPOBJECTS) $(COBJECTS)
	$(CPP) $(CFLAGS) -o $(BPMAX) $(BPMAXCPPOBJECTS) $(COBJECTS) -lm

$(BPPART): $(BPPARTCPPOBJECTS) $(COBJECTS)
	$(CPP) $(CFLAGS) -o $(BPPART) $(BPPARTCPPOBJECTS) $(COBJECTS) -lm

.c.o: $(CSRC)
	$(C) $(CFLAGS) -c $< -o $@

.C.o: $(PCPPSRC) $(ICPPSRC) $(I2CPPSRC) $(APPCPPSRC) $(BPMAXCPPSRC) $(BPPARTCPPSRC)
	$(CPP) $(CFLAGS) -c $< -o $@

clean :
	rm -f *.o *~ $(PART) $(INTERACT) $(INTERACT2) $(APPART) $(BPMAX) $(BPPART)
 





