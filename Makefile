# Makefile

GCC=gcc
SRC=deblur.c
OUT=deblur

PAPIDIR=/mnt/jc5/CS259/papi
INCLUDEDIR=-I${PAPIDIR}

CFLAGS=-g -pg
LDFLAGS=-L${PAPIDIR} -lm -lpapi -lutil_papi

all:
	${GCC} ${CFLAGS} ${INCLUDEDIR} ${SRC} -o ${OUT} ${LDFLAGS}

clean:
	rm ${OUT}
