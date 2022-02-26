CC=		clang13
#CC=		gcc10
#CC=		gcc
CFLAGS=		-g -I.
LDFLAGS=

PROG=		test

COMPONENTS=	rseq

all:		${PROG}

${PROG}:	${PROG}.c $(COMPONENTS:=.o) $(COMPONENTS:=.h)
		${CC} -Wall -O0 ${CFLAGS} ${LDFLAGS} $(COMPONENTS:=.o) ${PROG}.c -o $@

$(COMPONENTS:=.o):	$(COMPONENTS:=.c) $(COMPONENTS:=.h)
		${CC} -Wall -O0 ${CFLAGS} ${LDFLAGS} -c $(COMPONENTS:=.c)

.PHONY: clean
clean:
	rm -f ${PROG} *.o *.out
