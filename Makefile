CC=		cc
#CC=		clang
#CC=		gcc10
#CC=		gcc
CFLAGS=		-I. -ansi -Wall -DDEBUG -O0 -g
LDFLAGS=	-lm

COMPONENTS=	rseq ut2n_fn
PROG=		test

all:		${PROG}

${PROG}:	${PROG}.c $(COMPONENTS:=.o) $(COMPONENTS:=.h)
		${CC} ${CFLAGS} $(COMPONENTS:=.o) ${PROG}.c ${LDFLAGS} -o $@

$(COMPONENTS:=.o):	$(COMPONENTS:=.c) $(COMPONENTS:=.h)
		${CC} ${CFLAGS} -c $(COMPONENTS:=.c)

.PHONY: clean
clean:
	rm -f ${PROGS} *.o out/*
