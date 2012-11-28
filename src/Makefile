CC=			gcc
CXX=		g++
#CFLAGS=		`pkg-config --cflags --libs glib-2.0` -ggdb -g -Wall -O2 -pg
CFLAGS=		-ggdb -g -Wall -O2 -pg $(shell pkg-config --libs gthread-2.0 --cflags glib-2.0)
CXXFLAGS=	$(CFLAGS) $(shell pkg-config --cflags glib-2.0)
DFLAGS=		-DHAVE_PTHREAD #-D_FILE_OFFSET_BITS=64
OBJS=		rand.o utils.o roadmap.o bwt.o eva.o edge.o clean.o rnaseq.o \
			pepath.o bwtio.o bwtaln.o pet.o pealn.o edgelist.o pehash.o pelib.o \
			pechar.o pool.o bwtgap.o is.o bntseq.o peseq.o simu.o \
			bwtmisc.o bwtindex.o stdaln.o bwaseqio.o bamlite.o \
			bwase.o kstring.o cs2nt.o ass.o readrm.o scaffolding.o
PROG=		peta
INCLUDES=	
LIBS=		-lm -lz -lpthread -Lbwt_gen -lbwtgen -pg $(shell pkg-config --libs glib-2.0)
SUBDIRS=	. bwt_gen

export G_SLICE=always-malloc

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lib-recur all-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" CXX="$(CXX)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" $$target || exit 1; \
			cd $$wdir; \
		done;

lib:

peta:lib-recur $(OBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)
		
bwt.o:bwt.h
bwtio.o:bwt.h
bwtaln.o:bwt.h bwtaln.h kseq.h
bntseq.o:bntseq.h
bwtgap.o:bwtgap.h bwtaln.h bwt.h

#petasimu.o:petasimu.h

cleanlocal:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

clean:cleanlocal-recur
