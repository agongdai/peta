VERSION=	0.1
CC=			gcc
CXX=		g++
CFLAGS=		-ggdb -g -Wall -O2 $(shell pkg-config --libs gthread-2.0 --cflags glib-2.0)
CXXFLAGS=	-std=c++0x $(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD #-D_FILE_OFFSET_BITS=64
OBJS=		rand.o utils.o roadmap.o edge.o clean.o rnaseq.o correct.o mate.o \
			pepath.o bwtaln.o pealn.o edgelist.o pehash.o pelib.o \
			pechar.o pool.o peseq.o simu.o hits.o \
			bwaseqio.o list.o kmer.o kmers.o \
			kstring.o readrm.o scaffolding.o merge.o
PROG=		peta
SDARRAY=	bitarray/*.o
INCLUDES=	
LIBS=		-lm -lrt -lz -lpthread $(shell pkg-config --libs glib-2.0)
SUBDIRS=	bitarray

ROOTDIR=	peta_dev
SRCDIR=		$(ROOTDIR)/src
DOCDIR=		$(ROOTDIR)/doc
SCRIPTS=	$(ROOTDIR)/scripts

export G_SLICE=always-malloc

.SUFFIXES:.c .o .cpp

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lib-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			echo $$subdir; cd $$subdir; \
			$(MAKE) CC="$(CC)" CXX="$(CXX)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" || exit 1; \
			cd $$wdir; \
		done;

peta:lib-recur $(OBJS) main.o
		$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS) $(SDARRAY) main.o -o $@ $(LIBS)
		
cleanlocal: # $(RM) is the platform-independent command to remove a file.
		@$(RM) gmon.out *.o a.out $(PROG) *~ *.a core

clean:cleanlocal

dist:cleanlocal
		@cd ../../; \
		tar czf peta-$(VERSION).tgz $(SRCDIR) $(DOCDIR) $(SCRIPTS)