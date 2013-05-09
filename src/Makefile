VERSION=	0.1
CC=			gcc
CXX=		g++
CFLAGS=		-ggdb -g -Wall -O2 $(shell pkg-config --libs gthread-2.0 --cflags glib-2.0)
CXXFLAGS=	-std=c++0x $(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD #-D_FILE_OFFSET_BITS=64
OBJS=		utils.o edge.o rnaseq.o \
			bwtaln.o hash.o \
			pechar.o peseq.o \
			bwaseqio.o kmers.o ass.o \
			kstring.o 
PROG=		peta
INCLUDES=	
LIBS=		-lm -lrt -lz -lpthread $(shell pkg-config --libs glib-2.0)
SUBDIRS=	

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
		$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)
		
clean: # $(RM) is the platform-independent command to remove a file.
		@$(RM) gmon.out *.o a.out $(PROG) *~ *.a core 
#		wdir=`pwd`; \
#		for subdir in $(SUBDIRS); do \
#			cd $$subdir; \
#			$(MAKE) clean || exit 1; \
#			cd $$wdir; \
#		done;

dist:clean
		@cd ../../; \
		tar czf peta-$(VERSION).tgz $(SRCDIR) $(DOCDIR) $(SCRIPTS)