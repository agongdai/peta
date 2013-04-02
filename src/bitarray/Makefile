VERSION=	0.1
CC=			gcc
CXX=		g++
CFLAGS=		-ggdb -g -Wall -O2 
CXXFLAGS=	-std=c++0x $(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD #-D_FILE_OFFSET_BITS=64
OBJS=		bitop.o rank6p.o debug.o bitarray.o sdarray.o
INCLUDES=	

.SUFFIXES:.c .o .cpp

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(OBJS)
		
cleanlocal: # $(RM) is the platform-independent command to remove a file.
		@$(RM) gmon.out *.o a.out $(PROG) *~ *.a core

clean:cleanlocal
