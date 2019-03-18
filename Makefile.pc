CC = g++
CFLAGS = -Wall -O2 -std=c++11 -fopenmp
MKLINC = -I${MKLROOT}/include
CPLUSINC =
LFGSL = -lgsl -lgslcblas -fopenmp
LFMKL = #-L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl -liomp5
LFARP =
LFLAGS = -lpthread

SOURCES = main.cpp \
		  LDefinit.cpp \
		  rIntegral.cpp \
		  kzIntegral.cpp \

OBJS = ${SOURCES:.cpp=.o}


all: $(OBJS)
	$(CC) $(OBJS) -o vortexUR0 $(LFGSL) $(LFMKL) $(LFARP)
	rm -f *.o


$(OBJS): %.o: %.cpp
	$(CC) -c $(CFLAGS) $(MKLINC) $(CPLUSINC) $< -o $@

clean:
	rm -f *.o vortexUR*

