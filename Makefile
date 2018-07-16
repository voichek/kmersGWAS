#Good read for making a makefile: https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile"
CXX=g++
RM=rm -f
CPPFLAGS = -std=c++14 -Wall -O3 -I /ebio/abt6/yvoichek/smallproj/prefix/include/ -I ./include/

SRCEXT = cpp

INCLUDEDIR := include
SRCDIR := src
BUILDIR := build
BINDIR := bin

KMC_API = $(INCLUDEDIR)/KMC/kmc_api/
OBJ =   $(BUILDIR)/kmc_file.o $(BUILDIR)/kmer_api.o $(BUILDIR)/mmer.o
#SOURCES_OBJECTS := $(shell find $(KMC_API) -type f -name *.$(SRCEXT))
3OBJECTS := $(patsubst $(KMC_API)/%,$(BUILDDIR)/%,$(SOURCES_OBJECTS:.$(SRCEXT)=.o))


all: kmer_count_shareness kmer_hist kmer_DB

kmer_DB: $(SRCDIR)/kmer_DB.cpp $(SRCDIR)/kmer_general.h $(SRCDIR)/kmer_DB.h  $(OBJ) 
	$(CXX) $(OBJ) $(SRCDIR)/kmer_DB.cpp -o $(BINDIR)/kmer_DB $(CPPFLAGS) 

kmer_count_shareness: $(SRCDIR)/kmer_count_shareness.cpp $(SRCDIR)/kmer_general.h $(SRCDIR)/kmer_count_shareness.h $(OBJ) 
	$(CXX) $(OBJ) $(SRCDIR)/kmer_count_shareness.cpp -o $(BINDIR)/kmer_count_shareness $(CPPFLAGS) 

kmer_hist: $(SRCDIR)/kmer_hist.cpp $(SRCDIR)/kmer_general.h $(OBJ) 
	$(CXX) $(OBJ) $(SRCDIR)/kmer_hist.cpp -o $(BINDIR)/kmer_hist $(CPPFLAGS) 

$(BUILDIR)/mmer.o: 
	$(CXX) -c  $(KMC_API)/mmer.cpp -o $(BUILDIR)/mmer.o $(CPPFLAGS)  

$(BUILDIR)/kmer_api.o: 
	$(CXX) -c  $(KMC_API)/kmer_api.cpp -o $(BUILDIR)/kmer_api.o $(CPPFLAGS)

$(BUILDIR)/kmc_file.o:
	$(CXX) -c  $(KMC_API)/kmc_file.cpp -o $(BUILDIR)/kmc_file.o $(CPPFLAGS)


clean:
	rm $(BINDIR)/kmer_hist $(OBJ)
