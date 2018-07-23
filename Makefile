#Good read for making a makefile: https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile"
CXX=g++
RM=rm -f
CPPFLAGS = -std=c++14 -Wall -O3 -I /ebio/abt6/yvoichek/smallproj/prefix/include/ -I ./include/
LDFLAGS :=  -L/ebio/abt6/yvoichek/smallproj/prefix/lib -lstdc++ -lboost_program_options

SRCEXT = cpp

INCLUDEDIR := include
SRCDIR := src
BUILDIR := build
BINDIR := bin

KMC_API = $(INCLUDEDIR)/KMC/kmc_api/
OBJ_KMC =   $(BUILDIR)/kmc_file.o $(BUILDIR)/kmer_api.o $(BUILDIR)/mmer.o 
OBJ_YV = $(BUILDIR)/kmer_DB.o $(BUILDIR)/kmer_multipleDB.o
#SOURCES_OBJECTS := $(shell find $(KMC_API) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(KMC_API)/%,$(BUILDDIR)/%,$(SOURCES_OBJECTS:.$(SRCEXT)=.o))


all: kmer_count_shareness YV_kmerhist YV_correlate_kmers_to_phenotype


YV_correlate_kmers_to_phenotype: $(SRCDIR)/YV_correlate_kmers_to_phenotype.cpp $(SRCDIR)/kmer_DB.h $(SRCDIR)/kmer_DB.cpp $(OBJ_KMC) $(OBJ_YV)
	$(CXX) $(OBJ_KMC) $(OBJ_YV) $(SRCDIR)/YV_correlate_kmers_to_phenotype.cpp  -o $(BINDIR)/YV_correlate_kmers_to_phenotype $(CPPFLAGS) $(LDFLAGS)
	
YV_kmer_intersect_and_sort: $(SRCDIR)/YV_kmer_intersect_and_sort.cpp $(SRCDIR)/kmer_DB.h $(SRCDIR)/kmer_DB.cpp $(OBJ_KMC) 
	$(CXX) $(OBJ_KMC) $(SRCDIR)/YV_kmer_intersect_and_sort.cpp $(SRCDIR)/kmer_DB.cpp -o $(BINDIR)/YV_kmer_intersect_and_sort $(CPPFLAGS) 

YV_kmerhist: $(SRCDIR)/YV_kmerhist.cpp $(SRCDIR)/kmer_DB.h $(SRCDIR)/kmer_DB.cpp $(OBJ_KMC) 
	$(CXX) $(OBJ_KMC) $(SRCDIR)/YV_kmerhist.cpp $(SRCDIR)/kmer_DB.cpp -o $(BINDIR)/YV_kmerhist $(CPPFLAGS) 

kmer_count_shareness: $(SRCDIR)/kmer_count_shareness.cpp $(SRCDIR)/kmer_general.h  $(OBJ_KMC) 
	$(CXX) $(OBJ_KMC) $(SRCDIR)/kmer_count_shareness.cpp -o $(BINDIR)/kmer_count_shareness $(CPPFLAGS) 

$(BUILDIR)/kmer_multipleDB.o: $(SRCDIR)/kmer_multipleDB.cpp $(SRCDIR)/kmer_general.h $(SRCDIR)/kmer_multipleDB.h  $(OBJ_KMC) $(OBJ_YV) 
	$(CXX) -c $(SRCDIR)/kmer_multipleDB.cpp -o $(BUILDIR)/kmer_multipleDB.o $(CPPFLAGS) 

$(BUILDIR)/mmer.o: 
	$(CXX) -c  $(KMC_API)/mmer.cpp -o $(BUILDIR)/mmer.o $(CPPFLAGS)  

$(BUILDIR)/kmer_api.o: 
	$(CXX) -c  $(KMC_API)/kmer_api.cpp -o $(BUILDIR)/kmer_api.o $(CPPFLAGS)

$(BUILDIR)/kmc_file.o:
	$(CXX) -c  $(KMC_API)/kmc_file.cpp -o $(BUILDIR)/kmc_file.o $(CPPFLAGS)

$(BUILDIR)/kmer_DB.o: $(SRCDIR)/kmer_DB.cpp $(SRCDIR)/kmer_general.h $(SRCDIR)/kmer_DB.h  $(OBJ_KMC) 
	$(CXX) -c $(SRCDIR)/kmer_DB.cpp -o $(BUILDIR)/kmer_DB.o $(CPPFLAGS) 


clean:
	rm $(BINDIR)/* $(OBJ_KMC)
