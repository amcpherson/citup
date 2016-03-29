
#----------------------------------------------------------------------------------------
#  Makefile for CITUP
#
#    Set the following environment variables:
# 		CPLEX_DIRECTORY: cplex install
#		BOOST_DIRECTORY: unpacked boost distribution
#       PREFIX: installation prefix, binaries copied to $(PREFIX)/bin
#----------------------------------------------------------------------------------------

CC=g++

citupqip: citupqip.cpp
	$(CC) -ggdb citupqip.cpp -o citupqip -O3 -DIL_STD -DILOSTRICTPOD -I../External/ \
		-lconcert -l ilocplex -l cplex -l m -l pthread -I$(BOOST_DIRECTORY)


citupiter: citupiter.cpp
	$(CC) -ggdb citupiter.cpp -o citupiter -O3 -DIL_STD -DILOSTRICTPOD -I../External/ \
		-lconcert -l ilocplex -l cplex -l m -l pthread -I$(BOOST_DIRECTORY)

.PHONY : all
all : citupqip citupiter

.PHONY : install
install: all
	cp citupqip citupiter $(PREFIX)/bin
