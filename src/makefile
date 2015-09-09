
#------------------------------------------------
#  makefile for CITUP
#------------------------------------------------

# change the following two lines according to your system
BOOST_DIRECTORY=/Users/amcpherson/Build/boost_1_51_0
CPLEX_DIRECTORY=/Users/amcpherson/Applications/IBM/ILOG/CPLEX_Studio125
CPLEX_BUILD=x86-64_darwin

CC=g++-mp-4.8

citupqip: citupqip.cpp
	$(CC) -ggdb citupqip.cpp -o citupqip -O3 -DIL_STD -DILOSTRICTPOD -I../External/ -I$(CPLEX_DIRECTORY)/concert/include -I$(CPLEX_DIRECTORY)/cplex/include -L$(CPLEX_DIRECTORY)/concert/lib/$(CPLEX_BUILD)/static_pic/ -lconcert -L$(CPLEX_DIRECTORY)/cplex/lib/$(CPLEX_BUILD)/static_pic -l ilocplex -l cplex -l m -l pthread

	
