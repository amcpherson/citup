
#------------------------------------------------
#  makefile for CITUP
#------------------------------------------------

#CC=g++-mp-4.8
CC=gcc

citupqip: citupqip.cpp
	$(CC) -ggdb citupqip.cpp -o citupqip -O3 -DIL_STD -DILOSTRICTPOD -I../External/ -I$(CPLEX_DIRECTORY)/concert/include -I$(CPLEX_DIRECTORY)/cplex/include -L$(CPLEX_DIRECTORY)/concert/lib/$(CPLEX_BUILD)/static_pic/ -lconcert -L$(CPLEX_DIRECTORY)/cplex/lib/$(CPLEX_BUILD)/static_pic -l ilocplex -l cplex -l m -l pthread -I$(BOOST_DIRECTORY)


citupiter: citupiter.cpp
	$(CC) -ggdb citupiter.cpp -o citupiter -O3 -DIL_STD -DILOSTRICTPOD -I../External/ -I$(CPLEX_DIRECTORY)/concert/include -I$(CPLEX_DIRECTORY)/cplex/include -L$(CPLEX_DIRECTORY)/concert/lib/$(CPLEX_BUILD)/static_pic/ -lconcert -L$(CPLEX_DIRECTORY)/cplex/lib/$(CPLEX_BUILD)/static_pic -l ilocplex -l cplex -l m -l pthread -I$(BOOST_DIRECTORY)

	
