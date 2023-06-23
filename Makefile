#make file for parameter choosing programs fLPSparameters and SEGparameters

all: 
	gcc -O2 -o fLPSparameters fLPSparameters.c -lm 
	gcc -O2 -o SEGparameters SEGparameters.c -lm 

	
