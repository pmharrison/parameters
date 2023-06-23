/**** 
 **** fLPSparameters.c 
 **** 
 ****/ 
/****  Copyright 2023. Paul Martin Harrison. ****/ 
/****
 ****  Licensed under the 3-clause BSD license. See LICENSE.txt bundled with this program. 
 ****/ 
/****  
 ****  This program calculates recommended parameters for a given target length of low-complexity 
 ****  or compositionally-biased region in proteins.  
 **** 
 ****  to compile: 
 ****   gcc -O2 -o fLPSparameters fLPSparameters.c -lm 
 ****    OR 
 ****   make 
 ****
 ****  to run and get help: 
 ****   ./fLPSparameters -h 
 ****
 **** 
 ****  This code is bundled with a similar one that chooses parameters for the program SEG. 
 ****
 ****  The latest version of this code is available at:
 ****    http://biology.mcgill.ca/faculty/harrison/flps.html 
 ****      OR 
 ****    https://github.com/pmharrison/flps 
 **** 
 ****  Citation: 
 ****    Harrison, PM. "Optimal strategies for discovery of low-complexity or compositionally-biased regions
 ****     in proteins", submitted. 
 **** 
 ****   Also: 
 ****    Harrison, PM. 'fLPS: fast discovery of compositional biases for the protein universe', 
 ****    (2017) BMC Bioinformatics, 18: 476.  
 ****    Harrison, PM. 'fLPS 2.0: rapid annotation of compositional biases in 
 ****    biological sequences', (2021) PeerJ, 9: e12363. 
 **** 
 ****/ 
/*****************************************************************************************/

#include <stdio.h> 
#include <string.h> 
#include <stdlib.h> 
#include <math.h> 
#include <ctype.h> 
#include <unistd.h> 

enum calculation_type {DIVERSE, NARROW} focus;  
char focus_name[3][10] = {"DIVERSE", "NARROW"}; 

int small_m, big_m, max_m; 
int target_length=-1; 
int not_valid; 
double threshold; 

void print_help()
{
fprintf(stderr, "\nParameter choosing program for finding low-complexity or compositionally-biased regions using fLPS in proteins of a given target length\n");    
fprintf(stderr,   "=======================================================================================================================================\n\n" 
"The program options are:\n"
" -h   prints help\n"
" -f   focus of the parameters\n"
"      values: 'diverse' or 'narrow'; \n"
"      diverse = more diversity or variance of length is allowed (DEFAULT)\n"
"      narrow  = narrowest focus on a particular target length\n"
" -l   target length.\n" 
"      This must be in the range 5-300 inclusive.\n\n"
" The program outputs lists of suitable parameters for a given target length for low-complexity or compositionally-biased regions.\n"
" There are sets of parameters output for estimated protein coverage of approximately 2%%, 5%%, 10%%, 25%%, and 40%%.\n"
" The protein coverage is simply the proportion of proteins that are expected to be annotated or 'covered' when you choose\n"
" a certain set of parameters.\n"
" For some combinations of coverage level and target lengths, sets of parameters cannot be output because they are out of bounds.\n"
" This is an example of running the program:\n"
"        ./fLPSparameters -f diverse -l 15 > parameters.out\n\n"
" Here, diverse focus is specified with a target region length of 15 residues.\n\n"
"CITATION:\n"
" Harrison, PM. 'Optimal strategies for discovery of low-complexity or compositionally-biased regions',\n"
" submitted. \n" 
"URLs:\n http://biology.mcgill.ca/faculty/harrison/flps.html\n OR \nhttps://github.com/pmharrison/flps\n"); 
} /* end of print_help() */ 


void output_parameters(int upper_bound, int coverage)
{
if(target_length<5 || target_length>upper_bound) { not_valid=1; }
if(threshold>-3.0) { not_valid=1; } 
if(small_m<5) { not_valid=1; }

if(target_length<50 && focus==NARROW && coverage==25) { not_valid=1; } 
if(target_length<100 && focus==NARROW && coverage==40) { not_valid=1; } 
if(target_length<=15 && focus==DIVERSE && coverage==40) { not_valid=1; } 
if(target_length<=10 && focus==NARROW) { not_valid=1; } 

if(!not_valid) 
  { fprintf(stdout, "\t~%d%%\t\t\t%d\t%d\t%.1le\n", coverage, small_m, big_m, (double) pow(10.0, threshold) ); } 
else { /*not valid*/ fprintf(stdout, "\t~%d%%\t\t\tNA [ target length <5 OR >%d, OR t>0.001]\n", coverage, upper_bound); } 
} /* end of output_parameters() */ 


int main(int argc, char **argv)
{
int i, c, errflg=0; 
int cov=-1;

extern char *optarg;
extern int optind, optopt; 

/*  *  *  *  PROCESS COMMAND-LINE OPTIONS  *  *  *  */ 
while((c = getopt(argc, argv, "hf:l:")) != -1) { 
     switch(c) {
     case 'h': print_help(); exit(0);   
     case 'f': if(!strcmp(optarg, "narrow")) { focus=NARROW; } 
               else { focus=DIVERSE; } 
               break; 
     case 'l': sscanf(optarg,"%d", &target_length); 
               if(target_length<5 || target_length>300) 
                 { fprintf(stderr, " -l value is out of bounds, re-setting to a DEFAULT VALUE = 15\n"); 
                   target_length=15; } 
               break; 
     case ':': fprintf(stderr, "option -%c requires a value...\n", optopt); errflg++; break; 
     case '?': fprintf(stderr, "unrecognized option: -%c ...\n", optopt); errflg++; 
} /* end of switch */ 
} /* end of while() getopt */ 
if (errflg) { print_help(); exit(1); } 


/*  *  *  * HEADER OF OUTPUT *  *  *  */ 

fprintf(stdout, "\n%s has chosen the following parameters for target length %d and focus %s:\n\n", argv[0]+2, target_length, focus_name[focus]); 
if(focus==DIVERSE)
  { fprintf(stdout, "A DIVERSE focus means that a typical or average level of length variance for the annotated regions is allowed.\n"); } 
else { fprintf(stdout, "A NARROW focus means that length variance is minimized for the annotated regions.\n"); } 
fprintf(stdout, "\tEstimated_coverage\tm\tM\tt:\n"); 
fprintf(stdout, "\t------------------\t-\t-\t--\n"); 


/*  *  *  * CALCULATE THE PARAMETERS *  *  *  */ 
if(focus==DIVERSE)
  { 
  /* 2% */ 
  big_m = round(2.534 * pow((double) target_length, 0.506)); 
  small_m = big_m - 2; 
  threshold = -0.153 * (double) target_length - 3.994; 
  max_m=100; cov=2; 
  output_parameters(max_m, cov); 

  /* 5% */ 
  not_valid=0; 
  big_m = round(3.46 * pow((double) target_length, 0.508)); 
  small_m = big_m - 4; 
  threshold = -0.098 * (double) target_length - 3.305; 
  max_m=200; cov=5; 
  output_parameters(max_m, cov); 

  /* 10% */ 
  not_valid=0; 
  big_m = round(3.912 * pow((double) target_length, 0.543)); 
  small_m = big_m - 10; 
  threshold = -0.055 * (double) target_length - 3.635; 
  max_m=250; cov=10; 
  output_parameters(max_m, cov); 

  /* 25% */ 
  not_valid=0; 

  if(target_length<=105)
    {
    big_m = round(5.647 * pow((double) target_length, 0.56)); 
    small_m = round(0.872 * pow((double) target_length, 0.797)); 
    threshold = -0.039 * (double) target_length - 2.381;
    }
  else { 
       big_m = round(6.096 * pow((double) target_length, 0.552)); 
       small_m = big_m - 50; 
       threshold = -0.031 * (double) target_length - 2.93;
       }

   max_m=300; cov=25; 
   output_parameters(max_m, cov);   

  /* 40% */ 
  not_valid=0; 

  if(target_length<=105)
    {
    big_m = round(9.82 * pow((double) target_length, 0.522)); 
    small_m = round(0.481 * pow((double) target_length, 0.876)); 
    threshold = -0.022 * (double) target_length - 2.709;
    }
  else { 
       big_m = round(11.126 * pow((double) target_length, 0.484)); 
       small_m = big_m - 80; 
       threshold = -0.025 * (double) target_length - 2.762;
       }

   max_m=300; cov=40; 
   output_parameters(max_m, cov); 

  } /* end of focus==DIVERSE */ 

else { /* focus == NARROW */ 
     /* 2% */ 
     big_m = round(2.324 * pow((double) target_length, 0.539)); 
     small_m = big_m; 
     threshold = -0.149 * (double) target_length - 3.883; 
     max_m=100; cov=2; 
     output_parameters(max_m, cov); 

     /* 5% */ 
     not_valid=0; 
     big_m = round(2.976 * pow((double) target_length, 0.556)); 
     small_m = big_m; 
     if(target_length<=28)
       { threshold = -0.127 * (double) target_length - 2.183; }
     else if(target_length>=33) { threshold = -0.09 * (double) target_length - 3.173; }
     else { threshold = ((-0.127 * (double) target_length - 2.183) + (-0.09 * (double) target_length - 3.173))/2.0; } 
     max_m=200; cov=5; 
     output_parameters(max_m, cov); 

     /* 10% */ 
     not_valid=0; 
     big_m = round(3.493 * pow((double) target_length, 0.572)); 
     small_m = big_m; 
     threshold = -0.058 * (double) target_length - 2.731; 
     max_m=200; cov=10; 
     output_parameters(max_m, cov); 

     /* 25% */ 
     not_valid=0; 
     big_m = round(3.394 * pow((double) target_length, 0.672)); 
     small_m = big_m; 
     if(target_length<=90) { threshold = -4.0; } 
     else { threshold = -0.028 * (double) target_length - 1.695; }
     max_m=300; cov=25; 
     output_parameters(max_m, cov); 

     /* 40% */ 
     not_valid=0; 
     big_m = round(0.889 * pow((double) target_length, 0.977)); 
     small_m = big_m; 
     threshold = -4.0; 
     max_m=300; cov=40; 
     output_parameters(max_m, cov); 
     } /* end of focus==NARROW */ 


/*  *  *  * FOOTER OF OUTPUT *  *  *  */ 
fprintf(stdout, "\n\nCoverage is the proportion of protein sequences expected to be labelled by these parameter sets.\n"); 
fprintf(stdout, "\nIt is recommended to use all of the parameters progressively in separate runs of the fLPS program,\n"); 
fprintf(stdout, " and compare the outputs.\n"); 
fprintf(stdout, "If the calculated parameters are listed as 'NA', it means that at least one of them was out of bounds.\n\n"); 

exit(0); 
} /* end of main() */ 

/******** END OF CODE FILE ********/ 

