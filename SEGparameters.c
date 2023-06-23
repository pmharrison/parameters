/**** 
 **** SEGparameters.c 
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
 ****   gcc -O2 -o SEGparameters SEGparameters.c -lm 
 ****    OR 
 ****   make 
 ****
 ****  to run and get help: 
 ****   ./SEGparameters -h 
 ****
 ****  This code is bundled with a similar one that chooses parameters for the program fLPS. 
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

double K1, K2, max_l; 
int L; 
int target_length=-1; 
int not_valid; 

void print_help()
{
fprintf(stderr, "\nParameter choosing program for finding low-complexity or compositionally-biased regions using SEG in proteins of a given target length\n");    
fprintf(stderr,   "======================================================================================================================================\n\n" 
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
"        ./SEGparameters -f diverse -l 15 > parameters.out\n\n"
" Here, diverse focus is specified with a target region length of 15 residues.\n\n"
"CITATIONS:\n"
" Harrison, PM. 'Optimal strategies for discovery of low-complexity or compositionally-biased regions',\n"
" submitted. \n" 
"URLs:\n http://biology.mcgill.ca/faculty/harrison/flps.html\n OR \nhttps://github.com/pmharrison/flps\n"); 
} /* end of print_help() */ 


void output_parameters(int upper_bound, int coverage)
{
int a=0; 
if(target_length<5 || target_length>upper_bound) { not_valid=1; }
if(K2>4.2) { not_valid=1; } 

if(coverage==40 && focus==DIVERSE && target_length<10) 
  { a=1; not_valid=1; fprintf(stdout, "\t~%d%%\t\t\tNA [ target length <10 OR >%d, OR K2>4.2]\n", coverage, upper_bound);}

if(!not_valid) 
  { fprintf(stdout, "\t~%d%%\t\t\t%d\t%.2lf\t%.2lf\n", coverage, L, K1, K2); } 
else if(a==0){ /*not valid*/ fprintf(stdout, "\t~%d%%\t\t\tNA [ target length <5 OR >%d, OR K2>4.2]\n", coverage, upper_bound); } 
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
fprintf(stdout, "\tEstimated_coverage\tL\tK1\tK2:\n"); 
fprintf(stdout, "\t------------------\t-\t--\t---\n"); 


/*  *  *  * CALCULATE THE PARAMETERS *  *  *  */ 
if(focus==DIVERSE)
  { 
  /* 2% */ 
  if(target_length<=35)
    { 
    L = round(1.274 * pow((double) target_length, 0.823)); 
    K2 = 0.701*log(target_length) + 0.155; 
    }
  else if(target_length>45)
         { 
         L = round(1.004 * pow((double) target_length, 0.891)); 
         K2 = 0.447 * log(target_length) + 1.038; 
         }
  else { 
       L = round(((1.274 * pow((double) target_length, 0.823)) + (1.004 * pow((double) target_length, 0.891)))/2.0); 
       K2 = ((0.701*log(target_length) + 0.155) + (0.447 * log(target_length) + 1.038))/2.0; 
       }
  K1 = K2 - 0.3; 
  max_l=200; cov=2; 
  output_parameters(max_l, cov); 

  /* 5% */ 
  not_valid=0; 
  if(target_length<=50)
    {
    L = round(1.385 * pow((double) target_length, 0.801)); 
    K2 = 0.716 * log(target_length) + 0.381; 
    K1 = K2 - 0.3; 
    }
  else { 
       L = round(0.747 * pow((double) target_length, 0.912)); 
       K2 = 0.337 * log(target_length) + 1.883; 
       K1 = K2 - 0.4; 
       }
  max_l=300; cov=5; 
  output_parameters(max_l, cov); 

  /* 10% */ 
  not_valid=0; 
  if(target_length<=45)
    { 
    L = round(1.376 * pow((double) target_length, 0.799)); 
    K2 = 0.69 * log(target_length) + 0.625; 
    }
  else if(target_length>55) 
         { 
         L = round(1.298 * pow((double) target_length, 0.809)); 
         K2 = 0.347 * log(target_length) + 1.93; 
         }
  else { 
       L = round(((1.376 * pow((double) target_length, 0.799)) + (1.298 * pow((double) target_length, 0.809)))/2.0); 
       K2 = ((0.69 * log(target_length) + 0.625) + (0.347 * log(target_length) + 1.93))/2.0; 
       }
  K1 = K2 - 0.3; 
  max_l=300; cov=10; 
  output_parameters(max_l, cov); 

  /* 25% */ 
  not_valid=0; 
  L = round(1.507 * pow((double) target_length, 0.762)); 
  if(target_length<=45)
    { K2 = 0.476 * log(target_length) + 1.566; } 
  else if(target_length>55)
         { K2 = 0.314 * log(target_length) + 2.221; }
  else { K2 =((0.476 * log(target_length) + 1.566) + (0.314 * log(target_length) + 2.221))/2.0; }
  K1 = K2 - 0.3; 
  max_l=300; cov=25; 
  output_parameters(max_l, cov);   

  /* 40% */ 
  not_valid=0; 
  if(target_length<=55)
    {  
    L = round(1.491 * pow((double) target_length, 0.793)); 
    K2 = 0.581 * log(target_length) + 1.316;
    }
  else if(target_length>65) 
         {  
         L = round(1.138 * pow((double) target_length, 0.86)); 
         K2 = 0.28 * log(target_length) + 2.442;
         }
  else { 
       L = round(((1.491 * pow((double) target_length, 0.793)) + (1.138 * pow((double) target_length, 0.86)))/2.0); 
       K2 = ((0.581 * log(target_length) + 1.316) + (0.28 * log(target_length) + 2.442))/2.0; 
       }
  K1 = K2 - 0.2; 
  max_l=300; cov=40; 
  output_parameters(max_l, cov); 
  } /* end of focus==DIVERSE */ 

else { /* focus == NARROW */ 
     /* 2% */ 
     L=target_length; 

     if(target_length<=45) 
       { K2 = 0.818 * log(L) - 0.245; }
     else if(target_length>55){ K2 = 0.418 * log(L) + 1.206; }
     else { K2 = ((0.818 * log(L) - 0.245) + (0.418 * log(L) + 1.206))/2.0; }
     K1 = K2; 
     max_l=250; cov=2; 
     output_parameters(max_l, cov); 

     /* 5% */ 
     not_valid=0; 
     if(target_length<=45) 
       { K2 = 0.824 * log(L) - 0.003; }
     else if(target_length>55){ K2 = 0.355 * log(L) + 1.731; }
     else { K2 = ((0.824 * log(L) - 0.003) + (0.355 * log(L) + 1.731))/2.0; }
     K1 = K2; 
     max_l=300; cov=5; 
     output_parameters(max_l, cov); 

     /* 10% */ 
     not_valid=0; 
     if(target_length<=45) 
       { K2 = 0.803 * log(L) + 0.251; }
     else if(target_length>55){ K2 = 0.3 * log(L) + 2.135; }
     else { K2 = ((0.803 * log(L) + 0.251) + (0.3 * log(L) + 2.135))/2.0; }
     K1 = K2; 
     max_l=300; cov=10; 
     output_parameters(max_l, cov); 

     /* 25% */ 
     not_valid=0; 
     if(target_length<=45) 
       { K2 = 0.788 * log(L) + 0.499; }
     else if(target_length>55){ K2 = 0.278 * log(L) + 2.405; }
     else { K2 = ((0.788 * log(L) + 0.499) + (0.278 * log(L) + 2.405))/2.0; }
     K1 = K2; 
     max_l=300; cov=25; 
     output_parameters(max_l, cov); 

     /* 40% */ 
     not_valid=0; 
     if(target_length<=45) 
       { K2 = 0.705 * log(L) + 0.887; }
     else if(target_length>55){ K2 = 0.257 * log(L) + 2.596; }
     else { K2 = ((0.705 * log(L) + 0.887) + (0.257 * log(L) + 2.596))/2.0; }
     K1 = K2; 
     max_l=250; cov=40; 
     output_parameters(max_l, cov); 
     } /* end of focus==NARROW */ 


/*  *  *  * FOOTER OF OUTPUT *  *  *  */ 
fprintf(stdout, "\n\nCoverage is the proportion of protein sequences expected to be labelled by these parameter sets.\n"); 
fprintf(stdout, "\nIt is recommended to use all of the parameters progressively in separate runs of the SEG algorithm,\n"); 
fprintf(stdout, " and compare the outputs.\n"); 
fprintf(stdout, "If the calculated parameters are listed as 'NA', it means that at least one of them was out of bounds.\n\n"); 

exit(0); 
} /* end of main() */ 

/******** END OF CODE FILE ********/ 

