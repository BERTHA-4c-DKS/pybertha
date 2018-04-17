#include <string.h>
#include <stdio.h>

void __bertha_wrapper_MOD_bertha_init (char *, int *, int);
void __bertha_wrapper_MOD_bertha_main (char *, char *, char *, char *, 
    double *, int, int, int, int);
void __bertha_wrapper_MOD_bertha_finalize();

extern int __spec_MOD_ndim, __spec_MOD_nshift, __spec_MOD_nocc, 
       __opensh_MOD_nopen;
extern double __shiftr_MOD_sfact, __energy_MOD_etotal, 
       __bertha_wrapper_MOD_erep;

int init (char * filename, int verbosity)
{
  __bertha_wrapper_MOD_bertha_init(filename, &verbosity, 
      strlen(filename));

  return 0;
}

int mainrun(char * fittcoefffname, char * vctfilename, 
    char * ovapfilename, char * fittfname, double * a)
{
  /*
  int i;

  for (i=0; i<10; i++)
    printf("%f \n", a[i]);
  */

  __bertha_wrapper_MOD_bertha_main(fittcoefffname, vctfilename, 
        ovapfilename, fittfname, a, strlen(fittcoefffname), 
        strlen(vctfilename), strlen(ovapfilename), strlen(fittfname));

  return 0;
}

int get_ndim ()
{
  int val = __spec_MOD_ndim;

  return val;
}

int get_nshift ()
{
  int val = __spec_MOD_nshift;

  return val;
}

int get_nopen ()
{
  int val = __opensh_MOD_nopen;

  return val;
}

int get_nocc ()
{
  int val = __spec_MOD_nocc;

  return val;
}

double get_sfact ()
{
  double val = __shiftr_MOD_sfact;

  return val;
}

double get_erep ()
{
  double val = __bertha_wrapper_MOD_erep;

  return val;
}



double get_etotal ()
{
  double val = __energy_MOD_etotal;

  return val;
}



int finalize ()
{
  __bertha_wrapper_MOD_bertha_finalize();

  return 0;
}
