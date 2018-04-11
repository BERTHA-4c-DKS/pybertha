#include <string.h>
#include <stdio.h>

void __bertha_wrapper_MOD_bertha_init (char *, int *, int);
void __bertha_wrapper_MOD_bertha_main (char *, char *, char *, char *, 
    double *, int, int, int, int);
void __bertha_wrapper_MOD_bertha_finalize();

extern int __spec_MOD_ndim;

int bertha_init (char * filename, int verbosity)
{
  __bertha_wrapper_MOD_bertha_init(filename, &verbosity, 
      strlen(filename));

  return 0;
}

int bertha_main(char * fittcoefffname, char * vctfilename, 
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

int bertha_finalize ()
{
  __bertha_wrapper_MOD_bertha_finalize();

  return 0;
}
