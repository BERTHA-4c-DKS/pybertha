#include <string.h>
#include <stdio.h>

#ifdef USEINTELCMP

#else

void __bertha_wrapper_MOD_bertha_init (char *, int *, int *, int);
void __bertha_wrapper_MOD_bertha_realtime_init ();
void __bertha_wrapper_MOD_bertha_main (char *, char *, char *, char *, 
    double *, double *, double *, double *, int, int, int, int);
void __bertha_wrapper_MOD_bertha_finalize();
void __bertha_wrapper_MOD_bertha_realtime_finalize();
void __bertha_wrapper_MOD_bertha_realtime_dipolematrix(int *, 
      int *, double *);
void __bertha_wrapper_MOD_bertha_realtime_fock (double *, double *);

extern int _spec_MOD_ndim, __spec_MOD_nshift, __spec_MOD_nocc, 
       __opensh_MOD_nopen;
extern double __shiftr_MOD_sfact, __energy_MOD_etotal, 
       __bertha_wrapper_MOD_erep;

#define f_bertha_init __bertha_wrapper_MOD_bertha_init
#define f_bertha_realtime_init __bertha_wrapper_MOD_bertha_realtime_init
#define f_bertha_main __bertha_wrapper_MOD_bertha_main
#define f_bertha_finalize __bertha_wrapper_MOD_bertha_finalize
#define f_bertha_realtime_finalize __bertha_wrapper_MOD_bertha_realtime_finalize
#define f__bertha_realtime_dipolematrix __bertha_wrapper_MOD_bertha_realtime_dipolematrix
#define f__bertha_realtime_fock __bertha_wrapper_MOD_bertha_realtime_fock

#define f_ndim __spec_MOD_ndim
#define f_nshift __spec_MOD_nshift
#define f_nocc __spec_MOD_nocc
#define f_nopen __opensh_MOD_nopen

#define f_sfact __shiftr_MOD_sfact
/*
 __energy_MOD_etotal      
 __bertha_wrapper_MOD_erep
*/

#define 
#endif

// DATA METHODS

int get_ndim ()
{
  int val = f_ndim;

  return val;
}

int get_nshift ()
{
  int val = f_nshift;

  return val;
}

int get_nopen ()
{
  int val = f_nopen;

  return val;
}

int get_nocc ()
{
  int val = f_nocc;

  return val;
}

double get_sfact ()
{
  double val = f_sfact;

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

// METHODS

int init (char * filename, int verbosity, int dumpfiles)
{
  f_bertha_init(filename, &verbosity, 
      &dumpfiles, strlen(filename));

  return 0;
}

int realtime_fock (double * dens_ptr, double * fock_ptr)
{
  f_bertha_realtime_fock (dens_ptr, fock_ptr);

  return 0;
}

int realtime_dipolematrix (int direction, int norm, 
    double * vext_ptr)
{
  f_bertha_realtime_dipolematrix(&direction, 
      &norm, vext_ptr);

  return 0;
}

int realtime_init ()
{
  f_bertha_realtime_init();

  return 0;
}

int mainrun(char * fittcoefffname, char * vctfilename, 
    char * ovapfilename, char * fittfname, double * eigen,
    double * ovapin, double * eigenv, double * fockin)
{
  /*
  int i, j, ndim, counter;
  for (i=0; i<10; i++)
    printf("%f \n", a[i]);
  */

  f_bertha_main(fittcoefffname, vctfilename, 
        ovapfilename, fittfname, eigen, ovapin, eigenv, fockin,
        strlen(fittcoefffname), strlen(vctfilename), strlen(ovapfilename), 
        strlen(fittfname));

  /*
  ndim = get_ndim();
  counter = 0;
  for (i=0; i<ndim; ++i)
  {
    for (j=0; j<ndim; ++j)
    {
      //printf ("%10.5f , %10.5f i \n", ovap[counter], ovap[counter+1]);
      ovapin[counter] = ovap[counter];
      ovapin[counter+1] = ovap[counter+1]; 
      counter = counter + 2;
    }
  }

  ndim = get_ndim();
  counter = 0;
  for (i=0; i<ndim; ++i)
  {
    for (j=0; j<ndim; ++j)
    {
      printf ("%15.10f , %15.10fi \n", eigenv[counter], eigenv[counter+1]);
      counter = counter + 2;
    }
  }
  */

  /*
  printf ("%f %f \n", ovap[0], ovap[1]);
  printf ("%f %f \n", ovap[2], ovap[3]);
  */

  return 0;
}


int finalize ()
{
  f_bertha_finalize();

  return 0;
}

int realtime_finalize ()
{
  f_bertha_realtime_finalize();

  return 0;
}
