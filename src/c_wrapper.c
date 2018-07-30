#include <string.h>
#include <stdio.h>

#ifdef USEINTELCMP

#define f_bertha_init bertha_wrapper_mp_bertha_init_
#define f_bertha_realtime_init bertha_wrapper_mp_bertha_realtime_init_
#define f_bertha_main bertha_wrapper_mp_bertha_main_
#define f_bertha_finalize bertha_wrapper_mp_bertha_finalize_
#define f_bertha_realtime_finalize bertha_wrapper_mp_bertha_realtime_finalize_
#define f_bertha_realtime_dipolematrix bertha_wrapper_mp_bertha_realtime_dipolematrix_
#define f_bertha_realtime_fock bertha_wrapper_mp_bertha_realtime_fock_

#define f_ndim spec_mp_ndim_
#define f_nshift spec_mp_nshift_
#define f_nocc spec_mp_nocc_
#define f_nopen opensh_mp_nopen_
#define f_sfact shiftr_mp_sfact_
#define f_etotal energy_mp_etotal_
#define f_erep bertha_wrapper_mp_erep_
#define f_tresh bertha_wrapper_mp_tresh_

#else

#define f_bertha_init __bertha_wrapper_MOD_bertha_init
#define f_bertha_realtime_init __bertha_wrapper_MOD_bertha_realtime_init
#define f_bertha_main __bertha_wrapper_MOD_bertha_main
#define f_bertha_finalize __bertha_wrapper_MOD_bertha_finalize
#define f_bertha_realtime_finalize __bertha_wrapper_MOD_bertha_realtime_finalize
#define f_bertha_realtime_dipolematrix __bertha_wrapper_MOD_bertha_realtime_dipolematrix
#define f_bertha_realtime_fock __bertha_wrapper_MOD_bertha_realtime_fock

#define f_ndim __spec_MOD_ndim
#define f_nshift __spec_MOD_nshift
#define f_nocc __spec_MOD_nocc
#define f_nopen __opensh_MOD_nopen
#define f_sfact __shiftr_MOD_sfact
#define f_etotal __energy_MOD_etotal
#define f_erep __bertha_wrapper_MOD_erep
#define f_tresh __bertha_wrapper_MOD_tresh
#endif

void f_bertha_init (char *, int *, int *, int);
void f_bertha_realtime_init ();
void f_bertha_main (char *, char *, char *, char *, 
    double *, double *, double *, double *, int, int, int, int);
void f_bertha_finalize();
void f_bertha_realtime_finalize();
void f_bertha_realtime_dipolematrix(int *, 
      int *, double *);
void f_bertha_realtime_fock (double *, double *);

extern int f_ndim,f_nshift, f_nocc, f_nopen;
extern double f_sfact, f_etotal, f_erep, f_tresh;

// DATA METHODS

double get_tresh ()
{
  double val = f_tresh;

  return val;
}

void set_tresh (double val)
{
  f_tresh = val;
}

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
  double val = f_erep;

  return val;
}

double get_etotal ()
{
  double val = f_etotal;

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
