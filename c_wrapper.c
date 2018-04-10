#include <string.h>

void __bertha_wrapper_MOD_bertha_init (char *, int);
void __bertha_wrapper_MOD_bertha_main (char *, char *, char *, char *, 
    int, int, int, int);


int bertha_init (char * filename)
{
  __bertha_wrapper_MOD_bertha_init(filename, strlen(filename));

  return 0;
}

int bertha_main(char * fittcoefffname, char * vctfilename, 
    char * ovapfilename, char * fittfname)
{
  __bertha_wrapper_MOD_bertha_main(fittcoefffname, vctfilename, 
        ovapfilename, fittfname, strlen(fittcoefffname), 
        strlen(vctfilename), strlen(ovapfilename), strlen(fittfname));

  return 0;
}
