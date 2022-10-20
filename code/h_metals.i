# 0 "./code/h_metals.h"
# 1 "/cygdrive/c/Users/ry22aas/robyates/Astro/L-Galaxies/LGalaxies2020_code//"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "./code/h_metals.h"

/*struct elements_str
{
    float H;
    float He;
    float O;
    float Mg;
    float Fe;
};
#endif  //MAINELEMENTS*/
/* Define two views into the same element structure/array.
   For example define:
      union elements DiskMassElements;
   Then access as either
      DiskMassElements.str.H
   or
      DiskMassElements.arr[0]
*/
/*
union elements
{
    struct elements_str str;
    float arr[NUM_ELEMENTS];
};*/



//Number of chemical elements tracked:
