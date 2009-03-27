#include <stdio.h> 
#include <math.h>

int main (int argc, char** argv)
{
  int k = (int) (sqrt(4.) + .5);
  if (k == 2)
    printf ("ok\n");
  else
    printf ("not ok\n");
  return 0;
}
