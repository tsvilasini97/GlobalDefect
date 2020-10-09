
#include "global_defect.hpp"


int main(int argc, char **argv)
{

  int n,m;

  for (int i=1 ; i < argc ; i++ ){
    if ( argv[i][0] != '-' )continue;
    switch(argv[i][1])
    {
      case 'n':
        n = atoi(argv[++i]);
        break;
      case 'm':
        m =  atoi(argv[++i]);
        break;
    }
  }

  parallel.initialize(n,m);

  
  Global_defect sim(16,2);
  sim.evolve();



  return 0;
}
