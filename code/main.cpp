
#include "global_defect.hpp"


int main(int argc, char **argv)
{

  int n,m;
  string settings_filename;

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
      case 's':
        settings_filename = argv[++i];
        break;
    }
  }

  parallel.initialize(n,m);


  Global_defect sim(settings_filename);

  sim.evolve();



  return 0;
}
