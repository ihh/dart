#include "protpal/Placement.h"
#include "tree/phylogeny.h"
int main(int argc, char* argv[])
{
  try{
   // Create main placement object 
    Placement placement(argc, argv);
    // Run the analysis
    placement.Run();
  }
  
  catch(const Dart_exception& e)
    {
      // Exception during placement activities                                                  
      cerr<<e.what();
      exit(1);
    }
  return 0; 
}
