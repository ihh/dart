#include<iostream>
#include<fstream>
#include "protpal/Placement.h"

int main(int argc, char* argv[])
{
  try{
    
   // create main placement object and run the analysis
    Placement placement(argc, argv);
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
