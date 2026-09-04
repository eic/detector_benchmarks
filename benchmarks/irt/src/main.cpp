#include "Setup.h"
#include "RICHIRT.h"
int main(int argc, char* argv[]){
  if(argc!=4){
    std::cerr<<"Usage: "<<argv[0]<<" <RICHString> <infile> <outfile>"<<std::endl;
    return 1;
  }
  RICHIRT irt;

  irt.set_irtRICH(argv[1]);
  irt.set_inFile(argv[2]);
  irt.set_outFile(argv[3]);
  irt.HistoInit();
  irt.SetCollections();
  irt.Process();
}
