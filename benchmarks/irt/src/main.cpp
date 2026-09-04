#include "Setup.h"
#include "RICHIRT.h"
int main(int argc, char* argv[]){
  if(argc!=4){
    std::cerr<<"Usage: "<<argv[0]<<" <RICHString> <infile> <outfile>"<<std::endl;
    std::cerr<<"  <RICHString>  DRICH or PFRICH"<<std::endl;
    std::cerr<<"  <infile>      reconstructed edm4eic file"<<std::endl;
    std::cerr<<"  <outfile>     output file with histograms (written as given)"<<std::endl;
    return 1;
  }
  try{
    RICHIRT irt;

    irt.set_irtRICH(argv[1]);
    irt.set_inFile(argv[2]);
    irt.set_outFile(argv[3]);
    irt.HistoInit();
    irt.SetCollections();
    irt.Process();
  }catch(const std::exception& e){
    std::cerr<<"Error: "<<e.what()<<std::endl;
    return 1;
  }
  return 0;
}
