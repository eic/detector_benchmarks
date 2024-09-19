{
	// top-level include-dir
  gROOT->ProcessLine(".include include");

  // setup a local build directory so we don't polute our source code with
  // ROOT dictionaries etc. if desired
  const char* build_dir = gSystem->Getenv("ROOT_BUILD_DIR");
  if (build_dir) {
    gSystem->SetBuildDir(build_dir);
  }

  R__LOAD_LIBRARY(fmt)
  R__LOAD_LIBRARY(HepMC3)

  // Setting for Graphs
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1);
  gStyle->SetLineWidth(2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadLeftMargin(0.14);
}
