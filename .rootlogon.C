{
	// top-level include-dir
  gROOT->ProcessLine(".include include");

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
