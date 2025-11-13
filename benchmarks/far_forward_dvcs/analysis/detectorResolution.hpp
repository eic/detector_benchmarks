#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

TH1D* extractResolution(TString outputHistoName, TH2D* twoDHisto, bool saveOutput = false){

	int num_bins  = twoDHisto->GetNbinsX();
	double xBinWidth = twoDHisto->GetXaxis()->GetBinWidth(1); 
	double xMin = twoDHisto->GetXaxis()->GetBinCenter(1) - xBinWidth*0.5;
	double xMax = twoDHisto->GetXaxis()->GetBinCenter(num_bins) + xBinWidth*0.5;

	TH1D * finalResoHisto = new TH1D(outputHistoName, outputHistoName, num_bins, xMin, xMax);

	TH1D* tmp;
	double rmsReso = 0.0;
	double rmsErr = 0.0;
	for(int bin = 0; bin < num_bins+1; bin++){
	
		rmsReso = 0.0;
		tmp = (TH1D*)twoDHisto->ProjectionY("NEIN", bin, bin+1);
		TF1 * func = new TF1("fitFunc", "gaus", -1.0, 1.0);
		tmp->Fit(func, "Q0");
		
		//TCanvas * canRes = new TCanvas("canRes1", "canRes1", 500, 500);
		//tmp->Draw();
		//if(saveOutput){
		//	if(outputHistoName == "b0_extracted_pt_resolution"){
		//		canRes->SaveAs(Form("./resolutionFitOutput/pt_resolution_fits_bin_%d.pdf", bin));
		//	}
		//	else canRes->SaveAs(Form("./resolutionFitOutput/p_resolution_fits_bin_%d.pdf", bin));
		//}
		//rmsReso = tmp->GetRMS();
		//rmsErr  = tmp->GetRMSError();

		rmsReso = func->GetParameter(2);
		rmsErr  = func->GetParError(2);

		cout << "resolution from Gaussian fit = " << rmsReso << endl;
		cout << "error on resolution = " << rmsErr << endl;
	
		if(rmsErr/rmsReso > 0.5){ continue; }
	
		if(rmsErr > rmsReso){
			rmsReso = tmp->GetRMS();
			rmsErr  = tmp->GetRMSError();
		}
		
		finalResoHisto->SetBinContent(bin, rmsReso);
		finalResoHisto->SetBinError(bin, rmsErr);
		
		delete func;
	} 


	return finalResoHisto;

}
