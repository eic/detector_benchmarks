#include <TCanvas.h>
#include <TList.h>
#include <TString.h>

void SaveHistogramsFromList(TList &l, TString prefix) {
    TObject* obj = nullptr;
    TIter next(&l);

    while ((obj = next())) {
        if (obj->InheritsFrom(TH1D::Class()) || obj->InheritsFrom(TH2D::Class())) {
            TH1* hist = (TH1*)obj;

            TString name = hist->GetName();
            TString canvasName = name + "_canvas";

            TCanvas* c = new TCanvas(canvasName, canvasName, 800, 600);
            c->cd();

            if (hist->InheritsFrom(TH2::Class()))
                hist->Draw("COLZ");
            else
                hist->Draw();

            TString filename = prefix + name + ".png";
            c->SaveAs(filename);

            delete c;
        }
    }
}

