#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TTree.h>

void readFileAndDraw(){

    TFile *glauberFile = new TFile("../glauber_PbPb_default_v1p5_10k.root","read");
    TTree *tree = (TTree*)glauberFile->Get("nt_Pb_Pb");

    TH1D *hNpart = new TH1D("hNpart","Number of Participants Distribution;Npart;Counts",500,-0.5,499.5);
    tree->Draw("Npart>>hNpart");
    //cout<<"<Npart> = "<<hNpart->GetMean()<<endl;

    TH2D *hNpartB = new TH2D("hNpartB","Npart vs b; b; Npart",200,0,20.0,500,0,500);
    //tree->Draw("Npart:B>>hNpartB", "", "colz");  
    //cout<<"# of entries = "<<hNpartB->GetEntries()<<endl;

}

