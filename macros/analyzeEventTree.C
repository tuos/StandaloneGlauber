#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TTree.h>

void analyzeEventTree(){

    TFile *glauberFile = new TFile("../glauber_PbPb_default_v1p5_10k.root","read");
    TTree *tree = (TTree*)glauberFile->Get("nt_Pb_Pb");

    TH1D *hNpart = new TH1D("hNpart","Number of Participants Distribution;Npart;Counts",500,-0.5,499.5);
    TH1D *hB = new TH1D("hB","Impact Parameter Distribution;b (fm);Counts",200,0,20);
    TH2D *hNpartB = new TH2D("hNpartB","Npart vs b; b (fm); Npart",200,0,20.0,500,0,500);
    Float_t Npart, Ncoll, B;
    tree->SetBranchAddress("Npart", &Npart);
    tree->SetBranchAddress("Ncoll", &Ncoll);
    tree->SetBranchAddress("B", &B);

    Int_t nEvent = tree->GetEntries();
    for(int i=0; i<nEvent; i++){

      if(i%1000==0)  cout<<"Have run "<<i<<" out of "<<nEvent<<" events; "<<endl;
      tree->GetEntry(i);
    
      hNpart->Fill(Npart);
      hB->Fill(B);
      hNpartB->Fill(B, Npart);
    }

    cout<<"<Npart> = "<<hNpart->GetMean()<<endl;

    cout<<"# of entries = "<<hNpartB->GetEntries()<<endl;

    hB->Draw();
    //hNpartB->Draw("colz");

}

