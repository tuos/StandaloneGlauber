#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TTree.h>
#include <algorithm> 
#include <vector>

void centralityDetermination(){

    TFile *glauberFile = new TFile("../glauber_PbPb_default_v1p5_10k.root","read");
    TTree *tree = (TTree*)glauberFile->Get("nt_Pb_Pb");

    TH1D *hNpart = new TH1D("hNpart","Number of Participants Distribution;Npart;Counts",500,-0.5,499.5);
    TH2D *hNpartB = new TH2D("hNpartB","Npart vs b; b; Npart",200,0,20.0,500,0,500);
    Float_t Npart, Ncoll, B;
    tree->SetBranchAddress("Npart", &Npart);
    tree->SetBranchAddress("Ncoll", &Ncoll);
    tree->SetBranchAddress("B", &B);

    const int nBins = 10;
    double binBoundaries[nBins+1];
    vector<float> impactParameters;
    
    Int_t nEvent = tree->GetEntries();
    for(int i=0; i<nEvent; i++){

      if(i%1000==0)  cout<<"Have run "<<i<<" out of "<<nEvent<<" events; "<<endl;
      tree->GetEntry(i);
    
      hNpart->Fill(Npart);
      hNpartB->Fill(B, Npart);
     
      impactParameters.push_back(B);
    }
    sort(impactParameters.begin(),impactParameters.end());

    int size = impactParameters.size(); 
    binBoundaries[nBins] = impactParameters[size-1]+0.1;

    cout <<endl << "Impact parameter boundaries for each centrality range:"<<endl;
    for(int i = 0; i < nBins; i++) {
      int entry = (int)(i*(size/nBins));
      if(entry < 0 || i == 0) binBoundaries[i] = 0;
      else binBoundaries[i] = impactParameters[entry];
      cout << binBoundaries[i] << ", ";
    }
    cout << binBoundaries[nBins] << endl << endl; 


    TH1D *hCentrality = new TH1D("hCentrality","Number of Event in each Centrality bin;Centrality;Counts",nBins,0,100);

    cout<<"Event loop 2, to fill centrality histogram"<<endl;
    for(int i=0; i<nEvent; i++){

      if(i%1000==0)  cout<<"Have run "<<i<<" out of "<<nEvent<<" events in the Event loop 2;"<<endl;
      tree->GetEntry(i);

      for(int j = 0; j < nBins; j++){
        if(B>=binBoundaries[j] && B<binBoundaries[j+1])
          hCentrality->Fill(j*100/nBins);
      }
    }
    cout<<endl<<"Bin       Number of events"<<endl;
    for(int i = 0; i < nBins; i++) {
      cout<<i<<"      "<<hCentrality->GetBinContent(i+1)<<endl;
    } 
    cout<<endl<<endl;

    hCentrality->SetMinimum(0);
    hCentrality->Draw();

    //cout<<"<Npart> = "<<hNpart->GetMean()<<endl;

    //cout<<"# of entries = "<<hNpartB->GetEntries()<<endl;
   

}

