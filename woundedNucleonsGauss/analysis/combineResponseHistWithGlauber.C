#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <TF1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TTree.h>
#include <iomanip> 

using namespace std;

//bool descend(float i,float j) { return (i<j); }

const int nBins = 14; // CMS paper centrality ranges
double centralityRanges [nBins+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100};
TH1D* hNpartBin[nBins];
TH1D* hNcollBin[nBins];
TH1D* hSpart4Bin[nBins];
TH1D* hEcc2partBin[nBins];

//------------------------------------------------------------------------
void combineResponseHistWithGlauber() {

  const char* outFileName = "out_glauber.root";

  for (int ibin = 0; ibin < nBins; ibin++) {
   hNpartBin[ibin] = new TH1D(Form("npartBin_%d",ibin),"",450,0,450);
   hNcollBin[ibin] = new TH1D(Form("ncollBin_%d",ibin),"",2500,0,2500);
   hSpart4Bin[ibin] = new TH1D(Form("spart4Bin_%d",ibin),"",150,0,150);
   hEcc2partBin[ibin] = new TH1D(Form("ecc2partBin_%d",ibin),"",200,0,1);
  }
 
const char * infilename;
char fdirname[200];

sprintf(fdirname,"/scratch/tuos/glauber/CMSSW_10_3_0/src/v3p2/standalone/StandaloneGlauber");
infilename = Form("%s/slurm/woundedGaussOnly/glauber_PbPb_woundedGauss_v1p5_100k.root",fdirname);

 TChain * t = new TChain("nt_Pb_Pb");
 t->Add(infilename);

 //output
 TFile * outf = new TFile(outFileName,"RECREATE");
 ofstream txtfile("text_output_gausslike.txt");

 //Setting up variables and branches
 double binboundaries[nBins+1];
 vector<float> values;

 float b, npart, ncoll, spart4, ecc2part;

 t->SetBranchAddress("B",&b);
 t->SetBranchAddress("Npart",&npart);
 t->SetBranchAddress("Ncoll",&ncoll);
 t->SetBranchAddress("SPART4",&spart4);
 t->SetBranchAddress("Ecc2PART",&ecc2part);

 //Event loop 1
 unsigned int Nevents = t->GetEntries();
// Nevents=10000;
 txtfile << "Number of events = " << Nevents << endl << endl;
 for(unsigned int iev = 0; iev < Nevents; iev++) {
   if(iev%50000 == 0) cout<<"Processing event: " << iev << endl;
   t->GetEntry(iev);

   values.push_back(b);
 }

 //sort(values.begin(),values.end(),descend);
 sort(values.begin(),values.end());

 //Finding the bin boundaries
 txtfile << "-------------------------------------" << endl;
 txtfile << "b" << " based cuts are: " << endl;
 txtfile << "(";
//for(unsigned int iev = 0; iev < Nevents; iev++) txtfile<<values[iev]<<endl;
//txtfile<<endl<<endl<<endl;

 int size = values.size(); 
 //cout<<"size="<<size<<endl;
 binboundaries[nBins] = values[size-1];

 for(int i = 0; i < nBins; i++) {
   int entry = (int)(centralityRanges[i]*0.01*(size));
   //cout<<"entry="<<entry<<",  entry/size="<<entry*1.0/size<<endl;
   if(entry < 0 || i == 0) binboundaries[i] = 0;
   else binboundaries[i] = values[entry];
   txtfile << binboundaries[i] << ", ";
 }
 txtfile << binboundaries[nBins] << ")" << endl << "-------------------------------------" << endl;



 TH2D* hNpart = new TH2D("hNpart","",nBins,binboundaries,45000,0,450); 
 TH2D* hNcoll = new TH2D("hNcoll","",nBins,binboundaries,50000,0,2500);

 //Event loop 2
 for(unsigned int iev = 0; iev < Nevents; iev++) {
   if( iev % 50000 == 0 ) cout<<"Processing event : " << iev << endl;
   t->GetEntry(iev);

    hNpart->Fill(b,npart);
    hNcoll->Fill(b,ncoll);

    int bin = hNpart->GetXaxis()->FindBin(b) - 1;
    if(bin < 0) bin = 0;
    if(bin >= nBins) bin = nBins - 1;
   
    hNpartBin[bin]->Fill(npart);
    hNcollBin[bin]->Fill(ncoll);
    hSpart4Bin[bin]->Fill(spart4);
    hEcc2partBin[bin]->Fill(ecc2part);

 }

 txtfile<<"-------------------------------------"<<endl;
 txtfile<<"CentLow CentHigh NpartMean NcollMean Spart4Mean Ecc2PartMean BinEdge"<<endl;
 for(int i = 0; i < nBins; i++){ 

   txtfile << centralityRanges[i]<<"       "<< centralityRanges[i+1] << "       " << hNpartBin[i]->GetMean() << "   " << hNcollBin[i]->GetMean() << "    " << hSpart4Bin[i]->GetMean() << "    " << hEcc2partBin[i]->GetMean() << "   " << binboundaries[i+1] <<endl;
 }
 txtfile.close();
 
 outf->cd();
  for (int ibin = 0; ibin < nBins; ibin++) {
   hNpartBin[ibin]->Write(); 
   hNcollBin[ibin]->Write();
   hSpart4Bin[ibin]->Write();
   hEcc2partBin[ibin]->Write();
  }
   hNpart->Write();
   hNcoll->Write();
   outf->Write();
   outf->Close();


   TFile * outfile = new TFile("hist_vs_centrality.root","update");
   outfile->cd();
   TH1D *npartVsCent = new TH1D("npartvscent","npartvscent",nBins,centralityRanges);
   TH1D *ncollVsCent = new TH1D("ncollvscent","ncollvscent",nBins,centralityRanges);
   TH1D *spart4VsCent = new TH1D("spart4vscent","spart4vscent",nBins,centralityRanges);
   TH1D *ecc2partVsCent = new TH1D("ecc2partvscent","ecc2partvscent",nBins,centralityRanges);
   for(int i = 0; i < nBins; i++){
      npartVsCent->SetBinContent(i+1, hNpartBin[i]->GetMean());
      ncollVsCent->SetBinContent(i+1, hNcollBin[i]->GetMean());
      spart4VsCent->SetBinContent(i+1, hSpart4Bin[i]->GetMean());
      ecc2partVsCent->SetBinContent(i+1, hEcc2partBin[i]->GetMean());
   }
   npartVsCent->Write();
   ncollVsCent->Write();
   spart4VsCent->Write();
   ecc2partVsCent->Write();

   outfile->Write();
   outfile->Close();


}

