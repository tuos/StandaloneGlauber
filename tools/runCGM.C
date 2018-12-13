/*
 $idl: runCGM.C 150 2018-02-28 21:55:12Z loizides $

 Version 1.1, copyright 2016 C. Loizides <constantin.loizides@cern.ch>
 
 To run the code, you need to download at least v2.6 of TGlauberMC from http://tglaubermc.hepforge.org/ 
 and place runglauber_v2.6.C in the same directory
 To compile do in the root prompt
 .L runglauber_X.Y.C+
 .L runCGM.C+

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
//#include "runglauber_v2.6.C"
#include "./runglauber.C"
#include <TProfile2D.h>
#endif

//#define calccore
//#define calcdens

void runCGM(const Int_t n,            /*number of events*/
	    const char *sysA,         /*system A*/
	    const char *sysB,         /*system B*/
	    const Double_t signn,     /*NN cross section (mb)*/
	    const Double_t mind,      /*minimum distance between nucleons*/
	    const Int_t nc,           /*number of constituents / degrees of freedom*/
	    const Double_t sigcc,     /*constituent cross section (mb)*/
	    const Int_t type,         /*defines how to distribute degrees of freedom:
				        =0 no recentering, =1 recentering, =2 recentering + exp. rescaling, =3 gaussian profile, 
                                        =4 double gaussian profile, =5 modfied (PHENIX), = 8 free no recentering, =9 free recentering*/
	    const char *fname   = 0,  /*output filename to store ntuple*/
	    const Double_t bmin = 0., /*minimum impact parameter*/
	    const Double_t bmax = 20, /*maximum impact parameter*/
            const char *iname   = 0)  /*input glauber filename*/
{
  cout << "runCQM with" << endl;
  TObjString str(Form("n=%d, %s-%s, sNN=%.2fmb, mind=%.2f, nc=%d, sigcc=%.2fmb, type=%d", 
		      n, sysA, sysB, signn, mind, nc, sigcc, type));
  cout << str.GetName() << endl;

  TGlauberMC *mcg = new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);
  mcg->SetBmin(bmin);
  mcg->SetBmax(bmax);

  const Bool_t calc = 1;
  if (iname) { // read precomputed events
    if (!mcg->ReadNextEvent(calc,iname)) {
      cout << "Could not setup up reading from file, aborting!" << endl;
      return;
    }
  }

#ifdef calccore
  TH1D *h1d = new TH1D("h1d","",100,0,10);
  TProfile2D *prog = new TProfile2D("hArea2",";B;r (fm)",200,0,20,100,0,10);
  prog->SetStats(0);
#endif
#ifdef calcdens
  TH1D *hra1 = new TH1D("hra1",";f(r);r (fm)",250,0,10);
  hra1->SetStats(0);
  TH1D *hra2 = new TH1D("hra2",";f(r);r (fm)",250,0,10);
  hra2->SetStats(0);
  TH1D *hrb1 = new TH1D("hrb1",";f(r);r (fm)",250,0,10);
  hrb1->SetStats(0);
  TH1D *hrb2 = new TH1D("hrb2",";f(r);r (fm)",250,0,10);
  hrb2->SetStats(0);
#endif

  TFile *out = 0;
  if (fname==0) {
    TString tmp(Form("res_%s%s_snn=%.2fmb_nc=%d_scc=%.2fmb_type=%d_cqm.root", 
		      sysA, sysB, signn, nc, sigcc, type));
    out = TFile::Open(tmp,"recreate",tmp,9);
  } else
    out = TFile::Open(fname,"recreate",fname,9);

  if (!out)
    return;
  TNtuple *nt = new TNtuple("nt","nt",
			    "Npart:Npart0:Ncoll:B:Ncpart:Ncpart0:Nccoll:Ap:Ac"
			    ":Ecc1P:Ecc2P:Ecc3P:Ecc4P:Ecc5P:Ecc6P:Ecc7P:Ecc8P:Ecc9P:Ecc1C:Ecc2C:Ecc3C:Ecc4C:Ecc5C:Ecc6C:Ecc7C:Ecc8C:Ecc9C"
#ifdef calccore
			    ":c05:c1:c15:c2:c3:c4:c5:c6:c7:c8:c9"
#endif
);
  nt->SetDirectory(out);

  const Bool_t accAnyNucGlauberEvent = kTRUE;
  const TGlauNucleus *nucA   = mcg->GetNucleusA();
  const Int_t AN             = nucA->GetN();
  const Int_t nAN = AN*nc;
  Double_t *xA = new Double_t[nAN];
  Double_t *yA = new Double_t[nAN];
  Double_t *zA = new Double_t[nAN];
  Int_t *nA = new Int_t[nAN];
  const TGlauNucleus *nucB   = mcg-> GetNucleusB();
  const Int_t BN             = nucB->GetN();
  const Int_t nBN = BN*nc;
  Double_t *xB = new Double_t[nBN];
  Double_t *yB = new Double_t[nBN];
  Double_t *zB = new Double_t[nBN];
  Int_t *nB = new Int_t[nBN];
  Int_t nnperNA = AN;
  Int_t nnperNB = BN;
  Int_t nqperNA = nc;
  Int_t nqperNB = nc;

  TF1 *rad = 0;
  TF1 *rbd = 0;
  Bool_t doRecenter = 0;
  Bool_t doExpScale = 0;
  if (type>=0&&type<2)
    rad = new TF1("rad","[1]*x^2*TMath::Exp([0]*x)",0,5);
  else if (type==5) {
    rad = new TF1("rad","[1]*x^2*TMath::Exp([0]*x)*(1.21466-1.888*x+2.03*x^2)*(1+1.0/x-0.03/x^2)*(1+0.15*x)",0,5);
    doRecenter = 1;
  } else if (type==8||type==9) {
#if _runglauber_ == 3
    rad=nucA->GetFunc1();
    rbd=nucB->GetFunc1();
#else
    rad=nucA->GetFunction();
    rbd=nucB->GetFunction();
#endif
    nnperNA = 1;
    nnperNB = 1;
    nqperNA = nAN;
    nqperNB = nBN;
    if (type==9)
      doRecenter=1;
  }
  cout << "Indices: " << nnperNA << " " << nnperNB << " " << nqperNA << " " << nqperNB << endl;

  if (rad==0) {
    cerr << " unknown type " << type << endl;
    return;
  }

  if (type==1)
    doRecenter = 1;
  else if (type==2) {
    doRecenter = 1;
    doExpScale = 1;
  }

  if (rbd==0) {
    rad->SetNpx(1000);
    rad->SetParameter(0,-4.27);
    rad->SetParameter(1,1);
    if (doRecenter&&doExpScale)
      rad->SetParameter(0,-4.27*TMath::Sqrt(2./3));
    rbd = rad;
  }

  const Double_t d2 = (Double_t)sigcc/(TMath::Pi()*10);
  cout << "d2: " << d2 << endl;

  Int_t totalEvs = 0;
  Int_t accEvs   = 0;
  Double_t x[99999]={0},y[99999]={0},z[99999]={0};
  Double_t xqqvals[99999]={0},yqqvals[99999]={0};
  while (accEvs<n) {
    if (!iname) {
      if (accAnyNucGlauberEvent)
	mcg->NextEvent();
      else
	while (!mcg->NextEvent()) {;}
    } else 
      mcg->ReadNextEvent(calc);

    const TObjArray *nucleonsA = nucA->GetNucleons();
    for (Int_t i = 0,q=0; i<nnperNA; ++i) {

      Double_t xc=-mcg->GetB()/2,yc=0,zc=0;
      if (type<8) {
	TGlauNucleon *nucleonA = (TGlauNucleon*)(nucleonsA->At(i));
	xc = nucleonA->GetX();
	yc = nucleonA->GetY();
	zc = nucleonA->GetZ();
      }

      Double_t sumx=0,sumy=0,sumz=0;
      for (Int_t j=0; j<nqperNA; ++j) { 
	Double_t sr = rad->GetRandom();
	Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	Double_t ctheta = 2*gRandom->Rndm() - 1 ;
	Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
	x[j] = sr*stheta*TMath::Cos(sp);
	y[j] = sr*stheta*TMath::Sin(sp);
	z[j] = sr*ctheta;
        sumx += x[j];
        sumy += y[j];
        sumz += z[j];
#ifdef calcdens
	hra1->Fill(sr);
#endif
      }
      sumx /= nqperNA;
      sumy /= nqperNA;
      sumz /= nqperNA;
      if (doRecenter) {
	for (Int_t j=0; j<nqperNA; ++j) { 
	  x[j] -= sumx;
	  y[j] -= sumy;
	  z[j] -= sumz;
	}
      }
      for (Int_t j=0; j<nqperNA; ++j) { 
	xA[q] = xc + x[j];
	yA[q] = yc + y[j];
	zA[q] = zc + z[j];
#ifdef calcdens
	Double_t r=TMath::Sqrt(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]);
	hra2->Fill(r);
#endif
	nA[q] = 0;
	++q;
      }
    }

    const TObjArray *nucleonsB = nucB->GetNucleons();
    for (Int_t i = 0,q=0; i<nnperNB; ++i) {
      Double_t xc=+mcg->GetB()/2,yc=0,zc=0;
      if (type<8) {
	TGlauNucleon *nucleonB = (TGlauNucleon*)(nucleonsB->At(i));
	xc = nucleonB->GetX();
	yc = nucleonB->GetY();
	zc = nucleonB->GetZ();
      }
      Double_t sumx=0,sumy=0,sumz=0;
      for (Int_t j=0; j<nqperNB; ++j) { 
	Double_t sr = rbd->GetRandom();
	Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	Double_t ctheta = 2*gRandom->Rndm() - 1 ;
	Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
	x[j] = sr*stheta*TMath::Cos(sp);
	y[j] = sr*stheta*TMath::Sin(sp);
	z[j] = sr*ctheta;
        sumx += x[j];
        sumy += y[j];
        sumz += z[j];
#ifdef calcdens
	hrb1->Fill(sr);
#endif
      }
      sumx /= nqperNB;
      sumy /= nqperNB;
      sumz /= nqperNB;
      if (doRecenter) {
	for (Int_t j=0; j<nqperNB; ++j) { 
	  x[j] -= sumx;
	  y[j] -= sumy;
	  z[j] -= sumz;
	}
      }
      for (Int_t j=0; j<nqperNB; ++j) { 
	xB[q] = xc + x[j];
	yB[q] = yc + y[j];
	zB[q] = zc + z[j];
#ifdef calcdens
	Double_t r=TMath::Sqrt(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]);
	hrb2->Fill(r);
#endif
        nB[q] = 0;
	++q;
      }
    }

    Int_t ncollqq=0;
    for (Int_t i = 0; i<nBN; ++i) {
      for (Int_t j = 0; j<nAN; ++j) {
	Double_t dx = xB[i]-xA[j];
	Double_t dy = yB[i]-yA[j];
	Double_t dij = dx*dx+dy*dy;
	if (dij < d2) {
	  ++nA[j];
	  ++nB[i];
	  ++ncollqq;
	}
      }
    }

    Int_t npartqq = 0;
    Int_t npart0qq = 0;
    for (Int_t i = 0; i<nAN; ++i) {
      if (nA[i]==0)
	continue;
      xqqvals[npartqq] = xA[i];
      yqqvals[npartqq] = yA[i];
      ++npartqq;
      if (nA[i]>1)
	++npart0qq;
    }
    for (Int_t i = 0; i<nBN; ++i) {
      if (nB[i]==0)
	continue;
      xqqvals[npartqq] = xB[i];
      yqqvals[npartqq] = yB[i];
      ++npartqq;
      if (nB[i]>1)
	++npart0qq;
    }

    // *** event counting ***
    ++totalEvs;
    if (npartqq==0)
      continue;
    ++accEvs;
    if (!(accEvs%50)) {
      Double_t bm = mcg->GetBmax();
      Double_t totxs = (1.*accEvs/totalEvs)*TMath::Pi()*bm*bm/100;
      Double_t totes = totxs/TMath::Sqrt(accEvs)*TMath::Sqrt(1.-accEvs/totalEvs);
      if (!iname) {
	Double_t ratio = totxs/mcg->GetTotXSect()*100;
	cout << "Event # " << accEvs << " c-sect = " << totxs << " +- " << totes << " b (" << ratio << "%)                \r" << flush;
      } else 
	cout << "Event # " << accEvs << "                 \r" << flush;
    }
    
    Double_t MeanXQQParts   = 0;
    Double_t MeanYQQParts   = 0;
    Double_t MeanXYQQParts  = 0;
    Double_t MeanXQQParts2  = 0;
    Double_t MeanYQQParts2  = 0;
    Double_t sx2qq          = 0;
    Double_t sy2qq          = 0;
    Double_t sxyqq          = 0;

   for (Int_t i = 0; i<npartqq; ++i) {
      MeanXQQParts   += xqqvals[i];
      MeanYQQParts   += yqqvals[i];
      MeanXYQQParts  += xqqvals[i]*yqqvals[i];
      MeanXQQParts2  += xqqvals[i]*xqqvals[i];
      MeanYQQParts2  += yqqvals[i]*yqqvals[i];
    }
    MeanXQQParts  /= npartqq;
    MeanYQQParts  /= npartqq;
    MeanXQQParts2 /= npartqq;
    MeanYQQParts2 /= npartqq;
    MeanXYQQParts /= npartqq;
    sx2qq += MeanXQQParts2-MeanXQQParts*MeanXQQParts;
    sy2qq += MeanYQQParts2-MeanYQQParts*MeanYQQParts;
    sxyqq += MeanXYQQParts-MeanXQQParts*MeanYQQParts;

    Double_t sinphi[10] = {0};
    Double_t cosphi[10] = {0};
    Double_t rn[10]     = {0};
    Double_t ecc[10]    = {0};
    Double_t psi[10]    = {0};
    for (Int_t j=1; j<10; ++j) {
      for (Int_t i = 0; i<npartqq; ++i) {
        Double_t x   = xqqvals[i] - MeanXQQParts;
        Double_t y   = yqqvals[i] - MeanYQQParts;
        Double_t r   = TMath::Sqrt(x*x+y*y);
        Double_t phi = TMath::ATan2(y,x);
        Double_t w = j;
        if (j==1) 
          w = 3; // use r^3 weighting for Ecc1/Psi1
        cosphi[j] += TMath::Power(r,w)*TMath::Cos(j*phi);
        sinphi[j] += TMath::Power(r,w)*TMath::Sin(j*phi);
        rn[j]     += TMath::Power(r,w);
      }
    }
    for (Int_t j=1; j<10; ++j) {
      psi[j] = (TMath::ATan2(sinphi[j],cosphi[j]) + TMath::Pi())/j;
      ecc[j] = TMath::Sqrt(sinphi[j]*sinphi[j] + cosphi[j]*cosphi[j]) / rn[j];
    }
#ifdef calccore
    h1d->Reset();
    for (Int_t i = 0; i<npartqq; ++i) {
      Double_t dx   = xqqvals[i] - MeanXQQParts;
      Double_t dy   = yqqvals[i] - MeanYQQParts;
      Double_t d2 = TMath::Sqrt(dx*dx+dy*dy);
      h1d->Fill(d2);
    }
    for (Int_t i=1;i<=h1d->GetNbinsX();i++) {
      Double_t entries = h1d->GetBinContent(i);
      if (entries>0) {
        prog->Fill(mcg->GetB(),h1d->GetBinCenter(i),entries,1);
      }
    }
#endif

    Double_t ap = mcg->GetSx2()*mcg->GetSy2()-mcg->GetSxy()*mcg->GetSxy();
    if (ap>=0)
      ap = TMath::Sqrt(ap);
    else 
      ap = -1;
    Double_t aq = sx2qq*sy2qq-sxyqq*sxyqq;
    if (aq>=0) 
      aq = TMath::Sqrt(aq);
    else
      aq = -1;

    Float_t v[99]; Int_t f=0;
    v[f++] = mcg->GetNpart();  //number of nucleon participants
    v[f++] = mcg->GetNpart0(); //number of nucleon core participants
    v[f++] = mcg->GetNcoll();  //number of nucleon collisions
    v[f++] = mcg->GetB();      //impact parameter
    v[f++] = npartqq;          //number of constituent participants
    v[f++] = npart0qq;         //number of constituent core participants
    v[f++] = ncollqq;          //number of constituent collisions
    v[f++] = ap;               //area as defined from participant (co-)variances
    v[f++] = aq;               //area as defined from constituent (co-)variances
    v[f++] = mcg->GetEcc(1);   //ecc1 nucleon participants
    v[f++] = mcg->GetEcc(2);   //ecc2 nucleon participants
    v[f++] = mcg->GetEcc(3);   //ecc3 nucleon participants
    v[f++] = mcg->GetEcc(4);   //ecc4 nucleon participants
    v[f++] = mcg->GetEcc(5);   //ecc5 nucleon participants
    v[f++] = mcg->GetEcc(6);   //ecc6 nucleon participants
    v[f++] = mcg->GetEcc(7);   //ecc7 nucleon participants
    v[f++] = mcg->GetEcc(8);   //ecc8 nucleon participants
    v[f++] = mcg->GetEcc(9);   //ecc9 nucleon participants
    v[f++] = ecc[1];           //ecc1 constituent participants
    v[f++] = ecc[2];           //ecc2 constituent participants
    v[f++] = ecc[3];           //ecc3 constituent participants
    v[f++] = ecc[4];           //ecc4 constituent participants
    v[f++] = ecc[5];           //ecc5 constituent participants
    v[f++] = ecc[6];           //ecc6 constituent participants
    v[f++] = ecc[7];           //ecc7 constituent participants
    v[f++] = ecc[8];           //ecc8 constituent participants
    v[f++] = ecc[9];           //ecc9 constituent participants
#ifdef calccore
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(0.5)-1);
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(1)-1);
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(1.5)-1);
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(2)-1);
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(3)-1);
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(4)-1);
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(5)-1);
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(6)-1);
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(7)-1);
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(8)-1);
    v[f++] = h1d->Integral(h1d->FindBin(0),h1d->FindBin(9)-1);
#endif
    nt->Fill(v);
  }

#ifdef calccore
  delete h1d;
  prog->Write();
#endif
#ifdef calcdens
  hra1->Scale(1./hra1->GetEntries());
  hra1->Write();
  hra2->Scale(1./hra2->GetEntries());
  hra2->Write();
  hrb1->Scale(1./hrb1->GetEntries());
  hrb1->Write();
  hrb2->Scale(1./hrb2->GetEntries());
  hrb2->Write();
#endif

  cout << endl;
  cout << str.GetName() << endl;
  str.Write();
  Double_t bm = mcg->GetBmax();
  if (!iname) {
    Double_t totxs = (1.*accEvs/totalEvs)*TMath::Pi()*bm*bm/100;
    Double_t totes = totxs/TMath::Sqrt(accEvs)*TMath::Sqrt(1.-accEvs/totalEvs);
    Double_t ratio = totxs/mcg->GetTotXSect()*100;
    TObjString s1(Form("c-sect = %.5f +- %.5f b", totxs, totes));
    TObjString s2(Form("x-sect = %.5f +- %.5f b", mcg->GetTotXSect(), mcg->GetTotXSectErr()));
    TObjString s3(Form("ratio = %.5f", ratio));
    cout << s1.GetName() << endl;
    cout << s2.GetName() << endl;
    cout << s3.GetName() << endl;
    s1.Write();
    s2.Write();
    s3.Write();
  }
  out->Write();
  out->Close();
  delete out;
}


