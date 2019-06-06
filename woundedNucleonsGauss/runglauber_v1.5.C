// http://www.hepforge.org/downloads/tglaubermc
//run:
//root [1] .L runglauber_v1.5.C+                                                  
//root [2] runAndSaveNtuple(1000000,"p","Pb",70.0,0.4,"Phob_Glau_pPb_sNN70mb_v15_1M_dmin04.root")
//
//------------
// Version 1.5: modifications by djh to include additional "eccentricity" terms 1,4,5,6
//     Ecc1PART = v[26] = are the other harmonics, calculated as done below
//     Psi1PART = v[27] 
//     Ecc4PART = v[28] 
//     Psi4PART = v[29] 
//     Ecc5PART = v[30] 
//     Psi5PART = v[31] 
//     Ecc6PART = v[32] 
//     Psi6PART = v[33] 
//
//------------
// Version 1.4: modification by djh to include triangularity calculation.  Added the following:
//     Ecc2PART = v[22] = participant eccentricity, but calculated with formalism for the triangularity as detailed in Phys.Rev.C81:054905(2010) Equation 4
//     Psi2PART = v[23] = eccentricity minor axis from Phys.Rev.C81:054905(2010) Equation 5
//     Ecc3PART = v[24] = participant triangularity from Phys.Rev.C81:054905(2010) Equation 7
//     Psi3PART = v[25] = triangularity minor axis from Phys.Rev.C81:054905(2010) Equation 9
//
//------------
// Version 1.3:  modification by djh to include two definitions of the transverse overlap area of two nuclei.  Added the following
//     SPART = v[19] = overlap area that is relevant if one uses the "participant" eccentricity EccPART.  See: PRC77, 014906 (2008)
//     SPART4 = v[20] = overlap area above multiplied by a factor of 4 to better match Valdimir's result (and because we think the factor of 4 better matches the actual geometrical overlap)
//
//     S12RP_Npart = v[21] = from Vladimir Korotkikh. Eur.Phys.J.C66:173,2010; arXiv:0910.3029 [hep-ph]; Phys.Rev.C76 (2007) 024905: but here it is the overlap area using participants (not collisions!) and calculated relative to the original reaction plane of the nuclei (i.e. not modified event-by-event to match the maximal axis of the participant eccentricity)
//
//------------
// Version 1.2:  modification by david j hofman to include eccentricity in output ntuple.  Added the following
//      EccRP = v[17] = (fSy2-fSx2)/(fSy2+fSx2);  //v1.2  
//              THIS is the original definition of the eccentricity from participating nucleons - typically not used any more in physics results for flow.
//
//      EccPART = v[18] = TMath::Sqrt((fSy2-fSx2)*(fSy2-fSx2)+4.*fSxy*fSxy)/(fSy2+fSx2);  //v1.2  
//              THIS is the "participant" eccentricity, calculated in a plane that is rotated event-by-event to maximize the eccentricity as defined by the participants.  This variable seems relevant as it unifies Cu+Cu and Pb+Pb flow data from RHIC.  See: PRC77, 014906 (2008)
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TNamed.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TString.h>
#include <TEllipse.h>
#endif

#ifndef _runglauber_
#if !defined(__CINT__) || defined(__MAKECINT__)
#define _runglauber_
#endif

//---------------------------------------------------------------------------------
void runAndSaveNtuple(Int_t n,
                      const char *sysA="Au",
                      const char *sysB="Au",
                      Double_t signn=42,
                      Double_t mind=0.4,
                      const char *fname="glau_auau_ntuple.root");

//---------------------------------------------------------------------------------
void runAndSaveNucleons(Int_t n,                    
                        const char *sysA="Au",           
                        const char *sysB="Au",           
                        Double_t signn=42,           
                        Double_t mind=0.4,
                        Bool_t verbose=0,
                        const char *fname="glau_auau_nucleons.root");

//---------------------------------------------------------------------------------
class TGlauNucleon : public TNamed
{
   private:
      Double32_t fX;            //Position of nucleon
      Double32_t fY;            //Position of nucleon
      Double32_t fZ;            //Position of nucleon
      Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
      Int_t      fNColl;        //Number of binary collisions

   public:
      TGlauNucleon() : fX(0), fY(0), fZ(0), fInNucleusA(0), fNColl(0) {}
      virtual   ~TGlauNucleon() {}

      void       Collide()            {fNColl++;}
      Int_t      GetNColl()     const {return fNColl;}
      Double_t   GetX()         const {return fX;}
      Double_t   GetY()         const {return fY;}
      Double_t   GetZ()         const {return fZ;}
      Bool_t     IsInNucleusA() const {return fInNucleusA;}
      Bool_t     IsInNucleusB() const {return !fInNucleusA;}
      Bool_t     IsSpectator()  const {return !fNColl;}
      Bool_t     IsWounded()    const {return fNColl;}
      void       Reset()              {fNColl=0;}
      void       SetInNucleusA()      {fInNucleusA=1;}
      void       SetInNucleusB()      {fInNucleusA=0;}
      void       SetXYZ(Double_t x, Double_t y, Double_t z) {fX=x; fY=y; fZ=z;}

      ClassDef(TGlauNucleon,1)
};

//---------------------------------------------------------------------------------
class TGlauNucleus : public TNamed
{
   private:
      Int_t      fN;          //Number of nucleons
      Double_t   fR;          //Parameters of function
      Double_t   fA;          //Parameters of function
      Double_t   fW;          //Parameters of function
      Double_t   fMinDist;    //Minimum separation distance
      Int_t      fF;          //Type of radial distribution
      Int_t      fTrials;     //Store trials needed to complete nucleus
      TF1*       fFunction;   //Probability density function rho(r)
      TObjArray* fNucleons;   //Array of nucleons

      void       Lookup(const char* name);

   public:
      TGlauNucleus(const char* iname="Au", Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0);
      virtual ~TGlauNucleus();

      using      TObject::Draw;
      void       Draw(Double_t xs, Int_t col);
      Int_t      GetN()             const {return fN;}
      Double_t   GetR()             const {return fR;}
      Double_t   GetA()             const {return fA;}
      Double_t   GetW()             const {return fW;}
      TObjArray *GetNucleons()      const {return fNucleons;}
      Int_t      GetTrials()        const {return fTrials;}
      void       SetN(Int_t in)           {fN=in;}
      void       SetR(Double_t ir);
      void       SetA(Double_t ia);
      void       SetW(Double_t iw);
      void       SetMinDist(Double_t min) {fMinDist=min;}
      void       ThrowNucleons(Double_t xshift=0.);

      ClassDef(TGlauNucleus,1)
};

//---------------------------------------------------------------------------------
class TGlauberMC : public TNamed
{
   private:
      TGlauNucleus fANucleus;       //Nucleus A
      TGlauNucleus fBNucleus;       //Nucleus B
      Double_t     fXSect;          //Nucleon-nucleon cross section
      TObjArray*   fNucleonsA;      //Array of nucleons in nucleus A
      TObjArray*   fNucleonsB;      //Array of nucleons in nucleus B
      Int_t        fAN;             //Number of nucleons in nucleus A
      Int_t        fBN;             //Number of nucleons in nucleus B
      TNtuple*     fnt;             //Ntuple for results (created, but not deleted)
      Double_t     fMeanX2;         //<x^2> of wounded nucleons
      Double_t     fMeanY2;         //<y^2> of wounded nucleons
      Double_t     fMeanXY;         //<xy> of wounded nucleons
      Double_t     fMeanXParts;     //<x> of wounded nucleons
      Double_t     fMeanYParts;     //<x> of wounded nucleons
      Double_t     fMeanXSystem;    //<x> of all nucleons
      Double_t     fMeanYSystem;    //<x> of all nucleons  
      Double_t     fMeanX_A;        //<x> of nucleons in nucleus A
      Double_t     fMeanY_A;        //<x> of nucleons in nucleus A
      Double_t     fMeanX_B;        //<x> of nucleons in nucleus B
      Double_t     fMeanY_B;        //<x> of nucleons in nucleus B
      Double_t     fB_MC;           //Impact parameter (b)
      Int_t        fEvents;         //Number of events with at least one collision
      Int_t        fTotalEvents;    //All events within selected impact parameter range
      Double_t     fBMin;           //Minimum impact parameter to be generated
      Double_t     fBMax;           //Maximum impact parameter to be generated
      Int_t        fMaxNpartFound;  //Largest value of Npart obtained
      Int_t        fNpart;          //Number of wounded (participating) nucleons in current event
      Int_t        fNcoll;          //Number of binary collisions in current event
      Double_t     fSx2;            //Variance of x of wounded nucleons
      Double_t     fSy2;            //Variance of y of wounded nucleons
      Double_t     fSxy;            //Covariance of x and y of wounded nucleons

  Double_t     fEcc2PART;       //v1.4 participant eccentricity from Eqn(4) from PRC81, 054905 (2010)
  Double_t     fPsi2PART;       //v1.4 minor eccentricity axis Eqn(5) from PRC81, 054905 (2010)
  Double_t     fEcc3PART;       //v1.4 participant triangularity from Eqn(7) from PRC81, 054905 (2010)
  Double_t     fPsi3PART;       //v1.4 minor triangularity axis Eqn(8) from PRC81, 054905 (2010)

  Double_t     fEcc1PART;       //v1.5 
  Double_t     fPsi1PART;       //v1.5 
  Double_t     fEcc4PART;       //v1.5 
  Double_t     fPsi4PART;       //v1.5 
  Double_t     fEcc5PART;       //v1.5 
  Double_t     fPsi5PART;       //v1.5 
  Double_t     fEcc6PART;       //v1.5 
  Double_t     fPsi6PART;       //v1.5 

  Double_t     fEcc1Cos;       //v1.5x
  Double_t     fEcc1Sin;       //v1.5x
  Double_t     fEcc2Cos;       //v1.5x
  Double_t     fEcc2Sin;       //v1.5x
  Double_t     fEcc3Cos;       //v1.5x
  Double_t     fEcc3Sin;       //v1.5x
  Double_t     fEcc4Cos;       //v1.5x
  Double_t     fEcc4Sin;       //v1.5x
  Double_t     fEcc5Cos;       //v1.5x
  Double_t     fEcc5Sin;       //v1.5x
  Double_t     fEcc6Cos;       //v1.5x
  Double_t     fEcc6Sin;       //v1.5x
  Double_t     fEcc42Cos;       //v1.5x
  Double_t     fEcc42Sin;       //v1.5x
  Double_t     fEcc62Cos;       //v1.5x
  Double_t     fEcc62Sin;       //v1.5x
  Double_t     fEcc63Cos;       //v1.5x
  Double_t     fEcc63Sin;       //v1.5x

  Double_t     fCorr24Cos;       //v1.5x
  Double_t     fCorr24Sin;       //v1.5x

      Bool_t       CalcResults(Double_t bgen);
      Bool_t       CalcEvent(Double_t bgen);

   public:
      TGlauberMC(const char* NA = "Au", const char* NB = "Au", Double_t xsect = 42);
      virtual     ~TGlauberMC() {Reset();}

      void         Draw(Option_t* option);
      Double_t     GetB()               const {return fB_MC;}
      Double_t     GetBMin()            const {return fBMin;}
      Double_t     GetBMax()            const {return fBMax;}
      Int_t        GetNcoll()           const {return fNcoll;}
      Int_t        GetNpart()           const {return fNpart;}
      Int_t        GetNpartFound()      const {return fMaxNpartFound;}
      TNtuple*     GetNtuple()          const {return fnt;}
      TObjArray   *GetNucleons();
      Double_t     GetTotXSect()        const;
      Double_t     GetTotXSectErr()     const;
      Bool_t       NextEvent(Double_t bgen=-1);
      void         Reset()                    {delete fnt; fnt=0; }
      void         Run(Int_t nevents);
      void         SetBmin(Double_t bmin)      {fBMin = bmin;}
      void         SetBmax(Double_t bmax)      {fBMax = bmax;}
      void         SetMinDistance(Double_t d)  {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}

      static 
        void       PrintVersion()             {cout << "TGlauberMC " << Version() << endl;}
      static 
        const char *Version()                 {return "v1.1";}

      ClassDef(TGlauberMC,1)
};

//---------------------------------------------------------------------------------
void runAndSaveNtuple(Int_t n,
                      const char *sysA,
                      const char *sysB,
                      Double_t signn,
                      Double_t mind,
                      const char *fname)
{
   TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn);
   mcg->SetMinDistance(mind);
   mcg->Run(n);
   TNtuple  *nt=mcg->GetNtuple();
   TFile out(fname,"recreate",fname,9);
   if(nt) nt->Write();
   out.Close();
}

//---------------------------------------------------------------------------------
void runAndSaveNucleons(Int_t n,                    
                        const char *sysA,           
                        const char *sysB,           
                        Double_t signn,
                        Double_t mind,
                        Bool_t verbose,
                        const char *fname)
{
   TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn);
   mcg->SetMinDistance(mind);
   TFile *out=0;
   if(fname) 
      out=new TFile(fname,"recreate",fname,9);

   for(Int_t ievent=0;ievent<n;ievent++){

      //get an event with at least one collision
      while(!mcg->NextEvent()) {}

      //access, save and (if wanted) print out nucleons
      TObjArray* nucleons=mcg->GetNucleons();
      if(!nucleons) continue;
      if(out)
         nucleons->Write(Form("nucleonarray%d",ievent),TObject::kSingleKey);

      if(verbose) {
         cout<<endl<<endl<<"EVENT NO: "<<ievent<<endl;
         cout<<"B = "<<mcg->GetB()<<"  Npart = "<<mcg->GetNpart()<<endl<<endl;
         printf("Nucleus\t X\t Y\t Z\tNcoll\n");
         Int_t nNucls=nucleons->GetEntries();
         for(Int_t iNucl=0;iNucl<nNucls;iNucl++) {
            TGlauNucleon *nucl=(TGlauNucleon *)nucleons->At(iNucl);
            Char_t nucleus='A';
            if(nucl->IsInNucleusB()) nucleus='B';
            Double_t x=nucl->GetX();
            Double_t y=nucl->GetY();
            Double_t z=nucl->GetZ();
            Int_t ncoll=nucl->GetNColl();
            printf("   %c\t%2.2f\t%2.2f\t%2.2f\t%3d\n",nucleus,x,y,z,ncoll);
         }
      }
   }
   if(out) delete out;
}

//---------------------------------------------------------------------------------
ClassImp(TGlauNucleus)
//---------------------------------------------------------------------------------

TGlauNucleus::TGlauNucleus(const char* iname, Int_t iN, Double_t iR, Double_t ia, Double_t iw, TF1* ifunc) : 
   fN(iN),fR(iR),fA(ia),fW(iw),fMinDist(-1),
   fF(0),fTrials(0),fFunction(ifunc),
   fNucleons(0)
{
   if (fN==0) {
      cout << "Setting up nucleus " << iname << endl;
      Lookup(iname);
   }
}

TGlauNucleus::~TGlauNucleus()
{
   if (fNucleons) {
      delete fNucleons;
   }
   delete fFunction;
}

void TGlauNucleus::Draw(Double_t xs, Int_t col)
{
   Double_t r = 0.5*sqrt(xs/TMath::Pi()/10.);
   TEllipse e;
   e.SetLineColor(col);
   e.SetFillColor(0);
   e.SetLineWidth(1);

   for (Int_t i = 0;i<fNucleons->GetEntries();++i) {
      TGlauNucleon* gn = (TGlauNucleon*) fNucleons->At(i);
      e.SetLineStyle(1);
      if (gn->IsSpectator()) e.SetLineStyle(3);
      e.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
   }
}

void TGlauNucleus::Lookup(const char* name)
{
   SetName(name);

   if      (TString(name) == "p")    {fN = 1;   fR = 0.6;   fA = 0;      fW =  0;      fF = 0;}
   else if (TString(name) == "d")    {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 1;}
   else if (TString(name) == "dh")   {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 3;}
   else if (TString(name) == "dhh")  {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 4;}
   else if (TString(name) == "O")    {fN = 16;  fR = 2.608; fA = 0.513;  fW = -0.051;  fF = 1;}
   else if (TString(name) == "Si")   {fN = 28;  fR = 3.34;  fA = 0.580;  fW = -0.233;  fF = 1;}
   else if (TString(name) == "S")    {fN = 32;  fR = 2.54;  fA = 2.191;  fW =  0.16;   fF = 2;}
   else if (TString(name) == "Ca")   {fN = 40;  fR = 3.766; fA = 0.586;  fW = -0.161;  fF = 1;}
   else if (TString(name) == "Ni")   {fN = 58;  fR = 4.309; fA = 0.517;  fW = -0.1308; fF = 1;}
   else if (TString(name) == "Cu")   {fN = 63;  fR = 4.2;   fA = 0.596;  fW =  0;      fF = 1;}
   else if (TString(name) == "W")    {fN = 186; fR = 6.58;  fA = 0.480;  fW =  0;      fF = 1;}
   else if (TString(name) == "Au")   {fN = 197; fR = 6.38;  fA = 0.535;  fW =  0;      fF = 1;}
   else if (TString(name) == "Pb")   {fN = 208; fR = 6.62;  fA = 0.546;  fW =  0;      fF = 1;}
   else if (TString(name) == "U")    {fN = 238; fR = 6.81;  fA = 0.6;    fW =  0;      fF = 1;}
   else {
      cout << "Could not find nucleus " << name << endl;
      return;
   }

   switch (fF)
   {
      case 0: // Proton
         fFunction = new TF1("prot","x*x*exp(-x/[0])",0,10);
         fFunction->SetParameter(0,fR);
         break;
      case 1: // 3pF
         fFunction = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,15);
         fFunction->SetParameters(fR,fA,fW);
         break;
      case 2: // 3pG
         fFunction = new TF1("3pg","x*x*(1+[2]*(x/[0])**2)/(1+exp((x**2-[0]**2)/[1]**2))",0,15);
         fFunction->SetParameters(fR,fA,fW);
         break;
      case 3: // Hulthen
         fFunction = new TF1("f3","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,10);
         fFunction->SetParameters(1/4.38,1/.85);
         break;
      case 4: // Hulthen HIJING
         fFunction = new TF1("f4","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,20);
         fFunction->SetParameters(2/4.38,2/.85);
         break;
      default:
         cerr << "Could not find function type " << fF << endl;
         return;
   }
}

void TGlauNucleus::SetR(Double_t ir)
{
   fR = ir;
   switch (fF)
   {
      case 0: // Proton
         fFunction->SetParameter(0,fR);
         break;
      case 1: // 3pF
         fFunction->SetParameter(0,fR);
         break;
      case 2: // 3pG
         fFunction->SetParameter(0,fR);
         break;
   }
}

void TGlauNucleus::SetA(Double_t ia)
{
   fA = ia;
   switch (fF)
   {
      case 0: // Proton
         break;
      case 1: // 3pF
         fFunction->SetParameter(1,fA);
         break;
      case 2: // 3pG
         fFunction->SetParameter(1,fA);
         break;
   }
}

void TGlauNucleus::SetW(Double_t iw)
{
   fW = iw;
   switch (fF)
   {
      case 0: // Proton
         break;
      case 1: // 3pF
         fFunction->SetParameter(2,fW);
         break;
      case 2: // 3pG
         fFunction->SetParameter(2,fW);
         break;
   }
}

void TGlauNucleus::ThrowNucleons(Double_t xshift)
{
   if (fNucleons==0) {
      fNucleons=new TObjArray(fN);
      fNucleons->SetOwner();
      for(Int_t i=0;i<fN;i++) {
	 TGlauNucleon *nucleon=new TGlauNucleon(); 
	 fNucleons->Add(nucleon); 
      }
   } 
   
   fTrials = 0;

   Double_t sumx=0;       
   Double_t sumy=0;       
   Double_t sumz=0;       

   Bool_t hulthen = (TString(GetName())=="dh");
   if (fN==2 && hulthen) { //special treatmeant for Hulten

      Double_t r = fFunction->GetRandom()/2;
      Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
      Double_t ctheta = 2*gRandom->Rndm() - 1 ;
      Double_t stheta = sqrt(1-ctheta*ctheta);
     
      TGlauNucleon *nucleon1=(TGlauNucleon*)(fNucleons->At(0));
      TGlauNucleon *nucleon2=(TGlauNucleon*)(fNucleons->At(1));
      nucleon1->Reset();
      nucleon1->SetXYZ(r * stheta * cos(phi) + xshift,
		       r * stheta * sin(phi),
		       r * ctheta);
      nucleon2->Reset();
      nucleon2->SetXYZ(-nucleon1->GetX() + 2*xshift,
		       -nucleon1->GetY(),
		       -nucleon1->GetZ());
      fTrials = 1;
      return;
   }

   for (Int_t i = 0; i<fN; i++) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      while(1) {
         fTrials++;
         Double_t r = fFunction->GetRandom();
         Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
         Double_t ctheta = 2*gRandom->Rndm() - 1 ;
         Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
         Double_t x = r * stheta * cos(phi) + xshift;
         Double_t y = r * stheta * sin(phi);      
         Double_t z = r * ctheta;      
         nucleon->SetXYZ(x,y,z);
         if(fMinDist<0) break;
         Bool_t test=1;
         for (Int_t j = 0; j<i; j++) {
            TGlauNucleon *other=(TGlauNucleon*)fNucleons->At(j);
            Double_t xo=other->GetX();
            Double_t yo=other->GetY();
            Double_t zo=other->GetZ();
            Double_t dist = TMath::Sqrt((x-xo)*(x-xo)+
                                       (y-yo)*(y-yo)+
                                       (z-zo)*(z-zo));
	       
            if(dist<fMinDist) {
               test=0;
               break;
            }
         }
         if (test) break; //found nucleuon outside of mindist
      }
           
      sumx += nucleon->GetX();
      sumy += nucleon->GetY();
      sumz += nucleon->GetZ();
   }
      
   if(1) { // set the centre-of-mass to be at zero (+xshift)
      sumx = sumx/fN;  
      sumy = sumy/fN;  
      sumz = sumz/fN;  
      for (Int_t i = 0; i<fN; i++) {
         TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
         nucleon->SetXYZ(nucleon->GetX()-sumx-xshift,
                         nucleon->GetY()-sumy,
                         nucleon->GetZ()-sumz);
      }
   }
}

//---------------------------------------------------------------------------------
ClassImp(TGlauberMC)
//---------------------------------------------------------------------------------

TGlauberMC::TGlauberMC(const char* NA, const char* NB, Double_t xsect) :
   fANucleus(NA),fBNucleus(NB),
   fXSect(0),fNucleonsA(0),fNucleonsB(0),fnt(0),
   fMeanX2(0),fMeanY2(0),fMeanXY(0),fMeanXParts(0),
   fMeanYParts(0),fMeanXSystem(0),fMeanYSystem(0),  
   fMeanX_A(0),fMeanY_A(0),fMeanX_B(0),fMeanY_B(0),fB_MC(0),
   fEvents(0),fTotalEvents(0),fBMin(0),fBMax(0),fMaxNpartFound(0),
   fNpart(0),fNcoll(0),fSx2(0),fSy2(0),fSxy(0),
   fEcc2PART(0),fPsi2PART(0),fEcc3PART(0),fPsi3PART(0),
   fEcc1PART(0),fPsi1PART(0),fEcc4PART(0),fPsi4PART(0),
   fEcc5PART(0),fPsi5PART(0),fEcc6PART(0),fPsi6PART(0),
   fEcc1Cos(0),fEcc1Sin(0),fEcc2Cos(0),fEcc2Sin(0),
   fEcc3Cos(0),fEcc3Sin(0),fEcc4Cos(0),fEcc4Sin(0),
   fEcc5Cos(0),fEcc5Sin(0),fEcc6Cos(0),fEcc6Sin(0),
   fEcc42Cos(0),fEcc42Sin(0),fEcc62Cos(0),fEcc62Sin(0),
   fEcc63Cos(0),fEcc63Sin(0),fCorr24Cos(0),fCorr24Sin(0)
{
   fBMin = 0;
   fBMax = 20;
   fXSect = xsect;
   
   TString name(Form("Glauber_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
   TString title(Form("Glauber %s+%s Version",fANucleus.GetName(),fBNucleus.GetName()));
   SetName(name);
   SetTitle(title);
}

Bool_t TGlauberMC::CalcEvent(Double_t bgen)
{
   // prepare event
   fANucleus.ThrowNucleons(-bgen/2.);
   fNucleonsA = fANucleus.GetNucleons();
   fAN = fANucleus.GetN();
   for (Int_t i = 0; i<fAN; i++) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      nucleonA->SetInNucleusA();
   }
   fBNucleus.ThrowNucleons(bgen/2.);
   fNucleonsB = fBNucleus.GetNucleons();
   fBN = fBNucleus.GetN();
   for (Int_t i = 0; i<fBN; i++) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      nucleonB->SetInNucleusB();
   }

   // "ball" diameter = distance at which two balls interact
   Double_t d2 = (Double_t)fXSect/(TMath::Pi()*10); // in fm^2

   // for each of the A nucleons in nucleus B
   for (Int_t i = 0; i<fBN; i++) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      for (Int_t j = 0 ; j < fAN ;j++) {
	 TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(j));
         Double_t dx = nucleonB->GetX()-nucleonA->GetX();
         Double_t dy = nucleonB->GetY()-nucleonA->GetY();
         Double_t dij = dx*dx+dy*dy;
         if (dij < d2) {
            nucleonB->Collide();
            nucleonA->Collide();
         }
      }
  }
   
  return CalcResults(bgen);
}

Bool_t TGlauberMC::CalcResults(Double_t bgen)
{
   // calc results for the given event
   fNpart=0;
   fNcoll=0;
   fMeanX2=0;
   fMeanY2=0;
   fMeanXY=0;
   fMeanXParts=0;
   fMeanYParts=0;
   fMeanXSystem=0;
   fMeanYSystem=0;
   fMeanX_A=0;
   fMeanY_A=0;
   fMeanX_B=0;
   fMeanY_B=0;
  
  double xAgaus[300] = {0};
  double yAgaus[300] = {0};
  double xBgaus[300] = {0};
  double yBgaus[300] = {0};
  bool useGaus = true;

   for (Int_t i = 0; i<fAN; i++) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      Double_t xA=nucleonA->GetX();
      Double_t yA=nucleonA->GetY();

      fMeanXSystem  += xA;
      fMeanYSystem  += yA;
      fMeanX_A  += xA;
      fMeanY_A  += yA;

      if(nucleonA->IsWounded()) {

// using Gauss
        TF2 *fGauss = new TF2("bigaus","bigaus",xA-2,xA+2,yA-2,yA+2); // 2d gauss-like
        fGauss->SetParameters(1,xA,0.4,yA,0.4);
        fGauss->GetRandom2(xAgaus[i],yAgaus[i]);

         fNpart++;
         fMeanXParts  += xAgaus[i];
         fMeanYParts  += yAgaus[i];
         fMeanX2 += xAgaus[i] * xAgaus[i];
         fMeanY2 += yAgaus[i] * yAgaus[i];
         fMeanXY += xAgaus[i] * yAgaus[i];

        delete fGauss;

      }


   }

   for (Int_t i = 0; i<fBN; i++) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      Double_t xB=nucleonB->GetX();
      Double_t yB=nucleonB->GetY();

      fMeanXSystem  += xBgaus[i];
      fMeanYSystem  += yBgaus[i];
      fMeanX_B  += xBgaus[i];
      fMeanY_B  += yBgaus[i];

      if(nucleonB->IsWounded()) {

// using Gauss
        TF2 *fGauss = new TF2("bigaus","bigaus",xB-2,xB+2,yB-2,yB+2); // 2d gauss-like
        fGauss->SetParameters(1,xB,0.4,yB,0.4);
        fGauss->GetRandom2(xBgaus[i],yBgaus[i]);

         fNpart++;
         fMeanXParts  += xBgaus[i];
         fMeanYParts  += yBgaus[i];
         fMeanX2 += xBgaus[i] * xBgaus[i];
         fMeanY2 += yBgaus[i] * yBgaus[i];
         fMeanXY += xBgaus[i] * yBgaus[i];
	 fNcoll += nucleonB->GetNColl();

        delete fGauss;

      }

   }

   if (fNpart>0) {
      fMeanXParts /= fNpart;
      fMeanYParts /= fNpart;
      fMeanX2 /= fNpart;
      fMeanY2 /= fNpart;
      fMeanXY /= fNpart;
   } else {
      fMeanXParts = 0;
      fMeanYParts = 0;
      fMeanX2 = 0;
      fMeanY2 = 0;
      fMeanXY = 0;
   }
   
   if(fAN+fBN>0) {
      fMeanXSystem /= (fAN + fBN);
      fMeanYSystem /= (fAN + fBN);
   } else {
      fMeanXSystem = 0;
      fMeanYSystem = 0;
   }
   if(fAN>0) {
      fMeanX_A /= fAN;
      fMeanY_A /= fAN;
   } else {
      fMeanX_A = 0;
      fMeanY_A = 0;
   }

   if(fBN>0) {
      fMeanX_B /= fBN;
      fMeanY_B /= fBN;
   } else {
      fMeanX_B = 0;
      fMeanY_B = 0;
   }

   fSx2=fMeanX2-(fMeanXParts*fMeanXParts);
   fSy2=fMeanY2-(fMeanYParts*fMeanYParts);
   fSxy=fMeanXY-fMeanXParts*fMeanYParts;

   fB_MC = bgen;

   //----------------------------v1.4------------------------------
   //v1.4 Put here calculations for triangularity


   if (fNpart>0) {

     Double_t rsquared=0;
     Double_t phipart=0;
     Double_t avgr1sqcos1phi=0;
     Double_t avgr1sqsin1phi=0;
     Double_t avgr2sqcos2phi=0;
     Double_t avgr2sqsin2phi=0;
     Double_t avgr3sqcos3phi=0;
     Double_t avgr3sqsin3phi=0;
     Double_t avgr4sqcos4phi=0;
     Double_t avgr4sqsin4phi=0;
     Double_t avgr5sqcos5phi=0;
     Double_t avgr5sqsin5phi=0;
     Double_t avgr6sqcos6phi=0;
     Double_t avgr6sqsin6phi=0;
     Double_t avgrsquared=0;
     Double_t avgr1st=0;
     Double_t avgr3rd=0;
     Double_t avgr4th=0;
     Double_t avgr5th=0;
     Double_t avgr6th=0;          

     for (Int_t i = 0; i<fAN; i++) {
       TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
       if(nucleonA->IsWounded()) {
	Double_t xAshift = xAgaus[i]-fMeanXParts;
	Double_t yAshift = yAgaus[i]-fMeanYParts;
	 ///Double_t xAshift = nucleonA->GetX()-fMeanXParts;
	 ///Double_t yAshift = nucleonA->GetY()-fMeanYParts;
	 //Double_t xAshift = nucleonA->GetX();
	 //Double_t yAshift = nucleonA->GetY();
	 rsquared = xAshift*xAshift + yAshift*yAshift;
	 phipart = TMath::ATan2(yAshift,xAshift);
	 avgr1sqcos1phi += pow(rsquared,1.0/2)*TMath::Cos(phipart);
	 avgr1sqsin1phi += pow(rsquared,1.0/2)*TMath::Sin(phipart);
	 avgr2sqcos2phi += rsquared*TMath::Cos(2*phipart);
	 avgr2sqsin2phi += rsquared*TMath::Sin(2*phipart);
	 avgr3sqcos3phi += pow(rsquared,3.0/2)*TMath::Cos(3*phipart);
	 avgr3sqsin3phi += pow(rsquared,3.0/2)*TMath::Sin(3*phipart);
	 avgr4sqcos4phi += pow(rsquared,4.0/2)*TMath::Cos(4*phipart);
	 avgr4sqsin4phi += pow(rsquared,4.0/2)*TMath::Sin(4*phipart);
	 avgr5sqcos5phi += pow(rsquared,5.0/2)*TMath::Cos(5*phipart);
	 avgr5sqsin5phi += pow(rsquared,5.0/2)*TMath::Sin(5*phipart);
	 avgr6sqcos6phi += pow(rsquared,6.0/2)*TMath::Cos(6*phipart);
	 avgr6sqsin6phi += pow(rsquared,6.0/2)*TMath::Sin(6*phipart);
	 avgrsquared += rsquared;
         avgr1st += pow(rsquared,1.0/2);
         avgr3rd += pow(rsquared,3.0/2);
         avgr4th += pow(rsquared,4.0/2);
         avgr5th += pow(rsquared,5.0/2);
         avgr6th += pow(rsquared,6.0/2);
       }
     }

     for (Int_t i = 0; i<fBN; i++) {
       TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
       if(nucleonB->IsWounded()) {
	Double_t xBshift = xBgaus[i]-fMeanXParts;
	Double_t yBshift = yBgaus[i]-fMeanYParts;
	 ///Double_t xBshift = nucleonB->GetX()-fMeanXParts;
	 ///Double_t yBshift = nucleonB->GetY()-fMeanYParts;
	 rsquared = xBshift*xBshift + yBshift*yBshift;
	 phipart = TMath::ATan2(yBshift,xBshift);
	 avgr1sqcos1phi += pow(rsquared,1.0/2)*TMath::Cos(phipart);
	 avgr1sqsin1phi += pow(rsquared,1.0/2)*TMath::Sin(phipart);
	 avgr2sqcos2phi += rsquared*TMath::Cos(2*phipart);
	 avgr2sqsin2phi += rsquared*TMath::Sin(2*phipart);
	 avgr3sqcos3phi += pow(rsquared,3.0/2)*TMath::Cos(3*phipart);
	 avgr3sqsin3phi += pow(rsquared,3.0/2)*TMath::Sin(3*phipart);
	 avgr4sqcos4phi += pow(rsquared,4.0/2)*TMath::Cos(4*phipart);
	 avgr4sqsin4phi += pow(rsquared,4.0/2)*TMath::Sin(4*phipart);
	 avgr5sqcos5phi += pow(rsquared,5.0/2)*TMath::Cos(5*phipart);
	 avgr5sqsin5phi += pow(rsquared,5.0/2)*TMath::Sin(5*phipart);
	 avgr6sqcos6phi += pow(rsquared,6.0/2)*TMath::Cos(6*phipart);
	 avgr6sqsin6phi += pow(rsquared,6.0/2)*TMath::Sin(6*phipart);
	 avgrsquared += rsquared;
         avgr1st += pow(rsquared,1.0/2);
         avgr3rd += pow(rsquared,3.0/2);
         avgr4th += pow(rsquared,4.0/2);
         avgr5th += pow(rsquared,5.0/2);
         avgr6th += pow(rsquared,6.0/2);
       }
     }

     avgr1sqcos1phi /= fNpart;
     avgr1sqsin1phi /= fNpart;
     avgr2sqcos2phi /= fNpart;
     avgr2sqsin2phi /= fNpart;
     avgr3sqcos3phi /= fNpart;
     avgr3sqsin3phi /= fNpart;
     avgr4sqcos4phi /= fNpart;
     avgr4sqsin4phi /= fNpart;
     avgr5sqcos5phi /= fNpart;
     avgr5sqsin5phi /= fNpart;
     avgr6sqcos6phi /= fNpart;
     avgr6sqsin6phi /= fNpart;
     avgrsquared /= fNpart;
     avgr1st /= fNpart;
     avgr3rd /= fNpart;
     avgr4th /= fNpart;
     avgr5th /= fNpart;
     avgr6th /= fNpart;

     fEcc1PART=TMath::Sqrt(avgr1sqcos1phi*avgr1sqcos1phi+avgr1sqsin1phi*avgr1sqsin1phi)/avgr1st;
     fPsi1PART=(TMath::ATan2(avgr1sqsin1phi,avgr1sqcos1phi)+TMath::Pi())/1.;
     fEcc2PART=TMath::Sqrt(avgr2sqcos2phi*avgr2sqcos2phi+avgr2sqsin2phi*avgr2sqsin2phi)/avgrsquared;
     fPsi2PART=(TMath::ATan2(avgr2sqsin2phi,avgr2sqcos2phi)+TMath::Pi())/2.;
     fEcc3PART=TMath::Sqrt(avgr3sqcos3phi*avgr3sqcos3phi+avgr3sqsin3phi*avgr3sqsin3phi)/avgr3rd;
     fPsi3PART=(TMath::ATan2(avgr3sqsin3phi,avgr3sqcos3phi)+TMath::Pi())/3.;
     fEcc4PART=TMath::Sqrt(avgr4sqcos4phi*avgr4sqcos4phi+avgr4sqsin4phi*avgr4sqsin4phi)/avgr4th;
     fPsi4PART=(TMath::ATan2(avgr4sqsin4phi,avgr4sqcos4phi)+TMath::Pi())/4.;
     fEcc5PART=TMath::Sqrt(avgr5sqcos5phi*avgr5sqcos5phi+avgr5sqsin5phi*avgr5sqsin5phi)/avgr5th;
     fPsi5PART=(TMath::ATan2(avgr5sqsin5phi,avgr5sqcos5phi)+TMath::Pi())/5.;
     fEcc6PART=TMath::Sqrt(avgr6sqcos6phi*avgr6sqcos6phi+avgr6sqsin6phi*avgr6sqsin6phi)/avgr6th;
     fPsi6PART=(TMath::ATan2(avgr6sqsin6phi,avgr6sqcos6phi)+TMath::Pi())/6.;

   } else {

     fEcc1PART=-1000;
     fPsi1PART=-1000;
     fEcc2PART=-1000;
     fPsi2PART=-1000;
     fEcc3PART=-1000;
     fPsi3PART=-1000;
     fEcc4PART=-1000;
     fPsi4PART=-1000;
     fEcc5PART=-1000;
     fPsi5PART=-1000;
     fEcc6PART=-1000;
     fPsi6PART=-1000;

   }


  //:changed: add fEccijCos;
   if (fNpart>0) {

     Double_t ecc1c=0;
     Double_t ecc1s=0;
     Double_t ecc2c=0;
     Double_t ecc2s=0;
     Double_t ecc3c=0;
     Double_t ecc3s=0;
     Double_t ecc4c=0;
     Double_t ecc4s=0;
     Double_t ecc5c=0;
     Double_t ecc5s=0;
     Double_t ecc6c=0;
     Double_t ecc6s=0;
     Double_t ecc42c=0;
     Double_t ecc42s=0;
     Double_t ecc62c=0;
     Double_t ecc62s=0;
     Double_t ecc63c=0;
     Double_t ecc63s=0;
     Double_t rsquared=0;
     Double_t phipart=0;
     Double_t avgrsquared=0;
     Double_t avgr1st=0;
     Double_t avgr3rd=0;
     Double_t avgr4th=0;
     Double_t avgr5th=0;
     Double_t avgr6th=0;

     for (Int_t i = 0; i<fAN; i++) {
       TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
       if(nucleonA->IsWounded()) {
         //Double_t xAshift = nucleonA->GetX();//-fMeanXParts;   // phipart should not shift;
         //Double_t yAshift = nucleonA->GetY();//-fMeanYParts;
         ///Double_t xAshift = nucleonA->GetX()-fMeanXParts;   // phipart should shift;
         ///Double_t yAshift = nucleonA->GetY()-fMeanYParts;
        Double_t xAshift = xAgaus[i]-fMeanXParts;   // phipart should shift;
        Double_t yAshift = yAgaus[i]-fMeanYParts;
         rsquared = xAshift*xAshift + yAshift*yAshift;
         phipart = TMath::ATan2(yAshift,xAshift);
         avgrsquared += rsquared;
         avgr1st += pow(rsquared,1.0/2);
         avgr3rd += pow(rsquared,3.0/2);
         avgr4th += pow(rsquared,4.0/2);
         avgr5th += pow(rsquared,5.0/2);
         avgr6th += pow(rsquared,6.0/2);
         ecc1c += pow(rsquared,1.0/2)*TMath::Cos(1.0*(phipart-fPsi1PART+TMath::Pi()));
         ecc1s += pow(rsquared,1.0/2)*TMath::Sin(1.0*(phipart-fPsi1PART+TMath::Pi()));
         ecc2c += rsquared*TMath::Cos(2.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc2s += rsquared*TMath::Sin(2.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc3c += pow(rsquared,3.0/2)*TMath::Cos(3.0*(phipart-fPsi3PART+TMath::Pi()/3));
         ecc3s += pow(rsquared,3.0/2)*TMath::Sin(3.0*(phipart-fPsi3PART+TMath::Pi()/3));
         ecc4c += pow(rsquared,4.0/2)*TMath::Cos(4.0*(phipart-fPsi4PART+TMath::Pi()/4));
         ecc4s += pow(rsquared,4.0/2)*TMath::Sin(4.0*(phipart-fPsi4PART+TMath::Pi()/4));
         ecc5c += pow(rsquared,5.0/2)*TMath::Cos(5.0*(phipart-fPsi5PART+TMath::Pi()/5));
         ecc5s += pow(rsquared,5.0/2)*TMath::Sin(5.0*(phipart-fPsi5PART+TMath::Pi()/5));
         ecc6c += pow(rsquared,6.0/2)*TMath::Cos(6.0*(phipart-fPsi6PART+TMath::Pi()/6));
         ecc6s += pow(rsquared,6.0/2)*TMath::Sin(6.0*(phipart-fPsi6PART+TMath::Pi()/6));
         ecc42c += pow(rsquared,4.0/2)*TMath::Cos(4.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc42s += pow(rsquared,4.0/2)*TMath::Sin(4.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc62c += pow(rsquared,6.0/2)*TMath::Cos(6.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc62s += pow(rsquared,6.0/2)*TMath::Sin(6.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc63c += pow(rsquared,6.0/2)*TMath::Cos(6.0*(phipart-fPsi3PART+TMath::Pi()/3));
         ecc63s += pow(rsquared,6.0/2)*TMath::Sin(6.0*(phipart-fPsi3PART+TMath::Pi()/3));
       }
     }

     for (Int_t i = 0; i<fBN; i++) {
       TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
       if(nucleonB->IsWounded()) {
         ///Double_t xBshift = nucleonB->GetX()-fMeanXParts;
         ///Double_t yBshift = nucleonB->GetY()-fMeanYParts;
        Double_t xBshift = xBgaus[i]-fMeanXParts;
        Double_t yBshift = yBgaus[i]-fMeanYParts;
         rsquared = xBshift*xBshift + yBshift*yBshift;
         phipart = TMath::ATan2(yBshift,xBshift);
         avgrsquared += rsquared;
         avgr1st += pow(rsquared,1.0/2);
         avgr3rd += pow(rsquared,3.0/2);
         avgr4th += pow(rsquared,4.0/2);
         avgr5th += pow(rsquared,5.0/2);
         avgr6th += pow(rsquared,6.0/2);
         ecc1c += pow(rsquared,1.0/2)*TMath::Cos(1.0*(phipart-fPsi1PART+TMath::Pi()));
         ecc1s += pow(rsquared,1.0/2)*TMath::Sin(1.0*(phipart-fPsi1PART+TMath::Pi()));
         ecc2c += rsquared*TMath::Cos(2.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc2s += rsquared*TMath::Sin(2.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc3c += pow(rsquared,3.0/2)*TMath::Cos(3.0*(phipart-fPsi3PART+TMath::Pi()/3));
         ecc3s += pow(rsquared,3.0/2)*TMath::Sin(3.0*(phipart-fPsi3PART+TMath::Pi()/3));
         ecc4c += pow(rsquared,4.0/2)*TMath::Cos(4.0*(phipart-fPsi4PART+TMath::Pi()/4));
         ecc4s += pow(rsquared,4.0/2)*TMath::Sin(4.0*(phipart-fPsi4PART+TMath::Pi()/4));
         ecc5c += pow(rsquared,5.0/2)*TMath::Cos(5.0*(phipart-fPsi5PART+TMath::Pi()/5));
         ecc5s += pow(rsquared,5.0/2)*TMath::Sin(5.0*(phipart-fPsi5PART+TMath::Pi()/5));
         ecc6c += pow(rsquared,6.0/2)*TMath::Cos(6.0*(phipart-fPsi6PART+TMath::Pi()/6));
         ecc6s += pow(rsquared,6.0/2)*TMath::Sin(6.0*(phipart-fPsi6PART+TMath::Pi()/6));
         ecc42c += pow(rsquared,4.0/2)*TMath::Cos(4.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc42s += pow(rsquared,4.0/2)*TMath::Sin(4.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc62c += pow(rsquared,6.0/2)*TMath::Cos(6.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc62s += pow(rsquared,6.0/2)*TMath::Sin(6.0*(phipart-fPsi2PART+TMath::Pi()/2));
         ecc63c += pow(rsquared,6.0/2)*TMath::Cos(6.0*(phipart-fPsi3PART+TMath::Pi()/3));
         ecc63s += pow(rsquared,6.0/2)*TMath::Sin(6.0*(phipart-fPsi3PART+TMath::Pi()/3));
       }
     }

     avgrsquared /= fNpart;
     avgr1st /= fNpart;
     avgr3rd /= fNpart;
     avgr4th /= fNpart;
     avgr5th /= fNpart;
     avgr6th /= fNpart;
     fEcc1Cos = ecc1c/fNpart/avgr1st;
     fEcc1Sin = ecc1s/fNpart/avgr1st;
     fEcc2Cos = ecc2c/fNpart/avgrsquared;
     fEcc2Sin = ecc2s/fNpart/avgrsquared;
     fEcc3Cos = ecc3c/fNpart/avgr3rd;
     fEcc3Sin = ecc3s/fNpart/avgr3rd;
     fEcc4Cos = ecc4c/fNpart/avgr4th;
     fEcc4Sin = ecc4s/fNpart/avgr4th;
     fEcc5Cos = ecc5c/fNpart/avgr5th;
     fEcc5Sin = ecc5s/fNpart/avgr5th;
     fEcc6Cos = ecc6c/fNpart/avgr6th;
     fEcc6Sin = ecc6s/fNpart/avgr6th;
     fEcc42Cos = ecc42c/fNpart/avgr4th;
     fEcc42Sin = ecc42s/fNpart/avgr4th;
     fEcc62Cos = ecc62c/fNpart/avgr6th;
     fEcc62Sin = ecc62s/fNpart/avgr6th;
     fEcc63Cos = ecc63c/fNpart/avgr6th;
     fEcc63Sin = ecc63s/fNpart/avgr6th;

     fCorr24Cos = TMath::Cos(4.0*(fPsi2PART-fPsi4PART-TMath::Pi()/4)); //fPsi2PART is minor axis;
     fCorr24Sin = TMath::Sin(4.0*(fPsi2PART-fPsi4PART-TMath::Pi()/4)); //fPsi2PART-pi/2 is major axis;

   } else {

     fEcc1Cos=-1000;
     fEcc1Sin=-1000;
     fEcc2Cos=-1000;
     fEcc2Sin=-1000;
     fEcc3Cos=-1000;
     fEcc3Sin=-1000;
     fEcc4Cos=-1000;
     fEcc4Sin=-1000;
     fEcc5Cos=-1000;
     fEcc5Sin=-1000;
     fEcc6Cos=-1000;
     fEcc6Sin=-1000;
     fEcc42Cos=-1000;
     fEcc42Sin=-1000;
     fEcc62Cos=-1000;
     fEcc62Sin=-1000;
     fEcc63Cos=-1000;
     fEcc63Sin=-1000;

     fCorr24Cos=-1000;
     fCorr24Sin=-1000;

   }
   //----------------------------v1.4------------------------------

   fTotalEvents++;
   if (fNpart>0) fEvents++;

   if (fNpart==0) return kFALSE;
   if (fNpart > fMaxNpartFound) fMaxNpartFound = fNpart;

   return kTRUE;
}

void TGlauberMC::Draw(Option_t* /*option*/)
{
   fANucleus.Draw(fXSect, 2);
   fBNucleus.Draw(fXSect, 4);

   TEllipse e;
   e.SetFillColor(0);
   e.SetLineColor(1);
   e.SetLineStyle(2);
   e.SetLineWidth(1);
   e.DrawEllipse(GetB()/2,0,fBNucleus.GetR(),fBNucleus.GetR(),0,360,0);
   e.DrawEllipse(-GetB()/2,0,fANucleus.GetR(),fANucleus.GetR(),0,360,0);
}

Double_t TGlauberMC::GetTotXSect() const
{
   return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBMax*fBMax/100;
}

Double_t TGlauberMC::GetTotXSectErr() const
{
   return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) * 
      TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

TObjArray *TGlauberMC::GetNucleons() 
{
   if(!fNucleonsA || !fNucleonsB) return 0;
   fNucleonsA->SetOwner(0);
   fNucleonsB->SetOwner(0);
   TObjArray *allnucleons=new TObjArray(fAN+fBN);
   allnucleons->SetOwner();
   for (Int_t i = 0; i<fAN; i++) {
      allnucleons->Add(fNucleonsA->At(i));
   }
   for (Int_t i = 0; i<fBN; i++) {
      allnucleons->Add(fNucleonsB->At(i));
   }
   return allnucleons;
}

Bool_t TGlauberMC::NextEvent(Double_t bgen)
{
   if(bgen<0) 
      bgen = TMath::Sqrt((fBMax*fBMax-fBMin*fBMin)*gRandom->Rndm()+fBMin*fBMin);

   return CalcEvent(bgen);
}

void TGlauberMC::Run(Int_t nevents)
{
   TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
   TString title(Form("%s + %s (x-sect = %d mb)",fANucleus.GetName(),fBNucleus.GetName(),(Int_t) fXSect));
   if (fnt == 0) {
      fnt = new TNtuple(name,title,
			//			"Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB");
			//			"Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB:EccRP:EccPART"); //v1.2
			//			"Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB:EccRP:EccPART:SPART:SPART4:S12RP_Npart"); //v1.3
			//			"Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB:EccRP:EccPART:SPART:SPART4:S12RP_Npart:Ecc2PART:Psi2PART:Ecc3PART:Psi3PART"); //v1.4
			//"Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB:EccRP:EccPART:SPART:SPART4:S12RP_Npart:Ecc2PART:Psi2PART:Ecc3PART:Psi3PART:Ecc1PART:Psi1PART:Ecc4PART:Psi4PART:Ecc5PART:Psi5PART:Ecc6PART:Psi6PART"); //v1.5
			"Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB:EccRP:EccPART:SPART:SPART4:S12RP_Npart:Ecc2PART:Psi2PART:Ecc3PART:Psi3PART:Ecc1PART:Psi1PART:Ecc4PART:Psi4PART:Ecc5PART:Psi5PART:Ecc6PART:Psi6PART:Ecc1Cos:Ecc1Sin:Ecc2Cos:Ecc2Sin:Ecc3Cos:Ecc3Sin:Ecc4Cos:Ecc4Sin:Ecc5Cos:Ecc5Sin:Ecc6Cos:Ecc6Sin:Ecc42Cos:Ecc42Sin:Ecc62Cos:Ecc62Sin:Ecc63Cos:Ecc63Sin:Corr24Cos:Corr24Sin:XSect"); //v1.5x
      fnt->SetDirectory(0);
   }

   for (int i = 0;i<nevents;i++) {

      while(!NextEvent()) {}

      //      Float_t v[17];
      //      Float_t v[19];   //v1.2
      //      Float_t v[22];   //v1.3
      //      Float_t v[26];   //v1.4
      //Float_t v[34];   //v1.5
      Float_t v[55];   //v1.5x
      v[0]  = GetNpart();
      v[1]  = GetNcoll();
      v[2]  = fB_MC;
      v[3]  = fMeanXParts;
      v[4]  = fMeanYParts;
      v[5]  = fMeanX2;
      v[6]  = fMeanY2;
      v[7]  = fMeanXY;
      v[8]  = fSx2;
      v[9]  = fSy2;
      v[10] = fSxy;
      v[11] = fMeanXSystem;
      v[12] = fMeanYSystem;
      v[13] = fMeanX_A;
      v[14] = fMeanY_A;
      v[15] = fMeanX_B;
      v[16] = fMeanY_B;
      v[17] = (fSy2-fSx2)/(fSy2+fSx2);  //v1.2  Eqn A3 of PRC77, 014906 (2008)
      v[18] = TMath::Sqrt((fSy2-fSx2)*(fSy2-fSx2)+4.*fSxy*fSxy)/(fSy2+fSx2);  //v1.2  Eqn A6 of PRC77, 014906 (2008)
      v[19] = fSx2*fSy2-fSxy; //v1.3  Eqn A8 of PRC77, 014906 (2008)
      if(v[19]<0) {
	v[19]=0.0;
      }
      else {
	v[19] = TMath::Pi()*TMath::Sqrt(v[19]); //v1.3  Eqn A8 of PRC77, 014906 (2008)
      }
      v[20] = 4.0*v[19]; //v1.3  Four times Eqn A8 of PRC77, 014906 (2008)
      v[21] = 4.0*TMath::Pi()*TMath::Sqrt(fMeanX2)*TMath::Sqrt(fMeanY2); //v1.3  Eur.Phys.J.C66:173,2010 but with participants and relative to original Reaction Plane
      v[22] = fEcc2PART; //v1.4  PRC81, 054905 (2010) Equation 4
      v[23] = fPsi2PART; //v1.4  PRC81, 054905 (2010) Equation 5
      v[24] = fEcc3PART; //v1.4  PRC81, 054905 (2010) Equation 7
      v[25] = fPsi3PART; //v1.4  PRC81, 054905 (2010) Equation 9
      v[26] = fEcc1PART; //v1.5
      v[27] = fPsi1PART; //v1.5
      v[28] = fEcc4PART; //v1.5
      v[29] = fPsi4PART; //v1.5
      v[30] = fEcc5PART; //v1.5
      v[31] = fPsi5PART; //v1.5
      v[32] = fEcc6PART; //v1.5
      v[33] = fPsi6PART; //v1.5
      v[34] = fEcc1Cos; //v1.5x
      v[35] = fEcc1Sin; //v1.5x
      v[36] = fEcc2Cos; //v1.5x
      v[37] = fEcc2Sin; //v1.5x
      v[38] = fEcc3Cos; //v1.5x
      v[39] = fEcc3Sin; //v1.5x
      v[40] = fEcc4Cos; //v1.5x
      v[41] = fEcc4Sin; //v1.5x
      v[42] = fEcc5Cos; //v1.5x
      v[43] = fEcc5Sin; //v1.5x
      v[44] = fEcc6Cos; //v1.5x
      v[45] = fEcc6Sin; //v1.5x
      v[46] = fEcc42Cos; //v1.5x
      v[47] = fEcc42Sin; //v1.5x
      v[48] = fEcc62Cos; //v1.5x
      v[49] = fEcc62Sin; //v1.5x
      v[50] = fEcc63Cos; //v1.5x
      v[51] = fEcc63Sin; //v1.5x
      v[52] = fCorr24Cos; //v1.5x
      v[53] = fCorr24Sin; //v1.5x
      v[54] = GetTotXSect(); //v1.5x


      fnt->Fill(v);

      if (!(i%50)) cout << "Event # " << i << " x-sect = " << GetTotXSect() << " +- " << GetTotXSectErr() << " b        \r" << flush;  
   }
   cout << endl << "Done!" << endl;
}
#endif

