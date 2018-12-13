/*
 $Id: runglauber.C 173 2018-06-07 19:47:47Z loizides $
 -------------------------------------------------------------------------------------
 Latest documentation: https://arxiv.org/abs/1710.07098
 -------------------------------------------------------------------------------------
 To run the code, you need to have the ROOT (http://root.cern.ch/drupal/)
 environment. On the root prompt, then enter
 root [0] gSystem->Load("libMathMore")
 root [1] .L runglauber_X.Y.C+
 (where X.Y denotes the version number).
 If you do not have libMathMore comment out "#define HAVE_MATHMORE" below.
 See the documentation for more information.
 -------------------------------------------------------------------------------------
 v3.1:
  Fixes related to spherical nuclei, as well as consistent set of reweighted profiles 
  for Cu, Au and Xe, see https://arxiv.org/abs/1710.07098v2
 -------------------------------------------------------------------------------------
 v3.0:
  Major update to include separate profile for protons and neutrons, placement of nucleon 
  dof on lattice, as well as reweighted profiles for recentering, 
  see https://arxiv.org/abs/1710.07098v1
 -------------------------------------------------------------------------------------
 v2.6:
  Includes runAndCalcDens macro, as well as definition for Al, and fixes beta4 for Si2,
  see https://arxiv.org/abs/1408.2549v8
 -------------------------------------------------------------------------------------
 v2.5:
  Include core/corona determination in Npart, and if requested for area from mc and eccentricity,
  as well as various Xe parameterizations including deformation,
  see https://arxiv.org/abs/1408.2549v7
 -------------------------------------------------------------------------------------
 v2.4: 
  Minor update to include Xenon and fix of the TGlauberMC::Draw function, 
  see https://arxiv.org/abs/1408.2549v4
 -------------------------------------------------------------------------------------
 v2.3: 
  Small bugfixes, see https://arxiv.org/abs/1408.2549v3
 -------------------------------------------------------------------------------------
 v2.2:
  Minor update to provide higher harmonic eccentricities up to n=5, and the average
  nucleon--nucleon impact parameter (bNN) in tree output. 
 -------------------------------------------------------------------------------------
 v2.1: 
  Minor update to include more proton pdfs, see https://arxiv.org/abs/1408.2549v2
 -------------------------------------------------------------------------------------
 v2.0: 
  First major update with inclusion of Tritium, Helium-3, and Uranium, as well as the 
  treatment of deformed nuclei and Glauber-Gribov fluctuations of the proton in p+A 
  collisions, see https://arxiv.org/abs/1408.2549v1
 -------------------------------------------------------------------------------------
 v1.1: 
  First public release of the PHOBOS MC Glauber, see https://arxiv.org/abs/0805.4411
 -------------------------------------------------------------------------------------

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#define HAVE_MATHMORE

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TBits.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TNamed.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TRotation.h>
#include <TString.h>
#include <TSystem.h>
#include <TVector3.h>
#ifdef HAVE_MATHMORE
 #include <Math/SpecFuncMathMore.h>
#endif
using namespace std;
#endif

#ifndef _runglauber_
#if !defined(__CINT__) || defined(__MAKECINT__)
#define _runglauber_ 3
#endif

//---------------------------------------------------------------------------------
TF1 *getNNProf(Double_t snn=67.6, Double_t omega=0.4, Double_t G=1);

//---------------------------------------------------------------------------------
void runAndSaveNtuple(const Int_t n,
                      const char *sysA        = "Pb",
                      const char *sysB        = "Pb",
                      const Double_t signn    = 67.6,
                      const Double_t sigwidth = -1,
                      const Double_t mind     = 0.4,
		      const Double_t omega    = -1,
                      const Double_t noded    = -1,
                      const char *fname       = 0);

//---------------------------------------------------------------------------------
void runAndSaveNucleons(const Int_t n,                    
                        const char *sysA        = "Pb",           
                        const char *sysB        = "Pb",           
                        const Double_t signn    = 67.6,           
                        const Double_t sigwidth = -1,
                        const Double_t mind     = 0.4,
                        const Bool_t verbose    = 0,
                        const char *fname       = 0);

//---------------------------------------------------------------------------------
void runAndSmearNtuple(const Int_t n,
                       const Double_t sigs  = 0.4,
                       const char *sysA     = "p",
                       const char *sysB     = "Pb",
                       const Double_t signn = 67.6,
                       const Double_t mind  = 0.4,
                       const char *fname    = 0);

//---------------------------------------------------------------------------------
void runAndCalcDens(const Int_t n,
		    const Double_t alpha = 0.1,
		    const char *sysA     = "Pb",
		    const char *sysB     = "Pb",
		    const Double_t signn = 67.6,
		    const Double_t mind  = 0.4,
		    const char *fname    = "glau_dens_hists.root");

//---------------------------------------------------------------------------------
class TGlauNucleon : public TObject
{
  protected:
    Double32_t fX;            //Position of nucleon
    Double32_t fY;            //Position of nucleon
    Double32_t fZ;            //Position of nucleon
    Int_t      fType;         //0 = neutron, 1 = proton
    Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
    Int_t      fNColl;        //Number of binary collisions
    Double32_t fEn;           //Energy
  public:
    TGlauNucleon() : fX(0), fY(0), fZ(0), fInNucleusA(0), fNColl(0), fEn(0) {}
    virtual   ~TGlauNucleon() {}
    void       Collide()                                  {++fNColl;}
    Double_t   Get2CWeight(Double_t x) const              {return 2.*(0.5*(1-x)+0.5*x*fNColl);}
    Double_t   GetEnergy()             const              {return fEn;}
    Int_t      GetNColl()              const              {return fNColl;}
    Int_t      GetType()               const              {return fType;}
    Double_t   GetX()                  const              {return fX;}
    Double_t   GetY()                  const              {return fY;}
    Double_t   GetZ()                  const              {return fZ;}
    Bool_t     IsNeutron()             const              {return (fType==0);}
    Bool_t     IsInNucleusA()          const              {return fInNucleusA;}
    Bool_t     IsInNucleusB()          const              {return !fInNucleusA;}
    Bool_t     IsProton()              const              {return (fType==1);}
    Bool_t     IsSpectator()           const              {return !fNColl;}
    Bool_t     IsWounded()             const              {return fNColl>0;}
    void       Reset()                                    {fNColl=0;}
    void       RotateXYZ(Double_t phi, Double_t theta);
    void       RotateXYZ_3D(Double_t psiX, Double_t psiY, Double_t psiZ);
    void       SetEnergy(Double_t en)                     {fEn = en;}
    void       SetInNucleusA()                            {fInNucleusA=1;}
    void       SetInNucleusB()                            {fInNucleusA=0;}
    void       SetNColl(Int_t nc)                         {fNColl = nc;}
    void       SetType(Bool_t b)                          {fType = b;}
    void       SetXYZ(Double_t x, Double_t y, Double_t z) {fX=x; fY=y; fZ=z;}
    ClassDef(TGlauNucleon,4) // TGlauNucleon class
};

//---------------------------------------------------------------------------------
ClassImp(TGlauNucleon)
  //---------------------------------------------------------------------------------
void TGlauNucleon::RotateXYZ(Double_t phi, Double_t theta)
{
  TVector3 v(fX,fY,fZ);
  TVector3 vr;
  vr.SetMagThetaPhi(1,theta,phi);
  v.RotateUz(vr);
  fX = v.X();
  fY = v.Y();
  fZ = v.Z();
}

void TGlauNucleon::RotateXYZ_3D(Double_t psiX, Double_t psiY, Double_t psiZ)
{
  TVector3 v(fX,fY,fZ);
  v.RotateX(psiX);
  v.RotateY(psiY);
  v.RotateZ(psiZ);
  fX = v.X();
  fY = v.Y();
  fZ = v.Z();
}

//---------------------------------------------------------------------------------
class TGlauNucleus : public TNamed
{
  private:
    Int_t      fN;                   //Number of nucleons
    Int_t      fZ;                   //Number of protons
    Double_t   fR;                   //Parameters of function
    Double_t   fA;                   //Parameters of function (fA+fZ=fN)
    Double_t   fW;                   //Parameters of function
    Double_t   fR2;                  //Parameters of function (for p and n separately)
    Double_t   fA2;                  //Parameters of function (for p and n separately)
    Double_t   fW2;                  //Parameters of function (for p and n separately)
    Double_t   fBeta2;               //Beta2 (deformed nuclei) 
    Double_t   fBeta4;               //Beta4 (deformed nuclei) 
    Double_t   fMinDist;             //Minimum separation distance
    Double_t   fNodeDist;            //Average node distance (set to <=0 if you do not want the "crystal lattice")
    Double_t   fSmearing;            //Node smearing (relevant if fNodeDist>0)
    Int_t      fRecenter;            //=1 by default (0=no recentering, 1=recenter all, 2=recenter displacing only one nucleon, 3=recenter by rotation)
    Int_t      fLattice;             //=0 use HCP by default (1=PCS, 2=BCC, 3=FCC)
    Double_t   fSmax;                //Maximum magnitude of cms shift tolerated (99, ie all by default) 
    Int_t      fF;                   //Type of radial distribution
    Int_t      fTrials;              //Store trials needed to complete nucleus
    Int_t      fNonSmeared;          //Store number of non-smeared-node nucleons
    TF1*       fFunc1;               //!Probability density function rho(r)
    TF1*       fFunc2;               //!Probability density function rho(r) -> if set 1 is for p, 2 is for n
    TF2*       fFunc3;               //!Probability density function rho(r,theta) for deformed nuclei
    TObjArray* fNucleons;            //!Array of nucleons
    Double_t   fPhiRot;              //!Angle phi for nucleus
    Double_t   fThetaRot;            //!Angle theta for nucleus
    Double_t   fXRot;                //!Angle around X axis for nucleus
    Double_t   fYRot;                //!Angle around Y axis for nucleus
    Double_t   fZRot;                //!Angle around Z axis for nucleus
    Double_t   fHe3Arr[6000][3][3];  //!Array of events, 3 nucleons, 3 coordinates
    Int_t      fHe3Counter;          //!Event counter
    TBits     *fIsUsed;              //!Bits for lattice use  
    Double_t   fMaxR;                //!maximum radius (15fm)
    void       Lookup(const char* name);
    Bool_t     TestMinDist(Int_t n, Double_t x, Double_t y, Double_t z) const;

  public:
    TGlauNucleus(const char* iname="Pb", Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0);
    virtual ~TGlauNucleus();
    using      TObject::Draw;
    void       Draw(Double_t xs, Int_t colp, Int_t cols);
    Double_t   GetA()             const {return fA;}
    TF1*       GetFunc1()         const {return GetFuncP();}
    TF1*       GetFunc2()         const {return GetFuncN();}
    TF2*       GetFunc3()         const {return GetFuncDef();}
    TF1*       GetFuncP()         const {return fFunc1;}
    TF1*       GetFuncN()         const {return fFunc2;}
    TF2*       GetFuncDef()       const {return fFunc3;}
    Double_t   GetMinDist()       const {return fMinDist;}
    Int_t      GetN()             const {return fN;}
    Double_t   GetNodeDist()      const {return fNodeDist;}
    TObjArray *GetNucleons()      const {return fNucleons;}
    Int_t      GetRecenter()      const {return fRecenter;}
    Double_t   GetR()             const {return fR;}
    Double_t   GetPhiRot()        const {return fPhiRot;}
    Double_t   GetThetaRot()      const {return fThetaRot;}
    Int_t      GetTrials()        const {return fTrials;}
    Int_t      GetNonSmeared()    const {return fNonSmeared;}
    Double_t   GetShiftMax()      const {return fSmax;}
    Double_t   GetW()             const {return fW;}
    Double_t   GetXRot()          const {return fXRot;}
    Double_t   GetYRot()          const {return fYRot;}
    Double_t   GetZRot()          const {return fZRot;}
    void       SetA(Double_t ia, Double_t ia2=-1);
    void       SetBeta(Double_t b2, Double_t b4); 
    void       SetLattice(Int_t i)               {fLattice=i;}
    void       SetMinDist(Double_t min)          {fMinDist=min;}
    void       SetN(Int_t in)                    {fN=in;}
    void       SetNodeDist(Double_t nd)          {fNodeDist=nd;}
    void       SetR(Double_t ir, Double_t ir2=-1);
    void       SetRecenter(Int_t b)              {fRecenter=b;}
    void       SetShiftMax(Double_t s)           {fSmax=s;}
    void       SetSmearing(Double_t s)           {fSmearing=s;}
    void       SetW(Double_t iw);
    TVector3  &ThrowNucleons(Double_t xshift=0.);
    ClassDef(TGlauNucleus,6) // TGlauNucleus class
};

//---------------------------------------------------------------------------------
class TGlauberMC : public TNamed
{
  public:
    class Event {
      public:
        Float_t Npart;       //Number of wounded (participating) nucleons in current event
        Float_t Ncoll;       //Number of binary collisions in current event
        Float_t Nhard;       //Number of hard collisions in current event (based on fHardFrac)
        Float_t B;           //[0,0,16] Impact parameter (b)
        Float_t BNN;         //[0,0,16] Average NN impact parameter
        Float_t Ncollpp;     //Ncoll pp
        Float_t Ncollpn;     //Ncoll pn
        Float_t Ncollnn;     //Ncoll nn
        Float_t VarX;        //[0,0,16] Variance of x of wounded nucleons
        Float_t VarY;        //[0,0,16] Variance of y of wounded nucleons
        Float_t VarXY;       //[0,0,16] Covariance of x and y of wounded nucleons
        Float_t NpartA;      //Number of wounded (participating) nucleons in Nucleus A
        Float_t NpartB;      //Number of wounded (participating) nucleons in Nucleus B
        Float_t Npart0;      //Number of singly-wounded (participating) nucleons
        Float_t AreaW;       //[0,0,16] area defined by width of participants
        Float_t Psi1;        //[0,0,16] psi1
        Float_t Ecc1;        //[0,0,16] eps1
        Float_t Psi2;        //[0,0,16] psi2
        Float_t Ecc2;        //[0,0,16] eps2
        Float_t Psi3;        //[0,0,16] psi3
        Float_t Ecc3;        //[0,0,16] eps3
        Float_t Psi4;        //[0,0,16] psi4
        Float_t Ecc4;        //[0,0,16] eps4
        Float_t Psi5;        //[0,0,16] psi5
        Float_t Ecc5;        //[0,0,16] eps5
        Float_t AreaA;       //[0,0,16] area defined by "and" of participants
        Float_t AreaO;       //[0,0,16] area defined by "or" of participants
        Float_t X0;          //[0,0,16] production point in x
        Float_t Y0;          //[0,0,16] production point in y
        Float_t Phi0;        //[0,0,16] direction in phi
        Float_t Length;      //[0,0,16] length in phi0
        Float_t MeanX;       //[0,0,16] <x> of wounded nucleons
        Float_t MeanY;       //[0,0,16] <y> of wounded nucleons
        Float_t MeanX2;      //[0,0,16] <x^2> of wounded nucleons
        Float_t MeanY2;      //[0,0,16] <y^2> of wounded nucleons
        Float_t MeanXY;      //[0,0,16] <xy> of wounded nucleons
        Float_t MeanXSystem; //[0,0,16] <x> of all nucleons
        Float_t MeanYSystem; //[0,0,16] <y> of all nucleons  
        Float_t MeanXA;      //[0,0,16] <x> of nucleons in nucleus A
        Float_t MeanYA;      //[0,0,16] <y> of nucleons in nucleus A
        Float_t MeanXB;      //[0,0,16] <x> of nucleons in nucleus B
        Float_t MeanYB;      //[0,0,16] <y> of nucleons in nucleus B
        Float_t PhiA;        //[0,0,16] phi angle nucleus A
        Float_t ThetaA;      //[0,0,16] theta angle nucleus B
        Float_t PhiB;        //[0,0,16] phi angle nucleus B
        Float_t ThetaB;      //[0,0,16] theta angle nucleus B
        void    Reset()      {Npart=0;Ncoll=0;Nhard=0;B=0;BNN=0;Ncollpp=0;Ncollpn=0;Ncollnn=0;VarX=0;VarY=0;VarXY=0;NpartA=0;NpartB=0;Npart0=0;AreaW=0;
                              Psi1=0;Ecc1=0;Psi2=0;Ecc2=0;Psi3=0;Ecc3=0;Psi4=0;Ecc4=0;Psi5=0;Ecc5=0;
                              AreaA=0;AreaO=0;X0=0;Y0=0;Phi0=0;Length=0;
                              MeanX=0;MeanY=0;MeanX2=0;MeanY2=0;MeanXY=0;MeanXSystem=0;MeanYSystem=0;MeanXA=0;MeanYA=0;MeanXB=0;MeanYB=0;
                              PhiA=0;ThetaA=0;PhiB=0;ThetaB=0;} // order must match that given in vars below
        ClassDef(TGlauberMC::Event, 1)
    };

  protected:
    TGlauNucleus  fANucleus;       //Nucleus A
    TGlauNucleus  fBNucleus;       //Nucleus B
    Double_t      fXSect;          //Nucleon-nucleon cross section
    Double_t      fXSectOmega;     //StdDev of Nucleon-nucleon cross section
    Double_t      fXSectLambda;    //Jacobian from tot to inelastic (Strikman)
    Double_t      fXSectEvent;     //Event value of Nucleon-nucleon cross section
    TObjArray*    fNucleonsA;      //Array of nucleons in nucleus A
    TObjArray*    fNucleonsB;      //Array of nucleons in nucleus B
    TObjArray*    fNucleons;       //Array which joins Nucleus A & B
    Int_t         fAN;             //Number of nucleons in nucleus A
    Int_t         fBN;             //Number of nucleons in nucleus B
    TNtuple*      fNt;             //Ntuple for results (created, but not deleted)
    Int_t         fEvents;         //Number of events with at least one collision
    Int_t         fTotalEvents;    //All events within selected impact parameter range
    Double_t      fBmin;           //Minimum impact parameter to be generated
    Double_t      fBmax;           //Maximum impact parameter to be generated
    Double_t      fHardFrac;       //Fraction of cross section used for Nhard (def=0.65)
    Int_t         fDetail;         //Detail to store (99=all by default)
    Bool_t        fCalcArea;       //If true calculate overlap area via grid (slow, off by default)
    Bool_t        fCalcLength;     //If true calculate path length (slow, off by default)
    Bool_t        fDoCore;         //If true calculate area and eccentricy only for core participants (off by default)
    Int_t         fMaxNpartFound;  //Largest value of Npart obtained
    Double_t      fPsiN[10];       //Psi N
    Double_t      fEccN[10];       //Ecc N
    Double_t      f2Cx;            //Two-component x
    TF1          *fPTot;           //Cross section distribution
    TF1          *fNNProf;         //NN profile (hard-sphere == 0 by default)
    Event         fEv;             //Glauber event (results of calculation stored in tree)
    Bool_t        fBC[999][999];   //Array to record binary collision
    Bool_t        CalcResults(Double_t bgen);
    Bool_t        CalcEvent(Double_t bgen);

  public:
    TGlauberMC(const char* NA = "Pb", const char* NB = "Pb", Double_t xsect = 42, Double_t xsectsigma=0);
  virtual            ~TGlauberMC() {delete fNt; fNt=0;}

    Double_t            CalcDens(TF1 &prof, Double_t xval, Double_t yval) const;
    void                Draw(Option_t* option="");
    Double_t            GetB()                 const {return fEv.B;}
    Double_t            GetBNN()               const {return fEv.BNN;}
    Double_t            GetBmax()              const {return fBmax;}
    Double_t            GetBmin()              const {return fBmin;}
    Double_t            GetEcc(Int_t i=2)      const {return fEccN[i];}
    Double_t            GetHardFrac()          const {return fHardFrac;}
    Double_t            GetMeanX()             const {return fEv.MeanX;}
    Double_t            GetMeanXParts()        const {return fEv.MeanX;}
    Double_t            GetMeanXSystem()       const {return fEv.MeanXSystem;}
    Double_t            GetMeanY()             const {return fEv.MeanY;}
    Double_t            GetMeanYParts()        const {return fEv.MeanY;}
    Double_t            GetMeanYSystem()       const {return fEv.MeanYSystem;}
    Double_t            GetPsi(Int_t i=2)      const {return fPsiN[i];}
    Double_t            GetSx2()               const {return fEv.VarX;}    
    Double_t            GetSxy()               const {return fEv.VarXY;}    
    Double_t            GetSy2()               const {return fEv.VarY;}    
    Double_t            GetTotXSect()          const;
    Double_t            GetTotXSectErr()       const;
    Double_t            GetXSectEvent()        const {return fXSectEvent;}
    Int_t               GetNcoll()             const {return fEv.Ncoll;}
    Int_t               GetNcollnn()           const {return fEv.Ncollnn;}
    Int_t               GetNcollpn()           const {return fEv.Ncollpn;}
    Int_t               GetNcollpp()           const {return fEv.Ncollpp;}
    Int_t               GetNpart()             const {return fEv.Npart;}
    Int_t               GetNpart0()            const {return fEv.Npart0;}
    Int_t               GetNpartA()            const {return fEv.NpartA;}
    Int_t               GetNpartB()            const {return fEv.NpartB;}
    Int_t               GetNpartFound()        const {return fMaxNpartFound;}
    TF1*                GetXSectDist()         const {return fPTot;}
    TGlauNucleus*       GetNucleusA()                {return &fANucleus;}
    TGlauNucleus*       GetNucleusB()                {return &fBNucleus;}
    TNtuple*            GetNtuple()            const {return fNt;}
    TObjArray          *GetNucleons();
    const Event        &GetEvent()             const {return fEv;}
    const Event        *GetEvent()                   {return &fEv;}
    const TGlauNucleus* GetNucleusA()          const {return &fANucleus;}
    const TGlauNucleus* GetNucleusB()          const {return &fBNucleus;}
    Bool_t              IsBC(Int_t i, Int_t j) const {return fBC[i][j];}
    Bool_t              NextEvent(Double_t bgen=-1);
    void                Reset()                      {delete fNt; fNt=0; }
    Bool_t              ReadNextEvent(Bool_t calc=1, const char *fname=0);       
    void                Run(Int_t nevents,Double_t b=-1);
    void                Set2Cx(Double_t x)           {f2Cx = x;}
    void                SetBmax(Double_t bmax)       {fBmax = bmax;}
    void                SetBmin(Double_t bmin)       {fBmin = bmin;}
    void                SetCalcArea(Bool_t b)        {fCalcArea = b;}
    void                SetCalcCore(Bool_t b)        {fDoCore = b;}
    void                SetCalcLength(Bool_t b)      {fCalcLength = b;}
    void                SetDetail(Int_t d)           {fDetail = d;}
    void                SetHardFrac(Double_t f)      {fHardFrac=f;}
    void                SetLattice(Int_t i)          {fANucleus.SetLattice(i); fBNucleus.SetLattice(i);}
    void                SetMinDistance(Double_t d)   {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}
    void                SetNNProf(TF1 *f1)           {fNNProf = f1;}
    void                SetNodeDistance(Double_t d)  {fANucleus.SetNodeDist(d); fBNucleus.SetNodeDist(d);}
    void                SetRecenter(Int_t b)         {fANucleus.SetRecenter(b); fBNucleus.SetRecenter(b);}
    void                SetShiftMax(Double_t s)      {fANucleus.SetShiftMax(s); fBNucleus.SetShiftMax(s);}
    void                SetSmearing(Double_t s)      {fANucleus.SetSmearing(s); fBNucleus.SetSmearing(s);}
    const char         *Str()                  const {return Form("gmc-%s%s-snn%.1f-md%.1f-nd%.1f-rc%d-smax%.1f",fANucleus.GetName(),fBNucleus.GetName(),fXSect,fBNucleus.GetMinDist(),fBNucleus.GetNodeDist(),fBNucleus.GetRecenter(),fBNucleus.GetShiftMax());}
    static void         PrintVersion()               {cout << "TGlauberMC " << Version() << endl;}
    static const char  *Version()                    {return "v3.1a";}

    ClassDef(TGlauberMC,6) // TGlauberMC class
};

//---------------------------------------------------------------------------------
TF1 *getNNProf(Double_t snn, Double_t omega, Double_t G) 
{ // NN collisoin profile from https://arxiv.org/abs/1307.0636
  if ((omega<0) || (omega>1))
    return 0;
  Double_t R2 = snn/10./TMath::Pi();
  TF1 *nnprof = new TF1("nnprofgamma","[2]*(1-TMath::Gamma([0],[1]*x^2))",0,3);
  nnprof->SetParameters(1./omega,G/omega/R2,G);
  return nnprof;
}

//---------------------------------------------------------------------------------
void runAndSaveNtuple(const Int_t n,
                      const char *sysA,
                      const char *sysB,
                      const Double_t signn,
                      const Double_t sigwidth,
                      const Double_t mind,
		      const Double_t omega,
                      const Double_t noded,
                      const char *fname)
{
  TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn,sigwidth);
  mcg->SetMinDistance(mind);
  mcg->SetNodeDistance(noded);
  mcg->SetCalcLength(0);
  mcg->SetCalcArea(0);
  mcg->SetCalcCore(0);
  mcg->SetDetail(99);
  TString om;
  if ((omega>=0) && (omega<=1)) {
    TF1 *f1 = getNNProf(signn, omega);
    mcg->SetNNProf(f1);
    om=Form("-om%.1f",omega);
  }
  TString name;
  if (fname) 
    name = fname; 
  else {
    TString nd;
    if (noded>0) 
      nd=Form("-nd%.1f",noded);
    name = Form("%s%s%s.root",mcg->Str(),om.Data(),nd.Data());
  }
  mcg->Run(n);
  TFile out(name,"recreate",name,9);
  TNtuple  *nt=mcg->GetNtuple();
  nt->Write();
  out.Close();
}

//---------------------------------------------------------------------------------
void runAndSaveNucleons(const Int_t n,                    
                        const char *sysA,           
                        const char *sysB,           
                        const Double_t signn,
                        const Double_t sigwidth,
                        const Double_t mind,
                        const Bool_t verbose,
                        const char *fname)
{
  TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn,sigwidth);
  mcg->SetMinDistance(mind);
  TFile *out=0;
  if (fname) 
    out=new TFile(fname,"recreate",fname,9);

  for (Int_t ievent=0; ievent<n; ++ievent) {
    //get an event with at least one collision
    mcg->Run(1);
    if (ievent%100==0)
      cout << "\r" << 100.*ievent/n << "% done" << flush;

    //access, save and (if wanted) print out nucleons
    TObjArray* nucleons=mcg->GetNucleons();
    if (!nucleons) 
      continue;
    if (out)
      nucleons->Write(Form("nucleonarray%d",ievent),TObject::kSingleKey);

    if (verbose) {
      cout<<endl<<endl<<"EVENT NO: "<<ievent<<endl;
      cout<<"B = "<<mcg->GetB()<<"  Npart = "<<mcg->GetNpart()<<endl<<endl;
      printf("Nucleus\t X\t Y\t Z\tNcoll\n");
      Int_t nNucls=nucleons->GetEntries();
      for (Int_t iNucl=0; iNucl<nNucls; ++iNucl) {
        TGlauNucleon *nucl=(TGlauNucleon *)nucleons->At(iNucl);
        Char_t nucleus='A';
        if (nucl->IsInNucleusB()) 
	  nucleus='B';
        Double_t x=nucl->GetX();
        Double_t y=nucl->GetY();
        Double_t z=nucl->GetZ();
        Int_t ncoll=nucl->GetNColl();
        printf("   %c\t%2.2f\t%2.2f\t%2.2f\t%3d\n",nucleus,x,y,z,ncoll);
      }
    }
  }
  cout << endl << "Done!" << endl;
  if (out) {
    TNtuple *nt = mcg->GetNtuple();
    nt->Write();
    if (verbose)
      out->ls();
    delete out;
  }
}

//---------------------------------------------------------------------------------
void runAndSmearNtuple(const Int_t n,
                       const Double_t sigs,
                       const char *sysA,
                       const char *sysB,
                       const Double_t signn,
                       const Double_t mind,
                       const char *fname)
{
  // Run Glauber and store ntuple with smeared eccentricities in file.

  TGlauberMC *mcg = new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);
  TFile *out = TFile::Open(fname,"recreate",fname,9);
  if (!out)
    return;
  TNtuple *nt = new TNtuple("nt","nt",
      "Npart:Ncoll:B:Psi1P:Ecc1P:Psi2P:Ecc2P:Psi3P:Ecc3P:Psi4P:Ecc4P:Psi5P:Ecc5P:Psi1G:Ecc1G:Psi2G:Ecc2G:Psi3G:Ecc3G:Psi4G:Ecc4G:Psi5G:Ecc5G:Sx2P:Sy2P:Sx2G:Sy2G");
  nt->SetDirectory(out);

  const Int_t NSAMP = 100;
  TF1 *rad = new TF1("rad","x*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,3*sigs);
  rad->SetParameter(0,sigs);

  for (Int_t ievent=0; ievent<n; ++ievent) {
    while (!mcg->NextEvent()) {}

    const TGlauNucleus *nucA   = mcg->GetNucleusA();
    const TObjArray *nucleonsA = nucA->GetNucleons();
    const Int_t AN             = nucA->GetN();
    const TGlauNucleus *nucB   = mcg-> GetNucleusB();
    const TObjArray *nucleonsB = nucB->GetNucleons();
    const Int_t BN             = nucB->GetN();
    Double_t sinphi[10] = {0};
    Double_t cosphi[10] = {0};
    Double_t rn[10]     = {0};
    Double_t ecc[10]    = {0};
    Double_t psi[10]    = {0};
    Double_t sx2g       = 0;
    Double_t sy2g       = 0;

    for (Int_t s=0; s<NSAMP; ++s) {
      Int_t ni = 0;
      Double_t xvals[1000] = {0};
      Double_t yvals[1000] = {0};
      for (Int_t i = 0; i<AN; ++i) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
        if (!nucleonA->IsWounded())
          continue;
        Double_t sr = rad->GetRandom();
        Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
        xvals[ni]   = nucleonA->GetX() + sr*TMath::Cos(sp);
        yvals[ni]   = nucleonA->GetY() + sr*TMath::Sin(sp);
        ++ni;
      }
      for (Int_t i = 0; i<BN; ++i) {
        TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(i));
        if (!nucleonB->IsWounded())
          continue;
        Double_t sr = rad->GetRandom();
        Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
        xvals[ni]   = nucleonB->GetX() + sr*TMath::Cos(sp);
        yvals[ni]   = nucleonB->GetY() + sr*TMath::Sin(sp);
        ++ni;
      }

      Double_t MeanX  = 0;
      Double_t MeanY  = 0;
      Double_t MeanX2 = 0;
      Double_t MeanY2 = 0;
      for (Int_t i = 0; i<ni; ++i) {
        MeanX  += xvals[i];
        MeanY  += yvals[i];
        MeanX2 += xvals[i]*xvals[i];
        MeanY2 += yvals[i]*yvals[i];
      }
      MeanX  /= ni;
      MeanY  /= ni;
      MeanX2 /= ni;
      MeanY2 /= ni;
      sx2g        += MeanX2-MeanX*MeanX;
      sy2g        += MeanY2-MeanY*MeanY;

      for (Int_t j = 1; j<9; ++j) {
        for (Int_t i = 0; i<ni; ++i) {
          Double_t x   = xvals[i] - MeanX;
          Double_t y   = yvals[i] - MeanY;
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
    }
    for (Int_t j = 1; j<9; ++j) {
      psi[j] = (TMath::ATan2(sinphi[j],cosphi[j]) + TMath::Pi())/j;
      ecc[j] = TMath::Sqrt(sinphi[j]*sinphi[j] + cosphi[j]*cosphi[j]) / rn[j];
    }

    Float_t v[27]; Int_t i=0;
    v[i++] = mcg->GetNpart();
    v[i++] = mcg->GetNcoll();
    v[i++] = mcg->GetB();
    v[i++] = mcg->GetPsi(1); // point-like calculation values
    v[i++] = mcg->GetEcc(1);
    v[i++] = mcg->GetPsi(2);
    v[i++] = mcg->GetEcc(2);
    v[i++] = mcg->GetPsi(3);
    v[i++] = mcg->GetEcc(3);
    v[i++] = mcg->GetPsi(4);
    v[i++] = mcg->GetEcc(4);
    v[i++] = mcg->GetPsi(5);
    v[i++] = mcg->GetEcc(5);
    v[i++] = psi[1];         // Gaussian smeared values
    v[i++] = ecc[1];
    v[i++] = psi[2];
    v[i++] = ecc[2];
    v[i++] = psi[3];
    v[i++] = ecc[3];
    v[i++] = psi[4];
    v[i++] = ecc[4];
    v[i++] = psi[5];
    v[i++] = ecc[5];
    v[i++] = mcg->GetSx2();
    v[i++] = mcg->GetSy2();
    v[i++] = sx2g/NSAMP;
    v[i++] = sy2g/NSAMP;
    nt->Fill(v);
  }

  out->Write();
  out->Close();
  delete out;
}

//---------------------------------------------------------------------------------
void runAndCalcDens(const Int_t n,
		    const Double_t alpha,
		    const char *sysA,
		    const char *sysB,
		    const Double_t signn,
		    const Double_t mind,
		    const char *fname)
{
  // Run Glauber and store per event a density profile in x and y, calculated from participant and binary positions
  // with relative weight given by alpha.

  TGlauberMC *mcg = new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);

  TFile *out = 0;
  TCanvas *c1 = 0;
  if (fname) {
    out = TFile::Open(fname,"recreate",fname,9);
    if (!out)
      return;
  } else {
    c1 = new TCanvas;
  }

  const Int_t NSAMP = 100;
  const Double_t wp = (1-alpha)/2;
  const Double_t wb = alpha;
  const Double_t sigs = TMath::Sqrt(signn/20/TMath::Pi()); //from arXiv::0711.3724


  TF1 *rad = new TF1("rad","x*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,3*sigs);
  rad->SetParameter(0,sigs);

  TH2F *h2f = new TH2F("h2f",";x (fm);y (fm)",121,-15.5,15.5,121,-15.5,15.5);
  h2f->SetStats(0);
  
  for (Int_t ievent=0; ievent<n; ++ievent) {
    while (!mcg->NextEvent()) {}
    h2f->Reset();
    h2f->SetName(Form("Event_%d",ievent));
    h2f->SetTitle(Form("Npart=%d, Ncoll=%d",mcg->GetNpart(), mcg->GetNcoll()));

    const TGlauNucleus *nucA   = mcg->GetNucleusA();
    const TObjArray *nucleonsA = nucA->GetNucleons();
    const Int_t AN             = nucA->GetN();
    const TGlauNucleus *nucB   = mcg-> GetNucleusB();
    const TObjArray *nucleonsB = nucB->GetNucleons();
    const Int_t BN             = nucB->GetN();

    for (Int_t i = 0; i<AN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
      if (!nucleonA->IsWounded())
	continue;
      Double_t xA=nucleonA->GetX();
      Double_t yA=nucleonA->GetY();
      for (Int_t s=0; s<NSAMP; ++s) {
	Double_t sr = rad->GetRandom();
	Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	h2f->Fill(xA+sr*TMath::Cos(sp),yA+sr*TMath::Sin(sp),wp);
      }
    }

    for (Int_t j = 0; j<BN; ++j) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(j));
      if (!nucleonB->IsWounded())
	continue;
      Double_t xB=nucleonB->GetX();
      Double_t yB=nucleonB->GetY();
      for (Int_t s=0; s<NSAMP; ++s) {
	Double_t sr = rad->GetRandom();
	Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	h2f->Fill(xB+sr*TMath::Cos(sp),yB+sr*TMath::Sin(sp),wp);
      }
    }

    if (alpha>0) {
      for (Int_t i = 0; i<AN; ++i) {
	TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
	if (!nucleonA->IsWounded())
	  continue;
	Double_t xA=nucleonA->GetX();
	Double_t yA=nucleonA->GetY();
	for (Int_t j = 0; j<BN; ++j) {
	  TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(j));
	  if (!mcg->IsBC(i,j))
	    continue;
	  Double_t xB=nucleonB->GetX();
	  Double_t yB=nucleonB->GetY();
	  Double_t dX=(xA+xB)/2;
	  Double_t dY=(yA+yB)/2;
	  for (Int_t s=0; s<NSAMP; ++s) {
	    Double_t sr = rad->GetRandom();
	    Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	    h2f->Fill(dX+sr*TMath::Cos(sp),dY+sr*TMath::Sin(sp),wb);
	  }
	}
      }
    }
    h2f->Scale(1./h2f->Integral());
    if (out) {
      h2f->Write();
    } else {
      h2f->Draw("colz");
      c1->Update();
      gSystem->Sleep(1000);
    }
  }
  if (out) {
    out->Write();
    out->Close();
    delete out;
  }
}

//---------------------------------------------------------------------------------
ClassImp(TGlauNucleus)
//---------------------------------------------------------------------------------
TGlauNucleus::TGlauNucleus(const char* iname, Int_t iN, Double_t iR, Double_t ia, Double_t iw, TF1* ifunc) : 
  TNamed(iname,""),
  fN(iN),fR(iR),fA(ia),fW(iw),fR2(0),fA2(0),fW2(0),fBeta2(0),fBeta4(0),
  fMinDist(0.4),fNodeDist(0.0),fSmearing(0.0),fRecenter(1),fLattice(0),fSmax(99),
  fF(0),fTrials(0),fNonSmeared(0),fFunc1(ifunc),fFunc2(0),fFunc3(0),fNucleons(0),
  fPhiRot(0),fThetaRot(0),fHe3Counter(-1),fIsUsed(0),fMaxR(14)
{
  if (fN==0) {
    cout << "Setting up nucleus " << iname << endl;
    Lookup(iname);
  }
}

TGlauNucleus::~TGlauNucleus()
{
  if (fIsUsed)
    delete fIsUsed;
  if (fNucleons)
    delete fNucleons;
  delete fFunc1;
  delete fFunc2;
  delete fFunc3;
}

void TGlauNucleus::Draw(Double_t xs, Int_t colp, Int_t cols)
{
  Double_t r = 0.5*TMath::Sqrt(xs/TMath::Pi()/10.);
  TEllipse en;
  en.SetLineStyle(1);
  en.SetLineWidth(1);
  en.SetFillStyle(1001);
  for (Int_t i = 0; i<fNucleons->GetEntries(); ++i) {
    TGlauNucleon* gn = (TGlauNucleon*) fNucleons->At(i);
    if (!gn->IsSpectator()) {
      en.SetFillColor(colp);
      en.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
    } else {
      en.SetFillColor(cols);
      en.SetFillStyle(1001);
      en.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
    }
  }
}

void TGlauNucleus::Lookup(const char* name)
{
  SetName(name);

  TString tmp(name);

  Double_t r0 = 0;
  Double_t r1 = 0;
  Double_t r2 = 0;

  if      (TString(name) == "p")       {fN = 1;   fR = 0.234;      fA = 0;      fW =  0;       fF = 0;  fZ=1;}
  else if (TString(name) == "pg")      {fN = 1;   fR = 0.514;      fA = 0;      fW =  0;       fF = 9;  fZ=1;} 
  else if (TString(name) == "pdg")     {fN = 1;   fR = 1;          fA = 0;      fW =  0;       fF = 10; fZ=1;} // from arXiv:1101.5953
  else if (TString(name) == "d")       {fN = 2;   fR = 0.01;       fA = 0.5882; fW =  0;       fF = 1;  fZ=1;}
  else if (TString(name) == "dh")      {fN = 2;   fR = 0.01;       fA = 0.5882; fW =  0;       fF = 3;  fZ=1;}
  else if (TString(name) == "dhh")     {fN = 2;   fR = 0.01;       fA = 0.5882; fW =  0;       fF = 4;  fZ=1;}
  else if (TString(name) == "He3")     {fN = 3;   fR = 0.01;       fA = 0.5882; fW =  0;       fF = 6;  fZ=1;}
  else if (TString(name) == "H3")      {fN = 3;   fR = 0.01;       fA = 0.5882; fW =  0;       fF = 6;  fZ=2;}
  else if (TString(name) == "O")       {fN = 16;  fR = 2.608;      fA = 0.513;  fW = -0.051;   fF = 1;  fZ=8;}
  else if (TString(name) == "Si")      {fN = 28;  fR = 3.34;       fA = 0.580;  fW = -0.233;   fF = 1;  fZ=14;}
  else if (TString(name) == "Si2")     {fN = 28;  fR = 3.34;       fA = 0.580;  fW =  0;       fF = 8;  fZ=14; fBeta2=-0.478; fBeta4=0.250;}
  else if (TString(name) == "S")       {fN = 32;  fR = 2.54;       fA = 2.191;  fW =  0.16;    fF = 2;  fZ=16;}
  else if (TString(name) == "Ar")      {fN = 40;  fR = 3.53;       fA = 0.542;  fW =  0;       fF = 1;  fZ=18;}
  else if (TString(name) == "Ca")      {fN = 40;  fR = 3.766;      fA = 0.586;  fW = -0.161;   fF = 1;  fZ=20;}
  else if (TString(name) == "Ni")      {fN = 58;  fR = 4.309;      fA = 0.517;  fW = -0.1308;  fF = 1;  fZ=28;}
  else if (TString(name) == "Cu")      {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 1;  fZ=29;}
  else if (TString(name) == "Curw ")   {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 12; fZ=29; r0=1.00898; r1=-0.000790403; r2=-0.000389897;} 
  else if (TString(name) == "Cu2")     {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 8;  fZ=29; fBeta2=0.162; fBeta4=-0.006;}  
  else if (TString(name) == "Cu2rw")   {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 14; fZ=29; fBeta2=0.162; fBeta4=-0.006; r0=1.01269; r1=-0.00298083; r2=-9.97222e-05;}  
  else if (TString(name) == "CuHN")    {fN = 63;  fR = 4.28;       fA = 0.5;    fW =  0;       fF = 1;  fZ=29;} // from arXiv:0904.4080v1
  else if (TString(name) == "Xe")      {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 1;  fZ=54;} // adapted from arXiv:1703.04278
  else if (TString(name) == "Xes")     {fN = 129; fR = 5.42;       fA = 0.57;   fW =  0;       fF = 1;  fZ=54;} // scale from Sb (Antimony, A=122, r=5.32) by 1.019 = (129/122)**0.333
  else if (TString(name) == "Xe2")     {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 8;  fZ=54; fBeta2=0.161; fBeta4=-0.003;} // adapted from arXiv:1703.04278 and Z. Physik (1974) 270: 113
  else if (TString(name) == "Xe2a")    {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 8;  fZ=54; fBeta2=0.18; fBeta4=0;} // ALICE parameters (see public note from 2018 at https://cds.cern.ch/collection/ALICE%20Public%20Notes?ln=en)
  else if (TString(name) == "Xerw")    {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 12; fZ=54; r0=1.00911; r1=-0.000722999; r2=-0.0002663;}
  else if (TString(name) == "Xesrw")   {fN = 129; fR = 5.42;       fA = 0.57;   fW =  0;       fF = 12; fZ=54; r0=1.0096; r1=-0.000874123; r2=-0.000256708;}
  else if (TString(name) == "Xe2arw")  {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 14; fZ=54; fBeta2=0.18; fBeta4=0; r0=1.01246; r1=-0.0024851; r2=-5.72464e-05;} 
  else if (TString(name) == "W")       {fN = 186; fR = 6.58;       fA = 0.480;  fW =  0;       fF = 1;  fZ=74;}
  else if (TString(name) == "Au")      {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 1;  fZ=79;}
  else if (TString(name) == "Aurw")    {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 12; fZ=79; r0=1.00899; r1=-0.000590908; r2=-0.000210598;}
  else if (TString(name) == "Au2")     {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 8;  fZ=79; fBeta2=-0.131; fBeta4=-0.031; }
  else if (TString(name) == "Au2rw")   {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 14; fZ=79; fBeta2=-0.131; fBeta4=-0.031; r0=1.01261; r1=-0.00225517; r2=-3.71513e-05;}
  else if (TString(name) == "AuHN")    {fN = 197; fR = 6.42;       fA = 0.44;   fW =  0;       fF = 1;  fZ=79;} // from arXiv:0904.4080v1
  else if (TString(name) == "Pb")      {fN = 208; fR = 6.62;       fA = 0.546;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "Pbrw")    {fN = 208; fR = 6.62;       fA = 0.546;  fW =  0;       fF = 12; fZ=82; r0=1.00863; r1=-0.00044808; r2=-0.000205872;} //only Pb 207 was tested but should be the same for 208
  else if (TString(name) == "Pb*")     {fN = 208; fR = 6.624;      fA = 0.549;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "PbHN")    {fN = 208; fR = 6.65;       fA = 0.460;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "Pbpn")    {fN = 208; fR = 6.68;       fA = 0.447;  fW =  0;       fF = 11; fZ=82; fR2=6.69; fA2=0.56; fW2=0;}
  else if (TString(name) == "Pbpnrw")  {fN = 208; fR = 6.68;       fA = 0.447;  fW =  0;       fF = 13; fZ=82; fR2=6.69; fA2=0.56; fW2=0;}
  // Uranium description taken from Heinz & Kuhlman, nucl-th/0411054.  In this code, fR is defined as 6.8*0.91, fW=6.8*0.26
  else if (TString(name) == "U")       {fN = 238; fR = 6.188;      fA = 0.54;   fW =  1.77;    fF = 5;  fZ=92;}  
  else if (TString(name) == "U2")      {fN = 238; fR = 6.67;       fA = 0.44;   fW =  0;       fF = 8;  fZ=92; fBeta2=0.280; fBeta4=0.093;}
  else {
    cout << "Could not find nucleus " << name << endl;
    return;
  }

  switch (fF) {
    case 0: // Proton exp
      fFunc1 = new TF1("prot","x*x*exp(-x/[0])",0,5);
      fFunc1->SetParameter(0,fR);
      break;
    case 1: // 3pF
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      break;
    case 2: // 3pG
      fFunc1 = new TF1("3pg","x*x*(1+[2]*(x/[0])**2)/(1+exp((x**2-[0]**2)/[1]**2))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      break;
    case 3: // Hulthen (see nucl-ex/0603010)
      fFunc1 = new TF1("f3","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,fMaxR);
      fFunc1->SetParameters(1/4.38,1/.85);
      break;
    case 4: // Hulthen HIJING
      fFunc1 = new TF1("f4","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,fMaxR);
      fFunc1->SetParameters(2/4.38,2/.85);
      break;
    case 5: // Ellipsoid (Uranium)
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,0); // same as 3pF but setting W to zero
      break;
    case 6: // He3/H3
      fFunc1 = 0; // read in file instead
      break;
    case 7: // Deformed nuclei, box method
#ifndef HAVE_MATHMORE
      cerr << "Need libMathMore.so for deformed nuclei" << endl;
      gSystem->Exit(123);
#endif
     fFunc1 = 0; // no func: only need beta parameters and use uniform box distribution
      break;
    case 8: // Deformed nuclei, TF2 method
      fFunc3 = new TF2("f77","x*x*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))",0,fMaxR,0.0,TMath::Pi());
      fFunc3->SetNpx(120);
      fFunc3->SetNpy(120);
      fFunc3->SetParameters(fR,fA,fBeta2,fBeta4);
      break;
    case 9: // Proton gaus
      fFunc1 = new TF1("prot","x*x*exp(-x*x/[0]/[0]/2)",0,5);
      fFunc1->SetParameter(0,fR);
      break;
    case 10: // Proton dgaus
      fFunc1 = new TF1("prot","x*x*((1-[0])/[1]^3*exp(-x*x/[1]/[1])+[0]/(0.4*[1])^3*exp(-x*x/(0.4*[1])^2))",0,5);
      fFunc1->SetParameter(0,0.5);
      fFunc1->SetParameter(1,fR);
      break;
    case 11: // 3pF for proton and neutrons
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      fFunc2 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc2->SetParameters(fR2,fA2,fW2);
      break;
    case 12: // reweighted
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW,r0,r1,r2); 
      fRecenter=1;
      fSmax=0.1;
      break;
    case 13: // Pb for proton and neutrons reweighted
      fFunc1 = new TF1(Form("%s_prot",name),"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW,1.00866,-0.000461484,-0.000203571);
      fFunc2 = new TF1(Form("%s_neut",name),"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc2->SetParameters(fR2,fA2,fW2,1.00866,-0.000461484,-0.000203571);
      fRecenter=1;
      fSmax=0.1;
      break;
    case 14: // Deformed nuclei, TF2 method, reweighted
      fFunc3 = new TF2("f77","x*x*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))/([4]+[5]*x+[6]*x^2)",0,fMaxR,0.0,TMath::Pi());
      fFunc3->SetNpx(120);
      fFunc3->SetNpy(120);
      fFunc3->SetParameters(fR,fA,fBeta2,fBeta4,r0,r1,r2);
      fRecenter=1;
      fSmax=0.1;
      break;
    default:
      cerr << "Could not find function type " << fF << endl;
  }
  return;
}

void TGlauNucleus::SetA(Double_t ia, Double_t ia2)
{
  fA  = ia;
  fA2 = ia2;
  switch (fF) {
    case 1:  // 3pF
    case 12: // 3pF with pol2 normalization
    case 2:  // 3pG
    case 5:  // Ellipsoid (Uranium)
      fFunc1->SetParameter(1,fA);
      break;
    case 8:
      fFunc3->SetParameter(1,fA);
      break;
    case 11: //p&n
      fFunc1->SetParameter(1,fA);//proton
      fFunc2->SetParameter(1,fA2);//neutron
      break;
    default:
      cout << "Error: fA not needed for function " << fF <<endl;
  }
}

void TGlauNucleus::SetBeta(Double_t b2, Double_t b4) 
{
  fBeta2=b2; 
  fBeta4=b4;      
  if (fFunc3) {
    fFunc3->SetParameter(2,fBeta2);
    fFunc3->SetParameter(3,fBeta4);
  }
}

void TGlauNucleus::SetR(Double_t ir, Double_t ir2)
{
  fR  = ir;
  fR2 = ir2;
  switch (fF) {
    case 0:  // Proton exp
    case 9:  // Proton gaus
    case 1:  // 3pF
    case 12: // 3pF with pol2 normalization
    case 2:  // 3pG
    case 5:  // Ellipsoid (Uranium)
      fFunc1->SetParameter(0,fR);
      break;
    case 8:
      fFunc3->SetParameter(0,fR);
      break;
    case 10: // Proton
      fFunc1->SetParameter(1,fR);
      break;
    case 11: // p&n
      fFunc1->SetParameter(0,fR);//proton
      fFunc2->SetParameter(0,fR2);//neutron
      break;
    default:
      cout << "Error: fR not needed for function " << fF <<endl;
  }
}

void TGlauNucleus::SetW(Double_t iw)
{
  fW = iw;
  switch (fF) {
    case 1: // 3pF
    case 2: // 3pG
      fFunc1->SetParameter(2,fW);
      break;
    default:
      cout << "Error: fW not needed for function " << fF <<endl;
  }
}

Bool_t TGlauNucleus::TestMinDist(Int_t n, Double_t x, Double_t y, Double_t z) const
{
  if (fMinDist<=0)
    return kTRUE;
  const Double_t md2 = fMinDist*fMinDist; 
  for (Int_t j = 0; j<n; ++j) {
    TGlauNucleon *other=(TGlauNucleon*)fNucleons->At(j);
    Double_t xo=other->GetX();
    Double_t yo=other->GetY();
    Double_t zo=other->GetZ();
    Double_t dist2 = (x-xo)*(x-xo)+
		                 (y-yo)*(y-yo)+
		                 (z-zo)*(z-zo);
    if (dist2<md2) {
      return kFALSE;
    }
  }
  return kTRUE;
}

TVector3 &TGlauNucleus::ThrowNucleons(Double_t xshift)
{
  if (fNucleons==0) {
    fNucleons=new TObjArray(fN);
    fNucleons->SetOwner();
    for (Int_t i=0; i<fN; ++i) {
      TGlauNucleon *nucleon=new TGlauNucleon(); 
      nucleon->SetType(0);
      if (i<fZ) 
        nucleon->SetType(1);
      fNucleons->Add(nucleon); 
    }
  } 
  if (1) { //randomize p and n in nucleus
    for (Int_t i=0,iz=0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)fNucleons->At(i);
      Double_t frac=double(fZ-iz)/double(fN-i);
      Double_t rn=gRandom->Uniform(0,1);
      if (rn<frac) {
        nucleon->SetType(1);
        ++iz;
      } else {
        nucleon->SetType(0);
      }
    }
  }
 cmscheck: /* start over here in case shift was too large */
  fTrials = 0;
  fNonSmeared = 0;
  fPhiRot = gRandom->Rndm()*2*TMath::Pi();
  const Double_t cosThetaRot = 2*gRandom->Rndm()-1;
  fThetaRot = TMath::ACos(cosThetaRot);
  fXRot = gRandom->Rndm()*2*TMath::Pi();
  fYRot = gRandom->Rndm()*2*TMath::Pi();
  fZRot = gRandom->Rndm()*2*TMath::Pi();
  const Bool_t hulthen = (TString(GetName())=="dh" || TString(GetName())=="dhh");
  const Bool_t helium3 = (TString(GetName())=="He3") || (TString(GetName())=="H3");
  if (fN==1) { //special treatment for proton
    Double_t r = fFunc1->GetRandom();
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t ctheta = 2*gRandom->Rndm() - 1;
    Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(0));
    nucleon->Reset();
    nucleon->SetXYZ(r * stheta * TMath::Cos(phi),
		    r * stheta * TMath::Sin(phi),
		    r * ctheta);
    fTrials = 1;
  } else if (fN==2 && hulthen) { //special treatment for Hulten
    Double_t r = fFunc1->GetRandom()/2;
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t ctheta = 2*gRandom->Rndm() - 1;
    Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);

    TGlauNucleon *nucleon1=(TGlauNucleon*)(fNucleons->At(0));
    TGlauNucleon *nucleon2=(TGlauNucleon*)(fNucleons->At(1));
    nucleon1->Reset();
    nucleon1->SetXYZ(r * stheta * TMath::Cos(phi),
                     r * stheta * TMath::Sin(phi),
                     r * ctheta);
    nucleon2->Reset();
    nucleon2->SetXYZ(-nucleon1->GetX(),
                     -nucleon1->GetY(),
                     -nucleon1->GetZ());
    fTrials = 1;
  } else if (helium3) { 
    if (fHe3Counter == -1) {
      // read in the ascii file into the array and step through the counter
      char filename[100] = "he3_plaintext.dat";
      if ((TString(GetName())=="H3")) {
        sprintf(filename,"h3_plaintext.dat");
      }
      cout << "Reading in the " << filename << " file upon initialization" << endl;
      ifstream myfile;
      myfile.open(filename);
      if (!myfile) {
        cout << "ERROR:  no file for He3/H3 found with name = " << filename << endl;
        gSystem->Exit(123);
      }
      cout << "Reading file for He3/H3 found with name = " << filename << endl;
      Int_t inputcounter = 0;
      while (myfile) {
        if (inputcounter > 6000) break;
        Double_t foo;
        myfile >> fHe3Arr[inputcounter][0][0] >> fHe3Arr[inputcounter][0][1] >> fHe3Arr[inputcounter][0][2]
               >> fHe3Arr[inputcounter][1][0] >> fHe3Arr[inputcounter][1][1] >> fHe3Arr[inputcounter][1][2]
               >> fHe3Arr[inputcounter][2][0] >> fHe3Arr[inputcounter][2][1] >> fHe3Arr[inputcounter][2][2]
               >> foo >> foo >> foo >> foo;
        ++inputcounter;
      }  
      myfile.close();
      fHe3Counter=0;
    } // done reading in the file the first time
    if (fHe3Counter > 5999) 
      fHe3Counter = 0;
    TGlauNucleon *nucleon1=(TGlauNucleon*)(fNucleons->At(0));
    TGlauNucleon *nucleon2=(TGlauNucleon*)(fNucleons->At(1));
    TGlauNucleon *nucleon3=(TGlauNucleon*)(fNucleons->At(2));
    nucleon1->Reset();
    nucleon1->SetXYZ(fHe3Arr[fHe3Counter][0][0],
                     fHe3Arr[fHe3Counter][0][1],
                     fHe3Arr[fHe3Counter][0][2]);
    nucleon1->RotateXYZ(fPhiRot,fThetaRot);
    nucleon2->Reset();
    nucleon2->SetXYZ(fHe3Arr[fHe3Counter][1][0],
                     fHe3Arr[fHe3Counter][1][1],
                     fHe3Arr[fHe3Counter][1][2]);
    nucleon2->RotateXYZ(fPhiRot,fThetaRot);
    nucleon3->Reset();
    nucleon3->SetXYZ(fHe3Arr[fHe3Counter][2][0],
                     fHe3Arr[fHe3Counter][2][1],
                     fHe3Arr[fHe3Counter][2][2]);
    nucleon3->RotateXYZ(fPhiRot,fThetaRot);
    ++fHe3Counter;
    fTrials = 1;
  } else { // all other nuclei 
    const Double_t startingEdge  = 20; // throw nucleons within a cube of this size (fm)
    const Double_t startingEdgeX = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Double_t startingEdgeY = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Double_t startingEdgeZ = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Int_t nslots = 2*startingEdge/fNodeDist+1;
    if (fNodeDist>0) {
      if (fMinDist>fNodeDist) {
        cout << "Minimum distance (nucleon hard core diameter) [" 
          << fMinDist << "] cannot be larger than the nodal spacing of the grid [" 
          << fNodeDist << "]." << endl;
        cout << "Quitting...." << endl;
        gSystem->Exit(123);
      }
      if (!fIsUsed)
        fIsUsed = new TBits(nslots*nslots*nslots);
      else
        fIsUsed->ResetAllBits();
    }
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      while (1) {
        ++fTrials;
        Bool_t nucleon_inside = 0;
        Double_t x=999, xsmeared=999;
        Double_t y=999, ysmeared=999;
        Double_t z=999, zsmeared=999;
        if (fF==5||fF==7) { // the extended way, throw in a box and test the weight
          while (!nucleon_inside) {
            x = (fR*2)*(gRandom->Rndm() * 2 - 1);
            y = (fR*2)*(gRandom->Rndm() * 2 - 1);
            z = (fR*2)*(gRandom->Rndm() * 2 - 1);
            Double_t r = TMath::Sqrt(x*x+y*y);
            Double_t theta = TMath::ATan2(r,z);
            Double_t R = TMath::Sqrt(x*x+y*y+z*z);
            Double_t Rtheta = fR;
            if (fF==5)
              Rtheta= fR + fW*TMath::Cos(theta)*TMath::Cos(theta);
            if (fF==7)
#ifdef HAVE_MATHMORE
              Rtheta = fR*(1+fBeta2*ROOT::Math::sph_legendre(2,0,theta)+fBeta4*ROOT::Math::sph_legendre(4,0,theta));
#else
            cerr << "Should not end here because you do not have libMathMore" << endl;
#endif
            Double_t prob = 1/(1+TMath::Exp((R-Rtheta)/fA));
            if (gRandom->Rndm()<prob) 
              nucleon_inside=1;
          }
        } else if ((fF==8) || (fF==14)) { // use TF2
          Double_t r;
          Double_t theta;
          fFunc3->GetRandom2(r,theta);
          Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
          x = r * TMath::Sin(phi) * TMath::Sin(theta);
          y = r * TMath::Cos(phi) * TMath::Sin(theta);
          z = r *                   TMath::Cos(theta);
        } else { // all other types
          TF1 *ff = fFunc1;
          if ((fFunc2) && (nucleon->GetType()==0))
            ff = fFunc2;
          if (fNodeDist<=0) { // "continuous" mode
            Double_t r = ff->GetRandom();
            Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
            Double_t ctheta = 2*gRandom->Rndm() - 1 ;
            Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
            x = r * stheta * TMath::Cos(phi);
            y = r * stheta * TMath::Sin(phi);
            z = r * ctheta;
          } else { // "grid/lattice" mode
            Int_t iNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t jNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t kNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t index=iNode*nslots*nslots+jNode*nslots+kNode;
            if (fIsUsed->TestBitNumber(index))
              continue;
            if (fLattice==1) {       // Primitive cubic system (PCS) -> https://en.wikipedia.org/wiki/Cubic_crystal_system
              x = fNodeDist*(iNode) - startingEdgeX;
              y = fNodeDist*(jNode) - startingEdgeY;
              z = fNodeDist*(kNode) - startingEdgeZ;
            } else if (fLattice==2) { //Body centered cubic (BCC) -> http://mathworld.wolfram.com/CubicClosePacking.html
              x = 0.5*fNodeDist*(-iNode+jNode+kNode) - 0.5*startingEdgeX;
              y = 0.5*fNodeDist*(+iNode-jNode+kNode) - 0.5*startingEdgeY;
              z = 0.5*fNodeDist*(+iNode+jNode-kNode) - 0.5*startingEdgeZ;
            } else if (fLattice==3) { //Face Centered Cubic (FCC) -> http://mathworld.wolfram.com/CubicClosePacking.html
              x = 0.5*fNodeDist*(jNode+kNode) - startingEdgeX;
              y = 0.5*fNodeDist*(iNode+kNode) - startingEdgeY;
              z = 0.5*fNodeDist*(iNode+jNode) - startingEdgeZ;
            } else {                  //Hexagonal close packing (HCP) -> https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres
              x = 0.5*fNodeDist*(2*iNode+((jNode+kNode)%2))          - startingEdgeX;
              y = 0.5*fNodeDist*(TMath::Sqrt(3)*(jNode+(kNode%2)/3)) - startingEdgeY;
              z = 0.5*fNodeDist*(kNode*2*TMath::Sqrt(6)/3)           - startingEdgeZ;
            }
            const Double_t r2 = x*x + y*y + z*z;
            const Double_t r  = TMath::Sqrt(r2);
	    if ((r>fMaxR)||(r2*gRandom->Rndm()>ff->Eval(r)))
	      continue;
            if (fSmearing>0.0) {
              Int_t nAttemptsToSmear = 0;
              while (1) {
                xsmeared = x*gRandom->Gaus(1.0,fSmearing);
                ysmeared = y*gRandom->Gaus(1.0,fSmearing);
                zsmeared = z*gRandom->Gaus(1.0,fSmearing);
                nAttemptsToSmear++;
                if (TestMinDist(i,xsmeared,ysmeared,zsmeared)) {
                  x = xsmeared;
                  y = ysmeared;
                  z = zsmeared;
                  break;
                }
                if (nAttemptsToSmear>=99) {
                  cerr << "Could not place on this node :: [" << x <<","<< y <<","<< z <<"] r = " << TMath::Sqrt(x*x+y*y+z*z) << " fm; "
                    << "Node (" << iNode << "," << jNode << "," << kNode << ") not smeared !!!" << endl;
                  ++fNonSmeared;
                  break;
                }
              }
            }
            fIsUsed->SetBitNumber(index);
          } /* end "grid/lattice mode" */
	}
	nucleon->SetXYZ(x,y,z);
	if (fF==5||fF==7||fF==8||fF==14) 
	  nucleon->RotateXYZ(fPhiRot,fThetaRot); // Uranium etc.
	if (fNodeDist>0) {
	  nucleon->RotateXYZ_3D(fXRot,fYRot,fZRot);
	  break;
	}
	if (TestMinDist(i,x,y,z))
	  break;
      }
    }
  }    

  // calculate center of mass
  Double_t sumx=0;       
  Double_t sumy=0;       
  Double_t sumz=0;       
  for (Int_t i = 0; i<fN; ++i) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
    sumx += nucleon->GetX();
    sumy += nucleon->GetY();
    sumz += nucleon->GetZ();
  }
  sumx = sumx/fN;
  sumy = sumy/fN;
  sumz = sumz/fN;
  static TVector3 finalShift;
  finalShift.SetXYZ(sumx,sumy,sumz);
  if (finalShift.Mag()>fSmax)
    goto cmscheck;
  Double_t fsumx = 0;
  Double_t fsumy = 0;
  Double_t fsumz = 0;
  if (fRecenter==1) {
    fsumx = sumx;
    fsumy = sumy;
    fsumz = sumz;
  } else if (fRecenter==2) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(fN-1));
    Double_t x = nucleon->GetX() - fN*sumx;
    Double_t y = nucleon->GetY() - fN*sumy;
    Double_t z = nucleon->GetZ() - fN*sumz;
    nucleon->SetXYZ(x,y,z);
  } else if ((fRecenter==3)||(fRecenter==4)) {
    TVector3 zVec;
    zVec.SetXYZ(0,0,1);
    TVector3 shiftVec;
    shiftVec.SetXYZ(sumx,sumy,sumz);
    TVector3 orthVec;
    orthVec = shiftVec.Cross(zVec);
    TRotation myRot;
    myRot.Rotate(shiftVec.Angle(zVec),orthVec);
    TVector3 myNuc;
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      myNuc.SetXYZ(nucleon->GetX(),nucleon->GetY(),nucleon->GetZ());
      myNuc.Transform(myRot);
      nucleon->SetXYZ(myNuc.X(), myNuc.Y(), myNuc.Z());
    }
    if (fRecenter==3)
      fsumz = shiftVec.Mag();
  }

  // recenter and shift
  for (Int_t i = 0; i<fN; ++i) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
    nucleon->SetXYZ(nucleon->GetX()-fsumx + xshift,
		    nucleon->GetY()-fsumy,
		    nucleon->GetZ()-fsumz);
  }

  return finalShift;
}

//---------------------------------------------------------------------------------
  ClassImp(TGlauberMC)
  ClassImp(TGlauberMC::Event)
//---------------------------------------------------------------------------------

TGlauberMC::TGlauberMC(const char* NA, const char* NB, Double_t xsect, Double_t xsectsigma) :
  fANucleus(NA),fBNucleus(NB),
  fXSect(xsect),fXSectOmega(0),fXSectLambda(0),fXSectEvent(0),
  fNucleonsA(0),fNucleonsB(0),fNucleons(0),
  fAN(0),fBN(0),fNt(0),
  fEvents(0),fTotalEvents(0),fBmin(0),fBmax(20),fHardFrac(0.65),
  fDetail(99),fCalcArea(0),fCalcLength(0), fDoCore(0),
  fMaxNpartFound(0),f2Cx(0),fPTot(0),fNNProf(0),
  fEv()
{
  if (xsectsigma>0) {
    fXSectOmega = xsectsigma;
    fXSectLambda = 1;
    fPTot = new TF1("fPTot","((x/[2])/(x/[2]+[0]))*exp(-(((x/[2])/[0]-1 )**2)/([1]*[1]))/[2]",0,300);
    fPTot->SetParameters(fXSect,fXSectOmega,fXSectLambda);
    fPTot->SetNpx(1000);
    fXSectLambda = fXSect/fPTot->GetHistogram()->GetMean();
    cout << "final lambda=" << fXSectLambda << endl;
    fPTot->SetParameters(fXSect,fXSectOmega,fXSectLambda);
    cout << "final <sigma>=" << fPTot->GetHistogram()->GetMean() << endl;
  }

  TString name(Form("Glauber_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
  TString title(Form("Glauber %s+%s Version",fANucleus.GetName(),fBNucleus.GetName()));
  SetName(name);
  SetTitle(title);
}

Bool_t TGlauberMC::CalcEvent(Double_t bgen)
{
  // calc next event
  if (!fNucleonsA) {
    fNucleonsA = fANucleus.GetNucleons();
    fAN = fANucleus.GetN();
    for (Int_t i = 0; i<fAN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      nucleonA->SetInNucleusA();
    }
  }
  if (!fNucleonsB) {
    fNucleonsB = fBNucleus.GetNucleons();
    fBN = fBNucleus.GetN();
    for (Int_t i = 0; i<fBN; ++i) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      nucleonB->SetInNucleusB();
    }
  }

  if (fPTot)
    fXSectEvent = fPTot->GetRandom();
  else 
    fXSectEvent = fXSect;

  // "ball" diameter = distance at which two balls interact
  Double_t d2 = (Double_t)fXSectEvent/(TMath::Pi()*10); // in fm^2
  Double_t bh = TMath::Sqrt(d2*fHardFrac);
  if (fNNProf) {
    Double_t xmin=0,xmax=0;
    fNNProf->GetRange(xmin,xmax);
    d2 = xmax*xmax;
  }

  fEv.Reset();
  memset(fBC,0,sizeof(Bool_t)*999*999);
  Int_t nc=0,nh=0;
  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Bool_t tB=nucleonB->GetType();
    for (Int_t j = 0; j<fAN; ++j) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(j));
      Double_t dx = nucleonB->GetX()-nucleonA->GetX();
      Double_t dy = nucleonB->GetY()-nucleonA->GetY();
      Double_t dij = dx*dx+dy*dy;
      if (dij>d2) 
        continue;
      Double_t bij = TMath::Sqrt(dij);
      if (fNNProf) {
        Double_t val = fNNProf->Eval(bij);
        Double_t ran = gRandom->Uniform();
        if (ran>val)
          continue;
      }
      nucleonB->Collide();
      nucleonA->Collide();
      fBC[i][j] = 1;
      fEv.BNN  += bij;
      ++nc;
      if (bij<bh)
        ++nh;
      Bool_t tA=nucleonA->GetType();
      if (tA!=tB)
        ++fEv.Ncollpn;
      else if (tA==1)
        ++fEv.Ncollpp;
      else
        ++fEv.Ncollnn;
      if (nc==1) {
        fEv.X0 = (nucleonA->GetX()+nucleonB->GetX())/2;
        fEv.Y0 = (nucleonA->GetY()+nucleonB->GetY())/2;
      }
    }
  }
  fEv.B = bgen;
  ++fTotalEvents;
  if (nc>0) {
    ++fEvents;
    fEv.Ncoll     = nc;
    fEv.Nhard     = nh;
    fEv.BNN      /= nc;
    return CalcResults(bgen);
  }
  return kFALSE;
}

Bool_t TGlauberMC::CalcResults(Double_t bgen)
{
  // calc results for the given event
  Double_t sumW=0;
  Double_t sumWA=0;
  Double_t sumWB=0;

  Double_t sinphi[10] = {0};
  Double_t cosphi[10] = {0};
  Double_t rn[10]     = {0};

  const Int_t kNc = fDoCore; // used later for core/corona

  for (Int_t i = 0; i<fAN; ++i) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    Double_t xA=nucleonA->GetX();
    Double_t yA=nucleonA->GetY();
    fEv.MeanXSystem  += xA;
    fEv.MeanYSystem  += yA;
    fEv.MeanXA  += xA;
    fEv.MeanYA  += yA;
    if (nucleonA->IsWounded()) {
      Double_t w = nucleonA->Get2CWeight(f2Cx);
      ++fEv.Npart;
      if (nucleonA->GetNColl()==1)
	++fEv.Npart0;
      ++fEv.NpartA;
      sumW   += w;
      sumWA  += w;
      fEv.MeanX  += xA * w;
      fEv.MeanY  += yA * w;
      fEv.MeanX2 += xA * xA * w;
      fEv.MeanY2 += yA * yA * w;
      fEv.MeanXY += xA * yA * w;
    }
  }

  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Double_t xB=nucleonB->GetX();
    Double_t yB=nucleonB->GetY();
    fEv.MeanXSystem  += xB;
    fEv.MeanYSystem  += yB;
    fEv.MeanXB  += xB;
    fEv.MeanYB  += yB;
    if (nucleonB->IsWounded()) {
      Double_t w = nucleonB->Get2CWeight(f2Cx);
      ++fEv.Npart;
      if (nucleonB->GetNColl()==1)
	++fEv.Npart0;
      ++fEv.NpartB;
      sumW   += w;
      sumWB  += w;
      fEv.MeanX  += xB * w;
      fEv.MeanY  += yB * w;
      fEv.MeanX2 += xB * xB * w;
      fEv.MeanY2 += yB * yB * w;
      fEv.MeanXY += xB * yB * w;
    }
  }
  if (fEv.Npart>0) {
    fEv.MeanX  /= sumW;
    fEv.MeanY  /= sumW;
    fEv.MeanX2 /= sumW;
    fEv.MeanY2 /= sumW;
    fEv.MeanXY /= sumW;
  } else {
    fEv.MeanX = 0;
    fEv.MeanY  = 0;
    fEv.MeanX2 = 0;
    fEv.MeanY2 = 0;
    fEv.MeanXY = 0;
  }

  if (fAN+fBN>0) {
    fEv.MeanXSystem /= (fAN + fBN);
    fEv.MeanYSystem /= (fAN + fBN);
  } else {
    fEv.MeanXSystem = 0;
    fEv.MeanYSystem = 0;
  }
  if (fAN>0) {
    fEv.MeanXA /= fAN;
    fEv.MeanYA /= fAN;
  } else {
    fEv.MeanXA = 0;
    fEv.MeanYA = 0;
  }
  if (fBN>0) {
    fEv.MeanXB /= fBN;
    fEv.MeanYB /= fBN;
  } else {
    fEv.MeanXB = 0;
    fEv.MeanYB = 0;
  }

  fEv.VarX  = fEv.MeanX2-(fEv.MeanX*fEv.MeanX);
  fEv.VarY  = fEv.MeanY2-(fEv.MeanY*fEv.MeanY);
  fEv.VarXY = fEv.MeanXY-fEv.MeanX*fEv.MeanY;
  Double_t tmpa = fEv.VarX*fEv.VarY-fEv.VarXY*fEv.VarXY;
  if (tmpa<0) 
    fEv.AreaW = -1;
  else 
    fEv.AreaW = TMath::Sqrt(tmpa);

  if (fEv.Npart>0) {
    // do full moments relative to meanX and meanY
    for (Int_t n = 1; n<10; ++n) {
      for (Int_t ia = 0; ia<fAN; ++ia) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(ia));
	if (nucleonA->GetNColl()<=kNc) 
	  continue;
	Double_t xA=nucleonA->GetX() - fEv.MeanX;
	Double_t yA=nucleonA->GetY() - fEv.MeanY;
	Double_t r = TMath::Sqrt(xA*xA+yA*yA);
	Double_t phi = TMath::ATan2(yA,xA);
	Double_t w = n;
	if (n==1) 
	  w = 3; // use r^3 weighting for Ecc1/Psi1
	Double_t rw = TMath::Power(r,w);
	cosphi[n] += rw*TMath::Cos(n*phi);
	sinphi[n] += rw*TMath::Sin(n*phi);
	rn[n] += rw;
      }
      for (Int_t ib = 0; ib<fBN; ++ib) {
        TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(ib));
	if (nucleonB->GetNColl()<=kNc) 
	  continue;
	Double_t xB=nucleonB->GetX() - fEv.MeanX;
	Double_t yB=nucleonB->GetY() - fEv.MeanY;
	Double_t r = TMath::Sqrt(xB*xB+yB*yB);
	Double_t phi = TMath::ATan2(yB,xB);
	Double_t w = n;
	if (n==1)
	  w = 3; // use r^3 weighting for Ecc1/Psi1
	Double_t rw = TMath::Power(r,w);
	cosphi[n] += rw*TMath::Cos(n*phi);
	sinphi[n] += rw*TMath::Sin(n*phi);
	rn[n] += rw;
      }
      cosphi[n] /= fEv.Npart;
      sinphi[n] /= fEv.Npart;
      rn[n] /= fEv.Npart;
      if (rn[n]>0) {
	fPsiN[n] = (TMath::ATan2(sinphi[n],cosphi[n]) + TMath::Pi())/n;
	fEccN[n] = TMath::Sqrt(sinphi[n]*sinphi[n]+cosphi[n]*cosphi[n])/rn[n];
      } else {
	fPsiN[n] = -1;
	fEccN[n] = -1;
      }
    }
    if (!kNc) { //silly test but useful to catch errors 
      Double_t t=TMath::Sqrt(TMath::Power(fEv.VarY-fEv.VarX,2)+4.*fEv.VarXY*fEv.VarXY)/(fEv.VarY+fEv.VarX)/fEccN[2];
      if (t<0.99||t>1.01)
        cout << "Error: Expected t=1 but found t=" << t << endl;
    }
  }

  fEv.B      = bgen;
  fEv.PhiA   = fANucleus.GetPhiRot();
  fEv.ThetaA = fANucleus.GetThetaRot();
  fEv.PhiB   = fBNucleus.GetPhiRot();
  fEv.ThetaB = fBNucleus.GetThetaRot();
  fEv.Psi1   = fPsiN[1];
  fEv.Ecc1   = fEccN[1];
  fEv.Psi2   = fPsiN[2];
  fEv.Ecc2   = fEccN[2];
  fEv.Psi3   = fPsiN[3];
  fEv.Ecc3   = fEccN[3];
  fEv.Psi4   = fPsiN[4];
  fEv.Ecc4   = fEccN[4];
  fEv.Psi5   = fPsiN[5];
  fEv.Ecc5   = fEccN[5];

  if (fCalcArea) {
    const Int_t nbins=200;
    const Double_t ell=10;
    const Double_t da=2*ell*2*ell/nbins/nbins;
    const Double_t d2 = (Double_t)fXSectEvent/(TMath::Pi()*10); // in fm^2
    const Double_t r2 = d2/4.;
    const Double_t mx = fEv.MeanX;
    const Double_t my = fEv.MeanY;
    TH2D areaA("hAreaA",";x (fm);y (fm)",nbins,-ell,ell,nbins,-ell,ell);
    TH2D areaB("hAreaB",";x (fm);y (fm)",nbins,-ell,ell,nbins,-ell,ell);
    for (Int_t i = 0; i<fAN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      if (!nucleonA->IsWounded())
        continue;
      if (nucleonA->GetNColl()==kNc)
        continue;
      Double_t x = nucleonA->GetX()-mx;
      Double_t y = nucleonA->GetY()-my;
      for (Int_t xi=1; xi<=nbins; ++xi) {
        for (Int_t yi=1; yi<=nbins; ++yi) {
          Int_t bin = areaA.GetBin(xi,yi);
          Double_t val=areaA.GetBinContent(bin);
          if (val>0)
            continue;
          Double_t dx=x-areaA.GetXaxis()->GetBinCenter(xi);
          Double_t dy=y-areaA.GetYaxis()->GetBinCenter(yi);
          if (dx*dx+dy*dy<r2)
            areaA.SetBinContent(bin,1);
        }
      }
    }
    for (Int_t i = 0; i<fBN; ++i) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      if (!nucleonB->IsWounded())
        continue;
      if (nucleonB->GetNColl()==kNc)
        continue;
      Double_t x = nucleonB->GetX()-mx;
      Double_t y = nucleonB->GetY()-my;
      for (Int_t xi=1; xi<=nbins; ++xi) {
        for (Int_t yi=1; yi<=nbins; ++yi) {
          Int_t bin = areaB.GetBin(xi,yi);
          Double_t val=areaB.GetBinContent(bin);
          if (val>0)
            continue;
          Double_t dx=x-areaB.GetXaxis()->GetBinCenter(xi);
          Double_t dy=y-areaB.GetYaxis()->GetBinCenter(yi);
          if (dx*dx+dy*dy<r2)
            areaB.SetBinContent(bin,1);
        }
      }
    }
    Double_t overlap1=0;
    Double_t overlap2=0;
    for (Int_t xi=1; xi<=nbins; ++xi) {
      for (Int_t yi=1; yi<=nbins; ++yi) {
        Int_t bin = areaA.GetBin(xi,yi);
        Double_t vA=areaA.GetBinContent(bin);
        Double_t vB=areaB.GetBinContent(bin);
        if (vA>0&&vB>0)
          ++overlap1;
        if (vA>0||vB>0)
          ++overlap2;
      }
    }
    fEv.AreaO = overlap1*da;
    fEv.AreaA = overlap2*da;
  }

  if (fCalcLength) {
    const Double_t krhs = TMath::Sqrt(fXSectEvent/40./TMath::Pi());
    const Double_t ksg  = krhs/TMath::Sqrt(5);
    const Double_t kDL  = 0.1;
    TF1 rad("rad","2*pi/[0]/[0]*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,5*ksg); 
    rad.SetParameter(0,ksg);
    const Double_t minval = rad.Eval(5*ksg);
    fEv.Phi0         = gRandom->Uniform(0,TMath::TwoPi());
    Double_t kcphi0  = TMath::Cos(fEv.Phi0);
    Double_t ksphi0  = TMath::Sin(fEv.Phi0);
    Double_t x       = fEv.X0;
    Double_t y       = fEv.Y0;
    Double_t i0a     = 0;
    Double_t i1a     = 0;
    Double_t l       = 0;
    Double_t val     = CalcDens(rad,x,y);
    while (val>minval) {
      x     += kDL * kcphi0;
      y     += kDL * ksphi0;
      i0a   += val;
      i1a   += l*val;
      l+=kDL;
      val    = CalcDens(rad,x,y);
    }
    fEv.Length = 2*i1a/i0a;
  }

  if (fEv.Npart > fMaxNpartFound) 
    fMaxNpartFound = fEv.Npart;

  return kTRUE;
}

Double_t TGlauberMC::CalcDens(TF1 &prof, Double_t xval, Double_t yval) const
{
  Double_t rmin=0,rmax=0;
  prof.GetRange(rmin,rmax);
  Double_t r2max = rmax*rmax;
  Double_t ret = 0;
  for (Int_t i = 0; i<fAN; ++i) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    if (!nucleonA->IsWounded())
      continue;
    Double_t x = nucleonA->GetX();
    Double_t y = nucleonA->GetY();
    Double_t r2=(xval-x)*(xval-x)+(yval-y)*(yval-y);
    if (r2>r2max)
      continue;
    ret += prof.Eval(TMath::Sqrt(r2));
  }
  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    if (!nucleonB->IsWounded())
      continue;
    Double_t x = nucleonB->GetX();
    Double_t y = nucleonB->GetY();
    Double_t r2=(xval-x)*(xval-x)+(yval-y)*(yval-y);
    if (r2>r2max)
      continue;
    ret += prof.Eval(TMath::Sqrt(r2));
  }
  return ret;
}

void TGlauberMC::Draw(Option_t* option)
{
  static TH2F *h2f = new TH2F("hGlauberMC",";x (fm);y(fm)",1,-18,18,1,-12,12);

  h2f->Reset();
  h2f->SetStats(0);
  h2f->Draw();

  TEllipse e;
  e.SetFillColor(0);
  e.SetFillStyle(0);
  e.SetLineColor(1);
  e.SetLineStyle(2);
  e.SetLineWidth(1);
  e.DrawEllipse(GetB()/2,0,fBNucleus.GetR(),fBNucleus.GetR(),0,360,0);
  e.DrawEllipse(-GetB()/2,0,fANucleus.GetR(),fANucleus.GetR(),0,360,0);
  fANucleus.Draw(fXSect, kMagenta, kYellow);
  fBNucleus.Draw(fXSect, kMagenta, kOrange);

  TString opt(option);
  if (opt.IsNull())
    return;

  Double_t sy2 = GetSy2();
  Double_t sx2 = GetSx2();
  Double_t phase = 0;
  if (sy2<sx2) {
    Double_t d = sx2;
    sx2 = sy2;
    sy2 = d;
    phase = TMath::Pi()/2.;
  }
  Double_t x1 = (0.5*(sy2-sx2)+TMath::Sqrt(TMath::Power(sy2-sx2,2.)-4*TMath::Power(GetSxy(),2)));
  Double_t ang = TMath::ATan2(-GetSxy(),x1)+phase;
  TLine l;
  l.SetLineWidth(3);
  l.DrawLine(-10*TMath::Cos(ang),-10*TMath::Sin(ang),10*TMath::Cos(ang),10*TMath::Sin(ang));
}

Double_t TGlauberMC::GetTotXSect() const
{
  return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBmax*fBmax/100;
}

Double_t TGlauberMC::GetTotXSectErr() const
{
  return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) * TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

TObjArray *TGlauberMC::GetNucleons() 
{
  if (!fNucleonsA || !fNucleonsB) return 0;
  if (fNucleons) return fNucleons;

  fNucleonsA->SetOwner(0);
  fNucleonsB->SetOwner(0);
  TObjArray *allnucleons=new TObjArray(fAN+fBN);
  allnucleons->SetOwner();
  for (Int_t i = 0; i<fAN; ++i) {
    allnucleons->Add(fNucleonsA->At(i));
  }
  for (Int_t i = 0; i<fBN; ++i) {
    allnucleons->Add(fNucleonsB->At(i));
  }
  fNucleons = allnucleons;
  return allnucleons;
}

Bool_t TGlauberMC::NextEvent(Double_t bgen)
{
  if (bgen<0) 
    bgen = TMath::Sqrt((fBmax*fBmax-fBmin*fBmin)*gRandom->Rndm()+fBmin*fBmin);

  fANucleus.ThrowNucleons(-bgen/2.);
  fBNucleus.ThrowNucleons(bgen/2.);

  return CalcEvent(bgen);
}

Bool_t TGlauberMC::ReadNextEvent(Bool_t calc, const char *fname)
{
  static TFile *inf = 0;
  static Int_t iev  = 0;
  if (fname) {
    cout << "ReadNextEvent: Setting up file " << fname << endl;
    delete inf;
    inf = TFile::Open(fname);
    if (!inf) 
      return 0;
    if (!fNucleonsA) {
      fANucleus.ThrowNucleons();
      fNucleonsA = fANucleus.GetNucleons();
      fAN = fANucleus.GetN();
      for (Int_t i = 0; i<fAN; ++i) {
	TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
	nucleonA->SetInNucleusA();
      }
    }
    if (!fNucleonsB) {
      fBNucleus.ThrowNucleons();
      fNucleonsB = fBNucleus.GetNucleons();
      fBN = fBNucleus.GetN();
      for (Int_t i = 0; i<fBN; ++i) {
	TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
	nucleonB->SetInNucleusB();
      }
    }
    if (calc)
      return 1;
    fNt = dynamic_cast<TNtuple*>(inf->Get(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName())));
    if (!fNt) {
      cerr << "ReadNextEvent: Could not find ntuple!" << endl;
      inf->ls();
      return 0;
    }
    fNt->SetBranchAddress("Npart",&fEv.Npart);
    fNt->SetBranchAddress("Ncoll",&fEv.Ncoll);
    fNt->SetBranchAddress("B",&fEv.B);
    fNt->SetBranchAddress("BNN",&fEv.BNN);
    fNt->SetBranchAddress("VarX",&fEv.VarX);
    fNt->SetBranchAddress("VarY",&fEv.VarY);
    fNt->SetBranchAddress("VarXY",&fEv.VarXY);
    fNt->SetBranchAddress("NpartA",&fEv.NpartA);
    fNt->SetBranchAddress("NpartB",&fEv.NpartB);
    fNt->SetBranchAddress("Npart0",&fEv.Npart0);
    fNt->SetBranchAddress("Psi1",&fEv.Psi1);
    fNt->SetBranchAddress("Ecc1",&fEv.Ecc1);
    fNt->SetBranchAddress("Psi2",&fEv.Psi2);
    fNt->SetBranchAddress("Ecc2",&fEv.Ecc2);
    fNt->SetBranchAddress("Psi3",&fEv.Psi3);
    fNt->SetBranchAddress("Ecc3",&fEv.Ecc3);
    fNt->SetBranchAddress("Psi4",&fEv.Psi4);
    fNt->SetBranchAddress("Ecc4",&fEv.Ecc4);
    fNt->SetBranchAddress("Psi5",&fEv.Psi5);
    fNt->SetBranchAddress("Ecc5",&fEv.Ecc5);
    return 1;
  }
  if ((!inf)||(!fNt&&!calc)) {
    cerr << "ReadNextEvent was not initialized" <<endl;
    return 0;
  }
  TObjArray *arr = dynamic_cast<TObjArray*>(inf->Get(Form("nucleonarray%d",iev)));
  if (!arr) {
    if (iev==0) {
      cerr << "ReadNextEvent could not read nucleon array for event " << iev << endl;
      return 0;
    }
    iev = 0;
    cerr << "ReadNextEvent resetting to first event" << endl;
    arr = dynamic_cast<TObjArray*>(inf->Get(Form("nucleonarray%d",iev)));
  }

  Double_t bgenA=0, bgenB=0;
  Int_t inA=0, inB=0;
  const Int_t nNucls = arr->GetEntries();
  for (Int_t iNucl=0; iNucl<nNucls; ++iNucl) {
    TGlauNucleon *nuclinp = static_cast<TGlauNucleon*>(arr->At(iNucl));
    TGlauNucleon *nuclout = 0;
    if (nuclinp->IsInNucleusB()) { 
      nuclout = static_cast<TGlauNucleon*>(fNucleonsB->At(inB));
      bgenB += nuclinp->GetX();
      ++inB;
    } else {
      nuclout = static_cast<TGlauNucleon*>(fNucleonsA->At(inA));
      bgenA += nuclinp->GetX();
      ++inA;
    }
    nuclout->Reset();
    nuclout->SetXYZ(nuclinp->GetX(),nuclinp->GetY(),nuclinp->GetZ());
    nuclout->SetType(nuclinp->GetType());
    nuclout->SetEnergy(nuclinp->GetEnergy());
    if (!calc)
      nuclout->SetNColl(nuclinp->GetNColl());
  }
  delete arr;
  Double_t bgen = bgenB/inB-bgenA/inA;
  if (calc) {
    Bool_t ret = CalcEvent(bgen);
    if (0) 
      cout << iev << ": " << fEv.B << " " << fEv.Npart << " " << fEv.Ncoll << " " << fEv.Npart0 << endl;
    ++iev;
    return ret;
  }
  Int_t ret = fNt->GetEntry(iev);
  if (ret<=0) 
    return 0;
  fEccN[1]=fEv.Ecc1;
  fEccN[2]=fEv.Ecc2;
  fEccN[3]=fEv.Ecc3;
  fEccN[4]=fEv.Ecc4;
  fEccN[5]=fEv.Ecc5;
  if (0) 
    cout << iev << ": " << fEv.B << " " << fEv.Npart << " " << fEv.Ncoll << " " << fEv.Npart0 << endl;
  if (0) { // test ntuple values vs re-calculated values
    Double_t npart = fEv.Npart;
    Double_t ncoll = fEv.Ncoll;
    Double_t ecc2  = fEv.Ecc2;
    CalcEvent(bgen);
    if (npart!=fEv.Npart) 
      cout << iev << " differ in npart " << npart << " " << fEv.Npart << endl;
    if (ncoll!=fEv.Ncoll) 
      cout << iev << " differ in ncoll " << ncoll << " " << fEv.Ncoll << endl;
    if (TMath::Abs(ecc2-fEv.Ecc2)>0.001) 
      cout << iev << " differ in ecc2 " << ecc2 << " " << fEv.Ecc2 << endl;
  }
  ++iev;
  return 1;
}

void TGlauberMC::Run(Int_t nevents, Double_t b)
{
  if (fNt == 0) {
    TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
    TString title(Form("%s + %s (x-sect = %.1f mb) str %s",fANucleus.GetName(),fBNucleus.GetName(),fXSect,Str()));
    TString vars("Npart:Ncoll:Nhard:B:BNN:Ncollpp:Ncollpn:Ncollnn:VarX:VarY:VarXY:NpartA:NpartB:Npart0:AreaW");
    if (fDetail>1)
      vars+=":Psi1:Ecc1:Psi2:Ecc2:Psi3:Ecc3:Psi4:Ecc4:Psi5:Ecc5";
    if (fDetail>2)
      vars+=":AreaO:AreaA:X0:Y0:Phi0:Length";
    if (fDetail>3)
      vars+=":MeanX:MeanY:MeanX2:MeanY2:MeanXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB";
    if (fDetail>4)
      vars+=":PhiA:ThetaA:PhiB:ThetaB";
    fNt = new TNtuple(name,title,vars);
    fNt->SetDirectory(0);
    TObjArray *l = fNt->GetListOfBranches();
    for (Int_t i=0; i<l->GetEntries(); ++i) {
      TBranch *br = dynamic_cast<TBranch*>(l->At(i));
      if (br)
        br->SetCompressionLevel(9);
    }
  }

  for (Int_t i = 0; i<nevents; ++i) {
    while (!NextEvent(b)) {}
    fNt->Fill((Float_t*)(&fEv.Npart));
    if ((i>0)&&(i%100)==0) 
      cout << "Event # " << i << " x-sect = " << GetTotXSect() << " +- " << GetTotXSectErr() << " b        \r" << flush;
  }
  if (nevents>99)
    cout << endl << "Done!" << endl;
}
#endif



