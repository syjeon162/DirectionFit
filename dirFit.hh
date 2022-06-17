/*
 *  Eos directional fit header file.
 *
 *  Author: Guang Yang
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <sstream>
#include <TList.h>

#include <TROOT.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TH1.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TRint.h>
#include <TH2.h>
#include <TFormula.h>
#include <TF1.h>
#include <TSpline.h>

#include <TF2.h>
#include <TMath.h>
#include <Math/DistFunc.h>
#include <TLine.h>
#include <TTree.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TVirtualFFT.h>
#include <TFoamIntegrand.h>
#include <TMatrixD.h>
#include <TVectorT.h>
#include <TDecompChol.h>

#include <RooFit.h>
#include "RooGlobalFunc.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooRandom.h"
#include <RooMsgService.h>
#include <RooHist.h>
#include <RooTrace.h>
#include <RooCategory.h>
#include "RooConstVar.h"
#include "RooBinning.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TGraph2D.h"

#include "RooFit.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "TMinuit.h"
#include <RooRealVar.h>
#include "wbEvent.hh"

using namespace RooFit;

  class dirFit : public RooAbsReal {

  public:

    dirFit (const char* name);

    dirFit (const dirFit & other, const char* name = 0): RooAbsReal(other,name) {};
    virtual TObject* clone(const char* newname) const {return new dirFit (*this, newname);};
    virtual ~dirFit () ;
    dirFit (const dirFit & dirFit );

    std::vector<std::vector<double> > angleSplit(double px, double py, double pz, double variation) const;
    std::vector<std::vector<double> > landingLocation(double x, double y, double z, std::vector<std::vector<double> > dirList, double detR, double detZ) const;
    Double_t directionMatching ( RooListProxy* _pulls, bool addtime) const;
    Double_t matchingPMT(std::vector<std::vector<double> > pmtList, std::vector<std::vector<double> > predList, double x, double y, double z, bool addtime) const;
    Double_t calLikelihood( RooListProxy* _pulls, std::vector<wbEvent::wbHit> hitlist, std::vector<double> vertex, wbPDF* pdf ) const;
    Double_t calLikelihood( RooListProxy* _pulls, std::vector<wbEvent::wbHit> hitlist, std::vector<double> vertex, wbPDF* pdf, bool twodpdf ) const;
    Double_t calLikelihood( RooListProxy* _pulls, std::vector<wbEvent::wbHit> hitlist, std::vector<double> vertex, wbPDF* pdf, bool twodpdf, bool doCharge, bool doCos ) const;
    Double_t calPMTLikelihood( RooListProxy* _pulls, std::vector<wbEvent::wbHit> hitlist, std::vector<double> vertex, std::vector<wbPDF*> pdf, bool twodpdf, bool doCharge, bool doCos ) const;

    Double_t CalPromptCut(std::vector<std::vector<double> > list); 

    void ReadingEvents(std::vector<std::vector<double> > pmtlist);
    wbPDF* Reading_Processing_Events(std::vector<std::vector<double> > pmtlist, std::string mode, std::vector<double> vertex, bool twoDpdf, bool doCharge, bool doCos ); 
    std::vector<wbPDF*> Reading_Processing_Events_PerPMT(std::vector<std::vector<double> > pmtlist, std::string mode, std::vector<double> vertex, bool twoDpdf, bool doCharge, bool doCos );
    std::vector<wbPDF*> Reading_external_pdfs(TString pdf_filename);
    void SetPDFs(wbPDF* pdf) {_pdf = pdf;}
    void SetPDFs(std::vector<wbPDF*> pdf) {_pmtpdf = pdf;}
    std::vector<wbPDF*> GetPMTPDFs() {return _pmtpdf;}

    std::vector<wbEvent::wbHit> SetEvent(std::vector<std::vector<double> > list){
      std::vector<wbEvent::wbHit> event;
      for (std::vector<double> i: list){
	wbEvent::wbHit hit;
        hit.x = i.at(0);
        hit.y = i.at(1);
	if (_perPMT)  hit.pmtid = i.at(2);
	else	hit.z = i.at(2);
        hit.t = i.at(3);
	hit.charge = i.at(10);
	event.push_back(hit);
      }
      return event;
    }
    void SetPromptCut(double t){ _promptCut = t;}
    void SetPromptCutList(std::vector<double> t){ _promptCutList = t;}

    void randomGo(int newIsotope, RooListProxy* _pulls);

    dirFit & operator=(const dirFit & rhs);

    RooFormulaVar* Chi2() ;

    Double_t FillEv(RooListProxy* _pulls) const;

    void Fillup(RooListProxy* _pulls);

    Double_t ExtraPull(RooListProxy* _pulls) const;

    void setPull(TH1D* pullvecCV) ;
    void setPullUnc(TH1D* pullvecUnc) ;
    Double_t getPullUnc(Int_t pN) ;

    RooRealVar Par1 ;
    RooRealVar Par2 ;
    RooRealVar Par3 ;
    RooRealVar Par4 ;
    RooRealVar Par5 ;
    RooRealVar Par6 ;
    RooRealVar Par7 ;
    RooRealVar Par8 ;
    RooRealVar Par9 ;

    Double_t _par1;
    Double_t _par2;
    Double_t _par3;
    Double_t _par4;
    Double_t _par5;
    Double_t _par6;
    Double_t _par7;
    Double_t _par8;
    Double_t _par9;

    Double_t getPar(int i) ;
    void DataSwitch(Bool_t dataSwitch=true) const;

    Bool_t getDataSwitch() const;

    RooRealVar* getParVar(int i) ;
    RooListProxy* getParVar() ;

    RooArgList _parlist;
    RooListProxy* _pulls;

    TVectorD* pullCV;
    TVectorD* pullUnc;
    wbPDF* _wbpdf;
    int _nbins;
    bool _do2dpdf;
    bool _doCharge;
    bool _doCos;
    bool _ifScan;
    bool _fitVertex;
    bool _perPMT;
    void SetIfDo2dpdf (bool aaa) { _do2dpdf = aaa;}
    void SetIfDoCharge (bool aaa) { _doCharge = aaa;}
    void SetIfDoCos (bool aaa) { _doCos = aaa;}
    void SetNbins (int nbins) {_nbins = nbins;}
    void SetPerPMT (bool aaa) {_perPMT = aaa;}
    void ifScan (bool aaa) {_ifScan = aaa;} 

    Double_t _x;
    Double_t _y;
    Double_t _z;
    Double_t _detR;
    Double_t _detZ;
    std::vector<std::vector<double> > _pmtList;
    bool _addTime;
    std::vector<wbEvent::wbHit> _eventHitList;
    std::vector<double> _vertex;
    wbPDF* _pdf ;
    std::vector<wbPDF*> _pmtpdf;
    double _promptCut;
    std::vector<double> _promptCutList;

    void SetVertex(double x, double y, double z){ _x = x; _y = y; _z = z;}
    void SetDetectorSize (double detR, double detZ) { _detR = detR; _detZ = detZ;}
    void SetPMTList (std::vector<std::vector<double> > pmtList) {_pmtList = pmtList;}
    void SetAddTime (bool addtime) {_addTime = addtime;}
    void SetFitVertex (bool fitvertex) {_fitVertex = fitvertex;}

    std::vector<double> GetVertex() { std::vector<double> pmt; pmt.push_back(_x); pmt.push_back(_y); pmt.push_back(_z); return pmt; }    
    std::vector<double> GetDetectorSize() {std::vector<double> detSize; detSize.push_back(_detR); detSize.push_back(_detZ); return detSize;}
    std::vector<std::vector<double> > GetPMTList() {return _pmtList;}

    RooListProxy* getPullList() const { return _pulls; }

    TString fileName;
    double fissionFraction[100];
    double binEdge[100];
    Int_t  _nBins;
    TString singleExp;

    bool inSyst;
    bool equalIso;
    TString fileLocation;

    std::vector<TString> modelList;

    std::vector<TH1D*> dataList;

    //  private:

   virtual  Double_t evaluate() const;

   protected:
  };

