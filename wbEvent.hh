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

#include "RooFit.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "TMinuit.h"
#include <RooRealVar.h>

using namespace std;

  class wbPDF {

  public:

    wbPDF(const char* name);
    TH1F* _timePDF;
    TH1F* _thetaPDF;
    TH2F* _timeThetaPDF;
    TH2F* _PMTPDF;
    TH2F* _dirPDF;

    void SetTimePDFBinning (int nbin, double low, double high){ cout<<"set pdf binning "<< nbin<<endl; cout<<"original bin n. "<<_timePDF->TH1F::GetNbinsX()<<endl; _timePDF->TH1F::SetBins(nbin, low, high); cout<<"binning set"<<endl; }
    void SetThetaPDFBinning (int nbin, double low, double high){ _thetaPDF->TH1F::SetBins(nbin, low, high); }
    void SetTimeThetaPDFBinning (int nbin, double low, double high, int nbin2, double low2, double high2){ _timeThetaPDF->TH2F::SetBins(nbin, low, high, nbin2, low2, high2); }
    void SetPMTPDFBinning (int nbin, double low, double high, int nbin2, double low2, double high2){ _PMTPDF->TH2F::SetBins(nbin, low, high, nbin2, low2, high2); }
    void SetDirPDFBinning (int nbin, double low, double high, int nbin2, double low2, double high2){ _dirPDF->TH2F::SetBins(nbin, low, high, nbin2, low2, high2); }

    void SetPMTPDF(std::vector<std::vector<double>> a) {
	    //cout<<3456<<" "<<a.size()<<" "<<a.at(0).size()<<endl;
            if (a.size() != _PMTPDF->GetNbinsX() || a.at(0).size() != _PMTPDF->GetNbinsY()) {cout<<"setting a pdf with incorrect binning! "<<endl; exit(1);}
            for (int i=0;i<a.size();i++){
                   for (int j=0;j<a.at(0).size();j++)
                   {
                           //cout<<"setting bin i j "<<i<<" "<<j<<" "<<a.at(i).at(j)<<endl;
			   //cout<<_PMTPDF->GetNbinsX()<<endl;
                           _PMTPDF ->SetBinContent( i+1,j+1, a.at(i).at(j));
			   //cout<<3457<<endl;
                   }
            }
    }
    void SetDirPDF(std::vector<std::vector<double>> a) {
            if (a.size() != _dirPDF->GetNbinsX() || a.at(0).size() != _dirPDF->GetNbinsY()) {cout<<"setting a pdf with incorrect binning! "<<endl; cout<<"binning from get() "<<a.at(0).size()<<" and to-be-set binning "<<_dirPDF->GetNbinsY()<<endl; exit(1);}
            for (int i=0;i<a.size();i++){
                   for (int j=0;j<a.at(0).size();j++)
                   {
                           _dirPDF ->SetBinContent( i+1,j+1, a.at(i).at(j));
                   }
            }
    }

    void SetPMTPDF(TH2F* a) { 
            if (a->GetNbinsX() != _PMTPDF->GetNbinsX() || a->GetNbinsY() != _PMTPDF->GetNbinsY()) {cout<<"setting a pdf with incorrect binning! "<<endl; exit(1);}
            for (int i=0;i<a->GetNbinsX();i++){
                   for (int j=0;j<a->GetNbinsY();j++)
                   {
                           _PMTPDF ->SetBinContent( i+1,j+1, a->GetBinContent(i+1,j+1));
                   }
            }    
    }
    void SetDirPDF(TH2F* a) {
            if (a->GetNbinsX() != _dirPDF->GetNbinsX() || a->GetNbinsY() != _dirPDF->GetNbinsY()) {cout<<"setting a pdf with incorrect binning! "<<endl; cout<<" y binning from get() "<<a->GetNbinsY()<<" and to-be-set binning "<<_dirPDF->GetNbinsY()<<" "<<" x binning from get() "<<a->GetNbinsX()<<" and to-be-set binning "<<_dirPDF->GetNbinsX()<<endl; exit(1);}
            for (int i=0;i<a->GetNbinsX();i++){
                   for (int j=0;j<a->GetNbinsY();j++)
                   {
                           _dirPDF ->SetBinContent( i+1,j+1, a->GetBinContent(i+1,j+1));
                   }
            }
    }

    TH1F* GetTimePDF(){ return _timePDF; }
    TH1F* GetThetaPDF(){ return _thetaPDF; }
    TH2F* GetTimeThetaPDF(){ return _timeThetaPDF; }
    TH2F* GetPMTPDF(){ return _PMTPDF; }
    TH2F* GetDirPDF(){ return _dirPDF; }
  };


  class wbEvent {

  public:

    wbEvent (const char* name);

    //wbEvent (const wbEvent & other, const char* name = 0): RooAbsReal(other,name) {};
    //virtual TObject* clone(const char* newname) const {return new wbEvent (*this, newname);};
    virtual ~wbEvent () ;
    //wbEvent (const wbEvent & wbEvent );

    struct wbHit {
      double x; double y; double z; double t; double wavelength; double type;
      double px; double py; double pz; double trTime;
      double trX; double trY; double trZ; double trTheta; double trPhi;
      double charge; int pmtid;
      std::string something;
    };

    std::vector<wbHit> _hitlist;

    void SetHitList(std::vector<std::vector<double> > num);
    void SetHitListPMT(std::vector<std::vector<double> > num);
    std::vector<wbHit > GetHitList(){ return _hitlist; }

    std::vector<wbHit> sample_photons(  std::vector<wbHit> list, bool smearVtx = false, bool smearTime = false, bool smearQE = false, bool addNoise =false );

    wbPDF* createPDFs(std::vector<wbHit> list, std::vector<std::vector<double> > tPara, double prompt_cut, double wavelength_cut);
    wbPDF* createPDFs(std::vector<wbHit> list, double prompt_cut, double wavelength_cut, bool doCharge, bool doCos);
    wbPDF* createPDFs(std::vector<wbHit> list, double prompt_cut, double wavelength_cut, bool doCharge, bool doCos, int nbins);
    std::vector<wbPDF*> createPMTPDFs(std::vector<wbHit> list, double prompt_cut, double wavelength_cut, bool doCharge, bool doCos, int nbins);
    std::vector<wbPDF*> createPMTPDFs(std::vector<wbHit> list, double time1, double time2, double wavelength_cut, bool doCharge, bool doCos, int nbins);
    std::vector<wbPDF*> createDirPDFs(std::vector<wbHit> list, double time1, double time2, double wavelength_cut, bool doCharge, bool doCos, int nbin_pmt, int nbin_time);

    wbPDF* create2DPDFs(std::vector<wbHit> list, std::vector<std::vector<double> > tPara, double prompt_cut, double wavelength_cut);
    wbPDF* create2DPDFs(std::vector<wbHit> list, double prompt_cut, double wavelength_cut, bool doCharge, bool doCos);
    wbPDF* create2DPDFs(std::vector<wbHit> list, double prompt_cut, double wavelength_cut, bool doCharge, bool doCos, int nbins);
    private:
    protected:
  };


