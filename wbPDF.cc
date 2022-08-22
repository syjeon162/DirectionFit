#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

#include "TMath.h"
#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "TString.h"
#include <TVectorD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFrame.h>
#include <TFile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1D.h>
#include <TVector3.h>
#include "wbEvent.hh"

using namespace std;

wbPDF::wbPDF (const char* name)
{
    _timePDF = (TH1F*) new TH1F("","",30,-10,20);
    _thetaPDF = (TH1F*)new TH1F("","",30,-3.14,3.14);
    //_thetaPDF = (TH1F*)new TH1F("","",20,-1,1);

    _timeThetaPDF = (TH2F*) new TH2F("","",_timePDF->GetNbinsX(),_timePDF->GetBinLowEdge(1),_timePDF->GetBinLowEdge(_timePDF->GetNbinsX())+_timePDF->GetBinWidth(_timePDF->GetNbinsX()),          _thetaPDF->GetNbinsX(),_thetaPDF->GetBinLowEdge(1),_thetaPDF->GetBinLowEdge(_thetaPDF->GetNbinsX())+_thetaPDF->GetBinWidth(_thetaPDF->GetNbinsX()));

    _PMTPDF = (TH2F*) new TH2F("","",40,-3.14,3.14,40,-3.14,3.14);

    _dirPDF = (TH2F*) new TH2F("","", 500,0,500, 30, -2,4 );
}






