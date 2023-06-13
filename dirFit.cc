#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

#include "dirFit.hh"
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

using namespace std;

dirFit ::dirFit (const char* name)
  : RooAbsReal(name,name)
{
  _pulls     = new RooListProxy("_pulls","_pulls",this);
  RooRealVar* Par1 = new RooRealVar("par1","par1",1,-100,100);
  RooRealVar* Par2 = new RooRealVar("par2","par2",1,-100,100);
  RooRealVar* Par3 = new RooRealVar("par3","par3",1,-10,10);
  RooRealVar* Par4 = new RooRealVar("par4","par4",0,-1000,1000);
  RooRealVar* Par5 = new RooRealVar("par5","par5",0,-1000,1000);
  RooRealVar* Par6 = new RooRealVar("par6","par6",0,-1000,1000); // 0.00238,-10,10);
  RooRealVar* Par7 = new RooRealVar("par7","par7",0,0,10);
  RooRealVar* Par8 = new RooRealVar("par8","par8",1,0.,100);
  RooRealVar* Par9 = new RooRealVar("par9","par9",1,0.,100);

  Par1->setConstant(false);
  Par2->setConstant(false);
  Par3->setConstant(false);
  Par4->setConstant(true);
  Par5->setConstant(true);
  Par6->setConstant(true);
  Par7->setConstant(true);
  Par8->setConstant(true);
  Par9->setConstant(true);

  _parlist.add(*Par1);
  _parlist.add(*Par2);
  _parlist.add(*Par3);
  _parlist.add(*Par4);
  _parlist.add(*Par5);
  _parlist.add(*Par6);
  _parlist.add(*Par7);
  _parlist.add(*Par8);
  _parlist.add(*Par9);
  _pulls->add(_parlist);

  this->addServerList(*_pulls);

}

dirFit ::~dirFit ()
{;}

void dirFit::ReadingEvents(std::vector<std::vector<double> > pmtlist ){

  wbEvent _wbEvent("wbevent");
  _wbEvent.SetHitList(pmtlist);
  _wbEvent.sample_photons(_wbEvent.GetHitList());
}

wbPDF* dirFit::Reading_Processing_Events(std::vector<std::vector<double> > pmtlist, std::string mode, std::vector<double> vertex, bool twoDpdf, bool doCharge, bool doCos ){
  cout<<"mode .. "<<mode.c_str()<<endl;
  if (mode == "pdf"){
    wbEvent _wbEvent("wbevent");
    cout<<"setting pmtlist"<<endl;
    _wbEvent.SetHitList(pmtlist);
    cout<<"sampling photons"<<endl;
    std::vector<wbEvent::wbHit> list =  _wbEvent.sample_photons(_wbEvent.GetHitList());
    cout<<"creating pdfs"<<endl;
    if (twoDpdf) return _wbEvent.create2DPDFs(list, _promptCut, 100, doCharge, doCos, _nbins);
    return _wbEvent.createPDFs(list, _promptCut, 100, doCharge, doCos, _nbins);
  }
  if (mode == "event"){
    //cout<<"setting up event"<<endl;
    _eventHitList = this->SetEvent(pmtlist);
    //cout<<"setting up vertex"<<endl;
    _vertex.push_back(vertex[0]);
    _vertex.push_back(vertex[1]);
    _vertex.push_back(vertex[2]);
    return NULL;
  }
}

std::vector<wbPDF*> dirFit::Reading_Processing_Events_PerPMT(std::vector<std::vector<double> > pmtlist, std::string mode, std::vector<double> vertex, bool twoDpdf, bool doCharge, bool doCos ){
  cout<<"mode .. "<<mode.c_str()<<endl;
  if (mode == "pmtpdf"){
    wbEvent _wbEvent("wbevent");
    cout<<"setting pmtlist"<<endl;
    _wbEvent.SetHitListPMT(pmtlist);
    cout<<"sampling photons"<<endl;
    std::vector<wbEvent::wbHit> list =  _wbEvent.sample_photons(_wbEvent.GetHitList());
    cout<<"creating pdfs"<<endl;
    //if (twoDpdf) return _wbEvent.create2DPDFs(list, _promptCut, 100, doCharge, doCos, _nbins);
    return _wbEvent.createPMTPDFs(list, _promptCut, 100, doCharge, doCos, _nbins);
  }
}

std::vector<std::vector<wbPDF*>> dirFit::Reading_Processing_Events_PerPMT_timeSlice(std::vector<std::vector<double> > pmtlist, std::string mode, std::vector<double> vertex, bool twoDpdf, bool doCharge, bool doCos, double uppertime ){
  cout<<"mode .. "<<mode.c_str()<<endl;
  if (mode == "pmtpdf"){
    wbEvent _wbEvent("wbevent");
    cout<<"setting pmtlist"<<endl;
    _wbEvent.SetHitListPMT(pmtlist);
    cout<<"sampling photons"<<endl;
    std::vector<wbEvent::wbHit> list =  _wbEvent.sample_photons(_wbEvent.GetHitList());
    cout<<"creating pdfs"<<endl;
    //if (twoDpdf) return _wbEvent.create2DPDFs(list, _promptCut, 100, doCharge, doCos, _nbins);
    std::vector<std::vector<wbPDF*>> pdfs;
    for (double iloop=0 ; iloop<uppertime; iloop += _timeInterval){
      pdfs.push_back( _wbEvent.createPMTPDFs(list, iloop*_timeInterval, (iloop+1)*_timeInterval, 100, doCharge, doCos, _nbins));
    }
    return pdfs;
  }
}

std::vector<wbPDF*> dirFit::Reading_Processing_Events_PerDir(std::vector<std::vector<double> > pmtlist, std::string mode, std::vector<double> vertex, bool twoDpdf, bool doCharge, bool doCos){
  cout<<"mode .. "<<mode.c_str()<<endl;
  if (mode == "dirpdf"){
    wbEvent _wbEvent("wbevent");
    cout<<"setting pmtlist"<<endl;
    _wbEvent.SetHitListPMT(pmtlist);
    cout<<"sampling photons"<<endl;
    std::vector<wbEvent::wbHit> list =  _wbEvent.sample_photons(_wbEvent.GetHitList());
    cout<<"creating pdfs"<<endl;
    //if (twoDpdf) return _wbEvent.create2DPDFs(list, _promptCut, 100, doCharge, doCos, _nbins);
    std::vector<wbPDF*> backing = _wbEvent.createDirPDFs(list, -10, 10, 100, doCharge, doCos, _nbin_pmt, _nbin_time);
    return backing;
  }
}

std::vector<wbPDF*> dirFit::Reading_external_pdfs(TString pdf_filename){
  TFile f(pdf_filename.Data());

  std::vector<wbPDF*> pdfs(500);
  //std::vector<wbPDF*> pdfs;

  for (int ii=0;ii<500;ii++){
    pdfs[ii] = new wbPDF("_wbpdf");
    //cout<<_nbins<<endl;
    pdfs[ii]->SetPMTPDFBinning(_nbins, 0, 3.14, _nbins, 0, 6.28);
  }
  for (int ii=0;ii<500;ii++){
    TH2F* hpdf = (TH2F*)f.Get(Form("output_%d",ii));
    pdfs[ii]->SetPMTPDF(hpdf);
  }
  f.Close();
  return (std::vector<wbPDF*>)pdfs;
}

Double_t dirFit ::getPar(int i) {
(((RooAbsReal*)_pulls->at(i))->getVal());
}

RooRealVar* dirFit ::getParVar(int i) {
return ((RooRealVar*)_pulls->at(i));
}

void dirFit :: setPull(TH1D* pullvecCV){
  pullCV = new TVectorD(11);
  for(Int_t i=0;i<11;i++){
    (*pullCV)[i] =  pullvecCV->GetBinContent(i+1);
  }
}

void dirFit :: setPullUnc(TH1D* pullvecUnc){
  pullUnc = new TVectorD(11);
  for(Int_t i=0;i<11;i++){
    (*pullUnc)[i] = pullvecUnc->GetBinContent(i+1);
  }
}

Double_t dirFit::getPullUnc(Int_t pN){
  return (*pullUnc)[pN];
}

Double_t dirFit::ExtraPull (RooListProxy* _pulls) const
{
  return 0;
}

Double_t dirFit::calPMTLikelihood_timeSlice( RooListProxy* _pulls, std::vector<wbEvent::wbHit> hitlist, std::vector<double> vertex, vector<vector<wbPDF*>> pdff, bool twodpdf, bool doCharge, bool doCos, int hitlimit ) const {

  double totll = 1;
  double theta = -999;
  double phi = -999;
  int pmtcount = 0;

  theta = ((RooAbsReal*)_pulls->at(0))->getVal();
  phi   = ((RooAbsReal*)_pulls->at(1))->getVal();
  cout<<"theta phi "<<theta<<" "<<phi<<endl;

  for (wbEvent::wbHit i: hitlist){
    //cout<<"time and prompt cut "<<i.t<<" "<<_promptCut<<endl;
    if (i.t > _promptCut) continue;
    int pmtid = i.pmtid;
    //cout<<"size of pdf list : "<<pdf.size()<<" while pmt id "<<pmtid<<" size of hits "<<hitlist.size()<<" "<<pdff[itime].at(pmtid)->GetPMTPDF()->GetNbinsX()<<endl;

    if (i.t<-2 || i.t> 0) continue;
    if (pmtcount >= hitlimit ) break;
    //for (int itime = 0;itime< pdff.size(); itime++){
      double ctime = (i.t+2)/_timeInterval;
      //cout<<i.t<<" "<<_timeInterval<<endl;
      int itime = ctime;
      cout<<"size of pdf list : "<<pdff.size()<<" while pmt id and itime "<<pmtid<<" "<<itime<<" size of hits "<<hitlist.size()<<" "<<pdff[itime][pmtid]->GetPMTPDF()->GetNbinsX()<<endl;
      TH2F* h1 = pdff[itime][pmtid]->GetPMTPDF();
      //if (h1->Integral() == 0 || !h1->Integral() ){ cout<<"2D per pmt pdf has not been set ! Exit ! "<<endl; exit(0);}
      if (h1->Integral() == 0 || !h1->Integral() ) continue;
      h1->Scale(20./h1->Integral());
      TGraph2D pdf_pmt(h1);

      double currll;
      double temp;
      if (_ifScan){
        if (doCharge) { temp = h1->GetBinContent( (int)(theta)+1, (int)(phi)+1)* i.charge ;}
        else {temp = h1->GetBinContent(  (int)(theta)+1, (int)(phi)+1 ) ;}
      }
      else {
        if (doCharge) { temp = h1->GetBinContent( (int)((theta)/(3.14/_nbins))+1, (int)(phi/(6.28/_nbins))+1)* i.charge ;}
        else {temp = h1->GetBinContent(  (int)((theta)/(3.14/_nbins))+1, (int)(phi/(6.28/_nbins))+1 ) ;}
      }
      if (temp< 1e-3 || temp> 1e2 || temp< 1e-3 || temp> 1e2) currll = 1e-3;
      else currll = temp ;
      totll *= currll;
      pmtcount ++;
    //}
  }
  return totll;
}

Double_t dirFit::calPMTLikelihood( RooListProxy* _pulls, std::vector<wbEvent::wbHit> hitlist, std::vector<double> vertex, vector<wbPDF*> pdf, bool twodpdf, bool doCharge, bool doCos, int hitlimit ) const{

  double totll = 1; 
  double theta = -999;
  double phi = -999;
  int pmtcount = 0;
/*  
  double dirx = ((RooAbsReal*)_pulls->at(0))->getVal();
  double diry = ((RooAbsReal*)_pulls->at(1))->getVal();
  double dirz = ((RooAbsReal*)_pulls->at(2))->getVal();

  //cout<<"proposed direction: "<<dirx<<" "<<diry<<" "<<dirz<<endl;  
  if (dirz>0){
    theta = TMath::ATan(sqrt(dirx*dirx+diry*diry)/abs(dirz));
  }
  else {
    theta = TMath::Pi() - TMath::ATan(sqrt(dirx*dirx+diry*diry)/abs(dirz));
  }
  if (dirx>=0 && diry>0) {phi = TMath::ATan(diry/dirx);}
  else if (dirx<0 && diry>=0) {phi = TMath::ATan(diry/dirx)+TMath::Pi();}
  else if (dirx<0 && diry<0) {phi = TMath::ATan(diry/dirx)+ TMath::Pi();}
  else if (dirx>=0 && diry<=0) {phi = TMath::Pi()*2 -TMath::ATan(abs(diry)/dirx); }
*/
  theta = ((RooAbsReal*)_pulls->at(0))->getVal();
  phi   = ((RooAbsReal*)_pulls->at(1))->getVal();
  cout<<"theta phi "<<theta<<" "<<phi<<endl;

  //cout<<"point 1 "<<pdf.at(3)->GetPMTPDF()->GetNbinsX()<<endl;
  for (wbEvent::wbHit i: hitlist){
    //cout<<"time and prompt cut "<<i.t<<" "<<_promptCut<<endl;
    if (i.t > _promptCut) continue;	  
    if (pmtcount >= hitlimit ) break;
    int pmtid = i.pmtid;

    //cout<<"size of pdf list : "<<pdf.size()<<" while pmt id "<<pmtid<<" size of hits "<<hitlist.size()<<" "<<pdf.at(pmtid)->GetPMTPDF()->GetNbinsX()<<endl;
    TH2F* h1 = pdf[pmtid]->GetPMTPDF();
    //cout<<4<<endl;
    //for (int iii=0;iii<h1->GetNbinsX();iii++){
      //for (int jjj=0; jjj<h1->GetNbinsY(); jjj++){
	//cout<<h1->GetBinContent(iii+1,jjj+1)<<endl;
      //}
    //}
    if (h1->Integral() == 0 || !h1->Integral() ){ cout<<"2D per pmt pdf has not been set ! Exit ! "<<endl; exit(0);}
    //cout<<h1->Integral()<<endl;
    h1->Scale(200./h1->Integral());
    //cout<<6<<endl;
    TGraph2D pdf_pmt(h1);

    double currll;
    double temp;
    if (_ifScan){
      if (doCharge) { temp = h1->GetBinContent( (int)(theta)+1, (int)(phi)+1)* i.charge ;}
      else {temp = h1->GetBinContent(  (int)(theta)+1, (int)(phi)+1 ) ;}
    }
    else {
      if (doCharge) { temp = h1->GetBinContent( (int)((theta)/(3.14/_nbins))+1, (int)(phi/(6.28/_nbins))+1)* i.charge ;}
      else {temp = h1->GetBinContent(  (int)((theta)/(3.14/_nbins))+1, (int)(phi/(6.28/_nbins))+1 ) ;}    
    }

    //cout<<2<<endl;
    //if (doCharge) { temp = pdf_pmt.Interpolate(i.t, theta)* i.charge ;}
    //else {temp = pdf_pmt.Interpolate( theta, phi ) ;}
    //for (int ii=0;ii<20;ii++){
      //for (int jj=0;jj<20;jj++){
        //cout<<ii<<" "<<jj<<" "<<h1->ProjectionX()->GetBinCenter(ii+1)<<" "<<h1->ProjectionY()->GetBinCenter(jj+1)<<" "<< h1->GetBinContent(ii+1,jj+1)<<endl;
      //}
    //}
    if (temp< 1e-3 || temp> 1e2 || temp< 1e-3 || temp> 1e2) currll = 1e-3;
    else currll = temp ;
    //cout<<"pmt id and likelihood "<<pmtid<<" "<<currll<<endl;
    totll *= currll;
    pmtcount ++;
  }  
  return totll;
}

Double_t dirFit::calDirLikelihood( RooListProxy* _pulls, std::vector<wbEvent::wbHit> hitlist, std::vector<double> vertex, vector<wbPDF*> pdf, bool twodpdf, bool doCharge, bool doCos ) const{

  double totll = 1;
  double theta = -999;
  double phi = -999;
  int pmtcount = 0;
  theta = ((RooAbsReal*)_pulls->at(0))->getVal();
  phi   = ((RooAbsReal*)_pulls->at(1))->getVal();
  //cout<<"Dir likelihood theta phi "<<theta<<" "<<phi<<endl;

  for (wbEvent::wbHit i: hitlist){
    //cout<<"time and prompt cut "<<i.t<<" "<<_promptCut<<endl;
    if (i.t > _promptCut) continue;
    int pmtid = i.pmtid;
    int hdir = 20*((int)(phi)) + (int)(theta) ;
    int itime = (i.t+10)/0.2;

    //cout<<"hdir, itime and pmt id : "<<hdir<<" "<<itime<<" "<<pmtid<<endl;
    TH2F* h1 = pdf[hdir]->GetDirPDF();
    //cout<<"DirPDF x/y sizes and content: "<<h1->GetNbinsX()<<" "<<h1->GetNbinsY()<<" "<<h1->GetBinContent((int)(pmtid)+1, (int)(itime)+1)<<endl;
    //if (h1->Integral() == 0 || !h1->Integral() ){ cout<<"2D per pmt pdf has not been set ! Exit ! "<<endl; exit(0);}
    if (h1->Integral() == 0 || !h1->Integral()) {cout<<"there is an empty pdf at (theta, phi) "<<theta<<" "<<phi<<endl; return 1e12;}
    h1->Scale(200./h1->Integral());
    TGraph2D pdf_pmt(h1);

    double user_weight = i.user_weight;
    if (pmtid>=231) user_weight = 0;

    double currll;
    double temp;
    if (_ifScan){
      if (doCharge) { temp = h1->GetBinContent( (int)(pmtid)+1, (int)(itime)+1)* i.charge * user_weight ;}
      else {temp = h1->GetBinContent(  (int)(pmtid)+1, (int)(itime)+1, user_weight ) ;}
    }
    else {
      if (doCharge) { temp = h1->GetBinContent( (int)(pmtid)+1, (int)(itime)+1)* i.charge * user_weight ;}
      else {temp = h1->GetBinContent(  (int)(pmtid)+1, (int)(itime)+1, user_weight ) ;}
    }

    if (temp< 1e-3 || temp> 1e2 || temp< 1e-3 || temp> 1e2) currll = 1e-3;
    else currll = temp ;
    //cout<<"pmt id and likelihood "<<pmtid<<" "<<currll<<endl;
    totll *= currll;
    pmtcount ++;
  }
  return totll;
}

Double_t dirFit::calLikelihood( RooListProxy* _pulls, std::vector<wbEvent::wbHit> hitlist, std::vector<double> vertex, wbPDF* pdf, bool twodpdf, bool doCharge, bool doCos ) const{
  //cout<<"variable values: "<<((RooAbsReal*)_pulls->at(0))->getVal()<<" "<<((RooAbsReal*)_pulls->at(1))->getVal()<<" "<<((RooAbsReal*)_pulls->at(2))->getVal()<<endl;

  if (_fitVertex){
    vertex[0] = vertex[0] + ((RooAbsReal*)_pulls->at(3))->getVal();
    vertex[1] = vertex[1] + ((RooAbsReal*)_pulls->at(4))->getVal();
    vertex[2] = vertex[2] + ((RooAbsReal*)_pulls->at(5))->getVal();
  }

  TVector3 dir (((RooAbsReal*)_pulls->at(0))->getVal(), ((RooAbsReal*)_pulls->at(1))->getVal(), ((RooAbsReal*)_pulls->at(2))->getVal());
  TVector3 dirPred = dir.Unit();
  double totll= 1;
  double c = 300/1.333; //mm/ns
  //cout<<"event hit list size "<<hitlist.size()<<endl;
  for (wbEvent::wbHit i: hitlist){
    double D = TMath::Sqrt(TMath::Power(i.x - vertex[0], 2) + TMath::Power(i.y - vertex[1], 2) + TMath::Power(i.z - vertex[2], 2));
    double tt = i.t - D/c;
    if (i.t > _promptCut) continue;
    //if (tt > _promptCut) continue;

    TVector3 pm (i.x - vertex[0], i.y - vertex[1], i.z - vertex[2]); 
    TVector3 pmt = pm.Unit();
    double theta;
    if (doCos)
      theta = TMath::Cos(dirPred.Angle(pmt));
    else
      theta = dirPred.Angle(pmt);

    if (! twodpdf){
      TH1F* h1 = pdf->GetTimePDF();
      TH1F* h2 = pdf->GetThetaPDF();
      h1->Scale(1./h1->Integral());
      h2->Scale(1./h2->Integral());
      TSpline3 pdf_time(h1);
      TSpline3 pdf_theta(h2);

      double currll;
      if (pdf_time.Eval(i.t)< 1e-2 || pdf_time.Eval(i.t)> 1e6 || pdf_theta.Eval(theta)< 1e-2 || pdf_theta.Eval(theta)> 1e6) currll = 1e-2; 
      else { 
	if (doCharge){
          //cout<<i.charge<<endl;
	  currll = pdf_time.Eval(i.t) * pdf_theta.Eval(theta) * i.charge * 10;
	}
	else{
	  //currll = pdf_time.Eval(i.t) * pdf_theta.Eval(theta);
	  currll = pdf_theta.Eval(theta);
	}
      }
      if (currll < 1e-2) currll = 1e-2;
      //else currll = pdf_time.Eval(tt) * pdf_theta.Eval(theta);
      //cout<<"pmt location "<<pmt.x()<<" "<<pmt.y()<<" "<<pmt.z()<<endl;
      //cout<<dirPred.Angle(pmt)<<endl;
      //cout<<"one pmt hit .. "<<" time and theta "<<i.t<<" "<<theta<<" likelihood of time "<<pdf_time.Eval(i.t)<<"  and that of theta "<<pdf_theta.Eval(theta)<<endl;
      totll *= currll;
    }
    else {
      TH2F* h1 = pdf->GetTimeThetaPDF();
      if (h1->Integral() == 0 ){ cout<<"2D pdf has not been set ! Exit ! "<<endl; exit(0);}
      h1->Scale(1./h1->Integral());
      TGraph2D pdf_timeTheta(h1);

      double currll;
      double temp;
      if (doCharge) { temp = pdf_timeTheta.Interpolate(i.t, theta)* i.charge ;}
      else {temp = pdf_timeTheta.Interpolate(i.t, theta) ;}
      //double temp = pdf_timeTheta.Interpolate(tt, theta) ;
      //double temp = h1->GetBinContent((int)(i.t/0.25)+1, (int)(theta/0.1)+1);
      if (temp< 1e-2 || temp> 1e4 || temp< 1e-2 || temp> 1e4) currll = 1e-2;
      //if (temp< 1e-9 || temp> 1e9 || temp< 1e-9 || temp> 1e9) currll = 1;
      else currll = temp * temp;
      //cout<<"pmt location "<<pmt.x()<<" "<<pmt.y()<<" "<<pmt.z()<<endl;
      //cout<<dirPred.Angle(pmt)<<endl;
      //cout<<"one pmt hit .. "<<" time and theta "<<i.t<<" "<<theta<<" likelihood of time "<<pdf_time.Eval(i.t)<<"  and that of theta "<<pdf_theta.Eval(theta)<<endl;
      totll *= currll;    
    }
  }
  return totll;
}

Double_t dirFit::evaluate() const
{
  cout<<"evaluating .. "<<endl;
  //double result =  this->directionMatching(_pulls, _addTime);
  double result;
  if (_perDir) {
    cout<<"evaluated result  "<<this->calDirLikelihood(_pulls, _eventHitList, _vertex, _dirpdf, _do2dpdf, _doCharge, _doCos)<<"  "<<-TMath::Log(this->calDirLikelihood(_pulls, _eventHitList, _vertex, _dirpdf, _do2dpdf, _doCharge, _doCos))<<endl;	  
    return -TMath::Log(this->calDirLikelihood(_pulls, _eventHitList, _vertex, _dirpdf, _do2dpdf, _doCharge, _doCos));
  }
  if (_perPMT){ 
    if (!_do3D)
      result = this->calPMTLikelihood(_pulls, _eventHitList, _vertex, _pmtpdf, _do2dpdf, _doCharge, _doCos, _hitNumLimit);
    else
      result = this->calPMTLikelihood_timeSlice(_pulls, _eventHitList, _vertex, _pmtpdfs, _do2dpdf, _doCharge, _doCos, _hitNumLimit);
  }
  else {	  
    result = this->calLikelihood(_pulls, _eventHitList, _vertex, _pdf, _do2dpdf, _doCharge, _doCos);
  }
  cout<<"evaluated result  "<<result<<" "<<-TMath::Log(result)<<endl;
  return -TMath::Log(result);
}

Double_t dirFit::CalPromptCut(std::vector<std::vector<double> > list){
  double temp1 = 1e9;
  double temp2;
  for (int i =0;i<list.size(); i++){ 
    if( list[i].at(3)< temp1){
      temp1 = list[i].at(3);
    } 
  } 
  return temp1;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   Below is not the way using PDFs. Obselete.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t dirFit::directionMatching ( RooListProxy* _pulls, bool addtime) const 
{
  //cout<<"directionMatching .. "<<endl;
  //cout<<((RooAbsReal*)_pulls->at(0))->getVal()<<" "<<((RooAbsReal*)_pulls->at(1))->getVal()<<" "<<((RooAbsReal*)_pulls->at(2))->getVal()<<endl;
  std::vector<std::vector<double> > angleList = angleSplit( ((RooAbsReal*)_pulls->at(0))->getVal(), ((RooAbsReal*)_pulls->at(1))->getVal(), ((RooAbsReal*)_pulls->at(2))->getVal(), 0.02); 
  std::vector<std::vector<double> > landingList = landingLocation(_x, _y, _z, angleList, _detR, _detZ);
  double result = matchingPMT(_pmtList, landingList , _x, _y, _z, addtime);
  return result;

}

Double_t dirFit::matchingPMT(std::vector<std::vector<double> > pmtList, std::vector<std::vector<double> > predList, double x, double y, double z, bool addTime) const{

  //cout<<"matching PMT .. "<<"size of predList "<<predList.size()<<endl;
  double totalSum = 0;
  for (int j=0; j< pmtList.size(); j++){
    double saveSum = 1e9;
    int bestPred = -1;
    //cout<<"pmt list number "<<j<<endl;
    for (int i=0;i< predList.size(); i++){
      //cout<<"data list "<<pmtList[j].at(0)<<" "<<pmtList[j].at(1)<<" "<<pmtList[j].at(2)<<endl;
      //cout<<"pred list "<<predList[i]. at(0)<<" "<<predList[i]. at(1)<<" "<<predList[i]. at(2)<<endl;      
      double sum = sqrt((predList[i]. at(0) - pmtList[j].at(0)) * (predList[i]. at(0) - pmtList[j].at(0))
		      + (predList[i]. at(1) - pmtList[j].at(1)) * (predList[i]. at(1) - pmtList[j].at(1))
		      + (predList[i]. at(2) - pmtList[j].at(2)) * (predList[i]. at(2) - pmtList[j].at(2)) );
      //cout<<"current sum "<<sum<<"  at "<<i<<endl;
      double timeToVertex = sqrt((predList[i]. at(0) - x) * (predList[i]. at(0) - x)
                      + (predList[i]. at(1) - y) * (predList[i]. at(1) - y)
                      + (predList[i]. at(2) - z) * (predList[i]. at(2) - z) ) /200. ;
      if (addTime){
        sum *= abs(timeToVertex - pmtList[j].at(3));
      }
      if (sum < saveSum ){
        saveSum = sum;
	bestPred = i;
      }
    }
    if (bestPred == -1 ) saveSum = 1e9;
    else{
      //cout<<"  ::best fit sum::  "<<saveSum<<"  at "<<bestPred<<endl;
      //cout<<"data list      "<<pmtList[j].at(0)<<" "<<pmtList[j].at(1)<<" "<<pmtList[j].at(2)<<endl;
      //cout<<"best pred list "<<predList[bestPred]. at(0)<<" "<<predList[bestPred]. at(1)<<" "<<predList[bestPred]. at(2)<<endl;
    }
    totalSum += saveSum;
  }
  for (int j=0; j< pmtList.size(); j++){
    //cout<<"data list      "<<pmtList[j].at(0)<<" "<<pmtList[j].at(1)<<" "<<pmtList[j].at(2)<<endl;
  }
  for (int i=0;i< predList.size(); i++){
    //cout<<"best pred list "<<predList[i]. at(0)<<" "<<predList[i]. at(1)<<" "<<predList[i]. at(2)<<endl;
  }


  //cout<<" --- best fit sum ---  "<<totalSum<<endl;
  return totalSum*0.0001; //* (1./predList.size());
}

std::vector<std::vector<double> > dirFit::angleSplit(double px, double py, double pz, double variation) const{

  //cout<<"angle splitting .."<<endl; 
  double ppx = px/sqrt(px*px + py*py + pz*pz);
  double ppy = py/sqrt(px*px + py*py + pz*pz);
  double ppz = pz/sqrt(px*px + py*py + pz*pz);

  //cout<<" ppx ppy ppz :  "<<ppx<<" "<<ppy<<" "<<ppz<<endl;
  std::vector<std::vector<double> > hh;
  for (int i=0;i<40;i++){
    for (int j=0;j<40;j++){
      double xlist = -i/20.+1. ;
      double ylist = -j/20.+1. ;
      if (sqrt(xlist*xlist + ylist*ylist)>1) continue;
      double zlist = 1 - xlist*xlist - ylist*ylist;

      // 0.755 is 41 degree Cherenkov cone
      if ( ppx * xlist + ppy * ylist + ppz * zlist < 0.755- variation || ppx * xlist + ppy * ylist + ppz * zlist > 0.755+ variation) continue;
      //cout<<"adding one pred. direction "<<xlist<<" "<<ylist<<" "<<zlist<<endl;
      std::vector<double> h;
      h.push_back(xlist);
      h.push_back(ylist);
      h.push_back(zlist);
      hh.push_back(h);
      h.clear();
    }
  }

  //std::vector<double> h;
  //h.push_back(0);
  //std::vector<std::vector<double> > hh;
  //hh.push_back(h);
  return hh;
}

std::vector<std::vector<double> > dirFit::landingLocation(double x, double y, double z, std::vector<std::vector<double> > dirList, double detR, double detZ) const{

  double stepSize = 2;
  double currx = x;
  double curry = y;
  double currz = z;
  std::vector<std::vector<double> > hh;
  for (size_t i = 0; i< dirList.size(); i++){

    std::vector<double> h;
    for (int j=0;j< sqrt((detR*2 )*(detR*2 ) + (detZ*2)*(detZ*2)); j++ ){
      currx = currx + dirList[i].at(0)/sqrt(dirList[i].at(0)*dirList[i].at(0)+dirList[i].at(1)*dirList[i].at(1)+dirList[i].at(2)*dirList[i].at(2))* stepSize;
      curry = curry + dirList[i].at(1)/sqrt(dirList[i].at(0)*dirList[i].at(0)+dirList[i].at(1)*dirList[i].at(1)+dirList[i].at(2)*dirList[i].at(2))* stepSize;
      currz = currz + dirList[i].at(2)/sqrt(dirList[i].at(0)*dirList[i].at(0)+dirList[i].at(1)*dirList[i].at(1)+dirList[i].at(2)*dirList[i].at(2))* stepSize;
      if (sqrt(currx* currx + curry*curry)> detR || abs(currz) > detZ){
        h.push_back(currx);
        h.push_back(curry);
        h.push_back(currz);
        hh.push_back(h);
	h.clear();
	currx = x; curry = y; currz = z;
        break;
      }
    }
  }
  //std::vector<double> h;
  //h.push_back(0);
  //std::vector<std::vector<double> > hh;
  //hh.push_back(h);
  return hh;

}


