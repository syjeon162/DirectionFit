#include "dirFit.hh"
#include "TMath.h"
#include "TH3F.h"

#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "TString.h"
#include "TVector3.h"
#include <TVectorD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFrame.h>

using namespace std;

TString pmtLocation;
TString inputTxt;
TString outputTxt;
TString outputRoot;
TString pdf_filename; 
int detMass;
bool ifOutputRoot = false;
double trueOriginX = 0, trueOriginY = 0, trueOriginZ = 0;
double trueDirX = 0, trueDirY = 0, trueDirZ = -1;
int evtNum = 10000;
int branchSelection = 0;
bool all_light = false;
bool true_light = false;
float vertexSmear = 0;
int scinPct = 0;
bool digitize = false;
bool oTxt = true;
float sometime = 0;
bool do2Dpdf = false;
bool doCharge = false;
bool doCos = false;
bool doScan = false;
bool doFitVertex = false;
int upperEvt = evtNum;
int lowerEvt = 0;
int orient=0;
int specialConfig = 0;
int nbins = 10;
double zlowercut = 1e6;
bool perPMT = false;
int evtBase = 0;
bool externalPDF = false;
bool do3D = false;
bool perDir = false;
double timeInterval = 0.5;
double uppertime = 5;
int hitNumLimit = 1e6;
int npmt_bin = 252; 
int ntime_bin = 30;
int ndir_bin = 400;
float timeCorrection = 0;

void parseArguments(int argc, char**argv)
{
  //cout<<"inside parse"<<endl;
  for (int iarg=0; iarg<argc; iarg++)
  {
    //cout<<iarg<<endl;
    if (string( argv[iarg])=="-h" || string( argv[iarg])=="--help" )
    {
      cout << "**************************************" << endl;
      cout << "Macros run options:" << endl;
      cout << "   -h || --help      Print help info." << endl;
      cout << "   -i || --input     Input event file in txt format" << endl;
      cout << "   -p || --pmt       Input pmt location file in txt format" << endl;
      cout << "   -o || --outTxt    Output txt file" << endl;
      cout << "   -r || --outRoot   Output root file" << endl;
      cout << "   -m || --detMass   Detector mass" << endl;
      cout << "   -n || --evtNum    Event number" << endl;
      cout << "   -b || --branch    Branch Selection (event branch 1; mc branch 0)" << endl;
      cout << "   -s || --smearV    vertex smearing" << endl;
      cout << "   -a || --allLight  if all light included" << endl;
      cout << "   -t || --trueLight if only true light included" << endl;
      cout << "   -c || --scinPct   Scintillator concentration" << endl;
      cout << "   --upper   	    event number upper limit" << endl;
      cout << "   --lower   	    event number lower limit" << endl;
      cout << "   --sometime        Adding some time to prompt time cut" << endl;
      cout << "   --digitize        Using digitized simulation" <<endl;
      cout << "   --2Dpdf           Doing 2D PDFs" <<endl;
      cout << "   --orient          PDF building orientation, -999 means isotropic" <<endl;
      cout << "   --charge          Doing charge ?" <<endl;
      cout << "   --cos             PDF in cos ?" <<endl;
      cout << "   --scanning        Using space scanning instead of MINUIT?" <<endl;
      cout << "   --originZ         Z location of the source" <<endl;
      cout << "   --specialConfig   Special configuration? Any numbers not 1-4 means no special config." <<endl;
      cout << "   --nbins           Number of bins in PDF "<<endl;
      cout << "   --zlowercut       Z hit location lower cut? Ignore if don't want such a cut" <<endl;
      cout << "   --fitVertex       Fitting vertex at the same time?" <<endl;
      cout << "   --perPMT          Run per-PMT likelihood?" <<endl;
      cout << "   --externalPDF     Input PDF for fitting" <<endl;
      cout << "   --do3D            Do we want time slice? Always set this!" <<endl;
      cout << "   --perDir          Per-direction fit? This is the best so far." <<endl;
      cout << "   --timeInterval    What is the time step?" <<endl;
      cout << "   --uppertime       Upper limit for the residual time?" <<endl;
      cout << "   --hitNumLimit     Upper limit for the number of hits" <<endl;
      cout << "   --npmt            How many pmts in the detector? " <<endl;
      cout << "   --ntime           Number of bins for the time dimension" <<endl;
      cout << "   --ndir            Number of direction options" <<endl;
      cout << "   --timeCorrection     Time correction for all the residual time" <<endl;
      cout << "**************************************" << endl;
      exit(0);
    }


    else if (string( argv[iarg])=="-i" || string( argv[iarg])=="--input" )
    {
      iarg++;
      inputTxt = argv[iarg];
      if (inputTxt.Contains(".txt"))
      {
        cout << "Input txt file ok" << endl;
      }
      else {
        cerr << "input file must be a txt file!" << endl;
        exit(1);
      }
    }

    else if (string( argv[iarg])=="-p" || string( argv[iarg])=="--pmt" )
    {
      iarg++;
      pmtLocation = argv[iarg];
    }
    else if (string( argv[iarg])=="-o" || string( argv[iarg])=="--outTxt" )
    {
      iarg++;
      oTxt = true;
      outputTxt = argv[iarg];
    }
    else if (string( argv[iarg])=="-r" || string( argv[iarg])=="--outRoot" )
    {
      iarg++;
      outputRoot = argv[iarg];
      ifOutputRoot = true;
    }
    else if (string( argv[iarg])=="-m" || string( argv[iarg])=="--mass" )
    {
      iarg++;
      detMass = atoi(argv[iarg]);
    }
    else if (string( argv[iarg])=="-n" || string( argv[iarg])=="--evtNum" )
    {
      iarg++;
      evtNum = atoi(argv[iarg]);
    }    
    else if ( string( argv[iarg])=="--evtBase" )
    {
      iarg++;
      evtBase = atoi(argv[iarg]);
    }    
    else if (string( argv[iarg])=="-b" || string( argv[iarg])=="--branch" )
    {
      iarg++;
      branchSelection = atoi(argv[iarg]);
    }
    else if (string( argv[iarg])=="-s" || string( argv[iarg])=="--smearV" )
    {
      iarg++;
      vertexSmear = atof(argv[iarg]);
    }    
    else if (string( argv[iarg])=="-a" || string( argv[iarg])=="--allLight" )
    {
      all_light = true;
    }
    else if (string( argv[iarg])=="-t" || string( argv[iarg])=="--trueLight" )
    {
      true_light = true;
    }
    else if (string( argv[iarg])=="--upper" )
    {
      iarg++;
      upperEvt = atoi(argv[iarg]);
    }
    else if (string( argv[iarg])=="--lower" )
    {
      iarg++;
      lowerEvt = atoi(argv[iarg]);
    }
    else if (string( argv[iarg])=="-c" || string( argv[iarg])=="--scinPct" )
    {
      iarg++;
      scinPct = atoi(argv[iarg]);
    }
    else if (string( argv[iarg])=="--sometime" )
    {
      iarg++;
      sometime = atof(argv[iarg]);
    }
    else if (string( argv[iarg])=="--orient" )
    {
      iarg++;
      orient = atoi(argv[iarg]);
    }
    else if (string( argv[iarg])=="--digitize" )
    {
      digitize = true;
    }    
    else if (string( argv[iarg])=="--2Dpdf" )
    {
      do2Dpdf = true;
    }   
    else if (string( argv[iarg])=="--charge" )
    {
      doCharge = true;
    }
    else if (string( argv[iarg])=="--cos" )
    {
      doCos = true;
    }
    else if (string( argv[iarg])=="--scanning" )
    {
      doScan = true;
    }
    else if (string( argv[iarg])=="--originX" )
    {
      iarg++;
      trueOriginX = atof(argv[iarg]);
    }
    else if (string( argv[iarg])=="--originY" )
    {
      iarg++;
      trueOriginY = atof(argv[iarg]);
    }    
    else if (string( argv[iarg])=="--originZ" )
    {
      iarg++;
      trueOriginZ = atof(argv[iarg]);
    }
    else if (string( argv[iarg])=="--specialConfig" )
    {
      iarg++;
      specialConfig = atof(argv[iarg]);
    }
    else if (string( argv[iarg])=="--nbins" )
    {
      iarg++;
      nbins = atoi(argv[iarg]);
    }
    else if (string( argv[iarg])=="--zlowercut" )
    {
      iarg++;
      zlowercut = atof(argv[iarg]);
    }
    else if (string( argv[iarg])=="--fitVertex" )
    {
      doFitVertex = true;
    }
    else if (string( argv[iarg])=="--perPMT" )
    {
      perPMT = true;
    }
    else if (string( argv[iarg])=="--externalPDF" )
    {
      externalPDF = true;
      iarg++;
      pdf_filename = argv[iarg];
    }
    else if (string( argv[iarg])=="--do3D" )
    {
      do3D = true;
    }
    else if (string( argv[iarg])=="--perDir" )
    {
      perDir = true;
    }    
    else if (string( argv[iarg])=="--timeInterval" )
    {
      iarg++;
      timeInterval = atof(argv[iarg]);
    }

    else if (string( argv[iarg])=="--uppertime" )
    {
      iarg++;
      uppertime = atof(argv[iarg]);
    }
    
    else if (string( argv[iarg])=="--hitNumLimit" )
    {
      iarg++;
      hitNumLimit = atoi(argv[iarg]);
    }
    else if (string( argv[iarg])=="--npmt" )
    {
      iarg++;
      npmt_bin = atoi(argv[iarg]);
    }    
    else if (string( argv[iarg])=="--ntime" )
    {
      iarg++;
      ntime_bin = atoi(argv[iarg]);
    }
    else if (string( argv[iarg])=="--ndir" )
    {
      iarg++;
      ndir_bin = atoi(argv[iarg]);
    }
    else if (string( argv[iarg])=="--timeCorrection" )
    {
      iarg++;
      timeCorrection = atof(argv[iarg]);
    }
  }
}


int main(int argc, char**argv){
//int main(){

  parseArguments( argc, argv);
  cout<<"intput and output txt : "<<inputTxt.Data()<<" "<<outputTxt.Data()<<endl;
  if (ifOutputRoot) cout<<"output root file : "<<outputRoot.Data()<<endl;
  cout<<"pmt location and detector mass : "<<pmtLocation.Data()<<" "<<detMass<<endl;

  TH3F* h3 = new TH3F("","", 100,-1200,1200, 100, -1200,1200,100,-1200,1200);
  TH3F* h4 = new TH3F("","", 100,-1200,1200, 100, -1200,1200,100,-1200,1200);

  TH1F* hc = new TH1F("","",50,0,50);  
  TH1F* hhx[1000];
  TH1F* hhz[1000];
  for (int i=0;i<1000;i++){
    hhx[i] = new TH1F("","",20,0,1);
    hhz[i] = new TH1F("","",20,-1200,1200);
  }

  ifstream in;
  in.open(pmtLocation.Data());
  double pmtx[264]={};
  double pmty[264]={};
  double pmtz[264]={};

  int aa, ee;
  double bb,cc,dd;
  while (!in.eof()){
    in>>aa>>bb>>cc>>dd>>ee;
    pmtx[aa] = bb;
    pmty[aa] = cc;
    pmtz[aa] = dd;
  }

  TFile* outfile;
  TTree tree("tree","tree");
  double trueDir[3], recoDir[3];
  int nPMT; double PMTid[264][3]; double PMTtime[264]; double PMTangle[264];

  if (ifOutputRoot){
    outfile = TFile::Open( outputRoot.Data(), "recreate");
    tree.Branch("trueDir",&trueDir,"trueDir[3]/D");
    tree.Branch("recoDir",&recoDir,"recoDir[3]/D");
    tree.Branch("nPMT",&nPMT,"nPMT/I");
    tree.Branch("PMTid",&PMTid,"PMTid[264][3]/D");
    tree.Branch("PMTtime",&PMTtime,"PMTtime[264]/D");
    tree.Branch("PMTangle",&PMTangle,"PMTangle[264]/D");
  }

  std::vector<std::vector<double> > pmtlist;
  std::vector<std::vector<std::vector<double> > > fitlist(10000);
  ifstream in2;
  in2.open(inputTxt.Data());
  cout<<"reading file "<<inputTxt.Data()<<endl;

  int aaa,bbb;
  double ddd,fff,eee;
  double dddd,eeee,ffff;
  int cccc;
  double ccc;
  double ggg[10000],hhh[10000],mmm[10000];
  int current=0;
  int counter[10000]={};
  double iniX[10000] = {};
  double iniY[10000] = {};
  double iniZ[10000] = {};
  int eventTaker[10000] = {};
  int hitNumber[300][30] = {};

  TH1F* h_times = new TH1F("","",100,-5,5);
  TH1F* h_timec = new TH1F("","",100,-5,5);
  TH1F* h_pmtzs = new TH1F("","",30,-1500,1500);
  TH1F* h_pmtzc = new TH1F("","",30,-1500,1500);

  if (digitize && branchSelection == 0 ){
    cout<<"  --- be careful! You can't do digitization and true branch at the same time !! ---"<<endl;
    exit(0);
  }

  int count_top = 0;
  int count_side = 0;
  int count_bot = 0;
  while (!in2.eof() ){
      in2>>aaa>>bbb>>ccc>>ddd>>eee>>fff;
      if (bbb == -999){
        ggg[aaa] = ddd;
	hhh[aaa] = eee;
	mmm[aaa] = fff;
	continue;
      }

    double pmtxloc; double pmtyloc;
    pmtxloc = pmtx[bbb];
    pmtyloc = pmty[bbb];

    ccc = ccc - (sqrt((pmtxloc-trueOriginX)*(pmtxloc-trueOriginX)+(pmtyloc-trueOriginY)*(pmtyloc-trueOriginY) + (pmtz[bbb]-trueOriginZ)*(pmtz[bbb]-trueOriginZ))/200.) + timeCorrection;
    if (ddd == 0 ){
      h_timec->Fill(ccc);
      h_pmtzc->Fill(pmtz[bbb]);
    }
    else {
      h_times->Fill(ccc);
      h_pmtzs->Fill(pmtz[bbb]);    
    }

    if(aaa > evtNum || aaa< evtBase) continue;
    if (true_light && ddd != 0) continue;
    if (eee != branchSelection ) continue;
    if (true_light && eee == 1) continue;

    double uRadius = sqrt(pmtxloc*pmtxloc+pmtyloc*pmtyloc);
    if (specialConfig == 1) {

       if ( ! ( (uRadius > 243.84-10 && uRadius < 243.84 + 10) ||
		(uRadius > 487.68-10 && uRadius < 487.68 + 10) ||
		(uRadius > 731.52-10 && uRadius < 731.52 + 10)
		)) continue;  
    }
    if (specialConfig == 2) {

       if ( ! ( (uRadius > 134.112-10 && uRadius < 134.112 + 10 && pmtyloc<0)  ||
                (uRadius > 268.224-10 && uRadius < 268.224 + 10 && pmtyloc>=0) ||
                (uRadius > 402.336-10 && uRadius < 402.336 + 10 && pmtyloc<0)  ||
		(uRadius > 536.448-10 && uRadius < 536.448 + 10 && pmtyloc>=0) ||
                (uRadius > 670.560-10 && uRadius < 670.560 + 10 && pmtyloc<0)  ||
		(uRadius > 804.672-10 && uRadius < 804.672 + 10 && pmtyloc>=0)
                )) continue;
    }
    if (specialConfig == 3) {
       if ( ! ( (uRadius > 134.112-10 && uRadius < 134.112 + 10 && pmtyloc<0 && pmtxloc<=0)  || (uRadius > 184.112-10 && uRadius < 184.112 + 10 && pmtyloc<0 && pmtxloc>0) ||
                (uRadius > 268.224-10 && uRadius < 268.224 + 10 && pmtyloc>=0 && pmtxloc<=0) || (uRadius > 318.224-10 && uRadius < 318.224 + 10 && pmtyloc>=0 && pmtxloc>0) ||
                (uRadius > 402.336-10 && uRadius < 402.336 + 10 && pmtyloc<0 && pmtxloc<=0)  || (uRadius > 452.336-10 && uRadius < 452.336 + 10 && pmtyloc<0 && pmtxloc>0) ||
                (uRadius > 536.448-10 && uRadius < 536.448 + 10 && pmtyloc>=0 && pmtxloc<=0) || (uRadius > 586.448-10 && uRadius < 586.448 + 10 && pmtyloc>=0 && pmtxloc>0) ||
                (uRadius > 670.560-10 && uRadius < 670.560 + 10 && pmtyloc<0 && pmtxloc<=0)  || (uRadius > 720.560-10 && uRadius < 720.560 + 10 && pmtyloc<0 && pmtxloc>0) ||
                (uRadius > 804.672-10 && uRadius < 804.672 + 10 && pmtyloc>=0 && pmtxloc<=0) || (uRadius > 854.672-10 && uRadius < 854.672 + 10 && pmtyloc>=0 && pmtxloc>0) 
                )) continue;

    }
    if (specialConfig == 4){
    
       double doingNothing = 1;
    }
    // event branch eee == 1;  mc branch eee == 0
    std::vector<double> pmtloc;
    pmtloc.push_back(pmtxloc);
    pmtloc.push_back(pmtyloc);
    if (!perPMT && !perDir)
      pmtloc.push_back(pmtz[bbb]);
    else
      pmtloc.push_back(bbb);
    pmtloc.push_back(ccc);
    double locX = trueOriginX+ gRandom->Gaus(0, vertexSmear); 
    double locY = trueOriginY+ gRandom->Gaus(0, vertexSmear); 
    double locZ = trueOriginZ+ gRandom->Gaus(0, vertexSmear);

    double theta = -1; double phi = -1;

    // true location x, y, z
    pmtloc.push_back(locX); pmtloc.push_back(locY); pmtloc.push_back(locZ);
    // true time, theta, phi
    // The direction vector: TVector3 dvec (TMath::Cos(phi)*TMath::Sin(theta), TMath::Sin(phi)*TMath::Sin(theta), TMath::Cos(theta));
    // for direction (1,0,0) theta = TMath::Pi()/2.; phi = 0;  
    // for direction (0,1,0) theta = TMath::Pi()/2.; phi =  TMath::Pi()/2.;
    // for direction (0,0,-1) theta = TMath::Pi(); phi = 0;
    if (orient == 999) {
      pmtloc.push_back(0); 
      if (mmm[aaa]>0){
        pmtloc.push_back(TMath::ATan(sqrt(ggg[aaa]*ggg[aaa]+hhh[aaa]*hhh[aaa])/abs(mmm[aaa]))); 
	theta = TMath::ATan(sqrt(ggg[aaa]*ggg[aaa]+hhh[aaa]*hhh[aaa])/abs(mmm[aaa]));
      }
      else {
	pmtloc.push_back(TMath::Pi() - TMath::ATan(sqrt(ggg[aaa]*ggg[aaa]+hhh[aaa]*hhh[aaa])/abs(mmm[aaa])));
	theta = TMath::Pi() - TMath::ATan(sqrt(ggg[aaa]*ggg[aaa]+hhh[aaa]*hhh[aaa])/abs(mmm[aaa]));
      }
      if (ggg[aaa]>=0 && hhh[aaa]>0) {pmtloc.push_back(TMath::ATan(hhh[aaa]/ggg[aaa])); phi = TMath::ATan(hhh[aaa]/ggg[aaa]);}
      else if (ggg[aaa]<0 && hhh[aaa]>=0) {pmtloc.push_back(TMath::ATan(hhh[aaa]/ggg[aaa])+TMath::Pi()); phi = TMath::ATan(hhh[aaa]/ggg[aaa])+TMath::Pi();}
      else if (ggg[aaa]<0 && hhh[aaa]<0) {pmtloc.push_back(TMath::ATan(hhh[aaa]/ggg[aaa])+ TMath::Pi()); phi = TMath::ATan(hhh[aaa]/ggg[aaa])+ TMath::Pi();}
      else if (ggg[aaa]>=0 && hhh[aaa]<=0) {pmtloc.push_back(TMath::Pi()*2 -TMath::ATan(abs(hhh[aaa])/ggg[aaa])); phi = TMath::Pi()*2 -TMath::ATan(abs(hhh[aaa])/ggg[aaa]); }
    }
    else if (orient == 1) {pmtloc.push_back(0); pmtloc.push_back(TMath::Pi()/2.); pmtloc.push_back(TMath::Pi()/2.);}
    else if (orient == 2) {pmtloc.push_back(0); pmtloc.push_back(TMath::Pi()/2.); pmtloc.push_back(0);}
    else if (orient == 0) {pmtloc.push_back(0); pmtloc.push_back(TMath::Pi()); pmtloc.push_back(0);}
    else { cout<<"need to set up an orient "<<endl; exit(1);}
    pmtloc.push_back(fff);

    pmtlist.push_back(pmtloc);
    if (aaa>= lowerEvt && aaa< upperEvt ){//&& eee == 1){
    if (theta > -999){
      fitlist[aaa].push_back(pmtloc);
      cout<<"pmt id "<<pmtloc.at(2)<<endl;
      eventTaker[aaa] = 1;
    }
    if (iniX[aaa] == 0 && iniY[aaa] == 0 && iniZ[aaa] == 0)
      iniX[aaa] = pmtxloc; iniY[aaa] = pmtyloc; iniZ[aaa] = pmtz[bbb]; 
    }

    hitNumber[bbb][(int)(theta/0.2)] ++;
    pmtloc.clear();
    counter[aaa]++;
  }
  cout<<"event induced PMT loaded "<<endl;

  RooFitResult* res;
  char formula[10];
  dirFit * rep = new dirFit ("_rep");
  
  // 1 ton time cut 4.5, 2.4 ton time cut 5.5, 5 ton time cut 6.5
  cout<<"doing digitize ? "<<digitize<<"  adding sometime ? "<<sometime<<endl;

  if (!digitize){
    if (detMass == 1)
      //rep->SetPromptCut(4.5+sometime);
      rep->SetPromptCut(4.5+sometime/10.);
      //rep->SetPromptCut(5.5);
    if (detMass == 2)
      rep->SetPromptCut(5.5+sometime/10.);
    if (detMass == 3)
      rep->SetPromptCut(6+sometime/10.);
    if (detMass == 5 || detMass == 4)
      //rep->SetPromptCut(6.5+sometime);
      //rep->SetPromptCut(6.5+sometime/10.);
      rep->SetPromptCut(sometime);  
  }
  else{
      rep->SetPromptCut(1.01+sometime/10.);
  }

  rep->SetNbins(nbins);
  rep->SetIfDoCharge(doCharge);
  rep->SetIfDoCharge(doCos);
  rep->SetFitVertex(doFitVertex);
  rep->SetIfDo2dpdf(false);
  rep->SetDo3D(do3D);
  if (do2Dpdf) rep->SetIfDo2dpdf(true);
  if (perPMT) rep->SetPerPMT(true);
  else rep->SetPerPMT(false);
  if (perDir) rep->SetPerDir(true);
  else rep->SetPerDir(false);

  if (!perDir)
    ntime_bin = uppertime/timeInterval;

  cout<<"reading processing events"<<endl;
  std::vector<double> iniVertex(3);
  wbPDF* pdfs;
  std::vector<wbPDF*> pdfss(500);
  for (int ipmt = 0; ipmt< 500;ipmt++){
    pdfss[ipmt] = new wbPDF("_wbpdf");
    if (!perDir)  pdfss[ipmt]->SetPMTPDFBinning(nbins,0,3.14,nbins,0,6.28);
    else pdfss[ipmt]->SetDirPDFBinning(npmt_bin,0,npmt_bin,ntime_bin,-10,10);
  }
  cout<<"constructing 3D structure .. "<<endl;
  std::vector<std::vector<wbPDF*>> pdfsss(ntime_bin, std::vector<wbPDF*>(npmt_bin));
  if (perPMT && do3D && !perDir){
    for (int itime = 0; itime< pdfsss.size(); itime++){
      for (int ipmt = 0; ipmt< pdfsss[itime].size(); ipmt++){
        cout<<"time slice and pmt id "<<itime<<" "<<ipmt<<endl;
        pdfsss[itime][ipmt] = new wbPDF("_wbpdf");
        pdfsss[itime][ipmt]->SetPMTPDFBinning(nbins,0,3.14,nbins,0,6.28);
      }
    }
  }

  cout<<"extracting pdfs .. "<<endl;

  TH1F* timepdf;
  TH1F* thetapdf;
  TH2F* timeThetapdf;
  TH2F* pmtpdf[500];
  TH2F* pmtpdff[100][500];
  TH2F* dirpdf[500];
  if (!perPMT && !perDir){
    pdfs =  rep->Reading_Processing_Events(pmtlist,"pdf", iniVertex, do2Dpdf, doCharge, doCos );
    rep->SetPDFs(pdfs);
    cout<<"getting pdfs"<<endl;
    timepdf = pdfs->GetTimePDF();
    thetapdf = pdfs->GetThetaPDF();
    timeThetapdf = pdfs->GetTimeThetaPDF();
    cout<<"pdf obtained"<<endl;
  }
  else if(perDir){
    cout<<"reading perPMT PDFs .. "<<endl;
    rep->SetNbins_time(ntime_bin);
    rep->SetNdir(ndir_bin);
    rep->SetNbins_pmt(npmt_bin);
    if (!externalPDF){
      cout<<"doing a perDir pdf generation.."<<endl;
      pdfss = rep->Reading_Processing_Events_PerDir(pmtlist,"dirpdf", iniVertex, do2Dpdf, doCharge, doCos );
      cout<<"checking pdfss, size, and size of individual pdfss element:  "<<pdfss.size()<<" "<<pdfss[1]->GetDirPDF()->GetNbinsX()<<" "<<pdfss[1]->GetDirPDF()->GetNbinsY()<<endl;
      for (int iii=0;iii<pdfss[1]->GetDirPDF()->GetNbinsX(); iii++){
        for (int jjj=0;jjj<pdfss[1]->GetDirPDF()->GetNbinsY(); jjj++){
	   if (pdfss[368]->GetDirPDF()->GetBinContent(iii+1,jjj+1)>0) cout<<"in app2 iii, jjj, pdfss[1] element: "<<iii<<" "<<jjj<<" "<<pdfss[368]->GetDirPDF()->GetBinContent(iii+1,jjj+1)<<endl;
        }
      }
    }
    else{
        TFile fepdf(pdf_filename.Data());
        for (int iir=0;iir<ndir_bin;iir++){
          TH2F* hpdf = (TH2F*)fepdf.Get(Form("output_dir_%d",iir));
          pdfss[iir]->SetDirPDF(hpdf);
        }
        fepdf.Close();
	rep->SetDirPDFs(pdfss);
    }
    if (ifOutputRoot){
      for (int ipmt = 0; ipmt<pdfss.size() ; ipmt++){
        if (perDir){
          dirpdf[ipmt] = pdfss[ipmt]->GetDirPDF();
          dirpdf[ipmt]->Write(Form("output_dir_%d",ipmt));
        }
      }
    }
  } 
  else{
    cout<<"reading perPMT PDFs .. "<<endl;
    if (!externalPDF){
      if (!do3D){
        pdfss = rep->Reading_Processing_Events_PerPMT(pmtlist,"pmtpdf", iniVertex, do2Dpdf, doCharge, doCos );
      } 
      else{
	cout<<"doing a 3D pdf generation.."<<endl;
        rep->SetTimeInterval(timeInterval);
	pdfsss = rep->Reading_Processing_Events_PerPMT_timeSlice(pmtlist,"pmtpdf", iniVertex, do2Dpdf, doCharge, doCos, uppertime );
      }
    }
    else{
      if (!do3D){
        TFile fepdf(pdf_filename.Data());
        for (int iir=0;iir<npmt_bin;iir++){
          TH2F* hpdf = (TH2F*)fepdf.Get(Form("output_%d",iir));
          pdfss[iir]->SetPMTPDF(hpdf);
        }
        fepdf.Close();
      }
      else {
	rep->SetTimeInterval(timeInterval);      
        TFile fepdf(pdf_filename.Data());
        for (int iir=0;iir<npmt_bin;iir++){
	  for (int iit =0; iit< ntime_bin; iit++){
            TH2F* hpdf = (TH2F*)fepdf.Get(Form("output_%d_timeSlice_%d",iir,iit));
            pdfsss[iit][iir]->SetPMTPDF(hpdf);
	  }
        }
        fepdf.Close();      
      }
    }	
    if (!do3D && !perDir) rep->SetPDFs(pdfss);
    else rep->SetPDFS(pdfsss);
    std::vector<std::vector<wbPDF*>> set_pdfs = rep->GetPMTPDFS();
    cout<<"checking the setup pmt pdf info. "<<set_pdfs[0][0]->GetPMTPDF()->GetNbinsX()<<endl;
    cout<<"getting pmt pdfs"<<endl;
    for (int ipmt = 0; ipmt<pdfss.size() ; ipmt++){
      if (ifOutputRoot && !do3D && !perDir){ 
	pmtpdf[ipmt] = pdfss[ipmt]->GetPMTPDF();
	pmtpdf[ipmt]->Write(Form("output_%d",ipmt)); 
      }
      if (ifOutputRoot && do3D && !perDir){ 
	for (int itime= 0; itime< pdfsss.size();itime++){
	  pmtpdff[itime][ipmt] = pdfsss[itime][ipmt]->GetPMTPDF();
	  pmtpdff[itime][ipmt]->Write(Form("output_%d_timeSlice_%d",ipmt,itime)); 
  	}
      }
    }
    cout<<"pmt pdf obtained"<<endl;
  }
  if (ifOutputRoot){
    h_times->Write("h_times");
    h_pmtzs->Write("h_pmtzs");
    h_timec->Write("h_timec");
    h_pmtzc->Write("h_pmtzc");
    outfile->Write();
  }

  if (!externalPDF) { cout<<"just getting pdf, not performing fit .. "<<endl; exit(1);}

  for (Int_t aaa = lowerEvt;aaa< upperEvt; aaa++){

    if (eventTaker[aaa] == 0 ) continue;
    wbPDF* nothing =  rep->Reading_Processing_Events(fitlist[aaa],"event", iniVertex, do2Dpdf, doCharge, doCos );

    cout<<1<<endl;
    RooArgList list("list");
    list.add(*rep);
    sprintf(formula,"%s","@0");
    RooFormulaVar* fcn = new RooFormulaVar("fit","fit",formula,list);

    cout<<2<<endl;
    rep->getParVar(0)->setConstant(false);
    rep->getParVar(1)->setConstant(false);
    rep->getParVar(2)->setConstant(false);

    cout<<3<<endl;
    rep->getParVar(0)->setVal(iniX[aaa]/ sqrt(iniX[aaa]*iniX[aaa]+iniY[aaa]*iniY[aaa]+iniZ[aaa]*iniZ[aaa]));  
    rep->getParVar(1)->setVal(iniY[aaa]/ sqrt(iniX[aaa]*iniX[aaa]+iniY[aaa]*iniY[aaa]+iniZ[aaa]*iniZ[aaa]));
    rep->getParVar(2)->setVal(iniZ[aaa]/ sqrt(iniX[aaa]*iniX[aaa]+iniY[aaa]*iniY[aaa]+iniZ[aaa]*iniZ[aaa]));
    rep->ifScan(false);
    if (doScan){
      rep->ifScan(true);
      rep->getParVar(2)->setVal(-1);
      rep->getParVar(2)->setConstant(true);
    }
    rep->SetHitNumLimit(hitNumLimit);
    if (hitNumLimit < 1e5) {cout<<"--------------------- setting hit number limit !!! "<<endl; }

    cout<<4<<endl;
    ofstream out;
    //out.open(Form("output_result_timecut_finepdf_%dton_%dpct_full.txt",atoi(argv[1]),atoi(argv[2])),std::ofstream::out | std::ofstream::app);
    if (oTxt)
      out.open(outputTxt.Data(), std::ofstream::out | std::ofstream::app);
    cout<<"setting function"<<endl;
    double bestFit=1e9;
    double currX=1e9;
    double currY=1e9;
    double currZ=1e9;
    if (doScan){
      double currRes = 1e9;
      for (int xloop =1;xloop<20;xloop++){
        for (int yloop=0;yloop<20;yloop++){
          rep->getParVar(0)->setVal(xloop);
          rep->getParVar(1)->setVal(yloop);
	  rep->getParVar(0)->setConstant(true);
          rep->getParVar(1)->setConstant(true);
          double res = rep->evaluate();

          //cout<<"currX currY "<<xloop/50.-4<<" "<<yloop/50.-4<<"  res "<<res<<endl;
	  cout<<aaa<<" "<<xloop<<" "<<yloop<<" "<<res<<endl;
          if (oTxt)
	    out<<aaa<<" "<<xloop<<" "<<yloop<<" "<<res<<endl;
          if (res < currRes) {
            currRes = res;
            currX = xloop;
            currY = yloop;
          }
        }
      }
      bestFit = currRes;
    }

    else{
      RooMinuit m(*fcn);
      m.setStrategy(2);
      Double_t callsEDM[2] = {10500., 1.e-6};
      Int_t irf = 0;

      gMinuit->mnexcm("MIGRAD",callsEDM,2,irf);
      m.migrad();
      res = m.save();
      double bestFit = res->minNll();
      std::cout<<"fit status code is : "<< res->status()<<std::endl;
      cout<<"directional results: "<<rep->getParVar(0)->getVal()<<" "<<rep->getParVar(1)->getVal()<<" "<<rep->getParVar(2)->getVal()<<endl;
      currX = rep->getParVar(0)->getVal();
      currY = rep->getParVar(1)->getVal();
    }
    currZ = rep->getParVar(2)->getVal();

    if (ifOutputRoot){
      recoDir[0] = currX;
      recoDir[1] = currY;
      recoDir[2] = currZ;

      trueDir[0] = trueDirX;
      trueDir[1] = trueDirY;
      trueDir[2] = trueDirZ;
      nPMT = fitlist[aaa].size();
      for (int i = 0; i<fitlist[aaa].size(); i++){
        PMTid[i][0] = fitlist[aaa].at(i).at(0);
        PMTid[i][1] = fitlist[aaa].at(i).at(1);
        PMTid[i][2] = fitlist[aaa].at(i).at(2);
        PMTtime[i] = fitlist[aaa].at(i).at(3);
        TVector3 temp1  (PMTid[i][0] - trueOriginX, PMTid[i][1] - trueOriginY, PMTid[i][2] - trueOriginZ);
        TVector3 temp2  (trueDirX, trueDirY, trueDirZ);
        PMTangle[i] = temp1.Angle(temp2);
      }
      tree.Fill();
    }

    if (!doScan){
      if (abs(bestFit) > 0 && abs(rep->getParVar(0)->getVal())>0 ){
        cout<<"directional vector and bestFit: "<<currX<<" "<<currY<<" "<<currZ<<" "<<bestFit<<endl;
        if (oTxt)
          out<<currX<<" "<<currY<<" "<<currZ<<" "<<bestFit<<endl;
      }
    }
    else{
        cout<<"directional vector and bestFit: "<<currX<<" "<<currY<<" "<<currZ<<" "<<bestFit<<"   angle between true and reco. "<<endl;
        if (oTxt)
          out<<currX<<" "<<currY<<" "<<currZ<<" "<<bestFit<<endl;
    }
  }
  if (ifOutputRoot){
    timepdf->Write("timepdf");
    thetapdf->Write("thetapdf");
    tree.Write();
    outfile->Write();
  }
}
