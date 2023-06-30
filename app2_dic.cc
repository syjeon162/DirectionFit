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
int detMass = 4;
bool ifOutputRoot = false;
double trueOriginX = 0, trueOriginY = 0, trueOriginZ = 0;
double trueDirX = 0, trueDirY = 0, trueDirZ = -1;
int evtNum = 10000;
int branchSelection = 0;
bool inNtuple = true;
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
int orient=999;
int specialConfig = 0;
int nbins = 10;
double zlowercut = 1e6;
bool perPMT = false;
int evtBase = 0;
bool externalPDF = false;
bool do3D = false;
bool perDir = true;
double timeInterval = 0.5;
double uppertime = 5;
int hitNumLimit = 1e6;
int npmt_bin = 231; 
int ntime_bin = 30;
int ndir_bin = 400;
float timeCorrection = 0;
bool debug = false;
bool useDic = false;
bool useSource = false;

double mcx;
double mcy;
double mcz;
double mcu;
double mcv;
double mcw;
double mcke;
int evid;
int subev;
long nanotime;
int mcpcount;
int mcpecount;
int mcnhits;
int nentries;
std::vector<int>* pdgcodes;
std::vector<double>* mcKEnergies;
std::vector<double>* mcPosx;
std::vector<double>* mcPosy;
std::vector<double>* mcPosz;
std::vector<double>* mcDirx;
std::vector<double>* mcDiry;
std::vector<double>* mcDirz;
std::vector<int>*  hitPMTID;
std::vector<double>* hitPMTTime;
std::vector<double>* hitPMTDigitizedTime;
std::vector<double>* hitPMTCharge;
std::vector<double>* hitPMTDigitizedCharge;
std::vector<int>* mcPMTID;
std::vector<int>* mcPEIndex;
std::vector<double>* mcPETime;
std::vector<int>* mcPEProcess;
std::vector<int>* pmtId;
std::vector<double>* pmtX;
std::vector<double>* pmtY;
std::vector<double>* pmtZ;
TTree* ttree;

void parseArguments(int argc, char**argv)
{
  cout<<"inside parse"<<endl;
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
      cout << "   --input_ntuple    input format is ntuple" <<endl;
      cout << "   --dic             dichrocon" <<endl;
      cout << "**************************************" << endl;
      exit(0);
    }


    else if (string( argv[iarg])=="-i" || string( argv[iarg])=="--input" )
    {
      iarg++;
      inputTxt = argv[iarg];
      if (inputTxt.Contains(".txt") || inputTxt.Contains(".root"))
      {
        cout << "Input txt file ok" << endl;
      }
      else {
        cerr << "input file must be a txt or ntuple file!" << endl;
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
    else if (string( argv[iarg]) == "--debug" ){
      debug = true;
    }
    else if (string( argv[iarg]) == "--dic" ){
      useDic = true;
    }
    else if (string( argv[iarg]) == "--source" ){
      useSource = true;
    }
  }
}

void read_ntuple(std::string str){
  TFile ftemp(str.c_str());
  ttree = (TTree*)ftemp.Get("output");
  nentries = ttree->GetEntries();
  //ttree->Print();
  ttree->SetBranchAddress("mcx", &mcx);
  ttree->SetBranchAddress("mcy", &mcy);
  ttree->SetBranchAddress("mcz", &mcz);
  ttree->SetBranchAddress("mcu", &mcu);
  ttree->SetBranchAddress("mcv", &mcv);
  ttree->SetBranchAddress("mcw", &mcw);
  ttree->SetBranchAddress("mcke", &mcke);
  ttree->SetBranchAddress("evid", &evid);
  ttree->SetBranchAddress("subev", &subev);
  ttree->SetBranchAddress("nanotime", &nanotime);
  ttree->SetBranchAddress("mcpcount", &mcpcount);
  ttree->SetBranchAddress("mcpecount", &mcpecount);
  ttree->SetBranchAddress("mcnhits", &mcnhits);
  ttree->SetBranchAddress("pdgcodes", &pdgcodes);
  ttree->SetBranchAddress("mcKEnergies", &mcKEnergies);
  ttree->SetBranchAddress("mcPosx", &mcPosx);
  ttree->SetBranchAddress("mcPosy", &mcPosy);
  ttree->SetBranchAddress("mcPosz", &mcPosz);
  ttree->SetBranchAddress("mcDirx", &mcDirx);
  ttree->SetBranchAddress("mcDiry", &mcDiry);
  ttree->SetBranchAddress("mcDirz", &mcDirz);
  ttree->SetBranchAddress("hitPMTID", &hitPMTID);
  ttree->SetBranchAddress("hitPMTTime", &hitPMTTime);
  ttree->SetBranchAddress("hitPMTDigitizedTime", &hitPMTDigitizedTime);
  ttree->SetBranchAddress("hitPMTCharge", &hitPMTCharge);
  ttree->SetBranchAddress("hitPMTDigitizedCharge", &hitPMTDigitizedCharge);
  ttree->SetBranchAddress("mcPMTID", &mcPMTID);
  ttree->SetBranchAddress("mcPEIndex", &mcPEIndex);
  ttree->SetBranchAddress("mcPETime", &mcPETime);
  ttree->SetBranchAddress("mcPEProcess", &mcPEProcess);
}

int main(int argc, char**argv){

  cout<<"starting.."<<endl;
  parseArguments( argc, argv);
  cout<<"intput and output txt : "<<inputTxt.Data()<<" "<<outputTxt.Data()<<endl;
  if (ifOutputRoot) cout<<"output root file : "<<outputRoot.Data()<<endl;

  TH3F* h3 = new TH3F("","", 100,-1200,1200, 100, -1200,1200,100,-1200,1200);
  TH3F* h4 = new TH3F("","", 100,-1200,1200, 100, -1200,1200,100,-1200,1200);

  TH1F* hc = new TH1F("","",50,0,50);  
  TH1F* hhx[1000];
  TH1F* hhz[1000];
  for (int i=0;i<1000;i++){
    hhx[i] = new TH1F("","",20,0,1);
    hhz[i] = new TH1F("","",20,-1200,1200);
  }

  double pmtx[264]={};
  double pmty[264]={};
  double pmtz[264]={};

  TFile* outfile;
  TTree tree("tree","tree");
  double trueDir[3], recoDir[3];
  int nPMT; double PMTid[264][3]; double PMTtime[264]; double PMTangle[264];

  std::vector<std::vector<double> > pmtlist;
  std::vector<std::vector<std::vector<double> > > fitlist(10000);
  double ggg[10000],hhh[10000],mmm[10000];
  double iniX[10000] = {};
  double iniY[10000] = {};
  double iniZ[10000] = {};
  int eventTaker[10000] = {};
  int bbb;

  cout<<888<<endl;
  TH1F* h_times = new TH1F("","",100,-5,15);
  TH1F* h_timec = new TH1F("","",125,-10,15);
  TH1F* h_pmtzs = new TH1F("","",30,-1500,1500);
  TH1F* h_pmtzc = new TH1F("","",30,-1500,1500);
  TH1F* h_charges = new TH1F("","",200,0,1000);
  TH1F* h_chargec = new TH1F("","",200,0,1000);
  TH1F* h_mcpcount = new TH1F("","",300,0,300);
  TH1F* h_mcpecount = new TH1F("","",600,0,600);
  TH1F* h_mcnhits = new TH1F("","",300,0,300);
  TH1F* h_mcpecountc = new TH1F("","",600,0,600);
  TH1F* h_mcnhitsc = new TH1F("","",300,0,300);
  TH2F* h_timeVSprocess = new TH2F("","",100,-5,15,5,0,5);
  TH2F* h_timeVSpos = new TH2F("","",100,-10,15,30,-1500,1500);

  /////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////// using the ntuple input file
  // PDF:
  // 1. pmt locx 2. pmt locy 3. pmt locz 4. time 5. true locx 6. true locy 7. true locz
  //                         |  pmt id
  // 8. constant 0 9. theta 10. phi 11. charge 
  // fitting:
  // same thing

    //read_ntuple(inputTxt.Data());
      
    TFile ftemp(inputTxt.Data());
    ttree = (TTree*)ftemp.Get("output");
    nentries = ttree->GetEntries();
    //ttree->Print();
    ttree->SetBranchAddress("mcx", &mcx);
    ttree->SetBranchAddress("mcy", &mcy);
    ttree->SetBranchAddress("mcz", &mcz);
    ttree->SetBranchAddress("mcu", &mcu);
    ttree->SetBranchAddress("mcv", &mcv);
    ttree->SetBranchAddress("mcw", &mcw);
    ttree->SetBranchAddress("mcke", &mcke);
    ttree->SetBranchAddress("evid", &evid);
    ttree->SetBranchAddress("subev", &subev);
    ttree->SetBranchAddress("nanotime", &nanotime);
    ttree->SetBranchAddress("mcpcount", &mcpcount);
    ttree->SetBranchAddress("mcpecount", &mcpecount);
    ttree->SetBranchAddress("mcnhits", &mcnhits);
    ttree->SetBranchAddress("pdgcodes", &pdgcodes);
    ttree->SetBranchAddress("mcKEnergies", &mcKEnergies);
    ttree->SetBranchAddress("mcPosx", &mcPosx);
    ttree->SetBranchAddress("mcPosy", &mcPosy);
    ttree->SetBranchAddress("mcPosz", &mcPosz);
    ttree->SetBranchAddress("mcDirx", &mcDirx);
    ttree->SetBranchAddress("mcDiry", &mcDiry);
    ttree->SetBranchAddress("mcDirz", &mcDirz);
    ttree->SetBranchAddress("hitPMTID", &hitPMTID);
    ttree->SetBranchAddress("hitPMTTime", &hitPMTTime);
    ttree->SetBranchAddress("hitPMTDigitizedTime", &hitPMTDigitizedTime);
    ttree->SetBranchAddress("hitPMTCharge", &hitPMTCharge);
    ttree->SetBranchAddress("hitPMTDigitizedCharge", &hitPMTDigitizedCharge);
    ttree->SetBranchAddress("mcPMTID", &mcPMTID);
    ttree->SetBranchAddress("mcPEIndex", &mcPEIndex);
    ttree->SetBranchAddress("mcPETime", &mcPETime);
    ttree->SetBranchAddress("mcPEProcess", &mcPEProcess);
  
    //add the pmt mapping from ntuple here
    TTree* pmt_tree = (TTree*)ftemp.Get("meta");
    pmt_tree->SetBranchAddress("pmtId", &pmtId);
    pmt_tree->SetBranchAddress("pmtX", &pmtX);
    pmt_tree->SetBranchAddress("pmtY", &pmtY);
    pmt_tree->SetBranchAddress("pmtZ", &pmtZ);
    //pmt_tree->Print();

    pmt_tree->GetEntry(0);
    int spmt = -1;
    for (int i=0;i<pmtId->size(); i++){
      if (i != pmtId->at(i)) {
        cout<<"pmt id does not correspond to the order! exit!"<<endl;
	exit(0);
      }
      if (pmtId->at(i)>spmt) {
	spmt = pmtId->at(i);
      }
      pmtx[i] = pmtX->at(i);
      pmty[i] = pmtY->at(i);
      pmtz[i] = pmtZ->at(i);
    }
    npmt_bin = spmt+1;
    cout<<"ntuple entries "<<nentries<<endl;
    int counter[10000] = {};
    int pecounter[10000] = {};

    float sourceDir_x = -9;
    float sourceDir_y = -9;
    float sourceDir_z = -9;

    float xsign = 0, ysign=0, zsign = 0;
    int xadd = 0, yadd=0, zadd=0;

    if (useSource){

      stringstream ss;
      const char* stt;
      stt = inputTxt.Data();
      cout<<"input "<<stt<<endl;
      int num[1000];
      for (int i=0;i<strlen(stt); i++){
        //cout<<i<<" "<<stt[i]<<endl;
        if (stt[i]>='0' && stt[i]<='9'){
          num[i] = stt[i] - '0';
          //cout<<i<<" "<<stt[i]<<" "<<num[i]<<endl;
        }
      }
      if (stt[112] == '-'){
        xadd = 1;
        xsign = -1;
        if (stt[119] == '-'){
          yadd = 1;
          ysign = -1;
          if (stt[126] == '-'){
            zadd = 1;
            zsign = -1;
          }
          else{// x and y are negative, z is positive
            zadd = 0;
            zsign = 1;
          }
        }
        else{ // x is negative, y is positive
          yadd = 0;
          ysign = 1;
          if (stt[125] == '-'){
            zadd = 1;
            zsign = -1;
          }
          else{// x and y are negative, z is positive
            zadd = 0;
            zsign = 1;
          }
        }
      }
      else{
        xadd = 0;
        xsign = 1;
        if (stt[118] == '-'){
          yadd = 1;
          ysign = -1;
          if (stt[125] == '-'){
            zadd = 1;
            zsign = -1;
          }
          else{// x positive, y negative, z is positive
            zadd = 0;
            zsign = 1;
          }
        }
        else{ // x is positive, y is positive
          yadd = 0;
          ysign = 1;
          if (stt[124] == '-'){
            zadd = 1;
            zsign = -1;
          }
          else{// x and y are positive, z is positive
            zadd = 0;
            zsign = 1;
          }
        }
      }
      cout<<"add "<<xadd<<" "<<yadd<<" "<<zadd<<endl;
      cout<<"sign "<<xsign<<" "<<ysign<<" "<<zsign<<endl;
      sourceDir_x = xsign * (num[114+xadd]*0.1 + num[115+xadd]*0.01 + num[116+xadd]* 0.001);
      sourceDir_y = ysign * (num[120+xadd+yadd]*0.1 + num[121+xadd+yadd]*0.01 + num[122+xadd+yadd]* 0.001);
      sourceDir_z = zsign * (num[126+xadd+yadd+zadd]*0.1 + num[127+xadd+yadd+zadd]*0.01 + num[128+xadd+yadd+zadd]* 0.001);
      if (num[124+xadd+yadd+zadd] == 1 && stt[123+xadd+yadd+zadd] == '-') {sourceDir_z = -1;}
      cout<<"using source and source direction "<<sourceDir_x<<" "<<sourceDir_y<<" "<<sourceDir_z<<endl;
     
    }

    for (int i=0;i<nentries;i++){
      ttree->GetEntry(i);

/////////////////////////////////////////////////////////////// Looking at the true MC info. in order to get the dicrochon situation.
      // Dichrocon location with pmtloaction_201 id: 7,2,21,15,22,29,10,14,5,6 (a pictue provided by Sam on Slack)
      int registeredPMT[10] = {7,2,21,15,22,29,10,14,5,6};
      int hitDic[10] = {};

      if (useDic){
        double pmtxloc[mcPMTID->size()] ={};
        double pmtyloc[mcPMTID->size()] ={};
        double pmtzloc[mcPMTID->size()] ={};
        double pmtid[mcPMTID->size()] = {};

        for (int ihit = 0; ihit< mcPMTID->size(); ihit++){
          pmtxloc[ihit] = pmtx[mcPMTID->at(ihit)];
          pmtyloc[ihit] = pmty[mcPMTID->at(ihit)];
          pmtzloc[ihit] = pmtz[mcPMTID->at(ihit)];
          pmtid[ihit] = mcPMTID->at(ihit);
        }

        int dhit = 0;
        int hitcounter = 0;
        double locx = -999;
        double locy = -999;
        double locz = -999;
        int pmtidd = 0;

        for (int ihit = 0; ihit< mcPETime->size(); ihit++){

          if (mcPEIndex->at(ihit) == 0) {
            locx = pmtxloc[hitcounter];
            locy = pmtyloc[hitcounter];
            locz = pmtzloc[hitcounter];
            pmtidd =  pmtid[hitcounter];
            //1=Cherenkov, 0=Dark noise, 2=Scint., 3=Reem., 4=Unknown
            if (mcPEProcess->at(ihit) == 1 ){
              bool exists = std::find(std::begin(registeredPMT), std::end(registeredPMT), pmtidd) != std::end(registeredPMT);
              if (exists) {
                hitDic[dhit] = pmtidd;
                dhit ++;
              }
            }
            hitcounter ++;
            //cout<<"ievent, ihit, pmt id, x, y, z "<<i<<" "<<ihit<<" "<<pmtidd<<" "<<locx<<" "<<locy<<" "<<locz<<endl;    
          }
        }
     }
///////////////////////////////////////////////////////////

     
      if(i > evtNum || i< evtBase) continue;
      std::vector<double> pmtloc;
      double theta = -1; double phi = -1;
      int pmtgothit[300] = {};

      h_mcpcount->Fill(mcpcount);
      h_mcpecount->Fill(mcpecount);
      h_mcnhits->Fill(mcnhits);

      for (int ihit = 0; ihit< hitPMTTime->size(); ihit++){
      
        double userWeight = 1;
        if (useDic)
        {
          bool exists = std::find(std::begin(hitDic), std::end(hitDic), hitPMTID->at(ihit)) != std::end(hitDic);
          if (exists && hitPMTID->at(ihit) != 0){
            // above 450 nm, getting the intergral of all light at 5% 500 LY and Cherenkov light, using online plot digitizer.
            userWeight = 43774./20849.;
	    //userWeight = 1;
          }
        }

	pmtloc.push_back(pmtx[hitPMTID->at(ihit)]);
        pmtloc.push_back(pmty[hitPMTID->at(ihit)]);
	//cout<<"pmt id "<< hitPMTID->at(ihit)<<endl;
	if (perDir || perPMT)
	  pmtloc.push_back(hitPMTID->at(ihit));
	else
	  pmtloc.push_back(pmtz[hitPMTID->at(ihit)]);
	double pmtxloc = pmtx[hitPMTID->at(ihit)];
	double pmtyloc = pmty[hitPMTID->at(ihit)];
	double pmtzloc = pmtz[hitPMTID->at(ihit)];
	bbb = hitPMTID->at(ihit);
	double ccc = hitPMTTime->at(ihit) - (sqrt((pmtxloc-trueOriginX)*(pmtxloc-trueOriginX)+(pmtyloc-trueOriginY)*(pmtyloc-trueOriginY) + (pmtzloc-trueOriginZ)*(pmtzloc-trueOriginZ))/200.) + timeCorrection;
        double locX = mcx+ gRandom->Gaus(0, vertexSmear);
        double locY = mcy+ gRandom->Gaus(0, vertexSmear);
        double locZ = mcz+ gRandom->Gaus(0, vertexSmear);

        h_timec->Fill(ccc);
        h_pmtzc->Fill(pmtzloc );
        h_timeVSpos->Fill(ccc, pmtzloc);

	pmtloc.push_back(ccc);
	pmtloc.push_back(locX);
	pmtloc.push_back(locY);
	pmtloc.push_back(locZ);

        if (useSource){

          if (sourceDir_x == -9 || sourceDir_y == -9 || sourceDir_z == -9){
            cout<<"trying to use source but the source direciton is not set! exit!"<<endl;
            exit(0);
          }
          mcu = sourceDir_x;
          mcv = sourceDir_y;
          mcw = sourceDir_z;

	}

        ggg[i] = mcu / sqrt(mcu*mcu + mcv*mcv + mcw*mcw);	
	hhh[i] = mcv / sqrt(mcu*mcu + mcv*mcv + mcw*mcw);
	mmm[i] = mcw / sqrt(mcu*mcu + mcv*mcv + mcw*mcw);

        if (orient != 999) {
          std::cout<<"only isotropic PDF is supported as of Dec. 2022. Exiting!"<<std::endl;
          exit(0);
        }
        if (orient == 999) {
          pmtloc.push_back(0);
          if (mmm[i]>0){
            pmtloc.push_back(TMath::ATan(sqrt(ggg[i]*ggg[i]+hhh[i]*hhh[i])/abs(mmm[i])));
            theta = TMath::ATan(sqrt(ggg[i]*ggg[i]+hhh[i]*hhh[i])/abs(mmm[i]));
          }
          else {
            pmtloc.push_back(TMath::Pi() - TMath::ATan(sqrt(ggg[i]*ggg[i]+hhh[i]*hhh[i])/abs(mmm[i])));
            theta = TMath::Pi() - TMath::ATan(sqrt(ggg[i]*ggg[i]+hhh[i]*hhh[i])/abs(mmm[i]));
          }
          if (ggg[i]>=0 && hhh[i]>0) {pmtloc.push_back(TMath::ATan(hhh[i]/ggg[i])); phi = TMath::ATan(hhh[i]/ggg[i]);}
          else if (ggg[i]<0 && hhh[i]>=0) {pmtloc.push_back(TMath::ATan(hhh[i]/ggg[i])+TMath::Pi()); phi = TMath::ATan(hhh[i]/ggg[i])+TMath::Pi();}
          else if (ggg[i]<0 && hhh[i]<0) {pmtloc.push_back(TMath::ATan(hhh[i]/ggg[i])+ TMath::Pi()); phi = TMath::ATan(hhh[i]/ggg[i])+ TMath::Pi();}
          else if (ggg[i]>=0 && hhh[i]<=0) {pmtloc.push_back(TMath::Pi()*2 -TMath::ATan(abs(hhh[i])/ggg[i])); phi = TMath::Pi()*2 -TMath::ATan(abs(hhh[i])/ggg[i]); }
        }
        //cout<<"theta and phi "<<theta<<" "<<phi<<endl;
	pmtloc.push_back(hitPMTCharge->at(ihit));

        pmtloc.push_back(userWeight);

        pmtlist.push_back(pmtloc);
        if (i>= lowerEvt && i< upperEvt ){
          if (theta > -999){
            fitlist[i].push_back(pmtloc);
            eventTaker[i] = 1;
          }
          if (iniX[i] == 0 && iniY[i] == 0 && iniZ[i] == 0)
            iniX[i] = pmtx[bbb]; iniY[i] = pmty[bbb]; iniZ[i] = pmtz[bbb];
        }
        pmtloc.clear();
      }
      for (int iii=evtBase;iii<evtNum;iii++){
        h_mcnhitsc->Fill(counter[iii]);
        h_mcpecountc->Fill(pecounter[iii]);
      }
  }
  ////////////////////////////////////////////// end of using ntuple input file
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (ifOutputRoot){
    outfile = TFile::Open( outputRoot.Data(), "recreate");
    tree.Branch("trueDir",&trueDir,"trueDir[3]/D");
    tree.Branch("recoDir",&recoDir,"recoDir[3]/D");
    tree.Branch("nPMT",&nPMT,"nPMT/I");
    tree.Branch("PMTid",&PMTid,"PMTid[264][3]/D");
    tree.Branch("PMTtime",&PMTtime,"PMTtime[264]/D");
    tree.Branch("PMTangle",&PMTangle,"PMTangle[264]/D");
  }

  cout<<"event induced PMT loaded "<<endl;

  RooFitResult* res;
  char formula[10];
  dirFit * rep = new dirFit ("_rep");
  
  // 1 ton time cut 4.5, 2.4 ton time cut 5.5, 5 ton time cut 6.5
  cout<<"doing digitize ? "<<digitize<<"  adding sometime ? "<<sometime<<endl;

  rep->SetPromptCut(sometime);  

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
  //std::vector<wbPDF*> pdfss;
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
    cout<<"reading perDir PDFs .. "<<endl;
    rep->SetNbins_time(ntime_bin);
    rep->SetNdir(ndir_bin);
    rep->SetNbins_pmt(npmt_bin);
    if (!externalPDF){
      cout<<"doing a perDir pdf generation.."<<endl;
      pdfss = rep->Reading_Processing_Events_PerDir(pmtlist,"dirpdf", iniVertex, do2Dpdf, doCharge, doCos );
      cout<<"checking pdfss, size, and size of individual pdfss element:  "<<pdfss.size()<<" "<<pdfss[1]->GetDirPDF()->GetNbinsX()<<" "<<pdfss[1]->GetDirPDF()->GetNbinsY()<<endl;
    }
    else{
        TFile fepdf(pdf_filename.Data());
        for (int iir=0;iir<ndir_bin;iir++){
          TH2F* hpdf = (TH2F*)fepdf.Get(Form("output_dir_%d",iir));
          pdfss[iir]->SetDirPDF(hpdf);
        }
        fepdf.Close();
	rep->SetDirPDFs(pdfss);
	//cout<<"testing here (dir 1 integral) "<<rep->GetDirPDFs().at(1)->GetDirPDF()->Integral()<<endl;
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
    h_charges->Write("h_charges");
    h_chargec->Write("h_chargec");
    h_mcpcount->Write("h_mcpcount");
    h_mcpecount->Write("h_mcpecount");
    h_mcnhits->Write("h_mcnhits");
    h_mcpecountc->Write("h_mcpecountc");
    h_mcnhitsc->Write("h_mcnhitsc");
    h_timeVSprocess->Write("h_timeVSprocess");
    h_timeVSpos->Write("h_timeVSpos");
  }

  if (!externalPDF) { cout<<"just getting pdf, not performing fit .. "<<endl; outfile->Write(); outfile->Close(); exit(1);}

  for (Int_t aaa = lowerEvt;aaa< upperEvt; aaa++){

    //cout<<aaa<<" taker "<<eventTaker[aaa]<<endl;
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

    ofstream out;
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
          if (oTxt)
	    out<<aaa<<" "<<xloop<<" "<<yloop<<" "<<res<<endl;
          if (res < currRes && res > 0) {
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
    outfile->Close();
  }
}
