#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <iostream>
#include <cmath>

void proton_absorption_study(int runTypeInt = 1) {

  // =============================
  // 1. Open ROOT file
  // =============================
//   TFile *f = new TFile("deut_replay_prod_20856_singles.root", "READ");
//   TFile *outROOT = new TFile("Proton_abs_plots.root", "recreate");

//   if (!f || f->IsZombie()) {
//     std::cerr << "Error opening file\n";
//     return;
//   }

  // =============================
  // 2. Get trees
  // =============================
//   TTree *T = (TTree*)f->Get("T");
//   TTree *TSH = (TTree*)f->Get("TSH");

//   if (!T || !TSH) {
//     std::cerr << "Error: T or TSH tree not found\n";
//     return;
//   }

//   Long64_t nentries = T->GetEntries();
//   Long64_t nscal    = TSH->GetEntries();

// =============================
  // 1. Define run type
  // =============================
  enum RunType { kSingles = 0, kCoin = 1, kDummy = 2 };
  RunType runType = static_cast<RunType>(runTypeInt);

  std::cout << "=================================\n";
  if (runType == kSingles) std::cout << "Running SINGLE ARM analysis\n";
  if (runType == kCoin)    std::cout << "Running COINCIDENCE analysis\n";
  if (runType == kDummy)   std::cout << "Running DUMMY analysis\n";
  std::cout << "=================================\n";

  // =============================
  // 2. Choose input ROOT file
  // =============================
  TFile *infile = nullptr;

  if (runType == kSingles)
    infile = new TFile("deut_replay_prod_20859_singles.root", "READ");
  else if (runType == kCoin)
    infile = new TFile("deut_replay_prod_20858_-1.root", "READ");
  else if (runType == kDummy)
    infile = new TFile("deut_replay_prod_20860_dummy.root", "READ");

  if (!infile || infile->IsZombie()) {
    std::cerr << "ERROR: Cannot open input file\n";
    return;
  }

  TTree *T = (TTree*) infile->Get("T");

//OUTPUT FILE
 TDirectory *dir = nullptr;
 TFile *outROOT = new TFile("Proton_abs_plots_delta+4.root", "UPDATE");

if (runType == kCoin)
  dir = (TDirectory*) outROOT->Get("Coin");
else if (runType == kSingles)
  dir = (TDirectory*) outROOT->Get("Singles");
else if (runType == kDummy)
  dir = (TDirectory*) outROOT->Get("Dummy");

if (!dir) {
  if (runType == kCoin)      dir = outROOT->mkdir("Coin");
  else if (runType == kSingles) dir = outROOT->mkdir("Singles");
  else if (runType == kDummy)   dir = outROOT->mkdir("Dummy");
}

dir->cd();



  if (!T) {
    std::cerr << "ERROR: Tree T not found\n";
    return;
  }

  Long64_t nentries = T->GetEntries();

  // =============================
  // 3. Declare variables (T tree)
  // =============================
  Double_t h_delta, e_delta, W;
  Double_t ztar,hhod_beta_ntrk,hcal_etottracknorm,hhod_beta, hms_should, hms_did;
  Double_t e_xptar, e_yptar,hhod_GoodScinHit, hcer_npesum,hcal_etot,hcal_etotnorm ;
  Double_t h_xptar, h_yptar, evNum,ecal,hdc_ntrack, hTrkEff, hTrkEff_err, did_cut, should_cut,Em, Accp_cut;
  Double_t hTRIG1_tdc = 0;

  // =============================
  // 5. Set branch addresses (T)
  // =============================
  T->SetBranchAddress("P.react.z",        &ztar);
  T->SetBranchAddress("P.gtr.dp",          &e_delta);
  T->SetBranchAddress("P.kin.primary.W",  &W);
  T->SetBranchAddress("H.kin.secondary.emiss",  &Em);
  T->SetBranchAddress("P.cal.etotnorm", &ecal);
  T->SetBranchAddress("H.gtr.dp",          &h_delta);
  T->SetBranchAddress("H.gtr.th",          &h_xptar);
  T->SetBranchAddress("H.gtr.ph",          &h_yptar);
  T->SetBranchAddress("P.gtr.th",          &e_xptar);
  T->SetBranchAddress("P.gtr.ph",          &e_yptar);
  
  T->SetBranchAddress("T.coin.hTRIG1_ROC1_tdcTime", &hTRIG1_tdc);


  //T->SetBranchAddress("T.coin.hTRIG1_ROC1_tdcTime", &hTRIG1_tdc);
   //HMS DETECTORS
  T->SetBranchAddress("H.cer.npeSum",         &hcer_npesum);
  T->SetBranchAddress("H.cal.etot",           &hcal_etot);
  T->SetBranchAddress("H.cal.etotnorm",       &hcal_etotnorm);
  T->SetBranchAddress("H.cal.etottracknorm",  &hcal_etottracknorm);
  T->SetBranchAddress("H.hod.betanotrack",    &hhod_beta_ntrk);
  T->SetBranchAddress("H.hod.beta",           &hhod_beta);
  T->SetBranchAddress("H.hod.goodscinhit",    &hhod_GoodScinHit);    
  T->SetBranchAddress("H.dc.ntrack",          &hdc_ntrack);

 // =============================
  // Create Histograms
  // =============================


        // TH1F* H_W = new TH1F("W", "W", 100, 0.5, 1.7);
        // TH1F* H_W_did = new TH1F("W", "W_did", 100, 0.5, 1.7);
        // TH1F* H_W_should = new TH1F("W", "W_should", 100, 0.5, 1.7);
       // TH1F* H_Q2 = new TH1F("Q2", "Q2", 100, 1.3, 6);
        // TH1F* H_h_delta = new TH1F("h_delta", "h_delta", 100, -15, 15);
        // TH1F* H_e_delta = new TH1F("e_delta", "e_delta", 100, -12, 12);
        // TH1F* H_h_xptar = new TH1F("h_xptar", "h_xptar", 100, -0.2, 0.2);
        // TH1F* H_h_xptar = new TH1F("h_xptar", "h_xptar", 100, -0.2, 0.2);
        // TH1F* H_h_xptar = new TH1F("h_xptar", "h_xptar", 100, -0.2, 0.2);
        TH1F* H_ztar = new TH1F("ztar", "e_ztar", 100, -10, 10);
        TH1F* H_Em = new TH1F("Em", "Em", 100, -0.1, 0.5);
        TH1F* H_W          = new TH1F("H_W", "W", 100, 0.5, 1.7);
        TH1F* H_W_did      = new TH1F("H_W_did", "W_did", 100, 0.5, 1.7);
        TH1F* H_W_should   = new TH1F("H_W_should", "W_should", 100, 0.5, 1.7);

        TH1F* H_e_xptar          = new TH1F("H_e_xptar", "e_xptar", 100, -0.06, 0.06);
        TH1F* H_e_xptar_should   = new TH1F("H_e_xptar_should", "e_xptar_should", 100, -0.06, 0.06);
        TH1F* H_e_xptar_did      = new TH1F("H_e_xptar_did", "e_xptar_did", 100, -0.06, 0.06);

        TH1F* H_e_yptar          = new TH1F("H_e_yptar", "e_yptar", 100, -0.04, 0.04);
        TH1F* H_e_yptar_should   = new TH1F("H_e_yptar_should", "e_yptar_should", 100, -0.04, 0.04);
        TH1F* H_e_yptar_did      = new TH1F("H_e_yptar_did", "e_yptar_did", 100, -0.04, 0.04);
        // Momentum acceptance 2D Histograms
        TH2F* H_e_delta_vs_h_delta = new TH2F("e_delta_vs_h_delta", "e_delta vs h_delta; h_delta; e_delta",  100, -8, 8, 100, -5, 8 );
        TH2F* H_e_xptar_vs_e_yptar = new TH2F("e_xptar_vs_e_yptar", "e_xptar vs e_yptar; e_yptar ; e_xptar",  100, -0.03, 0.03, 100, -0.03, 0.03);
        TH2F* H_h_xptar_vs_h_yptar = new TH2F("h_xptar_vs_h_yptar", "h_xptar vs h_yptar; h_yptar ; h_xptar",  100, -0.1, 0.1, 100, -0.1, 0.1);
        

        //Tracking Efficiency check histos
        TH1F* H_hdc_ntrack     = new TH1F("H_hdc_ntrack", "H.dc.ntrack", 5, -0.5, 4.5);
        TH1F* H_hhod_beta_ntrk = new TH1F("H_hhod_beta_ntrk", "H.hod.beta (no track)", 100, 0, 2);
        TH1F* H_hhod_goodscin  = new TH1F("H_hhod_goodscin", "H.hod.goodscinhit", 5, -0.5, 4.5);


  // =============================
  // 7. Event-level loop (T tree)
  // =============================
  Long64_t e_should = 0;
  Long64_t e_should_hdelta = 0;
  Long64_t e_did = 0;
  Long64_t n_hTRIG1 =0;
  Long64_t h_did = 0;
  Long64_t h_should = 0;


  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);

    bool hTRIG1_tdc_cut = hTRIG1_tdc>0;
    bool edelta_cut  = (e_delta > 2 && e_delta < 6.5);
    bool W_cut       = (W > 0.9 && W < 1.0);
    bool ztar_cut    = (ztar > -2.5 && ztar < 2.5);
    bool hdelta_cut  = (std::abs(h_delta) < 7.0);
    bool eyptar_cut  = (std::abs(e_yptar) < 0.03);
    bool exptar_cut  = (std::abs(e_xptar) < 0.07);
    bool ecal_cut = ecal > 0.6;
    bool Em_cut = Em<0.03;

    //CUTS: HMS TRACKING EFFICIENCY
bool hdc_ntrk_cut =hdc_ntrack >= 1;
bool hScinGood_cut = hhod_GoodScinHit==1 ;
bool hcer_NPE_Sum_cut = hcer_npesum >= 0. && hcer_npesum <= 0.5;
bool hetotnorm_cut = hcal_etotnorm >= 0. && hcal_etotnorm <= 0.6;
bool hBeta_notrk_cut = hhod_beta_ntrk >= 0.5 && hhod_beta_ntrk <= 1.5;
//electrons (or hadrons) that 'SHOULD' have passed the cuts to form a track

hms_should = hScinGood_cut && hcer_NPE_Sum_cut && hetotnorm_cut && hBeta_notrk_cut;
//electrons (or hadrons) that 'DID' passed the cuts to form a track
hms_did = hdc_ntrk_cut && hms_should;

// For Proton absorption study
did_cut = edelta_cut && ztar_cut && eyptar_cut && W_cut && exptar_cut && ecal_cut && hdelta_cut && hTRIG1_tdc_cut;
should_cut = edelta_cut && ztar_cut && eyptar_cut && W_cut && exptar_cut && ecal_cut; 

//Acceptance cuts
Accp_cut = edelta_cut && hdelta_cut && Em_cut ;


    if (edelta_cut && ztar_cut && eyptar_cut && W_cut && exptar_cut&& ecal_cut ) {
      e_should++;
    }
    if (edelta_cut && ztar_cut && eyptar_cut && W_cut && exptar_cut && ecal_cut && hdelta_cut ) {
      e_should_hdelta++;
    }
     
    if (edelta_cut && ztar_cut && eyptar_cut && W_cut && exptar_cut && ecal_cut && hdelta_cut && hTRIG1_tdc_cut) {
       e_did++;
    }
    if (hTRIG1_tdc_cut) {
    n_hTRIG1++;
   }

    if(hms_did){ h_did++;}
	if(hms_should){ h_should++; }

     if (runType == kCoin && Accp_cut){
     H_Em->Fill(Em);
     H_hdc_ntrack->Fill(hdc_ntrack);
     H_hhod_beta_ntrk->Fill(hhod_beta_ntrk);
     H_hhod_goodscin->Fill(hhod_GoodScinHit);
     H_e_delta_vs_h_delta->Fill(h_delta, e_delta)  ;
     H_e_xptar_vs_e_yptar->Fill(e_yptar, e_xptar)  ;
     H_h_xptar_vs_h_yptar->Fill(h_yptar, h_xptar)  ;

    }
     if(runType == kSingles && should_cut){
    H_W_should->Fill(W);
    H_e_yptar_should->Fill(e_yptar);
    H_e_xptar_should->Fill(e_xptar);

    }
    if(runType == kSingles && did_cut){
    H_W_did->Fill(W);
    H_e_yptar_did->Fill(e_yptar);
    H_e_xptar_did->Fill(e_xptar);

    }
    if(runType == kDummy ){
        H_ztar->Fill(ztar);
    }

  }

  
//   for (Long64_t i = 0; i < nentries; i++) {
//     T->GetEntry(i);

   

//   }
  // =============================
  // 10. Results
  // =============================
   //Calculate HMS Tracking Efficiency                                                                                                                 
   if (h_should > 0) {
    hTrkEff = double(h_did)/double(h_should);
  }


   //Scale the "did" by HMS tracking efficiency
   if (runType == kSingles && hTrkEff > 0) {
  H_e_xptar_did->Scale(1./hTrkEff);
  H_e_yptar_did->Scale(1./hTrkEff);
  H_W_did->Scale(1./hTrkEff);
}

  double survival_fraction = (e_did/hTrkEff)/e_should;
  

  std::cout << "=============================\n";
  std::cout << " Proton Absorption Study\n";
  std::cout << "=============================\n";
  std::cout << "Protons SHOULD trigger : " << e_should << "\n";
  std::cout << "Protons SHOULD without ptrig : " << e_should_hdelta << "\n";
  std::cout << "Protons DID trigger    : " << e_did    << "\n";
  std::cout << "Events that pass hTRIG cut:" << n_hTRIG1 << "\n";
  std::cout << "HMS SHOULD  : " << h_should << "\n";
  std::cout << "HMS DID  : " << h_did << "\n";
  std::cout << "HMS Tracking Eff    : " << hTrkEff    << "\n";
  std::cout << "Survival fraction      : " << survival_fraction << "\n";
  std::cout << "Absorption loss        : " << 1.0 - survival_fraction << "\n";
  std::cout << "=============================\n";

  outROOT->Write();
  outROOT->Close();
}


