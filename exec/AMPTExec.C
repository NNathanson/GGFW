// root 

#include "TList.h"
#include "TProfile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TNamed.h"
#include <string>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TComplex.h"
#include "TBranch.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TGraph.h"

#include "../FlowContainer.h"
#include "../FlowPtContainer.h"
#include "../BootstrapProfile.h"
#include "../GFW.h"
#include "../GFWCumulant.h"


// STL
#include <iostream>
#include <vector>
#include <cmath>

// Costum classes here

// global var
TRandom* rndGenerator;
Int_t mMaxPtCorr = 8;
std::vector<std::vector<Double_t>> wp;
Double_t fEtaMpt = 0.4;
FlowPtContainer* fFCpt;
GFW* fGFW;
FlowContainer* fFC;
std::vector<GFW::CorrConfig> corrconfigs;
GFWCorrConfigs configs;
Int_t fNBootstrapProfiles = 10; 
UInt_t fConsistencyFlag = 3;
Bool_t fUseEventWeightOne = kTRUE;
Bool_t fUseNch = kFALSE;
Int_t nCentBins = 30;
Double_t* centBins = new Double_t[nCentBins+1];
Int_t nNchBins = 880;
Double_t* nchBins = new Double_t[nNchBins+1];
TRandom* fRndm;
TList* fptList;
TList *fCovList;
Int_t fPtMpar = 8;
Double_t ptMin = 0.2;
Double_t ptMax = 3;
Double_t fEtaAcceptance = 0.8;
Double_t fEtaV2Sep = 0.4;
unsigned int fEventWeight = 1;
TH1D *fMultiDist;


// DO NOT TOUCH VAR
Int_t ev, id, multiplicity;
Double_t fImpact;
std::vector<Double_t>* fPtTrack(0);
std::vector<Double_t>* fEtaTrack(0);
std::vector<Double_t>* fPhiTrack(0);
std::vector<Int_t>* fIdTrack(0);


void AMPTExec(TString deformation = "", Bool_t small = kFALSE, TString path = ""){ 
    gSystem->Load("../libGGFW");
    gROOT->ProcessLine("#include <vector>");
    printf("Processing simulations w.  %s...\n", deformation.Data());
    
    for(int i(0);i<=nNchBins;++i) nchBins[i] = i;
    for(int j(0);j<=nCentBins;++j) centBins[j] = j; 
    int nMultiBins = (fUseNch)?nNchBins:nCentBins;
    double* multiBins = (fUseNch)?nchBins:centBins;




    // set centrality estimator 
    TFile* getcentrestimator = TFile::Open(Form("/home/neo/Nina/GGFW/CentralityTask_b18_g00.root"), "READ");
    if(!getcentrestimator){ printf("Could not open centrality estimator"); return; };
    TGraph* getcentrality = (TGraph*)getcentrestimator->Get("fCentralityFromImpact");       // -->  multiplicity / vom / b
    getcentrestimator->Close();
    
    // setup chain of events
    TChain* chain = new TChain("T");
    TFile* fout;

    if(small==kTRUE){
      fout = new TFile(Form("output/AnalysisResult_%s_small.root", deformation.Data()), "RECREATE");
      chain->Add(path+deformation+"/central/ampt.root");
    }

    else{
      fout = new TFile(Form("output/AnalysisResult_%s.root", deformation.Data()), "RECREATE");
      chain->Add(path+deformation+"/central/ampt*.root");
    }
  

    // event var
    chain->SetBranchAddress("fNOfPart", &multiplicity);
    chain->SetBranchAddress("fImpact",  &fImpact);
    // track var
    chain->SetBranchAddress("fPtTrack",  &fPtTrack);
    chain->SetBranchAddress("fEtaTrack", &fEtaTrack);
    chain->SetBranchAddress("fPhiTrack", &fPhiTrack);
    chain->SetBranchAddress("fIdTrack",  &fIdTrack);


    printf("Total number of entries: %i \n",(int)chain->GetEntries());
    Int_t nOfEvents = chain->GetEntries();


    fRndm = new TRandom(0);
    fMultiDist = new TH1D("MultiDistribution","Multiplicity distribution; #it{N}_{ch}; N(events)",nMultiBins,multiBins);
    

    //Create GFW
    fGFW = new GFW();
    fGFW->AddRegion("refN",-fEtaAcceptance,-fEtaV2Sep,1,1);
    fGFW->AddRegion("refP",fEtaV2Sep,fEtaAcceptance,1,1);
    fGFW->AddRegion("mid",-fEtaAcceptance,fEtaAcceptance,1,2);
    

    //printf("setting configs\n");
    configs.SetCorrs({"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refP {2 2} refN {-2 -2}"});
    configs.SetHeads({"ChGap22", "ChGap32", "ChGap42", "ChGap24"});
    configs.SetpTDifs({0, 0, 0, 0});
    configs.SetpTCorrMasks({15, 1, 1, 3});

    for (auto i = 0; i < configs.GetSize(); ++i) {
      corrconfigs.push_back(fGFW->GetCorrelatorConfig(configs.GetCorrs()[i], configs.GetHeads()[i], configs.GetpTDifs()[i]));
    }
    //printf("configs set\n");

    fGFW->CreateRegions();

    //Setting up the FlowContainer
    //printf("Creating flow container\n");
    TObjArray *oba = new TObjArray();
    oba->Add(new TNamed("ChGap22","ChGap22")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap32","ChGap32")); //for gap (|eta|>0.4) case
    oba->Add(new TNamed("ChGap42","ChGap42")); //gap case
    oba->Add(new TNamed("ChGap24","ChGap24")); //for gap (|eta|>0.4) case


    fFC = new FlowContainer();
    TString fcname("FlowContainer");
    fFC->SetName(fcname.Data());
    fFC->Initialize(oba,nMultiBins,multiBins,fNBootstrapProfiles);
    delete oba;

    
    //printf("Flow container created\n");

    
    //printf("Creating flow-pt container\n");
    fFCpt = new FlowPtContainer("FlowPtContainer");
    fFCpt->setUseCentralMoments(kTRUE);
    fFCpt->setUseGapMethod(kTRUE);
    fFCpt->initialise(nMultiBins,multiBins,fPtMpar, configs, fNBootstrapProfiles);
    //printf("flow-pt container created\n");


    printf("Event loop entries : %i \n", nOfEvents);    
    //EVENT LOOP
    for(Int_t i(0); i < nOfEvents; i++){
        if (i % 1000 == 0) std::cout << "event no." << i << std::endl;   // Needed for event history           
        // get next event
        chain->GetEntry(i);   
	     // std::cout << "Deform: " << "		Event: " << i<<"/"<< nOfEvents << std::endl;        

      Double_t l_Random = fRndm->Rndm();
      Double_t l_Cent = getcentrality->Eval(fImpact); 
     // std::cout<<"Centrality: "<<l_Cent<<std::endl;

      
      Int_t lPosCount=0, lNegCount=0, lMidCount=0;
      Int_t nTotNoTracksMC=0;
      fGFW->Clear();
      fFCpt->clearVector();



        //TRACK LOOP
        //std::cout << "track loop starting" << std::endl;
        //std::cout << "Multiplicity: " << multiplicity << std::endl;
     
        for(Int_t ipart = 0; ipart < multiplicity; ipart++) {
          //std::cout << "ipart: " << ipart << std::endl;

          Double_t leta = fEtaTrack->at(ipart);
          if(TMath::Abs(leta) > 0.8) continue;
          if(std::isnan(leta)) continue;
          Double_t pt = fPtTrack->at(ipart);
          if(pt<ptMin || pt>ptMax) continue;
          if(leta<-fEtaV2Sep) lNegCount++;
          if(leta>fEtaV2Sep) lPosCount++;
          if(TMath::Abs(leta)<fEtaAcceptance) nTotNoTracksMC++; //Nch calculated in EtaNch region
          Double_t phi = fPhiTrack->at(ipart);
          //std::cout << "eta:  "<< leta << std::endl;
          //std::cout << "phi:  "<< phi << std::endl;
          fFCpt->fill(1,pt);
          //NB!!! ADD IF STATEMENT
          fGFW->Fill(leta,1,phi,1,1); //filling  gap (bit mask 1)? Since this is MC, weight is 1.
          //std::cout<<"gfw fill statement done"<<std::endl;

      
        }



      //std::cout<<"track loop complete"<<std::endl;

      //fEventCount->Fill("Tracks",1);
      //here in principle one could use the GFW output to check if the values are calculated, but this is more efficient
      // if(fConsistencyFlag&1) if(!lPosCount || !lNegCount) continue; // only events where v2{2, gap} could be calculated
      // if(fConsistencyFlag&2) if(nTotNoTracksMC<4) continue; //only events where v2{4} can be calculated (assuming same region as nch)
      // if(fConsistencyFlag&4) if(lPosCount < 2 && lNegCount < 2) continue; //only events where v2{4} gap can be calculated
     
      //Filling pT variance
      Double_t l_Multi = fUseNch?(1.0*nTotNoTracksMC):l_Cent;
      
      //A check in case l_Multi is completely off the charts (in MC, sometimes it ends up being... -Xe-310???)
      if(fUseNch && l_Multi<1) continue;
      //Double_t l_Random = fRndm->Rndm();
      //std::cout << wp[1][0] << std::endl;

      //std::cout<<"fFCpt functions below"<<std::endl;

      fFCpt->calculateCorrelations();
      //std::cout << "CalculateCorrelations executed" << std::endl;
      fFCpt->fillPtProfiles(l_Multi,l_Random);
      //std::cout << "FillProfiles executed" << std::endl;
      fFCpt->fillCMProfiles(l_Multi,l_Random);
      //std::cout << "FillCMProfiles executed" << std::endl;

      //std::cout<<"fProf bins: "<<fFC->fProf->GetNumberOfBins()<<std::endl;

      fMultiDist->Fill(l_Multi);
      
      for (uint l_ind = 0; l_ind < corrconfigs.size(); ++l_ind) {
        //std::cout<<"config: "<<corrconfigs.at(l_ind).Head.c_str()<<std::endl;
        //FillFCs(corrconfigs.at(l_ind),l_Multi,l_Random);
        if (!corrconfigs.at(l_ind).pTDif) {
          //std::cout<<"dnx calc:"<<std::endl;
          auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), 0, kTRUE).real();
          //std::cout<<"dnx: "<<dnx<<std::endl;
          if (dnx == 0)
            continue;
          
          //std::cout<<"val calc:"<<std::endl;
          auto val = fGFW->Calculate(corrconfigs.at(l_ind), 0, kFALSE).real() / dnx;
          //std::cout<<"???"<<std::endl;
          //std::cout<<"val: "<<val<<std::endl;
          if (std::abs(val) < 1) {
              fFC->FillProfile(corrconfigs.at(l_ind).Head.c_str(),l_Multi,val,(fUseEventWeightOne)?1:dnx,l_Random);
              

              fFCpt->fillVnDeltaPtProfiles(l_Multi, val, dnx, l_Random, configs.GetpTCorrMasks()[l_ind]);
              //std::cout<<"????"<<std::endl;
          }
          continue;
        }
      }

      
     // }
    }

    std::cout << "Loop finished" << std::endl;
    fout->cd();
    
    fMultiDist->Write("MultiDist", TObject::kSingleKey+TObject::kOverwrite);
    fFC->Write();
    fFCpt->Write();
 
    std::cout << "Write finished" << std::endl;
    fout->Close();
    std::cout << "fout closed" << std::endl;
    delete fout;
    return;
}
