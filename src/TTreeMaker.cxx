#include <memory>

// Include to get Event Loop.
#include <eventLoop.hxx>
#include <TReconTrack.hxx>
#include <TReconHit.hxx>
#include <TRealDatum.hxx>
#include "TTreeMaker.hxx"

#include <TEvent.hxx>

// Includes for ROOT classes

#include <HEPUnits.hxx>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TH1F.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>



CP::TTreeMakerLoop::TTreeMakerLoop() {

    run = 0;
    evt = 0;

    first_hit_X.clear();
    last_hit_X.clear();

    first_hit_Y.clear();
    last_hit_Y.clear();
	
    first_hit_Z.clear();
    last_hit_Z.clear();

    TPC_time = 0;
    PDS_RF_time.clear();
    PDS_delta_time.clear();
    PDS_trigger_type.clear();
    PDS_energy.clear();
    PDS_beam_trigger.clear();

    hfile= new TFile("tracks.root","RECREATE");
    tree = new TTree("tracks","");
    tree->Branch("run",&run,"run/I");
    tree->Branch("event",&evt,"event/I");

    tree->Branch("first_hit_X",&first_hit_X);
    tree->Branch("last_hit_X",&last_hit_X);
    tree->Branch("first_hit_Y",&first_hit_Y);
    tree->Branch("last_hit_Y",&last_hit_Y);
    tree->Branch("first_hit_Z",&first_hit_Z);
    tree->Branch("last_hit_Z",&last_hit_Z);

    tree->Branch("TPC_time",&TPC_time,"TPC_time/L");
    tree->Branch("PDS_RF_time","std::vector<Long64_t>",&PDS_RF_time);
    tree->Branch("PDS_delta_time","std::vector<Long64_t>",&PDS_delta_time);
    tree->Branch("PDS_trigger_type",&PDS_trigger_type);
    tree->Branch("PDS_energy",&PDS_energy);
    tree->Branch("PDS_beam_trigger",&PDS_beam_trigger);
	
}

CP::TTreeMakerLoop::~TTreeMakerLoop() {}

void CP::TTreeMakerLoop::Initialize(void) {}

bool CP::TTreeMakerLoop::operator () (CP::TEvent& event) {
		
    CP::THandle<CP::TReconObjectContainer> tracks = event.Get<CP::TReconObjectContainer>("~/fits/TCaptainRecon/final");
    CP::THandle<CP::TDataVector> dataPMT = event.Get<CP::TDataVector>("~/pmtData");

    run = event.GetContext().GetRun();
    evt = event.GetContext().GetEvent();

    TPC_time = event.GetTimeStamp();

    TString pdsEvent = "";
	if(dataPMT){ 
   for (int i=0; i<dataPMT->size(); i++) {
	pdsEvent.Form("~/pmtData/PDSEvent_%d",i);
	CP::THandle<CP::TEvent> eventPMT = event.Get<CP::TEvent>(pdsEvent);
	if (!eventPMT) {
	    std::cout<<"NO PMT EVENT"<<std::endl;
	}
	else {
	    //double TOF = (eventPMT->Get<CP::TRealDatum>("TOF_ns"))->GetValue();
	    PDS_RF_time.push_back((eventPMT->Get<CP::TRealDatum>("TimeFromFirstRF_ns"))->GetValue());
	    PDS_delta_time.push_back((eventPMT->Get<CP::TRealDatum>("DeltaT_ns"))->GetValue());
	    PDS_trigger_type.push_back((eventPMT->Get<CP::TRealDatum>("TriggerType"))->GetValue());
	    PDS_energy.push_back((eventPMT->Get<CP::TRealDatum>("Energy_MeV"))->GetValue());
	    PDS_beam_trigger.push_back((eventPMT->Get<CP::TRealDatum>("BeamTrig"))->GetValue());
	}
    }
}	
    if (tracks) {
	for (CP::TReconObjectContainer::const_iterator t = tracks->begin(); t != tracks->end(); ++t) {
	    TLorentzVector min_hit;
	    TLorentzVector max_hit;
		       CP::THandle<CP::TReconTrack> track = *t;
	    if(track){
	    	
	    if(track->GetFront()->GetPosition().X()>track->GetBack()->GetPosition().X()){
	      min_hit=track->GetFront()->GetPosition();
	      max_hit=track->GetBack()->GetPosition();
	    }else{
	      min_hit=track->GetBack()->GetPosition();
	      max_hit=track->GetFront()->GetPosition();	      
	    }
	     

	    first_hit_X.push_back(min_hit.X());
	    last_hit_X.push_back(max_hit.X());
	    first_hit_Y.push_back(min_hit.Y());
	    last_hit_Y.push_back(max_hit.Y());
	    first_hit_Z.push_back(min_hit.Z());
	    last_hit_Z.push_back(max_hit.Z());
	}
	}
    }
    else {
	std::cout<<"NO TRACKS"<<std::endl;
    }
    
    tree->Fill();
    run = 0;
    evt = 0;
    first_hit_X.clear();
    last_hit_X.clear();
    first_hit_Y.clear();
    last_hit_Y.clear();
    first_hit_Z.clear();
    last_hit_Z.clear();
    TPC_time = 0;
    PDS_RF_time.clear();
    PDS_delta_time.clear();
    PDS_trigger_type.clear();
    PDS_energy.clear();
    PDS_beam_trigger.clear();

    return true;
}
// Called at least once.  If multiple file are open, it will be called
// for each one.   Notice there are two forms...
void CP::TTreeMakerLoop::Finalize(CP::TRootOutput * const output) {
    hfile->Write();
    //gSystem->Exec("mv tracks.root "+fName);
}

