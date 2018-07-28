#include <memory>

// Include to get Event Loop.
#include <eventLoop.hxx>
#include <TReconTrack.hxx>
#include <TReconHit.hxx>
#include <TRealDatum.hxx>
#include <CaptGeomId.hxx>
#include <TChannelInfo.hxx>
#include "TTreeMaker.hxx"
#include "TGeometryInfo.hxx"
#include <TG4PrimaryVertex.hxx>
#include <TG4Trajectory.hxx>

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
    corrected_first_hit_Z.clear();
    corrected_last_hit_Z.clear();

    first_wire_X.clear();
    first_wire_U.clear();
    first_wire_V.clear();
    
    first_hit_charge_X.clear();
    first_hit_charge_U.clear();
    first_hit_charge_V.clear();

    track_length.clear();
    track_energy.clear();

    dQdx.clear();
    dEdx.clear();
    
    TPC_time = 0;
    PDS_RF_time.clear();
    PDS_delta_time.clear();
    PDS_trigger_type.clear();
    PDS_energy.clear();
    PDS_beam_trigger.clear();

    truth_vertex_X.clear();
    truth_vertex_Y.clear();
    truth_vertex_Z.clear();

    truth_particle_PDG.clear();
    truth_particle_E.clear();
    truth_particle_px.clear();
    truth_particle_py.clear();
    truth_particle_pz.clear();

    truth_trajectory_length.clear();

    truth_trajectory_first_X.clear();
    truth_trajectory_first_Y.clear();
    truth_trajectory_first_Z.clear();

    truth_trajectory_last_X.clear();
    truth_trajectory_last_Y.clear();
    truth_trajectory_last_Z.clear();

    tree = NULL;
}

CP::TTreeMakerLoop::~TTreeMakerLoop() {}

void CP::TTreeMakerLoop::Initialize(void) {
    tree = new TTree("tracks","");
    tree->Branch("run",&run,"run/I");
    tree->Branch("event",&evt,"event/I");

    tree->Branch("first_hit_X",&first_hit_X);
    tree->Branch("last_hit_X",&last_hit_X);
    tree->Branch("first_hit_Y",&first_hit_Y);
    tree->Branch("last_hit_Y",&last_hit_Y);
    tree->Branch("first_hit_Z",&first_hit_Z);
    tree->Branch("last_hit_Z",&last_hit_Z);
    tree->Branch("corrected_first_hit_Z",&corrected_first_hit_Z);
    tree->Branch("corrected_last_hit_Z",&corrected_last_hit_Z);

    tree->Branch("first_wire_X",&first_wire_X);
    tree->Branch("first_wire_U",&first_wire_U);
    tree->Branch("first_wire_V",&first_wire_V);
    
    tree->Branch("first_hit_charge_X",&first_hit_charge_X);
    tree->Branch("first_hit_charge_U",&first_hit_charge_U);
    tree->Branch("first_hit_charge_V",&first_hit_charge_V);
    
    tree->Branch("track_length",&track_length);
    tree->Branch("track_energy",&track_energy);

    tree->Branch("dQdx",&dQdx);
    tree->Branch("dEdx",&dEdx);
    
    tree->Branch("TPC_time",&TPC_time,"TPC_time/L");
    tree->Branch("PDS_RF_time","std::vector<Long64_t>",&PDS_RF_time);
    tree->Branch("PDS_delta_time","std::vector<double>",&PDS_delta_time);
    tree->Branch("PDS_trigger_type",&PDS_trigger_type);
    tree->Branch("PDS_energy",&PDS_energy);
    tree->Branch("PDS_beam_trigger",&PDS_beam_trigger);

    tree->Branch("truth_vertex_X",&truth_vertex_X);
    tree->Branch("truth_vertex_Y",&truth_vertex_Y);
    tree->Branch("truth_vertex_Z",&truth_vertex_Z);

    tree->Branch("truth_particle_PDG",&truth_particle_PDG);
    tree->Branch("truth_particle_E",&truth_particle_E);
    tree->Branch("truth_particle_px",&truth_particle_px);
    tree->Branch("truth_particle_py",&truth_particle_py);
    tree->Branch("truth_particle_pz",&truth_particle_pz);

    tree->Branch("truth_trajectory_length",&truth_trajectory_length);

    tree->Branch("truth_trajectory_first_X",&truth_trajectory_first_X);
    tree->Branch("truth_trajectory_first_Y",&truth_trajectory_first_Y);
    tree->Branch("truth_trajectory_first_Z",&truth_trajectory_first_Z);

    tree->Branch("truth_trajectory_last_X",&truth_trajectory_last_X);
    tree->Branch("truth_trajectory_last_Y",&truth_trajectory_last_Y);
    tree->Branch("truth_trajectory_last_Z",&truth_trajectory_last_Z);
}

bool CP::TTreeMakerLoop::operator () (CP::TEvent& event) {

    CP::THandle<CP::TDataVector> dataPMT =
	event.Get<CP::TDataVector>("~/pmtData");
    CP::THandle<CP::TReconObjectContainer> tracks =
	event.Get<CP::TReconObjectContainer>("~/fits/TCaptainRecon/final");
    CP::THandle<CP::TG4PrimaryVertexContainer> vertex =
	event.Get<CP::TG4PrimaryVertexContainer>("~/truth/G4PrimVertex00");
    CP::THandle<CP::TG4TrajectoryContainer> trajs =
	event.Get<CP::TG4TrajectoryContainer>("~/truth/G4Trajectories");

    run = event.GetContext().GetRun();
    evt = event.GetContext().GetEvent();

    TPC_time = event.GetTimeStamp();

    std::vector<double> PDS_deltaTs;

    CP::TChannelInfo& chanInfo = CP::TChannelInfo::Get();
    chanInfo.SetContext(event.GetContext());

    //CP::TGeometryInfo& geomInfo = CP::TGeometryInfo::Get();

    
    TString pdsEvent = "";
    if(dataPMT){ 
	for (u_int i=0; i<dataPMT->size(); i++) {
	    pdsEvent.Form("~/pmtData/PDSEvent_%d",i);
	    CP::THandle<CP::TEvent> eventPMT = event.Get<CP::TEvent>(pdsEvent);
	    if (!eventPMT) {
		std::cout<<"NO PMT EVENT"<<std::endl;
	    }
	    else {
		//double TOF = (eventPMT->Get<CP::TRealDatum>("TOF_ns"))->GetValue();
		PDS_RF_time.push_back((eventPMT->Get<CP::TRealDatum>("TimeFromFirstRF_ns"))->GetValue());
		PDS_delta_time.push_back((eventPMT->Get<CP::TRealDatum>("DeltaT_ns"))->GetValue());
		PDS_deltaTs.push_back((eventPMT->Get<CP::TRealDatum>("DeltaT_ns"))->GetValue());
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

		double corrected_Zposf = -999.;
		double corrected_Zposl = -999.;
		double Zposf;
		double Zposl;
		
		if (min_hit.Z() < -300 && min_hit.Z() > -1000) {
		    Zposf = min_hit.Z();
		    corrected_Zposf = min_hit.Z();
		    Zposl = max_hit.Z();
		    corrected_Zposl = max_hit.Z();
		
		    for (size_t k = 0; k<PDS_deltaTs.size(); k++) {
			if (Zposf + PDS_deltaTs[k]*1.6*1000000 < 10 && Zposf + PDS_deltaTs[k]*1.6*1000000 > -320) {
			    corrected_Zposf = Zposf + PDS_deltaTs[k]*1.6*1000000;
			    corrected_Zposl = Zposl + PDS_deltaTs[k]*1.6*1000000;
			    break;
			}		    
		    }	       		
		}

		corrected_first_hit_Z.push_back(corrected_Zposf);
		corrected_last_hit_Z.push_back(corrected_Zposl);

		float length = sqrt( (min_hit.X()-max_hit.X())*(min_hit.X()-max_hit.X()) + (min_hit.Y()-max_hit.Y())*(min_hit.Y()-max_hit.Y()) + (min_hit.Z()-max_hit.Z())*(min_hit.Z()-max_hit.Z()) );

		float energy = length*1.0;
		
		track_length.push_back(length);
		track_energy.push_back(energy);

		float sumCH = 0.0;
		float sumCorrCH = 0.0;
	
		int firstWireX = -999;
		int firstWireU = -999;
		int firstWireV = -999;

		double firstChargeX = 0.0;

		CP::THandle<CP::THitSelection> hitsCluster = track->GetHits();
		for (CP::THitSelection::const_iterator hi = hitsCluster->begin(); hi != hitsCluster->end(); ++hi) {
		    CP::THandle<CP::THit> h = *hi;
		    for (int ic = 0; ic < h->GetConstituentCount(); ic++) {			
		    	CP::TGeometryId geomId = h->GetConstituent(ic)->GetGeomId();
			if (CP::GeomId::Captain::IsXWire(geomId)) {

			    Int_t wireN = CP::GeomId::Captain::GetWireNumber(geomId);
			
			    if ( wireN > firstWireX) {
				firstWireX = wireN;
				firstChargeX = h->GetConstituent(ic)->GetCharge();
			    }
			    
			    float corrCH = h->GetCharge()*exp(( abs(h->GetPosition().Z()) )/(1.6*80));
			    sumCorrCH += corrCH;
			    sumCH += h->GetConstituent(ic)->GetCharge();
			}
		    }
		}

		first_wire_X.push_back(firstWireX);
		first_hit_charge_X.push_back(firstChargeX);
		
		float tdQdx = -1.0;
		float tdEdx = -1.0;
		float tdEdxBox = -1.0;

		if (length > 10) {
		    tdQdx = (sumCorrCH)/(length);
		    float A_B = 0.8;
		    float W_ion = 23.6e-6;// * unit::MeV;
		    float k_B = 0.0486;// * (unit::g/unit::cm)/(unit::MeV/unit::cm);
		    float epsilon = 0.5;

		    float alpha = 0.93;
		    float beta = 0.3; // *unit::cm/unit::MeV;
		    
		    tdEdx = tdQdx/(A_B/W_ion - k_B*tdQdx/epsilon);
		    tdEdxBox = (exp(beta*W_ion*(tdQdx)) - alpha)/beta;
		    
		    dQdx.push_back(tdQdx);
		    dEdx.push_back(tdEdx);		
		    
		 }
		
	    }
	}
    }
    else {
	std::cout<<"NO TRACKS"<<std::endl;
    }

    if (vertex) {
	for (CP::TG4PrimaryVertexContainer::const_iterator v = vertex->begin(); v != vertex->end(); ++v) {
	    truth_vertex_X.push_back(v->GetPosition().X());
	    truth_vertex_Y.push_back(v->GetPosition().Y());
	    truth_vertex_Z.push_back(v->GetPosition().Z());

	    for (CP::TG4PrimaryParticleContainer::const_iterator p = (v->GetPrimaryParticles()).begin(); p != (v->GetPrimaryParticles()).end(); ++p) {
		truth_particle_PDG.push_back(p->GetPDGCode());
		truth_particle_E.push_back(p->GetMomentum().E() - p->GetMomentum().M());
		truth_particle_px.push_back(p->GetMomentum().Px());
		truth_particle_py.push_back(p->GetMomentum().Py());
		truth_particle_pz.push_back(p->GetMomentum().Pz());
		if (trajs) {
		    CP::THandle<CP::TG4Trajectory> t = trajs->GetTrajectory(p->GetTrackId());
		    TLorentzVector t1 = t->GetInitialPosition();
		    TLorentzVector t2 = t->GetFinalPosition();
		    for (CP::TG4TrajectoryContainer::const_iterator ti = trajs->begin(); ti != trajs->end(); ++ti) {
		        CP::TG4Trajectory traj = (*ti).second;
			if (traj.GetPDGEncoding() == 2212) {
			    if (traj.GetInitialPosition().X() < -500) continue;
			    if (traj.GetInitialPosition() == t2) {
				t2 = traj.GetFinalPosition();
			    }
			}
		    }

		    // Define the edges of the detector.
		    // If trajectory escapes, cut it off.
		    float t1X = 600;
		    float t2X = -600;
		    float t1Y = 600;
		    float t2Y = -600;
		    float t1Z = 100;
		    float t2Z = -400;

		    if (t1.X() < 500 && t1.X() > -500)
			t1X = t1.X();
		    if (t2.X() < 500 && t2.X() > -500)
			t2X = t2.X();
		    if (t1.Y() < 500 && t1.Y() > -500)
			t1Y = t1.Y();
		    if (t2.Y() < 500 && t2.Y() > -500)
			t2Y = t2.Y();
		    if (t1.Z() < 0 && t1.Z() > -350)
			t1Z = t1.Z();
		    if (t2.Z() < 0 && t2.Z() > -350)
			t2Z = t2.Z();
		    
		    truth_trajectory_first_X.push_back(t1X);
		    truth_trajectory_first_Y.push_back(t1Y);
		    truth_trajectory_first_Z.push_back(t1Z);

		    truth_trajectory_last_X.push_back(t2X);
		    truth_trajectory_last_Y.push_back(t2Y);
		    truth_trajectory_last_Z.push_back(t2Z);

		    float truth_length = sqrt( (t1X-t2X)*(t1X-t2X) + (t1Y-t2Y)*(t1Y-t2Y) + (t1Z-t2Z)*(t1Z-t2Z) );
		    truth_trajectory_length.push_back(truth_length);
		}
	    }	
	}
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
    corrected_first_hit_Z.clear();
    corrected_last_hit_Z.clear();

    first_wire_X.clear();
    first_wire_U.clear();
    first_wire_V.clear();
    
    first_hit_charge_X.clear();
    first_hit_charge_U.clear();
    first_hit_charge_V.clear();

    track_length.clear();
    track_energy.clear();

    dQdx.clear();
    dEdx.clear();

    TPC_time = 0;
    PDS_RF_time.clear();
    PDS_delta_time.clear();
    PDS_trigger_type.clear();
    PDS_energy.clear();
    PDS_beam_trigger.clear();

    truth_vertex_X.clear();
    truth_vertex_Y.clear();
    truth_vertex_Z.clear();

    truth_particle_PDG.clear();
    truth_particle_E.clear();
    truth_particle_px.clear();
    truth_particle_py.clear();
    truth_particle_pz.clear();

    truth_trajectory_length.clear();

    truth_trajectory_first_X.clear();
    truth_trajectory_first_Y.clear();
    truth_trajectory_first_Z.clear();

    truth_trajectory_last_X.clear();
    truth_trajectory_last_Y.clear();
    truth_trajectory_last_Z.clear();

    return false;
}
// Called at least once.  If multiple file are open, it will be called
// for each one.   Notice there are two forms...
void CP::TTreeMakerLoop::Finalize(CP::TRootOutput * const output) {
    
}

