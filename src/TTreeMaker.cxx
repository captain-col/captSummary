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

bool inBeamXY(double x, double y) {
    double a = -0.1188;
    double b = -54.9894;
    double width = 54/2;
    double dist=fabs(a*x-y+b)/sqrt(a*a+1);
    if(dist<width ) return true;

    else return false;
}

bool inPeaksZ(double z) {

    if((z>-195 && z<-145) || (z>-510 && z<-460) || (z>-830 && z<-780)) return true;
    else return false;
}

std::vector<std::pair<double,double>> Recalculation(std::vector<double>& minHitX,std::vector<double>& minHitY,std::vector<double>& minHitZ,std::vector<double>& maxHitZ,const std::vector<double>& deltaT,std::vector<double>* pdsEn,std::vector<double>* pdsHQ,const std::vector<int>& coincNumber, const std::vector<double>& eventNumber, std::vector<double>& minEnergy, std::vector<double>& maxEnergy){
    std::vector<double> ene;
    std::vector<double> eneBest;
    std::vector<double> hQ;
    std::vector<double> sumhQ;
    std::vector<double> delta;

    std::vector<double>* minX(new std::vector<double>);
    std::vector<double>* minY(new std::vector<double>);
    std::vector<double>* minZ(new std::vector<double>);
    std::vector<double>* maxZ(new std::vector<double>);

    for(std::size_t i=0;i<minHitZ.size();++i){
	minX->push_back(minHitX[i]);
	minY->push_back(minHitY[i]);
	minZ->push_back(minHitZ[i]);
	maxZ->push_back(maxHitZ[i]);
    }
    
    for(std::size_t i=0;i<deltaT.size();++i){
	if (deltaT[i] > 500e3) continue;
        delta.push_back(1.0*deltaT[i]*1.6/1000.0);
	double sumQ = 0.0;//pdsHQ->at(i);
	double eBest = pdsEn->at(i);

	//std::cout<<"S="<<deltaT.size() << " "<<pdsHQ->size()<<std::endl;

	//std::cout<<"TOF="<<deltaT[i]<<" "<<coincNumber[i]<<" "<<sumQ<<std::endl;
	for (std::size_t j=0;j<deltaT.size();++j) {
	    if (eventNumber[i] == eventNumber[j] && coincNumber[i] == coincNumber[j] && pdsHQ->at(j) != -9999) {
		// std::cout<<"TOF="<<deltaT[i]<<" "<<coincNumber[i]<<" "<<sumQ<<std::endl;
		//std::cout<<"Ev="<<eventNumber[i]<<" "<<eventNumber[j]<<std::endl;
		sumQ += pdsHQ->at(j);
		if ( pdsEn->at(j) > pdsEn->at(i))
		    eBest = pdsEn->at(j);
	    }
	}

	eneBest.push_back(eBest);
	
	// std::cout<<"TOF2="<<deltaT[i]<<" "<<eventNumber[i]<< " " << coincNumber[i]<<" "<<sumQ<<" "<< pdsHQ->at(i)<<" "<<pdsEn->at(i)<<" "<<eBest<<" "<<sumQ<<std::endl;
	// std::cout<<"sumQ="<<coincNumber[i]<<" "<<pdsHQ->at(i)<<" "<<sumQ<<std::endl;
	// ene.push_back((*pdsEn)[i]);
	sumhQ.push_back(sumQ);
    }

    std::vector<std::pair<double,double>> newZ;
    int ntracks = minZ->size();
    int npds = delta.size();
    if(npds>0 ){//&& npds<10){

    
	for(int i=0;i<ntracks;++i){
	    if ((*minZ)[i]>0) {
		newZ.push_back(std::make_pair((*minZ)[i],(*maxZ)[i]));
		ene.push_back(-50);
		hQ.push_back(-50);
		minEnergy.push_back(-50);
		maxEnergy.push_back(-50);
		continue;
	    }
	    if ((*maxZ)[i]>0) {
		newZ.push_back(std::make_pair((*minZ)[i],(*maxZ)[i]));
		ene.push_back(-50);
		hQ.push_back(-50);
		minEnergy.push_back(-50);
		maxEnergy.push_back(-50);
		continue;
	    }
        
	    if ((*minZ)[i]<-1000) {
		newZ.push_back(std::make_pair((*minZ)[i],(*maxZ)[i]));
		ene.push_back(-50);
		hQ.push_back(-50);
		minEnergy.push_back(-50);
		maxEnergy.push_back(-50);
		continue;
	    }
	    if ((*maxZ)[i]<-1000) {
		newZ.push_back(std::make_pair((*minZ)[i],(*maxZ)[i]));
		ene.push_back(-50);
		hQ.push_back(-50);
		minEnergy.push_back(-50);
		maxEnergy.push_back(-50);
		continue;
	    }
        
	    if (fabs((*minZ)[i]-(*maxZ)[i])>9999){
		newZ.push_back(std::make_pair((*minZ)[i],(*maxZ)[i]));
		ene.push_back(-50);
		hQ.push_back(-50);
		minEnergy.push_back(-50);
		maxEnergy.push_back(-50);
		continue;
	    }

	    if (!inBeamXY((*minX)[i],(*minY)[i]) || !inPeaksZ((*minZ)[i])) {
		newZ.push_back(std::make_pair((*minZ)[i],(*maxZ)[i]));
		ene.push_back(-50);
		hQ.push_back(-50);
		minEnergy.push_back(-50);
		maxEnergy.push_back(-50);
		continue;
	    }
	    int npds_d = delta.size();
	    //std::cout<<"Z="<<(*minZ)[i]<<std::endl;
            if(npds_d>0){
		// && npds_d<10){
                double min=9999;
                double maxChar=0.0;
                double minEne=99999.0;
                double maxEne=-999.0;
                int k=0;
		double aveEne = 0.0;
		double aveHQ = 0.0;
		int nPass = 0;
                for(int j=0;j<npds_d;++j){
                    double minimum = 0.0;
		    // std::cout<<"corrZ="<<(*minZ)[i]+delta[j]<<std::endl;
		    // std::cout<<"dT="<<delta[j]*1000.0/1.6<<std::endl;
                    minimum=fabs((*minZ)[i]+delta[j]);

		    if(minimum+170<min && -1.0*(*minZ)[i]>delta[j] && eneBest[j] > 1) {
			k=j;
			min=minimum;
			maxChar = sumhQ[j];
			if (eneBest[j] > maxEne)
			    maxEne = eneBest[j];
			if (eneBest[j] < minEne)
			    minEne = eneBest[j];
			// std::cout<<"dT="<<delta[j]*1000.0/1.6<<std::endl;
			// std::cout<<"corrZ="<<(*minZ)[i]+delta[j]<<std::endl;
			// std::cout<<"sumhQ="<<sumhQ[j]<<std::endl;
			// std::cout<<"neutE="<<eneBest[j]<<std::endl;
			// std::cout<<"coincN="<<coincNumber[j]<<std::endl;
		    }
		    // if(minimum+170<25 && -1.0*(*minZ)[i]>delta[j]) {

		    // }
                    else if((*minZ)[i]+delta[j]>-195 && (*minZ)[i]+delta[j] < -145 && -1.0*(*minZ)[i]>delta[j] && sumhQ[j] > maxChar && eneBest[j] > 1) {
			//else if((*minZ)[i]+delta[j]>-195 && (*minZ)[i]+delta[j] < -145 && -1.0*(*minZ)[i]>delta[j]) {
		    	k=j;
		    	min=minimum;
		    	maxChar = sumhQ[j];
			if (eneBest[j] > maxEne)
			    maxEne = eneBest[j];
			if (eneBest[j] < minEne)
			    minEne = eneBest[j];
			// std::cout<<"Other"<<std::endl;
			// std::cout<<"dT="<<delta[j]*1000.0/1.6<<std::endl;
			// std::cout<<"corrZ="<<(*minZ)[i]+delta[j]<<std::endl;
			// std::cout<<"sumhQ="<<sumhQ[j]<<std::endl;
			// std::cout<<"neutE="<<eneBest[j]<<std::endl;
			// std::cout<<"coincN="<<coincNumber[j]<<std::endl;
			if (eneBest[j] > 0) {
			    aveHQ += sumhQ[j];
			    aveEne += eneBest[j];
			    nPass++;
			}

		    }
                }
                
                /*
		  if((*minZ)[i]>-350 && (*minZ)[i]<-0 && (*minZ)[i]+delta[k]>0){
		  std::cout<<"OUTOF1stPEAK"<<std::endl;
		  std::cout<<(*minZ)[i]<<" ; "<<(*minZ)[i]+delta[k]<<std::endl;}
		*/
                
                if((*minZ)[i]+delta[k]>0){
                    newZ.push_back(std::make_pair((*minZ)[i],(*maxZ)[i]));
		    ene.push_back(-50);
		    hQ.push_back(-50);
		    minEnergy.push_back(-50);
		    maxEnergy.push_back(-50);
                    continue;
		}
		//std::cout<<"corrZ="<<(*minZ)[i]+delta[k]<<std::endl;
             
                newZ.push_back(std::make_pair((*minZ)[i]+delta[k],(*maxZ)[i]+delta[k]));
		ene.push_back((eneBest)[k]);
		hQ.push_back((sumhQ)[k]);
		minEnergy.push_back(minEne);
		maxEnergy.push_back(maxEne);

		// ene.push_back(aveEne/float(nPass));
		// hQ.push_back(aveHQ/float(nPass));
		// std::cout<<"Best"<<std::endl;
		// std::cout<<"corrZ="<<(*minZ)[i]+delta[k]<<std::endl;
		// std::cout<<"sumhQ="<<sumhQ[k]<<std::endl;
		// std::cout<<"neutE="<<eneBest[k]<<std::endl;
		// std::cout<<"ave neutE="<<aveEne/float(nPass)<<std::endl;
		// std::cout<<"ave sumhQ="<<aveHQ/float(nPass)<<std::endl;
		
                for(int j=i+1;j<ntracks;j++){
                    if(fabs((*minZ)[i]-(*minZ)[j])<40){                        
			ene.push_back((eneBest)[k]);
			hQ.push_back((sumhQ)[k]);
                        newZ.push_back(std::make_pair((*minZ)[i+1]+delta[k],(*maxZ)[i+1]+delta[k]));
			minEnergy.push_back(minEne);
			maxEnergy.push_back(maxEne);
                        i++;
                    }
                    else{
			break;
		    }
                }
                
                delta.erase(delta.begin()+k);
                eneBest.erase(eneBest.begin()+k);
                sumhQ.erase(sumhQ.begin()+k);
	    }else{
  
		newZ.push_back(std::make_pair((*minZ)[i],(*maxZ)[i]));
		ene.push_back(-50);
		hQ.push_back(-50);
		minEnergy.push_back(-50);
		maxEnergy.push_back(-50);		
	    }
        
	}
	*pdsEn=ene;
	*pdsHQ=hQ;
	delete minZ;
	delete maxZ;
	return newZ;
    }else {
	delete minZ;
	delete maxZ;
	return newZ;
    }
}



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
    corrected_PDS_energy.clear();
    corrected_PDS_hit_charge.clear();

    first_wire_X.clear();
    first_wire_U.clear();
    first_wire_V.clear();

    hit2DWireNX.clear();
 hit2DChargeX.clear();
 hit2DTimeX.clear();
hit2DWireNU.clear();
 hit2DChargeU.clear();
 hit2DTimeU.clear();
hit2DWireNV.clear();
 hit2DChargeV.clear();
hit2DTimeV.clear();
    
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
    PDS_min_energy.clear(); 	
    PDS_max_energy.clear(); 	
    PDS_hit_charge.clear(); 	
    PDS_beam_trigger.clear();
    PDS_qsum.clear();
    PDS_qmax.clear();
    PDS_event.clear();

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
    tree->Branch("corrected_PDS_energy",&corrected_PDS_energy);
    tree->Branch("corrected_PDS_hit_charge",&corrected_PDS_hit_charge);

    tree->Branch("first_wire_X",&first_wire_X);
    tree->Branch("first_wire_U",&first_wire_U);
    tree->Branch("first_wire_V",&first_wire_V);

    tree->Branch("hits_2D_wireN_Xplane",&hit2DWireNX);
    tree->Branch("hits_2D_Charge_Xplane",&hit2DChargeX);
    tree->Branch("hits_2D_Time_Xplane",&hit2DTimeX);
    tree->Branch("hits_2D_wireN_Uplane",&hit2DWireNU);
    tree->Branch("hits_2D_Charge_Uplane",&hit2DChargeU);
    tree->Branch("hits_2D_Time_Uplane",&hit2DTimeU);
    tree->Branch("hits_2D_wireN_Vplane",&hit2DWireNV);
    tree->Branch("hits_2D_Charge_Vplane",&hit2DChargeV);
    tree->Branch("hits_2D_Time_Vplane",&hit2DTimeV);
    
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
    tree->Branch("PDS_min_energy",&PDS_min_energy);
    tree->Branch("PDS_max_energy",&PDS_max_energy);
    tree->Branch("PDS_beam_trigger",&PDS_beam_trigger);
    tree->Branch("PDS_qSum",&PDS_qsum);
    tree->Branch("PDS_qMax",&PDS_qmax);
    tree->Branch("PDS_event",&PDS_event);

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
    std::vector<int> PDS_coincNumber;

    CP::TChannelInfo& chanInfo = CP::TChannelInfo::Get();
    chanInfo.SetContext(event.GetContext());

    // CP::TGeometryInfo& geomInfo = CP::TGeometryInfo::Get();
    // std::cout<<event.GetContext()<<std::endl;
    // std::cout<<TPC_time<<std::endl;


    
    TString pdsEvent = "";
    if(dataPMT){
	//std::cout<<"PDS size="<<dataPMT->size()<<std::endl;

	for (int i=0; i< int(dataPMT->size()); i++) {
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
		//PDS_qsum.push_back((eventPMT->Get<CP::TRealDatum>("qSum"))->GetValue());	
		//PDS_qmax.push_back((eventPMT->Get<CP::TRealDatum>("qMax"))->GetValue());    
		PDS_event.push_back((eventPMT->Get<CP::TRealDatum>("eventNumber"))->GetValue());
		PDS_hit_charge.push_back((eventPMT->Get<CP::TRealDatum>("HitCharge"))->GetValue());
		PDS_tof.push_back((eventPMT->Get<CP::TRealDatum>("TOF_ns"))->GetValue());
		PDS_coincNumber.push_back((eventPMT->Get<CP::TRealDatum>("CoincNumber"))->GetValue());
		//std::cout<<"EN="<<(eventPMT->Get<CP::TRealDatum>("eventNumber"))->GetValue()<<std::endl;
		//std::cout<<"Co="<<(eventPMT->Get<CP::TRealDatum>("CoincNumber"))->GetValue()<<std::endl;
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

#ifdef correction1
		{
		    double corrected_Zposf = -999.;
		    double corrected_Zposl = -999.;
		    double Zposf;
		    double Zposl;

		    corrected_Zposf = min_hit.Z();
		    corrected_Zposl = max_hit.Z();
		
		    if (min_hit.Z() < -300 && min_hit.Z() > -1000) {
			Zposf = min_hit.Z();
			Zposl = max_hit.Z();
		
			for (size_t k = 0; k<PDS_deltaTs.size(); k++) {
			    if (Zposf + PDS_deltaTs[k]*1.6*1.000000e-3 < 10 && Zposf + PDS_deltaTs[k]*1.6*1.000000e-3 > -320) {
				corrected_Zposf = Zposf + PDS_deltaTs[k]*1.6*1.000000e-3;
				corrected_Zposl = Zposl + PDS_deltaTs[k]*1.6*1.000000e-3;
				break;
			    }		    
			}	       		
		    }

		    corrected_first_hit_Z.push_back(corrected_Zposf);
		    corrected_last_hit_Z.push_back(corrected_Zposl);
		}
#endif
		float length = sqrt( (min_hit.X()-max_hit.X())*(min_hit.X()-max_hit.X()) + (min_hit.Y()-max_hit.Y())*(min_hit.Y()-max_hit.Y()) + (min_hit.Z()-max_hit.Z())*(min_hit.Z()-max_hit.Z()) );

		float energy = length*1.0;
		
		track_length.push_back(length);
		track_energy.push_back(energy);

		float sumCH = 0.0;
		float sumCorrCH = 0.0;
	
		int firstWireX = -999;
		// int firstWireU = -999;
		// int firstWireV = -999;

		double firstChargeX = 0.0;

		CP::THandle<CP::THitSelection> hitsCluster = track->GetHits();
		std::set<CP::THandle<CP::THit>> xHits_set;
		std::set<CP::THandle<CP::THit>> uHits_set;
		std::set<CP::THandle<CP::THit>> vHits_set;
		for (CP::THitSelection::const_iterator hi = hitsCluster->begin(); hi != hitsCluster->end(); ++hi) {
		    CP::THandle<CP::THit> h = *hi;
		    for (int ic = 0; ic < h->GetConstituentCount(); ic++) {
		    	CP::TGeometryId geomId = h->GetConstituent(ic)->GetGeomId();
			if (CP::GeomId::Captain::IsXWire(geomId)) {
			  xHits_set.insert(h->GetConstituent(ic));
			    Int_t wireN = CP::GeomId::Captain::GetWireNumber(geomId);
			
			    if ( wireN > firstWireX) {
				firstWireX = wireN;
				firstChargeX = h->GetConstituent(ic)->GetCharge();
			    }
			    
			    float corrCH = h->GetCharge()*exp(( abs(h->GetPosition().Z()) )/(1.6*80));
			    sumCorrCH += corrCH;
			    sumCH += h->GetConstituent(ic)->GetCharge();
			}
			if (CP::GeomId::Captain::IsUWire(geomId)) {
			  uHits_set.insert(h->GetConstituent(ic));
			}
			if (CP::GeomId::Captain::IsVWire(geomId)) {
			  vHits_set.insert(h->GetConstituent(ic));
			}
		    }
		}

		for(std::set<CP::THandle<CP::THit>>::iterator it = xHits_set.begin();it!=xHits_set.end();++it){
		    CP::TGeometryId geomId = (*it)->GetGeomId();
		    int wireN = CP::GeomId::Captain::GetWireNumber(geomId);
		    double hitCharge = (*it)->GetCharge();
		    double hitTime = (*it)->GetTime();
		    hit2DWireNX.push_back(wireN);
		    hit2DChargeX.push_back(hitCharge);
		    hit2DTimeX.push_back(hitTime);
		  }
		for(std::set<CP::THandle<CP::THit>>::iterator it = uHits_set.begin();it!=uHits_set.end();++it){
		    CP::TGeometryId geomId = (*it)->GetGeomId();
		    int wireN = CP::GeomId::Captain::GetWireNumber(geomId);
		    double hitCharge = (*it)->GetCharge();
		    double hitTime = (*it)->GetTime();
		    hit2DWireNU.push_back(wireN);
		    hit2DChargeU.push_back(hitCharge);
		    hit2DTimeU.push_back(hitTime);
		  }
		for(std::set<CP::THandle<CP::THit>>::iterator it = vHits_set.begin();it!=vHits_set.end();++it){
		    CP::TGeometryId geomId = (*it)->GetGeomId();
		    int wireN = CP::GeomId::Captain::GetWireNumber(geomId);
		    double hitCharge = (*it)->GetCharge();
		    double hitTime = (*it)->GetTime();
		    hit2DWireNV.push_back(wireN);
		    hit2DChargeV.push_back(hitCharge);
		    hit2DTimeV.push_back(hitTime);
		  }

		first_wire_X.push_back(firstWireX);
		first_hit_charge_X.push_back(firstChargeX);
		
		float tdQdx = -1.0;
		float tdEdx = -1.0;
		// float tdEdxBox = -1.0;

		if (length > 10) {
		    tdQdx = (sumCorrCH)/(length);
		    float A_B = 0.8;
		    float W_ion = 23.6e-6;// * unit::MeV;
		    float k_B = 0.0486;// * (unit::g/unit::cm)/(unit::MeV/unit::cm);
		    float epsilon = 0.5;

		    // float alpha = 0.93;
		    // float beta = 0.3; // *unit::cm/unit::MeV;
		    
		    tdEdx = tdQdx/(A_B/W_ion - k_B*tdQdx/epsilon);
		    // tdEdxBox = (exp(beta*W_ion*(tdQdx)) - alpha)/beta;
		    
		    dQdx.push_back(tdQdx);
		    dEdx.push_back(tdEdx);		
		    
		}
		
	    }
	}
    }
    else {
	std::cout<<"NO TRACKS"<<std::endl;
    }
#define correction2
#ifdef correction2

    std::vector<std::pair<double,double>> newZ;
    int npds = PDS_deltaTs.size();

    if (npds>0) {
	// && npds<10){
	std::vector<double>* pdsEn(new std::vector<double>);
	*pdsEn = PDS_energy;
	std::vector<double>* pdsHQ(new std::vector<double>);
	*pdsHQ = PDS_hit_charge;
	newZ=Recalculation(first_hit_X,first_hit_Y,first_hit_Z,last_hit_Z,PDS_deltaTs,pdsEn,pdsHQ,PDS_coincNumber, PDS_event, PDS_min_energy, PDS_max_energy);
	
	for(std::size_t j=0; j<newZ.size();++j){
	    corrected_first_hit_Z.push_back(newZ[j].first);
	    corrected_last_hit_Z.push_back(newZ[j].second);
	    corrected_PDS_energy.push_back((*pdsEn)[j]);
	    corrected_PDS_hit_charge.push_back((*pdsHQ)[j]);
	    
	}
 	delete pdsEn;
    }
    else {
	for (std::size_t j=0; j<first_hit_Z.size();++j) {
	    corrected_first_hit_Z.push_back(first_hit_Z[j]);
	    corrected_last_hit_Z.push_back(last_hit_Z[j]);
	    corrected_PDS_energy.push_back(-50);
	}
    }    

#endif
    
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
    corrected_PDS_energy.clear();
    corrected_PDS_hit_charge.clear();

    first_wire_X.clear();
    first_wire_U.clear();
    first_wire_V.clear();

        hit2DWireNX.clear();
 hit2DChargeX.clear();
 hit2DTimeX.clear();
hit2DWireNU.clear();
 hit2DChargeU.clear();
 hit2DTimeU.clear();
hit2DWireNV.clear();
 hit2DChargeV.clear();
hit2DTimeV.clear();
    
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
    PDS_min_energy.clear();
    PDS_max_energy.clear();
    PDS_hit_charge.clear();
    PDS_coincNumber.clear();
    PDS_beam_trigger.clear();
    PDS_qsum.clear();
    PDS_qmax.clear();
    PDS_event.clear();

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

