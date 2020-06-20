#include <memory>

// Include to get Event Loop.
#include <eventLoop.hxx>
#include <TEventLoopFunction.hxx>
#include <TReconTrack.hxx>
#include <TReconHit.hxx>
#include <TRealDatum.hxx>

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
#include <TChain.h>
#include <TSystem.h>


namespace CP {
    class TTreeMakerLoop;
}

class CP::TTreeMakerLoop: public CP::TEventLoopFunction {
public:
    TTreeMakerLoop();
    virtual ~TTreeMakerLoop();
    bool operator () (CP::TEvent& event);
    virtual void Initialize(void);
    void Finalize(CP::TRootOutput * const output);
private:
    TTree* tree;

    Int_t run;
    Int_t evt;
  
    std::vector<double> first_hit_X;
    std::vector<double> last_hit_X;

    std::vector<double> first_hit_Y;
    std::vector<double> last_hit_Y;

    std::vector<double> first_hit_Z;
    std::vector<double> last_hit_Z;

    std::vector<double> corrected_first_hit_Z;
    std::vector<double> corrected_last_hit_Z;
    std::vector<double> corrected_PDS_energy;     
    std::vector<double> corrected_PDS_hit_charge;     

    std::vector<int> first_wire_X;
    std::vector<int> first_wire_U;
    std::vector<int> first_wire_V;

  
  std::vector<int> trackHitSyncX; // for each hit has index of starting point of correcponding track from first_wire_X array
   std::vector<int> trackHitSyncU;// for each hit has index of starting point of correcponding track from first_wire_U array
   std::vector<int> trackHitSyncV;// for each hit has index of starting point of correcponding track from first_wire_V array
  std::vector<int> hit2DWireNX;
  std::vector<double> hit2DChargeX;
  std::vector<double> hit2DTimeX;
   std::vector<int> hit2DWireNU;
  std::vector<double> hit2DChargeU;
  std::vector<double> hit2DTimeU;
   std::vector<int> hit2DWireNV;
  std::vector<double> hit2DChargeV;
  std::vector<double> hit2DTimeV;

    std::vector<double> first_hit_charge_X;
    std::vector<double> first_hit_charge_U;
    std::vector<double> first_hit_charge_V;

    std::vector<double> track_length;
    std::vector<double> track_energy;

    std::vector<double> dQdx;
    std::vector<double> dEdx;

    Long64_t TPC_time;
    std::vector<Long64_t> PDS_RF_time;
    std::vector<double> PDS_delta_time;
    std::vector<int> PDS_trigger_type;
    std::vector<double> PDS_energy;
    std::vector<double> PDS_min_energy;
    std::vector<double> PDS_max_energy;
    std::vector<int> PDS_beam_trigger;
    std::vector<double> PDS_qsum;
    std::vector<double> PDS_qmax;
    std::vector<double> PDS_event;    
    std::vector<double> PDS_hit_charge;    
    std::vector<double> PDS_tof;
  std::vector<int> PDS_coincidence;

  std::vector<int> true_traj_Id;
  std::vector<int> track_traj_Id;

    std::vector<double> truth_vertex_X;
    std::vector<double> truth_vertex_Y;
    std::vector<double> truth_vertex_Z;

    std::vector<int> truth_particle_PDG;
    std::vector<double> truth_particle_E;
    std::vector<double> truth_particle_px;
    std::vector<double> truth_particle_py;
    std::vector<double> truth_particle_pz;
    
    std::vector<double> truth_trajectory_length;

    std::vector<double> truth_trajectory_first_X;
    std::vector<double> truth_trajectory_first_Y;
    std::vector<double> truth_trajectory_first_Z;

    std::vector<double> truth_trajectory_last_X;
    std::vector<double> truth_trajectory_last_Y;
    std::vector<double> truth_trajectory_last_Z;

      std::vector<double> truth_trajectory_last_X_corr;
    std::vector<double> truth_trajectory_last_Y_corr;
    std::vector<double> truth_trajectory_last_Z_corr;

  std::vector<int> truth_particle_ParentId;
  std::vector<int> truth_StartInDriftVolume;
 

    TChain *chainz = new TChain("new_tree");   
    Int_t entryPDS = 0;
    Float_t kePDS = 0;

    
};
