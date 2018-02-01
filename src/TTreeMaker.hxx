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
  
    TFile* hfile;
    TTree* tree;
  
    std::vector<double> first_hit_X;
    std::vector<double> last_hit_X;

    std::vector<double> first_hit_Y;
    std::vector<double> last_hit_Y;

    std::vector<double> first_hit_Z;
    std::vector<double> last_hit_Z;

    Long64_t TPC_time;
    std::vector<Long64_t> PDS_RF_time;
    std::vector<Long64_t> PDS_delta_time;
    std::vector<int> PDS_trigger_type;
    std::vector<double> PDS_energy;
    
    
};
