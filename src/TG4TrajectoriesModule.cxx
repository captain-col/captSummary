#include "TG4TrajectoriesModule.hxx"
#include "dstUtils.hxx"

#include <HEPUnits.hxx>
#include <HEPConstants.hxx>

#include <TCaptLog.hxx>
#include <TG4Trajectory.hxx>

#include <TLorentzVector.h>

#include <cstdlib>
#include <set>


ClassImp(CP::TG4TrajectoriesModule);
ClassImp(CP::TG4TrajectoriesModule::TG4Trajectory);

//
// TG4TrajectoriesModule
//
CP::TG4TrajectoriesModule::TG4TrajectoriesModule(
    const char *name, const char *title)
    : fSaveParentChain(kTRUE) {
    SetNameTitle(name, title);
    SetEnabled();
}

CP::TG4TrajectoriesModule::~TG4TrajectoriesModule() {}

Bool_t CP::TG4TrajectoriesModule::Configure(const std::string& option) {
    if (option == "saveall" || option == "SaveAll") {
        CaptInfo("Setting TG4TrajectoriesModule to Save All Trajectories");
        return kTRUE;
    }
    return kFALSE;
}

void CP::TG4TrajectoriesModule::InitializeModule() {} 

void CP::TG4TrajectoriesModule::InitializeBranches() {
    GetOutputTree()->Branch(
        "Trajectory",
        "std::vector<CP::TG4TrajectoriesModule::TG4Trajectory>", 
        &fTrajectories, 
        GetBufferSize(), GetSplitLevel());
}

Bool_t CP::TG4TrajectoriesModule::ProcessFirstEvent(CP::TEvent& event) {
    bool IsDataFile = event.GetContext().IsDetector();
    // This only gets run on MC files.
    if (IsDataFile) throw CP::EDataFile();
    return true;
}

void CP::TG4TrajectoriesModule::FillSaveList(
    CP::THandle<CP::TG4TrajectoryContainer> trajectories) {
	
    // First go through the entire trajectory container and mark the
    // trajectories that should be saved.
    for (CP::TG4TrajectoryContainer::iterator t = trajectories->begin();
        t != trajectories->end(); ++t) {
        const CP::TG4Trajectory& traj = t->second;
        if (!SaveTraj(traj)) continue;
        int trajId = traj.GetTrackId();
        if (trajId < 1) continue;
        // Make sure the primary is saved.
        int primId = trajectories->GetPrimaryId(trajId);
        if (primId > 0) fSaveList.insert(primId);
        // Now save the chain of every trajectory.
        while (trajId > 0) {
            fSaveList.insert(trajId);
            if (!fSaveParentChain) break;
            trajId = trajectories->GetTrajectory(trajId)->GetParentId();
        }
    }
}

namespace {
    bool trajLessThan(const CP::TG4TrajectoriesModule::TG4Trajectory& a,
                      const CP::TG4TrajectoriesModule::TG4Trajectory& b) {
        return a.TrajId < b.TrajId;
    } 
}

bool CP::TG4TrajectoriesModule::FillTree(CP::TEvent& event) {
    if (!event.GetContext().IsMC()) {
        CaptInfo("Event Context reports event is non-MC.");
        return false; // disable module
    }

    if (! event.Has<CP::TG4TrajectoryContainer>("truth/G4Trajectories")) {
        CaptVerbose("No trajectories container found, skipping event.");
        return true; // this happens occasionally, so don't disable the module
    }
	
    CP::THandle<CP::TG4TrajectoryContainer> trajectories
        = event.Get<CP::TG4TrajectoryContainer>("truth/G4Trajectories");
	
    fSaveList.clear();
    FillSaveList(trajectories);

    // Reserve space for the trajectories.
    fTrajectories.resize(fSaveList.size());

    CP::TG4TrajectoryContainer::iterator traj_it = trajectories->begin();
    CP::TG4TrajectoryContainer::iterator traj_end = trajectories->end();
    CP::TG4Trajectory* traj;
    CP::TG4TrajectoriesModule::TG4Trajectory* trajToFill;
	
    int filledTraj = 0;
    for (; traj_it != traj_end; ++traj_it) {
        traj = &traj_it->second;
		
        // Check if this trajectory is worth saving
        if (fSaveList.count( traj->GetTrackId() ) == 0 ) continue;
		
        // Protect from empty trajectories (bug encountered while testing)
        if (traj->GetTrackId() == 0 ) {
            CaptError("Found TG4Trajectory with TrajId = 0");
            continue;
        }
		
        //
        // Do the saving
		
        trajToFill = new(&(fTrajectories[filledTraj])) TG4Trajectory();
        ++ filledTraj;
		
        trajToFill->TrajId = traj->GetTrackId();
        trajToFill->ParentId = traj->GetParentId();
        trajToFill->PrimaryId = trajectories->GetPrimaryId(trajToFill->TrajId);
        trajToFill->PDG = traj->GetPDGEncoding();
        trajToFill->Mass = CP::utils::FindPDGMass(traj->GetPDGEncoding());
        trajToFill->Charge = CP::utils::FindPDGCharge(traj->GetPDGEncoding());
        FillPoints(traj, trajToFill);
    }
		
    // Sort the vector array so trajectories are in order of ascending Id.
    std::sort(fTrajectories.begin(), fTrajectories.end(), trajLessThan);
	
    return true;
}

bool 
CP::TG4TrajectoriesModule::SaveTraj(const CP::TG4Trajectory& traj) const {
    // Never save a trajectory with a zero trajectory id.
    if (traj.GetTrackId() == 0) return false;
    
    //Always save primary trajectories
    if (traj.GetParentId() == 0) return true;
	
    // Non-primary trajectories, even those in detectors are only saved if
    // their length is sufficiently long
    TLorentzVector length 
        = traj.GetInitialPosition() - traj.GetFinalPosition();
    if (length.Vect().Mag() < GetMinimumTrajectoryLengthToSave()) return false;

    return true;
}

void CP::TG4TrajectoriesModule::FillPoints(
    CP::TG4Trajectory *const traj,
    CP::TG4TrajectoriesModule::TG4Trajectory* trajToFill) {

    trajToFill->Position.clear();
    trajToFill->Momentum.clear();
    trajToFill->Region.clear();
        
    // Get the points from the trajectory.
    const CP::TG4Trajectory::Points& points = traj->GetTrajectoryPoints();
    if (points.empty()) {
        CaptError("Found TG4Trajectory with no trajectory points");
        return;
    }

    for (CP::TG4Trajectory::Points::const_iterator p = points.begin();
         p != points.end(); ++p) {
        std::string volumeName = p->GetVolumeName();

        // Figure out the volume number for this.  Be careful with the
        // conditions since the first one that's met determine the volume
        // number.
        int thisVolume = 0;
        if (volumeName.find("Drift") != std::string::npos) {
            thisVolume = 1;
        }
        else if (volumeName.find("Cryostat") != std::string::npos) {
            thisVolume = 50;
        }
        else {
            thisVolume = 999;
        }

        // Find the next volume name.
        CP::TG4Trajectory::Points::const_iterator n = p;
        int nextVolume = 0;
        if (++n != points.end()) {
            volumeName = n->GetVolumeName();
            if (volumeName.find("Drift") != std::string::npos) {
                nextVolume = 1;
            }
            else if (volumeName.find("Cryostat") != std::string::npos) {
                nextVolume = 50;
            }
            else {
                nextVolume = 999;
            }
        }

        trajToFill->Position.push_back(p->GetPosition());
        trajToFill->Momentum.push_back(p->GetMomentum());
        trajToFill->Region.push_back(1000000 + 1000*nextVolume + thisVolume);
    }    
}

//
// G4 Trajectory
// =====================================================================
//
CP::TG4TrajectoriesModule::TG4Trajectory::TG4Trajectory()
    : TrajId(-1),
      ParentId(-1),
      PrimaryId(-1),
      PDG(-1),
      Mass(-999999),
      Charge(-999999) {
    Position.clear();
    Momentum.clear();
    Region.clear();
}
