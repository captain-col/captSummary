#include "TG4VerticesModule.hxx"
#include "dstUtils.hxx"

#include <HEPUnits.hxx>
#include <HEPConstants.hxx>

#include <TCaptLog.hxx>
#include <TG4PrimaryVertex.hxx>
#include <TG4PrimaryParticle.hxx>

#include <TLorentzVector.h>

#include <cstdlib>
#include <set>


ClassImp(CP::TG4VerticesModule);
ClassImp(CP::TG4VerticesModule::TG4Vertex);
ClassImp(CP::TG4VerticesModule::TG4Particle);

//
// TG4VerticesModule
//
CP::TG4VerticesModule::TG4VerticesModule(
    const char *name, const char *title) {
    SetNameTitle(name, title);
    SetEnabled();
}

CP::TG4VerticesModule::~TG4VerticesModule() {}

void CP::TG4VerticesModule::InitializeModule() {} 

void CP::TG4VerticesModule::InitializeBranches() {
    GetOutputTree()->Branch(
        "Vertex",
        "std::vector<CP::TG4VerticesModule::TG4Vertex>", 
        &fVertices, 
        GetBufferSize(), GetSplitLevel());
    GetOutputTree()->Branch(
        "Particle",
        "std::vector<CP::TG4VerticesModule::TG4Particle>", 
        &fParticles, 
        GetBufferSize(), GetSplitLevel());
}

Bool_t CP::TG4VerticesModule::ProcessFirstEvent(CP::TEvent& event) {
    bool IsDataFile = event.GetContext().IsDetector();
    // This only gets run on MC files.
    if (IsDataFile) throw CP::EDataFile();
    return true;
}

namespace {
    bool vtxLessThan(const CP::TG4VerticesModule::TG4Vertex& a,
                      const CP::TG4VerticesModule::TG4Vertex& b) {
        return a.VertexId < b.VertexId;
    } 

    bool trkLessThan(const CP::TG4VerticesModule::TG4Particle& a,
                      const CP::TG4VerticesModule::TG4Particle& b) {
        return a.TrackId < b.TrackId;
    } 
}

bool CP::TG4VerticesModule::FillTree(CP::TEvent& event) {
    if (!event.GetContext().IsMC()) {
        CaptInfo("Event Context reports event is non-MC.");
        return false; // disable module
    }

    CP::THandle<CP::TG4PrimaryVertexContainer> vertices
        = event.Get<CP::TG4PrimaryVertexContainer>("truth/G4PrimVertex00");

    fVertices.clear();
    int vertexCount = 0;
    for (CP::TG4PrimaryVertexContainer::iterator v = vertices->begin();
         v != vertices->end(); ++v) {
        CP::TG4VerticesModule::TG4Vertex vertex;
        vertex.VertexId = ++vertexCount;
        vertex.Position = v->GetPosition();
        vertex.Particles.clear();
        CP::TG4PrimaryParticleContainer& particles = v->GetPrimaryParticles();
        for (CP::TG4PrimaryParticleContainer::iterator p = particles.begin();
             p != particles.end(); ++p) {
            CP::TG4VerticesModule::TG4Particle particle;
            particle.TrackId = p->GetTrackId();
            particle.VertexId = vertex.VertexId;
            particle.PDG = p->GetPDGCode();
            particle.Mass = CP::utils::FindPDGMass(particle.PDG);
            particle.Charge = CP::utils::FindPDGCharge(particle.PDG);
            particle.Momentum = p->GetMomentum();
            vertex.Particles.push_back(particle.TrackId);
            fParticles.push_back(particle);
        }
        fVertices.push_back(vertex);
    }

    // Sort the vector array so trajectories are in order of ascending Id.
    std::sort(fVertices.begin(), fVertices.end(), vtxLessThan);
    std::sort(fParticles.begin(), fParticles.end(), trkLessThan);
	
    return true;
}

CP::TG4VerticesModule::TG4Vertex::TG4Vertex() 
    : VertexId(-1) {
    Particles.clear();
}

CP::TG4VerticesModule::TG4Particle::TG4Particle()
    : TrackId(-1), VertexId(-1), PDG(-1) {
}
