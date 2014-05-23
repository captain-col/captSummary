#include "dstUtils.hxx"

#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <HEPConstants.hxx>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TGeoElement.h>

double CP::utils::FindPDGMass(int pdg) {
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdg);
    if (particle)  {
        // In ROOT the TParticlePDG mass is given in GeV/c2, but we
        // want MeV/c2
        return particle->Mass() * unit::GeV;
    }

    int i = 0;
    int a = pdg/10 % 1000;
    int z = pdg/10000 % 1000;
    const TGeoElementRN* rn 
        = TGeoElement::GetElementTable()->GetElementRN(a,z,i);
    if (!rn) {
        // ROOT TGeoElementRN stores mass in atomic mass units, but we want
        // MeV/c2
        return rn->MassNo() * unit::amu_c2;
    }
    CaptError("PDG number " << pdg
              << " does not exist in database."
              << " Cannot fill Mass and Charge values.");
    return a * unit::amu_c2;
}

int CP::utils::FindPDGCharge(int pdg) {
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdg);
    if (particle)  {
        return static_cast<int>(particle->Charge());
    }

    int i = 0;
    int a = pdg/10 % 1000;
    int z = pdg/10000 % 1000;
    const TGeoElementRN* rn 
        = TGeoElement::GetElementTable()->GetElementRN(a,z,i);
    if (!rn) {
        return static_cast<int>( rn->AtomicNo() * 3);
    }
    return z * 3;
}
