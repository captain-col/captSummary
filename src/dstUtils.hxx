#ifndef dstUtils_hxx_seen
#define dstUtils_hxx_seen

namespace CP {
    namespace utils {
        
        /// Find the mass of a particle based on it's PDG code.  This handles
        /// both elementary particles and nuclei.
        double FindPDGMass(int pdg);

        /// Find the charge of a particle based on it's PDG code.
        int FindPDGCharge(int pdg);
    };
};

#endif
