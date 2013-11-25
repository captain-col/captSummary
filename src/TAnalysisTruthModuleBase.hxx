#ifndef TAnalysisTruthModuleBase_hxx_seen
#define TAnalysisTruthModuleBase_hxx_seen

#include "TAnalysisModuleBase.hxx"
namespace CP {
  class TAnalysisTruthModuleBase;
};

/// A base class for analysis output modules which contain event truth
/// information.  This class adds very little extra functionality, but
/// GetTreeType is defined to return kTruth.
class CP::TAnalysisTruthModuleBase : public TAnalysisModuleBase {
  public:

    // Define the tree type.  This sets the tree directory name.
    virtual EType GetTreeType() const { return kTruth;}

  protected:

    virtual ~TAnalysisTruthModuleBase();

  private:

    ClassDef(TAnalysisTruthModuleBase,1);

};
#endif
