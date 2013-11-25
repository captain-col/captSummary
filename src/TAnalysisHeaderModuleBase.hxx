#ifndef TAnalysisHeaderModuleBase_hxx_seen
#define TAnalysisHeaderModuleBase_hxx_seen

#include "TAnalysisModuleBase.hxx"
namespace CP {
  class TAnalysisHeaderModuleBase;
};

/// A base class for analysis output modules which contain event header
/// information.  This class adds very little extra functionality, but
/// GetTreeType is defined to return kHeader.
class CP::TAnalysisHeaderModuleBase : public TAnalysisModuleBase {
  public:
    // Define the tree type.  This sets the tree directory name.
    EType GetTreeType() const { return kHeader;}

  protected:
    virtual ~TAnalysisHeaderModuleBase() {}

  private:
    ClassDef(TAnalysisHeaderModuleBase,1);
};
#endif
