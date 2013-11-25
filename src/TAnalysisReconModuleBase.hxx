#ifndef TAnalysisReconModuleBase_hxx_seen
#define TAnalysisReconModuleBase_hxx_seen

#include "TAnalysisModuleBase.hxx"

namespace CP {
  class TAnalysisReconModuleBase;
};

/// A base class for analysis output modules which contain reconstructed event
/// information.  This class adds very little extra functionality, but
/// GetTreeType is defined to return kRecon.
class CP::TAnalysisReconModuleBase : public TAnalysisModuleBase
{
  public:
    // Define the tree type.  This sets the tree directory name.
    virtual EType GetTreeType() const { return kRecon;}

  protected:
    virtual ~TAnalysisReconModuleBase();

  private:

    ClassDef(TAnalysisReconModuleBase,1);

};
#endif
