#ifndef TAnalysisModuleBase_hxx_seen
#define TAnalysisModuleBase_hxx_seen
#include "ECaptainSummary.hxx"

#include <TEvent.hxx>

#include <TNamed.h>
#include <TTree.h>
#include <TGeoManager.h>

#include <string>

namespace CP {
    class TAnalysisModuleBase;
    EXCEPTION(EAnalysisModuleBase,ECaptainSummary);
    EXCEPTION(EUndefinedTreeType,EAnalysisModuleBase);
    EXCEPTION(EAnalysisFailure,EAnalysisModuleBase);
};

/// A base class for classes which specify how to set up an Analysis format
/// output tree, and fill it. These classes will be read in and run by the
/// oaAnalysis package software, namely by TAnalysisLoop.  The classes must be
/// accompanied by tests to check that they work, and runtime tests to see if
/// they are appropriate for the input root files being used.
///
/// A note on naming of data members in analysis modules classes and their
/// sub-classes:
/// 
///  * The names of members belonging to the modules themselves all start with
///    an 'f' to indicate that they are fields in the classes. These names are
///    used in coding up the modules, so it is better that they follow the
///    convention.
///   
///  * The members of storage sub-classes for the modules, which get saved in
///    TTrees within TClonesArrays, DO NOT get an 'f' at the beginning,
///    because they represent names that get used directly in output trees.
///   
///  * The modules themselves get saved in the output files, including their
///    data members (except those with a comment with an exclamation mark
///    "//!"  after them, which tells ROOT not to save them). These can be
///    accessed in the analysis macros, but since they are explicitly members
///    of the module classes, it is consistent to have them start with 'f'
///    too.
class CP::TAnalysisModuleBase : public TNamed {
public:
    TAnalysisModuleBase();
    virtual ~TAnalysisModuleBase();
    
protected:
    /// Define an enumerated list of the possible tree types.
    enum EType {
        /// A tree containing header information.
        kHeader = 0, 
        /// A tree containing truth information
        kTruth, 
        /// A tree containing reconstruction information.
        kRecon, 
        /// The number of tree types.
        kNTypes};
    
    /////////////////////////////////////////////////////////////////////
    // These are methods that must be overridden in the derived modules.
    /////////////////////////////////////////////////////////////////////
protected:
    /// Initialize Module.  This is where the internal fields of the modules
    /// should be set-up, but don't create the branches.  Modules might not
    /// need to do any thing in this method, but it must be defined.
    virtual void InitializeModule() = 0;
    
    /// Initialize Branches (don't do anything else in this function).
    virtual void InitializeBranches() = 0;
    
    /// Fill all the stuff that goes in the output tree. Return true if
    /// everything went well. Otherwise, the module may be disabled!
    virtual bool FillTree(CP::TEvent&) = 0;

    /////////////////////////////////////////////////////////////////////
    // These are methods might be overriden in the derived modules.  The
    // mostly control the behavior of the module.
    /////////////////////////////////////////////////////////////////////
public:
    /// Is the module is enabled by default. Default is to enable module. To
    /// set to disable override this method in the derived module.
    virtual Bool_t IsEnabledByDefault() const {return kTRUE;}
    
    /// Is called after the first event is loaded in.  This is a good time to
    /// save persistent quantities in the module's data members, which will be
    /// retrievable from the output file.  Not intended for filling in the
    /// tree with the first event, as Process() will also be called.  
    virtual Bool_t ProcessFirstEvent(CP::TEvent&);

    /// Whether the module thinks it is worth saving the entire event
    /// tree entry for this event.  The summary can be used for event
    /// pre-selection in this way.
    virtual bool IsFullEventWorthSaving(CP::TEvent& event);

    /// A function that allows the module to be configured from an external
    /// class without any dependencies.  Should be overridden with a function
    /// that responds to the string option, and returns kTRUE if configuration
    /// succeeded.  Used in TAnalysisLoop.cxx for options of the form: -O
    /// TTruthTrajectoriesModule=SaveAll
    virtual Bool_t Configure(std::string &option);

public:
    /////////////////////////////////////////////////////////////////////
    // The remaining methods won't need to be overridden in the derived
    // classes (There are some minor exceptions where extremely special case
    // modules might).
    /////////////////////////////////////////////////////////////////////
    
    /// Returns the type of tree, header, truth, or recon.  This is overridden
    /// in the derived base classes such as TAnalysisReconModuleBase, so users
    /// do not need to override it explicitly
    virtual EType GetTreeType() const = 0;
    
    /// Returns the name of the directory which the output of a particular
    /// module will be saved in.
    std::string const GetDirectoryName() const;
    

    /// Sets whether the module is enabled. This only refer to modules which
    /// have been included for consideration by being instantiated in
    /// TAnalysisLoop.cxx or similar.
    virtual void SetEnabled(Bool_t yesorno = kTRUE) {
        fIsEnabled = yesorno; 
    }

    /// Disables the module. Is called when an exception is thrown inside the
    /// module.
    void SetDisabled() {
        SetEnabled(kFALSE);
    }
    
    /// Sets whether the module should call IsFullEventWorthSaving() function
    /// for each event, to decide if the full event info for that event
    /// should be saved in the output
    void SetUsedForPreselection(Bool_t yesorno = kTRUE) { 
        fIsUsedForPreselection = yesorno; 
    }
    
    /// Whether the module is enable or not.
    Bool_t IsEnabled() const { return fIsEnabled; }
    
    /// Whether the module should call IsFullEventWorthSaving() function for
    /// each event, to decide if the full event info for that event should
    /// be saved in the output
    Bool_t IsUsedForPreselection() const { 
        return fIsUsedForPreselection;
    }
    
    /// Prints a simple message describing the module.  Should be overridden
    /// if more detail is needed (it's probably not needed).
    virtual void PrintMessage();
    
    /// Construct the fields for the tree.  This must not be overridden in the
    /// derived modules.  Derived modules will provide InitializeModule() and
    /// InitializeBranches().
    void Initialize(TTree *tree);
    
    /// Gets the run and event IDs and calls FillTree with the event, and then
    /// actually calls fOutputTree->Fill.  This doesn't need to be changed by
    /// most modules.
    virtual bool Process(CP::TEvent& event);

    /// ROOT output parameters, usually no need to touch
    Int_t GetBufferSize() {return fBufferSize;}

    /// ROOT output parameters, usually no need to touch
    void SetBufferSize(Int_t  buffersize) {fBufferSize =  buffersize;}

    /// ROOT output parameters, usually no need to touch
    Int_t GetSplitLevel() {return fSplitLevel;}

    /// ROOT output parameters, usually no need to touch
    void SetSplitLevel(Int_t splitlevel) {fSplitLevel = splitlevel;}

    /// The output tree
    TTree const * GetOutputTree() const {return fOutputTree;}

    /// Set the location of a directory where the module can look for extra
    /// files.
    void SetInputDirectory(const std::string& dir) {fInputDirectory = dir;}
    
    /// Get the location of a directory where the module can look for extra
    /// files.
    const std::string& GetInputDirectory() const {return fInputDirectory;}

    /// The derived module can override this if it needs access to the file
    /// pointer.  Examples of modules that will need the file pointer are the
    /// primary vertex trees where the vertex needs to be found in the tree.
    virtual void SetBeginFile(TFile* file) {}

private: 
    /// This is true if the module should be run.
    Bool_t fIsEnabled;

    /// This is true if the module might cause full events to be saved to the
    /// output.
    Bool_t fIsUsedForPreselection;

    /// Keep the root tree handy.
    TTree *fOutputTree;

    /// Buffer Size for TBranch. Has a default value that
    /// can be changed per module.
    Int_t fBufferSize;  // Buffer Size for TBranch.

    /// Split Level for TBranch.
    Int_t fSplitLevel; 

    ///////////////////////////////////////////////////////////////////
    /// Default Tree Entries
    Int_t fRunID;
    Int_t fSubrunID;
    Int_t fEventID;
    Int_t fPreselected;

    /// An input directory where analysis modules can search for extra files.
    std::string fInputDirectory;    

private:

    ClassDef(TAnalysisModuleBase,1);

};
#endif
