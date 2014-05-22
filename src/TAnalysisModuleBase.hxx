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
    /// Initialize Module.  This is where the internal fields of the modules
    /// should be set-up, but don't create the branches.  Modules might not
    /// need to do any thing in this method, but it must be defined.
    virtual void InitializeModule() = 0;
    
    /// Initialize Branches (don't do anything else in this function).
    virtual void InitializeBranches() = 0;
    
    /// Fill all the stuff that goes in the output tree. Return true if
    /// everything went well. Otherwise, the module may be disabled!
    virtual bool FillTree(CP::TEvent&) = 0;

public:
    /////////////////////////////////////////////////////////////////////
    // These are methods might be overriden in the derived modules.  The
    // mostly control the behavior of the module.
    /////////////////////////////////////////////////////////////////////

    /// Is the module is enabled by default. Default is to enable module. To
    /// set to disable override this method in the derived module.
    virtual bool IsEnabledByDefault() const {return kTRUE;}
    
    /// Is called after the first event is loaded in.  This is a good time to
    /// save persistent quantities in the module's data members, which will be
    /// retrievable from the output file.  Not intended for filling in the
    /// tree with the first event, as Process() will also be called.  
    virtual bool ProcessFirstEvent(CP::TEvent&);

    /// Return true with the full event informat should be saved.  A lot of
    /// the information in the full event isn't needed for an analysis, so
    /// analysis modules build data summary trees (dst) based on the full
    /// event information.  However sometimes it is useful to save a
    /// preselected sample of the full events to allow further reconstruction.
    /// If this method returns true, then the full event information will also
    /// be saved to the output file.  Note: This causes the output files to be
    /// much larger.
    virtual bool SaveFullEvent(CP::TEvent& event);

    /// A function that allows the module to be configured from an external
    /// class without any dependencies.  Should be overridden with a function
    /// that responds to the string option, and returns kTRUE if configuration
    /// succeeded.  Used in TAnalysisLoop.cxx for options of the form: -O
    /// TTruthTrajectoriesModule=SaveAll
    virtual bool Configure(const std::string &option);

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
    std::string GetDirectoryName() const;

    /// Sets whether the module is enabled. This only refer to modules which
    /// have been included for consideration by being instantiated in
    /// TAnalysisLoop.cxx or similar.
    virtual void SetEnabled(bool yesorno = kTRUE) {fIsEnabled = yesorno;}

    /// Disables the module. Is called when an exception is thrown inside the
    /// module.
    void SetDisabled() {SetEnabled(kFALSE);}
    
    /// Sets whether the analysis loop should call SaveFullEvent() function
    /// for each event, to decide if the full event info for that event
    /// should be saved in the output
    void SetUsedForPreselection(bool yesorno = kTRUE) { 
        fIsUsedForPreselection = yesorno; 
    }
    
    /// Whether the module is enable or not.
    bool IsEnabled() const {return fIsEnabled;}
    
    /// Whether the analysis loop should call SaveFullEvent() function for
    /// each event, to decide if the full event info for that event should
    /// be saved in the output
    bool IsUsedForPreselection() const {return fIsUsedForPreselection;}
    
    /// Prints a simple message describing the module.  Should be overridden
    /// if more detail is needed (it's probably not needed).
    virtual void PrintMessage();
    
    /// Construct the fields for the tree.  This must not be overridden in the
    /// derived modules.  Derived modules will provide InitializeModule() and
    /// InitializeBranches().
    void Initialize(TTree *tree);
    
    /// Gets the run and event Id, calls FillTree and SaveFullEvent with the
    /// event, and then calls fOutputTree->Fill.  This doesn't need to be
    /// changed by most modules.
    virtual bool Process(CP::TEvent& event);

    /// ROOT output parameters, usually no need to touch
    int GetBufferSize() const {return fBufferSize;}

    /// ROOT output parameters, usually no need to touch
    void SetBufferSize(int  buffersize) {fBufferSize =  buffersize;}

    /// ROOT output parameters, usually no need to touch
    int GetSplitLevel() const {return fSplitLevel;}

    /// ROOT output parameters, usually no need to touch
    void SetSplitLevel(int splitlevel) {fSplitLevel = splitlevel;}

    /// The output tree
    TTree * const GetOutputTree() const {return fOutputTree;}

    /// Set the location of a directory where the module can look for extra
    /// files.
    void SetInputDirectory(const std::string& dir) {fInputDirectory = dir;}
    
    /// Get the location of a directory where the module can look for extra
    /// files.
    std::string GetInputDirectory() const {return fInputDirectory;}

    /// The derived module can override this if it needs access to the file
    /// pointer.  Examples of modules that will need the file pointer are the
    /// primary vertex trees where the vertex needs to be found in the tree.
    virtual void SetBeginFile(TFile* file) {}

private: 
    /// The tree that will be filled.
    TTree* fOutputTree;

    /// This is true if the module should be used in this run.
    bool fIsEnabled;

    /// This is true if the module might cause full events to be saved to the
    /// output.
    bool fIsUsedForPreselection;

    /// Buffer Size for TBranch. Has a default value that
    /// can be changed per module.
    int fBufferSize;  // Buffer Size for TBranch.

    /// Split Level for TBranch.
    int fSplitLevel; 

    ///////////////////////////////////////////////////////////////////
    /// Default Tree Entries
    int fRunID;
    int fSubrunID;
    int fEventID;
    int fPreselected;

    /// An input directory where analysis modules can search for extra files.
    std::string fInputDirectory;    

    ClassDef(TAnalysisModuleBase,1);

};
#endif
