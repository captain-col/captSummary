#include "ECaptainSummary.hxx"

#include <TEventLoopFunction.hxx>

#include <cmath>
#include <iostream>
#include <exception>
#include <list>

namespace CP {
    class TAnalysisLoop;
    class TAnalysisModuleBase;
    EXCEPTION(EProductionException, ECaptainSummary);
}

/// Generate an analysis tree to summarize the event reconstruction.  This
/// should be replaced by the general oaAnalysis tree, but that's not ready.
class CP::TAnalysisLoop: public CP::TEventLoopFunction {
public:
    /// Initialize any class specific variables, but most of the work can be
    /// done in Initialize.  Don't create histograms here!
    TAnalysisLoop();
    virtual ~TAnalysisLoop();

    /// Print a usage message.  This is generally called when there is a
    /// command line input error.
    void Usage(void);

    /// Set an option and return true if it is valid.  This is called by the
    /// event loop command line argument processing code for each "-O
    /// [name]=[value]" option found on the command line.  If the command line
    /// has "-O [name]" without a value, then the value string will be equal
    /// to "".  This must return false if the option was not correctly
    /// processed.
    virtual bool SetOption(std::string, std::string);

    /// Called for each event inside the event loop, and returns true if the
    /// event should be saved to the output file.  If the remainder of the
    /// current file should be skipped, this should through the
    /// ENextEventLoopFile exception.
    bool operator () (CP::TEvent& event);

    /// Called after the arguments are processes by before reading the first
    /// event.  The output file is open so any histograms will be added to the
    /// output file.
    virtual void Initialize(void);

    /// Called before the first event of a file is read, but you should prefer
    /// Initialize() for general initialization.  This method will be called
    /// once for each input file.
    virtual void BeginFile(CP::TVInputFile *input);

    /// Called after the last event of a file is read, but you should prefer
    /// Finalize() for general finalization.  This method will be called once
    /// for each input file.
    virtual void EndFile(CP::TVInputFile *);

    /// Called after reading the last event.  The output file is still open,
    /// so you can add extra information.  Because of an idiosyncrasy in the
    /// way root handles histograms, objects created in Initialize() will
    /// already be stored in the output file.
    virtual void Finalize(CP::TRootOutput* const output);

private:
    /// Mehod to set the the default enable/disable values for all modules
    /// in fAnalysisModules list.
    void SetModuleDefaults();

    /// Keep track of the events that have been read.
    int fEventCount;

    /// Flag if the full event (not the analysis format) is also to be saved
    /// This overrides any -preselection=<modulename> settings
    bool fSaveOriginalFullEvent;

    /// Set automatically if the Event tree needs to be saved in output
    /// file.
    bool fSaveOutputEventTree;

    /// Internal list of all modules.
    typedef std::list<CP::TAnalysisModuleBase*> ModuleList;
    ModuleList fAnalysisModules;

    /// Store whether a module has been set to disabled/enables by user
    /// as well as whether the disable/enable all options have been set.
    /// fUserEnabled overrides fDisableAll and fDisableAll overrides 
    /// fEnableAll. This overrides fDisable(Enable)All.
    std::map<std::string, bool> fUserEnabled; 

    /// Store whether a module has been set to disabled/enables by user
    /// as well as whether the disable/enable all options have been set.
    /// fUserEnabled overrides fDisableAll and fDisableAll overrides 
    /// fEnableAll.  This overrides fEnableAll.
    bool fDisableAll;

    /// Store whether a module has been set to disabled/enables by user
    /// as well as whether the disable/enable all options have been set.
    /// fUserEnabled overrides fDisableAll and fDisableAll overrides 
    /// fEnableAll.  
    bool fEnableAll;

    /// Booleans to control production mode
    bool fProduction;

    /// Says to exit if there is a problem
    bool fExit;

    /// Boolean to control validation mode
    bool fValidation;
 
};
