#include "TAnalysisLoop.hxx"
#include "TAnalysisModuleBase.hxx"
#include "TG4TrajectoriesModule.hxx"
#include "TG4VerticesModule.hxx"

#include <TRootInput.hxx>
#include <TRootOutput.hxx>
#include <TManager.hxx>

#include <TString.h>
#include <TSystem.h>
#include <TFile.h>

#include <cmath>
#include <iostream>
#include <exception>
#include <list>

/////////////////////////////////////////////////////////////////
// Initialize any class specific variables, but most of the work can be
// done in Initialize.  Don't create histograms here!
CP::TAnalysisLoop::TAnalysisLoop()
    : fEventCount(0),
      fSaveOriginalFullEvent(false), fSaveOutputEventTree(false),
      fDisableAll(false), fEnableAll(false) {

    // Add the "must have" modules.
    // fAnalysisModules.push_back(new CP::TBasicHeaderModule);
    fAnalysisModules.push_back(new CP::TG4TrajectoriesModule);
    fAnalysisModules.push_back(new CP::TG4VerticesModule);

    // Save pointer to these modules; this has to be done, because we
    // no longer have access to the input file pointer, except in BeginFile.
    // fGRooTracker = new CP::TGRooTrackerVtxModule();
    // fAnalysisModules.push_back(fGRooTracker);
    // fNRooTracker = new CP::TNRooTrackerVtxModule();
    // fAnalysisModules.push_back(fNRooTracker);

    // Set defaults for enabling or disabling modules.
    SetModuleDefaults();
}

CP::TAnalysisLoop::~TAnalysisLoop() { }

/////////////////////////////////////////////////////////////////
// Print a usage message.  This is generally called when there is a command
// line input error.
void CP::TAnalysisLoop::Usage(void) { }

/////////////////////////////////////////////////////////////////
// Set an option and return true if it is valid.  This is called by the
// event loop command line argument processing code for each "-O
// [name]=[value]" option found on the command line.  If the command line
// has "-O [name]" without a value, then the value string will be equal to
// "".  This must return false if the option was not correctly processed.
// Any -O name=value options where name is not explicitly given here will
// be checked against all the modules, and the value sent to
// module::Configure() in the case of a match.
bool CP::TAnalysisLoop::SetOption(std::string option,std::string value) {
    if (option == "save") {
        fSaveOriginalFullEvent = true;
        fSaveOutputEventTree = true;
        return true;
    }
    return false;
}

/////////////////////////////////////////////////////////////////
// Called for each event inside the event loop, and returns true if the
// event should be saved to the output file.  If the remainder of the
// current file should be skipped, this should through the
// ENextEventLoopFile exception.
bool CP::TAnalysisLoop::operator () (CP::TEvent& event) {
    // Need to make sure we can get geoemtry before going through analysis
    // modules.
    CP::TManager::Get().Geometry();
    
    bool saveThisFullEvent = false;
    
    ModuleList::iterator rooModIter;
    for (ModuleList::iterator modIter = fAnalysisModules.begin();
        modIter != fAnalysisModules.end(); modIter++) {
        rooModIter = modIter;
        CP::TAnalysisModuleBase *mod = *modIter;
        if (mod->IsEnabled()) {
            try {
                if (fEventCount == 1) {
                    // As a special case, call ProcessFirstEvent() for the
                    // first event.  This is a good time to save persistent
                    // quantities in modules' data members, which will be
                    // retrievable from the oaAnalysis output file.  Not
                    // intended for filling in the tree with the first event,
                    // as Process() will also be called.
                    mod->ProcessFirstEvent(event);
                }     
                // Call the module's own Process, passing it the event. If
                // returns true and module is selected for preselection then
                // save original nd280 event as well as the output analysis
                // trees.
                saveThisFullEvent = mod->Process(event);
                if (!mod->IsUsedForPreselection()) saveThisFullEvent = false;
            }
            catch (CP::EDataFile &ex) {
                mod->SetDisabled();
                continue;
            }
            catch(std::exception &ex) {
                CaptError(ex.what());
                CaptError("#### Disabling module "<<mod->GetName());
                mod->SetDisabled();
                continue;
            }
            catch (...) {
                CaptError("#### Disabling module "<<mod->GetName());
                mod->SetDisabled();
                continue;
            }
        }
    }
    
    // Check if the original event should also be saved.
    if (saveThisFullEvent || fSaveOriginalFullEvent) {
        fSaveOutputEventTree = true;
        return true;
    }

    return false;
}

/////////////////////////////////////////////////////////////////
// Called after the arguments are processes by before reading the first
// event.  The output file is open so any histograms will be added to the
// output file.
void CP::TAnalysisLoop::Initialize(void) {

    // Check that the output file is writable first
    if (!gFile) {
        CaptError("#### Output file must be provided.");
        exit(1);
    }

    // Override default enable/disable based on usr cmd line options.
    for (ModuleList::iterator m = fAnalysisModules.begin();
        m != fAnalysisModules.end(); ++m) {
        std::string className((*m)->ClassName());
        // remove the namespace.
        className.replace(className.find("CP::"), 4, ""); 
        // Ignore enable/disable options if selected for preselection.
        if (!(*m)->IsUsedForPreselection()){
            if (fDisableAll) (*m)->SetDisabled();  // Overrides default values.
            if (fEnableAll) (*m)->SetEnabled();    // Overrides disable all.
            if (fUserEnabled.find(className) != fUserEnabled.end()){
                // Override enable/disable all. 
                if (fUserEnabled.find(className)->second == true) {
                    (*m)->SetEnabled(); 
                }
                else {
                    (*m)->SetDisabled();
                }
            }
        }
    }


    // Move modules with fIsUsedForPreselection set to the end of the 
    // fAnalysisModules list
    for (ModuleList::reverse_iterator m = fAnalysisModules.rbegin();
        m != fAnalysisModules.rend(); ++m) {
        if ((*m)->IsUsedForPreselection()) {
            fAnalysisModules.push_back(*m);
            fAnalysisModules.erase(--m.base());
        }
    }
    
    // Print list of modules
    for (ModuleList::iterator m = fAnalysisModules.begin();
        m != fAnalysisModules.end(); ++m) {
        (*m)->Print();
    }
    
    for (ModuleList::iterator modIter  = fAnalysisModules.begin();
        modIter != fAnalysisModules.end(); ++modIter) {
        CP::TAnalysisModuleBase *mod = *modIter;
        // Ask the modules to fill the appropriate directory
        std::string dirname = mod->GetDirectoryName();
        gFile->cd("/");
        gDirectory->cd("/");
        if (gFile && !gFile->GetDirectory(dirname.c_str())) {
            CaptVerbose("Making "<<dirname);
            gDirectory->mkdir(dirname.c_str());
            gDirectory->cd(dirname.c_str());
            gDirectory->ls();
            gDirectory->mkdir("Modules");
            gDirectory->cd("/");
        }
        if (mod->IsEnabled()) {
            gDirectory->cd(dirname.c_str());
            try {
                // Call the module's own Initialize, passing it its own tree
                mod->Initialize(new TTree(mod->GetName(), mod->GetTitle()));
            }
            catch (std::exception &ex) {
                CaptError(ex.what());
                CaptError("### Disabling module (unable to initialise) "
                          << mod->GetName());
                mod->SetDisabled();
                continue;
            }
            catch (...) {
                CaptError("### Disabling module (unable to initialise) " 
                          << mod->GetName());
                mod->SetDisabled();
                continue;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////
// Called before the first event of a file is read, but you should prefer
// Initialize() for general initialization.  This method will be called
// once for each input file.
void CP::TAnalysisLoop::BeginFile(CP::TVInputFile *vinput) {
    CP::TRootInput *input = dynamic_cast<CP::TRootInput *>(vinput);
    if (!input) return;
    TFile* file = input->GetFilePointer();
    if (!file) return;
    for (ModuleList::iterator m  = fAnalysisModules.begin();
        m != fAnalysisModules.end(); ++m) {
        (*m)->SetBeginFile(file);
    }
}

/////////////////////////////////////////////////////////////////
// Called after the last event of a file is read, but you should
// prefer Finalize() for general finalization.  This method will
// be called once for each input file.
void CP::TAnalysisLoop::EndFile(CP::TVInputFile *) {}

/////////////////////////////////////////////////////////////////
// Called after reading the last event.  The output file is still open, so
// you can add extra information.  Because of an idiosyncrasy in the way
// root handles histograms, objects created in Initialize() will already be
// stored in the output file.
void CP::TAnalysisLoop::Finalize(CP::TRootOutput* const output) {
    // Write geometry if needed

    // Write module itself inside [Header|Truth|Recon]Dir/Modules/
    for (std::list<CP::TAnalysisModuleBase*>::iterator modIter
            = fAnalysisModules.begin(); modIter != fAnalysisModules.end();
        ++modIter) {
        CP::TAnalysisModuleBase *mod = *modIter;

        // Write the modules themselves in the Module directory.
        std::string moduledir("/");
        moduledir += mod->GetDirectoryName();
        moduledir += "/Modules";
        CaptLog("ModulesDir: "<<moduledir);

        if (output->cd(moduledir.c_str())==false) {
            output->mkdir(moduledir.c_str());
        }

        // I'd rather not change the name, but the browser gets
        // confused if the module has the same name as the tree....
        std::string moduleName(mod->GetName());
        moduleName += "Mod";
        mod->Write(moduleName.c_str());
        CaptLog("Saved module as "<<moduledir<<"/"<<moduleName);

        if (mod->IsEnabled() == false) {
            // Delete the output tree if the module is disabled
            std::string dirname = mod->GetDirectoryName();
            output->cd("/");
            gDirectory->cd(dirname.c_str());
            gDirectory->Delete(TString(mod->GetName()) + ";*");
        }
        delete (*modIter);
    }
    CaptLog("Total number of events processed in Analysis Format: "
            << fEventCount);
    output->cd("//");

    if (!fSaveOutputEventTree){
        output->Delete("captainEventTree;*");
    }
}

/////////////////////////////////////////////////////////////////
// Called just before setting user options so that modules start off with
// default enabled/disabled option.
void CP::TAnalysisLoop::SetModuleDefaults() {
    for (ModuleList::iterator m=fAnalysisModules.begin();
        m != fAnalysisModules.end(); ++m) {
        (*m)->SetEnabled((*m)->IsEnabledByDefault());
    }
}

