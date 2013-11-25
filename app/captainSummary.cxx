#include "TAnalysisLoop.hxx"

#include <eventLoop.hxx>

int main(int argc, char **argv) {
    CP::TAnalysisLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
