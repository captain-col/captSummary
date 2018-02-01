#include "TTreeMaker.hxx"

#include <eventLoop.hxx>

int main(int argc, char **argv) {
    CP::TTreeMakerLoop userCode;
    CP::eventLoop(argc,argv,userCode);
}
