
#include <string>
#include <iostream>

#include "AnHiMaHGCAL/TriggerJetRate/interface/AnalysisJetRate.h"

using namespace std;
using namespace AnHiMa;

int main(int argc, char** argv)
{
    if(argc!=2)
    {
        cout<<"Usage: jetrate.exe parFile\n";
        return 1;
    }
    string parFile(argv[1]);

    AnalysisJetRate* analysis = new AnalysisJetRate();
    
    try
    {
        bool status = analysis->initialize(parFile);
        if(!status)
        {
            delete analysis;
            return 1;
        }
        analysis->loop();
    }
    catch(string s)
    {
        cout<<"ERROR: "<<s<<"\n";
        delete analysis;
        return 1;
    }
    
    delete analysis;
    return 0;
}