
#include <string>
#include <iostream>

#include "AnHiMaHGCAL/ProjectionMapping/interface/ProjectionMapProducer.h"

using namespace std;
using namespace AnHiMa;

int main(int argc, char** argv)
{
    if(argc!=1)
    {
        cout<<"Usage: projectionMapping.exe\n";
        return 1;
    }
    ProjectionMapProducer mapProducer;
    mapProducer.initialize();
    mapProducer.produce();

    
    return 0;
}
