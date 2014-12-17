
#include <string>
#include <iostream>

#include "AnHiMaHGCAL/TestNavigation/interface/NavigationTester.h"

using namespace std;
using namespace AnHiMa;

int main(int argc, char** argv)
{
    if(argc!=1)
    {
        cout<<"Usage: testNavigation.exe\n";
        return 1;
    }
    NavigationTester test;
    test.initialize();
    test.leftRightTest();

    
    return 0;
}
