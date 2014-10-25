/**
 *  @file  Analysis.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    19/10/2014
 *
 *  @internal
 *     Created :  19/10/2014
 * Last update :  19/10/2014 13:44:32
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#include <iostream>


#include "AnHiMaHGCAL/Test/interface/Analysis.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"



using namespace AnHiMa;
using namespace std;



/*****************************************************************/
Analysis::Analysis():IAnalysis()
/*****************************************************************/
{

}



/*****************************************************************/
Analysis::~Analysis()
/*****************************************************************/
{
}



/*****************************************************************/
bool Analysis::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);


    m_hgcalNavigator.initialize();

    // Copy-paste from https://github.com/PFCal-dev/cmssw/blob/CMSSW_6_2_X_SLHC/Geometry/CaloTopology/test/HGCalTopologyTester.cc 
    for (int izz=0; izz<=1; izz++) 
    {
        int iz = (2*izz-1);
        for (int subsec=0; subsec<=1; ++subsec) 
        {
            for (int sec=1; sec<=36; ++sec) 
            {
                for (int lay=1; lay<=10; ++lay) 
                {
                    for (int cell=0; cell<8000; ++cell) 
                    {
                        const HGCEEDetId id(HGCEE,iz,lay,sec,subsec,cell);
                        if (m_hgcalNavigator.valid(id)) 
                        {
                            std::cout << "Neighbours for Tower "  << id << std::endl;
                            std::vector<DetId> idE = m_hgcalNavigator.east(id);
                            std::vector<DetId> idW = m_hgcalNavigator.west(id);
                            std::vector<DetId> idN = m_hgcalNavigator.north(id);
                            std::vector<DetId> idS = m_hgcalNavigator.south(id);
                            std::vector<HGCEEDetId> idU = m_hgcalNavigator.up(id);
                            std::vector<HGCEEDetId> idD = m_hgcalNavigator.down(id);
                            std::cout << "          " << idE.size() << " sets along East:";
                            for (unsigned int i=0; i<idE.size(); ++i) std::cout << " " << (HGCEEDetId)(idE[i]());
                            std::cout << std::endl;
                            std::cout << "          " << idW.size() << " sets along West:";
                            for (unsigned int i=0; i<idW.size(); ++i) std::cout << " " << (HGCEEDetId)(idW[i]());
                            std::cout << std::endl;
                            std::cout << "          " << idN.size() << " sets along North:";
                            for (unsigned int i=0; i<idN.size(); ++i) std::cout << " " << (HGCEEDetId)(idN[i]());
                            std::cout << std::endl;
                            std::cout << "          " << idS.size() << " sets along South:";
                            for (unsigned int i=0; i<idS.size(); ++i) std::cout << " " << (HGCEEDetId)(idS[i]());
                            std::cout << std::endl;
                            std::cout << "          " << idU.size() << " sets along Up:";
                            for (unsigned int i=0; i<idU.size(); ++i) std::cout << " " << (HGCEEDetId)(idU[i]());
                            std::cout << std::endl;
                            std::cout << "          " << idD.size() << " sets along Down:";
                            for (unsigned int i=0; i<idD.size(); ++i) std::cout << " " << (HGCEEDetId)(idD[i]());
                            std::cout << std::endl;
                        }
                        cell += 100;
                    }
                }
            }
        }
    }


    return true;
}


/*****************************************************************/
void Analysis::execute()
/*****************************************************************/
{
    event().update();
    if(!event().passSelection()) return;

    fillHistos();
}

/*****************************************************************/
void Analysis::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;

    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events

}


