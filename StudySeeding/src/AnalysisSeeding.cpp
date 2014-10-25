/**
 *  @file  AnalysisSeeding.cpp
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


#include "AnHiMaHGCAL/StudySeeding/interface/AnalysisSeeding.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"



using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisSeeding::AnalysisSeeding():IAnalysis()
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisSeeding::~AnalysisSeeding()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisSeeding::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);


    m_hgcalNavigator.initialize();


    return true;
}


/*****************************************************************/
void AnalysisSeeding::execute()
/*****************************************************************/
{
    event().update();
    if(!event().passSelection()) return;

    fillHistos();
}

/*****************************************************************/
void AnalysisSeeding::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;

    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events

}


