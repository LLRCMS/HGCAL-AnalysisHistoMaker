/**
 *  @file  AnalysisPileup.cpp
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

#include "TVector2.h"

#include "AnHiMaHGCAL/StudyPileup/interface/AnalysisPileup.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"






using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisPileup::AnalysisPileup():IAnalysis()
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisPileup::~AnalysisPileup()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisPileup::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);


    return true;
}


/*****************************************************************/
void AnalysisPileup::execute()
/*****************************************************************/
{
    event().update();
    if(!event().passSelection()) return;

    //double mip = 0.000055;



}

/*****************************************************************/
void AnalysisPileup::fillHistos()
/*****************************************************************/
{
    //double mip = 0.000055;

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;


    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events

}


