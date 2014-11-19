/**
 *  @file  AnalysisClustering.cpp
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
#include <fstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/StudyClustering/interface/AnalysisClustering.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisClustering::AnalysisClustering():IAnalysis()
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisClustering::~AnalysisClustering()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisClustering::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    // Read parameters
    string pileupParamsFile = m_reader.params().GetValue("PileupParams", "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/src/AnHiMaHGCAL/ElectronClusterThreshold/data/thresholdParameters.txt");

    m_egammaAlgo.initialize(event(), pileupParamsFile);
    return true;
}


/*****************************************************************/
void AnalysisClustering::execute()
/*****************************************************************/
{
    event().update();
    if(!event().passSelection()) return;

    m_seeds.clear();
    m_egammaAlgo.seeding(event(), m_seeds);


    fillHistos();

}

/*****************************************************************/
void AnalysisClustering::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;


    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    //
    m_histos.FillHisto(10+hoffset, m_seeds.size(), weight, sysNum);

    for(const auto& tower : m_seeds)
    {
        m_histos.FillHisto(11+hoffset, fabs(tower.eta()), weight, sysNum);
        m_histos.FillHisto(12+hoffset, tower.phi(), weight, sysNum);
        m_histos.FillHisto(13+hoffset, fabs(tower.eta()), tower.phi(), weight, sysNum);
        m_histos.FillHisto(14+hoffset, tower.calibratedEt(), weight, sysNum);

        m_histos.Fill1BinHisto(100+hoffset, event().event(), fabs(tower.eta()), tower.phi(), weight, sysNum);
    }
}

