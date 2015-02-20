/**
 *  @file  AnalysisLayerCalibration.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    20/02/2015
 *
 *  @internal
 *     Created :  20/02/2015
 * Last update :  20/02/2015 16:29:59
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */






#include <iostream>
#include <fstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/TriggerLayerCalibration/interface/AnalysisLayerCalibration.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"




using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisLayerCalibration::AnalysisLayerCalibration():IAnalysis(),
    m_genParticle(0),
    m_matchedElectronCone(0)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisLayerCalibration::~AnalysisLayerCalibration()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisLayerCalibration::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    m_egammaAlgo.initialize(event(), m_reader.params());
    return true;
}


/*****************************************************************/
void AnalysisLayerCalibration::execute()
/*****************************************************************/
{
    cerr<<"Event "<<event().event()<<"\n";
    event().update();
    if(!event().passSelection()) return;


    for(int i=0;i<2;i++)
    {
        // get gen electron
        m_genParticle = &event().genparticles()[i];
        if(fabs(m_genParticle->Eta())<3.0 && fabs(m_genParticle->Eta())>1.5) 
        {
            m_electronCones.clear();
            m_egammaAlgo.idealClustering(event(), m_electronCones, m_genParticle->Eta(), m_genParticle->Phi()); // cone of 0.3 around gen electron
            m_matchedElectronCone = (m_electronCones.size()>0 ? &m_electronCones[0] : 0);
            //
            if(m_matchedElectronCone) fillHistos();
        }
    }
}

/*****************************************************************/
void AnalysisLayerCalibration::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;



    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events

}

