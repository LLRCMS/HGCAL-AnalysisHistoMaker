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
    event().update();
    if(!event().passSelection()) return;

    fillGenElectrons();

    for(const auto& ele : m_genElectrons)
    {
        // get gen electron
        m_genParticle = ele;
        if(fabs(m_genParticle->Eta())<3.0 && fabs(m_genParticle->Eta())>1.5) 
        {
            m_electronCones.clear();
            m_egammaAlgo.coneClustering(event(), m_electronCones, m_genParticle->Eta(), m_genParticle->Phi()); // cone of 0.3 around gen electron
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


    double Ehalf = 0.;
    for(unsigned l=2;l<=30;l+=2)
    {
        double Ei   = m_matchedElectronCone->layerCalibratedEnergy(l);
        double Eim1 = m_matchedElectronCone->layerCalibratedEnergy(l-1);
        if(Ei>0.) m_histos.Fill1BinHisto(10+hoffset, l, (Ei+Eim1)/Ei, weight, sysNum);
        Ehalf += Ei;
    }
    for(unsigned l=1;l<=30;l++)
    {
        double E = m_matchedElectronCone->layerCalibratedEnergy(l);
        m_histos.FillHisto(100+l+hoffset, E, weight, sysNum);
    }
    m_histos.FillHisto(99+hoffset, m_matchedElectronCone->eta(), weight, sysNum);
    m_histos.FillHisto(100+hoffset, Ehalf, weight, sysNum);
    m_histos.FillNtuple(500+hoffset, event().run(), event().event(), weight, sysNum);

}



/*****************************************************************/
void AnalysisLayerCalibration::fillGenElectrons()
/*****************************************************************/
{
    m_genElectrons.clear();
    for(const auto& gen : event().genparticles())
    {
        if(gen.status()!=3) continue;
        if(abs(gen.id())!=11) continue;
        if(gen.Pt()<10.) continue;
        if(fabs(gen.Eta())<1.5 || fabs(gen.Eta())>3.) continue;
        m_genElectrons.push_back(&gen);
    }
}

