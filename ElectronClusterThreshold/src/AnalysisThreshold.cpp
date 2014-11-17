/**
 *  @file  AnalysisThreshold.cpp
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

#include "AnHiMaHGCAL/ElectronClusterThreshold/interface/AnalysisThreshold.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"






using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisThreshold::AnalysisThreshold():IAnalysis(),
    m_chosenHit(0),
    m_genParticle(0),
    m_pileupParamsFile(""),
    m_electronHitEnergy(0.)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisThreshold::~AnalysisThreshold()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisThreshold::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    // Read parameters
    m_pileupParamsFile = m_reader.params().GetValue("PileupParams", "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/src/AnHiMaHGCAL/ElectronClusterThreshold/data/thresholdParameters.txt");
    fillPileupParams();


    // initialize vectors
    for(int s=0;s<12;s++)
    {
        for(int z=0;z<2;z++)
        {
            for(int u=0;u<2;u++)
            {
                int triggerRegion = s + 12*u + 24*z;
                for(int l=1;l<=30;l++)
                {
                    m_regionSimHits[make_pair(triggerRegion, l)] = vector<pair<HGCEEDetId,float> >(0);
                    m_regionHitsAboveTh[make_pair(triggerRegion, l)] = vector<const SimHit*>(0);
                }
            }
        }
    }
    // fill vectors
    for (int izz=0; izz<=1; izz++) 
    {
        int iz = (2*izz-1);
        for (int subsec=0; subsec<=1; ++subsec) 
        {
            int ssec = (2*subsec-1);
            for (int sec=1; sec<=18; ++sec) 
            {
                int sector30deg = (2*(sec-1) + (subsec==0 ? 1 : 0))/3;
                for (int lay=1; lay<=30; ++lay) 
                {
                    for (int cell=0; cell<4000; ++cell) 
                    {
                        const HGCEEDetId id(HGCEE,iz,lay,sec,ssec,cell);
                        if (event().hgcalNavigator().valid(id)) 
                        {
                            double eta = event().hgcalNavigator().geometry()->getPosition(id).eta();
                            //int up = (cell>=1200 ? 1 : 0);
                            int up = (fabs(eta)<=1.8 ? 1 : 0);
                            int triggerRegion = sector30deg + 12*up + 24*izz;
                            m_regionSimHits[make_pair(triggerRegion, lay)].push_back(make_pair(id,eta));
                        }
                    }
                }
            }
        }
    }

    return true;
}


/*****************************************************************/
void AnalysisThreshold::execute()
/*****************************************************************/
{
    event().update();
    if(!event().passSelection()) return;

    fillRegionHits();

    for(int i=0; i<2; i++)
    {
        findMaxHit(i);
        if(!m_chosenHit)
        {
            cout<<"WARNING: cannot find hit close to electron\n";
            continue;
        }
        findElectronHits();
        findElectronHitsBelowThreshold();

        fillHistos();
    }


}

/*****************************************************************/
void AnalysisThreshold::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;


    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    //
    double missEnergy = 0.;
    int missHits = m_electronHitsBelowThreshold.size();
    for(const auto hit : m_electronHitsBelowThreshold)
    {
        int layer = hit->layer();
        double deta = (hit->eta() - m_chosenHit->eta());
        double dphi = TVector2::Phi_mpi_pi(hit->phi() - m_chosenHit->phi());
        missEnergy += hit->energy();
        m_histos.Fill1BinHisto(200+hoffset, fabs(m_genParticle->Eta()), layer, weight, sysNum);
        m_histos.Fill1BinHisto(250+hoffset, fabs(m_genParticle->Eta()), deta, dphi, weight, sysNum);
    }
    m_histos.Fill1BinHisto(100+hoffset, fabs(m_genParticle->Eta()), missHits, weight, sysNum);
    m_histos.Fill1BinHisto(150+hoffset, fabs(m_genParticle->Eta()), missEnergy, weight, sysNum);
    //
    for(const auto& id : m_electronHits)
    {
        const SimHit* hit = event().simhit(id);
        double deta = (hit->eta() - m_chosenHit->eta());
        double dphi = TVector2::Phi_mpi_pi(hit->phi() - m_chosenHit->phi());
        m_histos.Fill1BinHisto(1000+hoffset, fabs(m_genParticle->Eta()),  deta, weight, sysNum);
        m_histos.Fill1BinHisto(1050+hoffset, fabs(m_genParticle->Eta()),  dphi, weight, sysNum);
        m_histos.Fill1BinHisto(1100+hoffset, fabs(m_genParticle->Eta()), deta, dphi, weight, sysNum);
    }

}


/*****************************************************************/
void AnalysisThreshold::findMaxHit(int i)
/*****************************************************************/
{
    m_chosenHit = NULL;

    m_genParticle = &event().genparticles()[i];
    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxE = 0.;
    for(auto itrHit=event().hardsimhits().cbegin(); itrHit!=event().hardsimhits().end(); itrHit++)
    {
        const SimHit& hit = *itrHit;
        // match only to hits in layers 10-20
        if(hit.layer()<10 || hit.layer()>20) continue;
        //if(hit.layer()!=15 || hit.energy()<5.*mip) continue;

        double deta = (hit.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(hit.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.1 && hit.energy()>maxE)
        {
            m_chosenHit = &hit;
            maxE = hit.energy();
        }
    }

}

/*****************************************************************/
void AnalysisThreshold::findElectronHits()
/*****************************************************************/
{
    m_electronHits.clear();
    m_electronHitEnergy = 0.;
    for(const auto& hit : event().hardsimhits())
    {
        double deta = (hit.eta() - m_chosenHit->eta());
        double dphi = TVector2::Phi_mpi_pi(hit.phi() - m_chosenHit->phi());       
        if(fabs(deta)>0.05 || fabs(dphi)>0.1) continue;
        m_electronHits.push_back( HGCEEDetId(hit.detid()) );
        m_electronHitEnergy += hit.energy();
    }
}

/*****************************************************************/
void AnalysisThreshold::findElectronHitsBelowThreshold()
/*****************************************************************/
{
    m_electronHitsBelowThreshold.clear();
    for(const auto& id : m_electronHits)
    {
        const SimHit* hit = event().simhit(id);
        int triggerRegion = triggerRegionIndex(hit->eta(), hit->zside(), hit->sector(), hit->subsector());
        int nhits = triggerRegionHits(triggerRegion, hit->layer());
        float threshold = pileupThreshold(hit->eta(), hit->layer(), nhits);
        if(hit->energy()<threshold*mip) m_electronHitsBelowThreshold.push_back(hit);

    }
}


/*****************************************************************/
void AnalysisThreshold::fillRegionHits()
/*****************************************************************/
{
    for(auto& region_hits : m_regionHitsAboveTh)
    {
        region_hits.second.clear();
    }

    for(const auto& region_id : m_regionSimHits)
    {
        pair<int,int> region_layer = region_id.first;
        for(const auto& id_eta : region_id.second)
        {
            const HGCEEDetId& id = id_eta.first;
            const SimHit* hit = event().simhit(id);
            if(hit && hit->energy()>mip) m_regionHitsAboveTh[region_layer].push_back(hit);
        }
    }
}



/*****************************************************************/
void AnalysisThreshold:: fillPileupParams()
/*****************************************************************/
{
    ifstream stream(m_pileupParamsFile);
    if(!stream.is_open())
    {
        cout<<"ERROR: Cannot open pileup parameters file '"<<m_pileupParamsFile<<"'\n";
        return;
    }
    while(!stream.eof())
    {
        int up, layer, subeta;
        float a, b;
        stream >> up >> layer >> subeta >> a >> b;
        int eta = subeta + 6*up;
        m_pileupParams[make_pair(eta,layer)] = make_pair(a,b);
    }
    stream.close();
}


/*****************************************************************/
float AnalysisThreshold::pileupThreshold(float eta, int layer, int nhits)
/*****************************************************************/
{
    int etaIndex = 0;
    if     (fabs(eta)>=1.8 && fabs(eta)<2.)  etaIndex = 0;
    else if(fabs(eta)>=2.  && fabs(eta)<2.2) etaIndex = 1;
    else if(fabs(eta)>=2.2 && fabs(eta)<2.4) etaIndex = 2;
    else if(fabs(eta)>=2.4 && fabs(eta)<2.6) etaIndex = 3;
    else if(fabs(eta)>=2.6 && fabs(eta)<2.8) etaIndex = 4;
    else if(fabs(eta)>=2.8 && fabs(eta)<3.)  etaIndex = 5;
    else if(fabs(eta)>=1.5 && fabs(eta)<1.6) etaIndex = 6;
    else if(fabs(eta)>=1.6 && fabs(eta)<1.7) etaIndex = 7;
    else if(fabs(eta)>=1.7 && fabs(eta)<1.8) etaIndex = 8;
    const auto itr = m_pileupParams.find(make_pair(etaIndex,layer));
    if(itr==m_pileupParams.end())
    {
        cout<<"ERROR: cannot find pileup parameters for eta="<<eta<<",layer="<<layer<<"\n";
        return 1.;
    }
    double a = itr->second.first;
    double b = itr->second.second;
    return a*(double)nhits + b;
}


/*****************************************************************/
int AnalysisThreshold::triggerRegionIndex(float eta, int zside, int sector, int subsector)
/*****************************************************************/
{
    int iz = (zside==-1 ? 0 : 1);
    int subsec = (subsector==-1 ? 1 : 0);
    int sector30deg = (2*(sector-1) + subsec)/3;
    int up = (abs(eta)<=1.8 ? 1 : 0);
    int triggerRegion = sector30deg + 12*up + 24*iz;
    return triggerRegion;
}

/*****************************************************************/
int AnalysisThreshold::triggerRegionHits(int triggerRegion, int layer)
/*****************************************************************/
{
    const auto itr = m_regionHitsAboveTh.find(make_pair(triggerRegion, layer));
    if(itr==m_regionHitsAboveTh.end())
    {
        cout<<"ERROR: cannot find number of hits in region "<<triggerRegion<<" layer "<<layer<<"\n";
        return 0;
    }
    return itr->second.size();
}
