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
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"






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

    // initialize vectors
    for(int s=0;s<12;s++)
    {
        for(int z=0;z<2;z++)
        {
            for(int u=0;u<2;u++)
            {
                int triggerRegion = s + 12*u + 24*z;
                // ECAL
                for(int l=1;l<=30;l++)
                {
                    m_regionSimHits[make_pair(triggerRegion, l)] = vector<pair<HGCEEDetId,float> >(0);
                }
                // HCAL
                for(int l=1;l<=12;l++)
                {
                    m_regionSimHitsHCal[make_pair(triggerRegion, l)] = vector<pair<HGCHEDetId,float> >(0);
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
                // ECAL
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
                // HCAL
                for (int lay=1; lay<=12; ++lay) 
                {
                    for (int cell=0; cell<4000; ++cell) 
                    {
                        const HGCHEDetId id(HGCHEF,iz,lay,sec,ssec,cell);
                        if (event().hgcalNavigator().valid(id,4)) 
                        {
                            double eta = event().hgcalNavigator().geometry(4)->getPosition(id).eta();
                            //int up = (cell>=1200 ? 1 : 0);
                            int up = (fabs(eta)<=1.8 ? 1 : 0);
                            int triggerRegion = sector30deg + 12*up + 24*izz;
                            m_regionSimHitsHCal[make_pair(triggerRegion, lay)].push_back(make_pair(id,eta));
                        }
                    }
                }
            }
        }
    }

    return true;
}


/*****************************************************************/
void AnalysisPileup::execute()
/*****************************************************************/
{
    event().update();
    if(!event().passSelection()) return;


    fillHistos();


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

    // ECAL
    for(auto& region_vector : m_regionSimHits)
    {
        int region = region_vector.first.first;
        int layer = region_vector.first.second;
        int up = (region/12)%2;
        int nHits = 0;
        for(const auto& id: region_vector.second)
        {
            const SimHit* hit = event().simhit(id.first);
            if(!hit) continue;
            if(hit->energy()>=1.*mip) nHits++;
        }
        m_histos.Fill2BinHisto(100+hoffset, up, layer, nHits, weight, sysNum);
        m_histos.FillHisto(50000+hoffset, up, weight, sysNum);
        m_histos.FillHisto(50001+hoffset, 3, weight, sysNum);
        m_histos.FillHisto(50002+hoffset, layer, weight, sysNum);
        m_histos.FillHisto(50003+hoffset, nHits, weight, sysNum);
        m_histos.FillNtuple(50020+hoffset, event().run(), event().event(), weight, sysNum);
        for(const auto& id: region_vector.second)
        {
            double energy = 0.;
            double eta = id.second;
            m_histos.FillHisto(10+up+hoffset, fabs(eta), weight, sysNum);
            const SimHit* hit = event().simhit(id.first);
            if(!hit) continue;
            energy = hit->energy();
            m_histos.Fill2BinHisto(1000+hoffset, fabs(eta), layer, energy/mip, weight, sysNum);
            //cout<<eta<<", "<<layer<<", "<<nHits<<", "<<energy/mip<<"\n";
            m_histos.Fill2BinHisto(2000+3000*abs(up-1)+100*(layer-1)+hoffset,fabs(eta), nHits, energy/mip, weight, sysNum);
        }
    }

    // HCAL
    for(auto& region_vector : m_regionSimHitsHCal)
    {
        int region = region_vector.first.first;
        int layer = region_vector.first.second;
        int up = (region/12)%2;
        int nHits = 0;
        for(const auto& id: region_vector.second)
        {
            const SimHit* hit = event().simhitFH(id.first);
            if(!hit) continue;
            if(hit->energy()>=1.*mip_HF) nHits++;
        }
        m_histos.Fill2BinHisto(200+hoffset, up, layer, nHits, weight, sysNum);
        m_histos.FillHisto(50000+hoffset, up, weight, sysNum);
        m_histos.FillHisto(50001+hoffset, 4, weight, sysNum);
        m_histos.FillHisto(50002+hoffset, layer, weight, sysNum);
        m_histos.FillHisto(50003+hoffset, nHits, weight, sysNum);
        m_histos.FillNtuple(50020+hoffset, event().run(), event().event(), weight, sysNum);
        for(const auto& id: region_vector.second)
        {
            double energy = 0.;
            double eta = id.second;
            //m_histos.FillHisto(10+up+hoffset, fabs(eta), weight, sysNum);
            const SimHit* hit = event().simhitFH(id.first);
            if(!hit) continue;
            energy = hit->energy();
            //m_histos.Fill2BinHisto(1000+hoffset, fabs(eta), layer, energy/mip, weight, sysNum);
            //cout<<eta<<", "<<layer<<", "<<nHits<<", "<<energy/mip<<"\n";
            m_histos.Fill2BinHisto(8000+3000*abs(up-1)+100*(layer-1)+hoffset,fabs(eta), nHits, energy/mip_HF, weight, sysNum);
        }
    }

}


