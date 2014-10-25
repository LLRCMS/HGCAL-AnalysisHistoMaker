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
#include <stdio.h>


#include "AnHiMaHGCAL/StudySeeding/interface/AnalysisSeeding.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"



using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisSeeding::AnalysisSeeding():IAnalysis(),
    m_random(12345)
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


    double eta = m_random.Uniform(1.5, 3.);
    double phi = m_random.Uniform(-M_PI, M_PI);
    const SimHit* choosenHit = NULL;
    double mindr = 9999.;
    for(auto itrHit=event().simhits().cbegin(); itrHit!=event().simhits().end(); itrHit++)
    {
        const SimHit& hit = *itrHit;
        //cout<<"detid="<<hit.detid()<<", subdet="<<hit.subdet()<<", zside="<<hit.zside()<<", layer="<<hit.layer()<<", sector="<<hit.sector()<<", subsector="<<hit.subsector()<<", cell="<<hit.cell()<<"\n";
        //uint32_t detid = hit.detid();
        //HGCEEDetId idtest(detid);
        //HGCEEDetId idtest2(ForwardSubdetector(hit.subdet()), hit.zside(), hit.layer(), hit.sector(), hit.subsector(), hit.cell());
        //if(idtest.cell()!=idtest2.cell()) cout<<"DIFF\n";
        //cout<<" Test1 detid="<<idtest.rawId()<<", subdet="<<idtest.subdet()<<", zside="<<idtest.zside()<<", layer="<<idtest.layer()<<", sector="<<idtest.sector()<<", subsector="<<idtest.subsector()<<", cell="<<idtest.cell()<<"\n";
        //cout<<" Test1 detid="<<idtest2.rawId()<<", subdet="<<idtest2.subdet()<<", zside="<<idtest2.zside()<<", layer="<<idtest2.layer()<<", sector="<<idtest2.sector()<<", subsector="<<idtest2.subsector()<<", cell="<<idtest2.cell()<<"\n";

        double deta = (hit.eta() - eta);
        double dphi = (hit.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<mindr)
        {
            choosenHit = &hit;
            mindr = dr;
        }
    }
    HGCEEDetId detid(ForwardSubdetector(choosenHit->subdet()), choosenHit->zside(), choosenHit->layer(), choosenHit->sector(), choosenHit->subsector(), choosenHit->cell());

    vector<DetId> northIds = m_hgcalNavigator.north(detid);
    vector<DetId> southIds = m_hgcalNavigator.south(detid);
    vector<DetId> eastIds  = m_hgcalNavigator.east(detid);
    vector<DetId> westIds = m_hgcalNavigator.west(detid);
    vector<HGCEEDetId> upIds    = m_hgcalNavigator.up(detid);
    vector<HGCEEDetId> downIds  = m_hgcalNavigator.down(detid);

    cout<<"Fired eta,phi = "<<eta<<","<<phi<<"\n";
    cout<<"  Chosen detid eta="<<choosenHit->eta()<<", phi="<<choosenHit->phi()<<", layer="<<choosenHit->layer()<<", cell="<<choosenHit->cell()<<"\n";
    cout<<"   detid="<<detid.rawId()<<", subdet="<<detid.subdet()<<", zsdetide="<<detid.zside()<<", layer="<<detid.layer()<<", sector="<<detid.sector()<<", subsector="<<detid.subsector()<<", cell="<<detid.cell()<<"\n";
    cout<<"   North Ids:";
    for(unsigned i=0;i<northIds.size();i++)
    {
        HGCEEDetId id(northIds[i]);
        cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    }
    cout<<"\n";
    cout<<"   South Ids:";
    for(unsigned i=0;i<southIds.size();i++)
    {
        HGCEEDetId id(southIds[i]);
        cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    }
    cout<<"\n";
    cout<<"   East Ids:";
    for(unsigned i=0;i<eastIds.size();i++)
    {
        HGCEEDetId id(eastIds[i]);
        cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    }
    cout<<"\n";
    cout<<"   West Ids:";
    for(unsigned i=0;i<westIds.size();i++)
    {
        HGCEEDetId id(westIds[i]);
        cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    }
    cout<<"\n";
    cout<<"   Up Ids:";
    for(unsigned i=0;i<upIds.size();i++)
    {
        HGCEEDetId id(upIds[i]);
        cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    }
    cout<<"\n";
    cout<<"   Down Ids:";
    for(unsigned i=0;i<downIds.size();i++)
    {
        HGCEEDetId id(downIds[i]);
        cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    }
    cout<<"\n";

    //std::pair<float,float> xy = m_hgcalNavigator.dddConstants()->locateCell(detid.cell(), detid.layer(), detid.subsector(), true);
    //cout<<"DDD ctes: x,y  = "<<xy.first<<","<<xy.second<<"\n";
    //if(detid.layer()+1<30)
    //{
    //    std::pair<int,int>     kcell = m_hgcalNavigator.dddConstants()->assignCell(xy.first, xy.second, detid.layer()+1, detid.subsector(), true);
    //    cout<<"DDD ctes: cell = "<<kcell.second<<"\n";
    //}


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

