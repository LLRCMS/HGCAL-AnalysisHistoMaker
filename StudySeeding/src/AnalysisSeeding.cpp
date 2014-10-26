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

    vector<DetId> northIds      = m_hgcalNavigator.north(detid);
    vector<DetId> southIds      = m_hgcalNavigator.south(detid);
    vector<DetId> eastIds       = m_hgcalNavigator.east(detid);
    vector<DetId> westIds       = m_hgcalNavigator.west(detid);
    vector<HGCEEDetId> upIds    = m_hgcalNavigator.up(detid);
    vector<HGCEEDetId> downIds  = m_hgcalNavigator.down(detid);

    cout<<"Fired eta,phi = "<<eta<<","<<phi<<"\n";
    cout<<"  Chosen detid eta="<<choosenHit->eta()<<", phi="<<choosenHit->phi()<<", layer="<<choosenHit->layer()<<", cell="<<choosenHit->cell()<<"\n";

    const SimHit* hitNorth = ( northIds.size()>0 ? event().simhit( (HGCEEDetId)northIds[0] ) : NULL);
    const SimHit* hitSouth = ( southIds.size()>0 ? event().simhit( (HGCEEDetId)southIds[0] ) : NULL);
    const SimHit* hitEast  = ( eastIds.size()>0  ? event().simhit( (HGCEEDetId)eastIds[0] )  : NULL);
    const SimHit* hitWest  = ( westIds.size()>0  ? event().simhit( (HGCEEDetId)westIds[0] )  : NULL);
    const SimHit* hitUp    = ( upIds.size()>0    ? event().simhit( (HGCEEDetId)upIds[0] )    : NULL);
    const SimHit* hitDown  = ( downIds.size()>0  ? event().simhit( (HGCEEDetId)downIds[0] )  : NULL);

    cout<<"  North hit: ";
    if(hitNorth) cout << "eta="<<hitNorth->eta()<<", phi="<<hitNorth->phi()<<", layer="<<hitNorth->layer()<<", cell="<<hitNorth->cell()<<", sec="<<hitNorth->sector()<<", subsec="<<hitNorth->subsector()<<"\n";
    else cout << "NULL\n";

    cout<<"  South hit: ";
    if(hitSouth) cout << "eta="<<hitSouth->eta()<<", phi="<<hitSouth->phi()<<", layer="<<hitSouth->layer()<<", cell="<<hitSouth->cell()<<", sec="<<hitSouth->sector()<<", subsec="<<hitSouth->subsector()<<"\n";
    else cout << "NULL\n";

    cout<<"  East hit: ";
    if(hitEast) cout << "eta="<<hitEast->eta()<<", phi="<<hitEast->phi()<<", layer="<<hitEast->layer()<<", cell="<<hitEast->cell()<<", sec="<<hitEast->sector()<<", subsec="<<hitEast->subsector()<<"\n";
    else cout << "NULL\n";

    cout<<"  West hit: ";
    if(hitWest) cout << "eta="<<hitWest->eta()<<", phi="<<hitWest->phi()<<", layer="<<hitWest->layer()<<", cell="<<hitWest->cell()<<", sec="<<hitWest->sector()<<", subsec="<<hitWest->subsector()<<"\n";
    else cout << "NULL\n";

    cout<<"  Up hit: ";
    if(hitUp) cout << "eta="<<hitUp->eta()<<", phi="<<hitUp->phi()<<", layer="<<hitUp->layer()<<", cell="<<hitUp->cell()<<", sec="<<hitUp->sector()<<", subsec="<<hitUp->subsector()<<"\n";
    else cout << "NULL\n";

    cout<<"  Down hit: ";
    if(hitDown) cout << "eta="<<hitDown->eta()<<", phi="<<hitDown->phi()<<", layer="<<hitDown->layer()<<", cell="<<hitDown->cell()<<", sec="<<hitDown->sector()<<", subsec="<<hitDown->subsector()<<"\n";
    else cout << "NULL\n";



    //cout<<"   detid="<<detid.rawId()<<", subdet="<<detid.subdet()<<", zsdetide="<<detid.zside()<<", layer="<<detid.layer()<<", sector="<<detid.sector()<<", subsector="<<detid.subsector()<<", cell="<<detid.cell()<<"\n";
    //cout<<"   North Ids:";
    //for(unsigned i=0;i<northIds.size();i++)
    //{
    //    HGCEEDetId id(northIds[i]);
    //    cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    //}
    //cout<<"\n";
    //cout<<"   South Ids:";
    //for(unsigned i=0;i<southIds.size();i++)
    //{
    //    HGCEEDetId id(southIds[i]);
    //    cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    //}
    //cout<<"\n";
    //cout<<"   East Ids:";
    //for(unsigned i=0;i<eastIds.size();i++)
    //{
    //    HGCEEDetId id(eastIds[i]);
    //    cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    //}
    //cout<<"\n";
    //cout<<"   West Ids:";
    //for(unsigned i=0;i<westIds.size();i++)
    //{
    //    HGCEEDetId id(westIds[i]);
    //    cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    //}
    //cout<<"\n";
    //cout<<"   Up Ids:";
    //for(unsigned i=0;i<upIds.size();i++)
    //{
    //    HGCEEDetId id(upIds[i]);
    //    cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    //}
    //cout<<"\n";
    //cout<<"   Down Ids:";
    //for(unsigned i=0;i<downIds.size();i++)
    //{
    //    HGCEEDetId id(downIds[i]);
    //    cout<<" detid="<<id.rawId()<<", subdet="<<id.subdet()<<", zside="<<id.zside()<<", layer="<<id.layer()<<", sector="<<id.sector()<<", subsector="<<id.subsector()<<", cell="<<id.cell();
    //}
    //cout<<"\n";


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


