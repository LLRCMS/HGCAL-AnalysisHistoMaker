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
#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"

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


    const SimHit* chosenHit = NULL;

    while(!chosenHit)
    {
        double eta = m_random.Uniform(1.5, 3.);
        double phi = m_random.Uniform(-M_PI, M_PI);
        double maxE = 0.;
        for(auto itrHit=event().simhits().cbegin(); itrHit!=event().simhits().end(); itrHit++)
        {
            const SimHit& hit = *itrHit;
            double deta = (hit.eta() - eta);
            double dphi = (hit.phi() - phi);
            double dr = sqrt(deta*deta + dphi*dphi);
            if(dr<0.2 && hit.energy()>maxE)
            {
                chosenHit = &hit;
                maxE = hit.energy();
            }
        }
    }

    HGCEEDetId detid(ForwardSubdetector(chosenHit->subdet()), chosenHit->zside(), chosenHit->layer(), chosenHit->sector(), chosenHit->subsector(), chosenHit->cell());

    // find cells in the tower at (eta,phi)
    int layer = detid.layer();
    vector<HGCEEDetId> idLayers(31);
    vector<HGCEEDetId> idLayer15 = m_hgcalNavigator.upProj(detid, 15-layer);
    if(idLayer15.size()>0) idLayers[15] = idLayer15[0];
    else
    {
        cout<<"[WARNING] Cannot find cell in layer 15\n";
        return;
    }
    for(int l=1;l<=30;l++)
    {
        if(l==15) continue;
        vector<HGCEEDetId> ids = m_hgcalNavigator.upProj(idLayers[15], l-15);
        if(ids.size()==0)
        {
            cout<<"[WARNING] Cannot find cell in layer "<<l<<"\n";
            idLayers[l] = (HGCEEDetId)DetId(0);
        }
        else
        {
            idLayers[l] = ids[0];
        }
    }

    // build tower from the list of cells
    Tower tower;
    tower.setEta(chosenHit->eta());
    tower.setPhi(chosenHit->phi());
    for(unsigned l=0; l<idLayers.size(); l++)
    {
        HGCEEDetId id = idLayers[l];
        if(id==DetId(0)) continue;
        const SimHit* hit = event().simhit(id);
        if(!hit) continue;
        tower.addHit(*hit);
    }

    cout<<"Tower (eta,phi)=("<<tower.eta()<<","<<tower.phi()<<"): #hits="<<tower.nHits()<<", energy="<<tower.energy()<<"\n";




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


