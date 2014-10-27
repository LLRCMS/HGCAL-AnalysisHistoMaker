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

#include "TVector2.h"

#include "AnHiMaHGCAL/StudySeeding/interface/AnalysisSeeding.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"


#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"



using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisSeeding::AnalysisSeeding():IAnalysis(),
    m_chosenHit(0),
    m_genParticle(0),
    m_random(12345),
    m_sample("")
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

    m_sample = m_reader.params().GetValue("Sample", "minbias");

    m_hgcalNavigator.initialize();


    return true;
}


/*****************************************************************/
void AnalysisSeeding::execute()
/*****************************************************************/
{
    event().update();
    if(!event().passSelection()) return;

    unsigned n = (m_sample=="minbias" ? 10 : 2);

    for(unsigned i=0;i<n;i++)
    {
        findMaxHit(i);
        if(!m_chosenHit)
        {
            cout<<"[WARNING] Cannot match a hit close to the chosen (eta,phi) direction\n";
            continue;
        }
        buildTower();
        if(m_tower.energy()==0) 
        {
            cout<<"[WARNING] Tower has 0 energy\n";
            continue;
        }

        fillHistos();
    }

}

/*****************************************************************/
void AnalysisSeeding::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;

    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events

    for(unsigned l=1; l<=30; l++)
    {
        double fraction = m_tower.layerEnergy(l)/m_tower.energy();
        m_histos.Fill2BinHisto(100+hoffset, m_tower.eta(), l, fraction, weight, sysNum);
    }
    double fraction1 = m_tower.layersEnergy(1,10)/m_tower.energy();
    double fraction2 = m_tower.layersEnergy(11,20)/m_tower.energy();
    double fraction3 = m_tower.layersEnergy(21,30)/m_tower.energy();
    m_histos.Fill2BinHisto(200+hoffset, m_tower.eta(), 1, fraction1, weight, sysNum);
    m_histos.Fill2BinHisto(200+hoffset, m_tower.eta(), 2, fraction2, weight, sysNum);
    m_histos.Fill2BinHisto(200+hoffset, m_tower.eta(), 3, fraction3, weight, sysNum);

    double middleOverFront = ( m_tower.layersEnergy(1,10)>0. ? m_tower.layersEnergy(11,20)/m_tower.layersEnergy(1,10) : 19.999);
    if(middleOverFront>=20.) middleOverFront = 19.999;
    m_histos.Fill1BinHisto(250+hoffset, m_tower.eta(), middleOverFront, weight, sysNum);

}


/*****************************************************************/
void AnalysisSeeding::buildTower()
/*****************************************************************/
{
    m_tower.clear();

    HGCEEDetId detid(ForwardSubdetector(m_chosenHit->subdet()), m_chosenHit->zside(), m_chosenHit->layer(), m_chosenHit->sector(),m_chosenHit->subsector(), m_chosenHit->cell());

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
    m_tower.setEta(m_chosenHit->eta());
    m_tower.setPhi(m_chosenHit->phi());
    for(unsigned l=0; l<idLayers.size(); l++)
    {
        HGCEEDetId id = idLayers[l];
        if(id==DetId(0)) continue;
        const SimHit* hit = event().simhit(id);
        if(!hit) continue;
        m_tower.addHit(*hit);
    }



    //cout<<"Tower (eta,phi)=("<<m_tower.eta()<<","<<m_tower.phi()<<"): #hits="<<m_tower.nHits()<<", energy="<<m_tower.energy()<<"\n";



}


/*****************************************************************/
void AnalysisSeeding::findMaxHit(int i)
/*****************************************************************/
{
    m_chosenHit = NULL;

    if(m_sample=="electron")
    {
        m_genParticle = &event().genparticles()[i];
        double eta = m_genParticle->Eta();
        double phi = m_genParticle->Phi();
        double maxE = 0.;
        for(auto itrHit=event().simhits().cbegin(); itrHit!=event().simhits().end(); itrHit++)
        {
            const SimHit& hit = *itrHit;
            // match only to hits in layers 10-20
            if(hit.layer()<10 || hit.layer()>20) continue;

            double deta = (hit.eta() - eta);
            double dphi = TVector2::Phi_mpi_pi(hit.phi() - phi);
            double dr = sqrt(deta*deta + dphi*dphi);
            if(dr<0.2 && hit.energy()>maxE)
            {
                m_chosenHit = &hit;
                maxE = hit.energy();
            }
        }
    }
    else if(m_sample=="electronPU")
    {
        m_genParticle = &event().genparticles()[i];
        double eta = m_genParticle->Eta();
        double phi = m_genParticle->Phi();
        double maxE = 0.;
        for(auto itrHit=event().hardsimhits().cbegin(); itrHit!=event().hardsimhits().end(); itrHit++)
        {
            const SimHit& hit = *itrHit;
            // match only to hits in layers 10-20
            if(hit.layer()<10 || hit.layer()>20) continue;

            double deta = (hit.eta() - eta);
            double dphi = TVector2::Phi_mpi_pi(hit.phi() - phi);
            double dr = sqrt(deta*deta + dphi*dphi);
            if(dr<0.2 && hit.energy()>maxE)
            {
                m_chosenHit = &hit;
                maxE = hit.energy();
            }
        }
    }
    else // MinBias. Choose random direction
    {
        while(!m_chosenHit)
        {
            double eta = m_random.Uniform(1.5, 3.);
            double phi = m_random.Uniform(-M_PI, M_PI);
            double maxE = 0.;
            for(auto itrHit=event().simhits().cbegin(); itrHit!=event().simhits().end(); itrHit++)
            {
                const SimHit& hit = *itrHit;
                double deta = (hit.eta() - eta);
                double dphi = TVector2::Phi_mpi_pi(hit.phi() - phi);
                double dr = sqrt(deta*deta + dphi*dphi);
                if(dr<0.2 && hit.energy()>maxE)
                {
                    m_chosenHit = &hit;
                    maxE = hit.energy();
                }
            }
        }
    }

}

