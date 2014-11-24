/**
 *  @file  Tower.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    26/10/2014
 *
 *  @internal
 *     Created :  26/10/2014
 * Last update :  26/10/2014 16:09:15
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"

#include <iostream>


using namespace AnHiMa;
using namespace std;


/*****************************************************************/
bool AnHiMa::towerSort(const Tower* tw1, const Tower* tw2)
/*****************************************************************/
{
    return tw1->calibratedEt() > tw2->calibratedEt();
};




/*****************************************************************/
Tower::Tower():
    m_eta(0.),
    m_phi(0.),
    m_energy(0.),
    m_calibratedEnergy(0.),
    m_layerEnergies(31),
    m_layerHits(31)
/*****************************************************************/
{
    for(unsigned l=0;l<m_layerEnergies.size();l++)
    {
        m_layerEnergies[l] = 0.;
        m_layerHits[l] = 0;
    }
}

/*****************************************************************/
Tower::Tower(const Tower& tower):
    m_eta(tower.eta()),
    m_phi(tower.phi()),
    m_energy(tower.energy()),
    m_calibratedEnergy(tower.calibratedEnergy()),
    m_layerEnergies(31),
    m_layerHits(31)
/*****************************************************************/
{
    //cerr<<"Copy tower eta="<<m_eta<<",phi="<<m_phi<<",Et="<<calibratedEt()<<"\n";
    for(unsigned l=0;l<m_layerEnergies.size();l++)
    {
        m_layerEnergies[l] = tower.layerEnergy(l);
        m_layerHits[l] = tower.layerNHits(l);
    }
    for(const auto hit : tower.hits())
    {
        m_hits.push_back(hit);
    }
}



/*****************************************************************/
Tower::~Tower()
/*****************************************************************/
{
    //cerr<<"Destroying tower eta="<<m_eta<<",phi="<<m_phi<<",Et="<<calibratedEt()<<"\n";
}


/*****************************************************************/
void Tower::addHit(const SimHit& hit)
/*****************************************************************/
{
    m_hits.push_back(&hit);
    m_eta *= m_energy;
    m_phi *= m_energy;

    m_energy += hit.energy();
    m_eta += hit.eta()*hit.energy();
    m_phi += hit.phi()*hit.energy();
    m_eta /= m_energy;
    m_phi /= m_energy;
    int layer = hit.layer();
    m_layerEnergies[layer] += hit.energy();
    m_layerHits[layer]++;
}


/*****************************************************************/
float Tower::layersEnergy(unsigned l1, unsigned l2) const
/*****************************************************************/
{

    double energy = 0.;
    for(unsigned l=l1; l<=l2; l++)
    {
        energy += (double)m_layerEnergies[l];
    }
    return (float)energy;
}

/*****************************************************************/
int Tower::nLayers() const
/*****************************************************************/
{

    int n = 0;
    for(unsigned l=0; l<m_layerEnergies.size(); l++)
    {
        if(m_layerEnergies[l]>0) n++;
    }
    return n;
}


/*****************************************************************/
void Tower::clear()
/*****************************************************************/
{
    m_eta = 0.;
    m_phi = 0.;
    m_energy = 0.;
    m_hits.clear();
    for(unsigned l=0;l<m_layerEnergies.size();l++)
    {
        m_layerEnergies[l] = 0.;
        m_layerHits[l] = 0;
    }
}



/*****************************************************************/
float Tower::longitudinalBarycenter()
/*****************************************************************/
{
    double meanLayer = 0.;
    for(unsigned l=1;l<m_layerEnergies.size();l++)
    {
        meanLayer += (double)l * m_layerEnergies[l];
    }
    meanLayer /= m_energy;
    return meanLayer;

}
