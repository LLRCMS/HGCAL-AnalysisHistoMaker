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

#include "TVector2.h"

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
bool AnHiMa::towerSortE(const Tower* tw1, const Tower* tw2)
/*****************************************************************/
{
    return tw1->calibratedEnergy() > tw2->calibratedEnergy();
};




/*****************************************************************/
Tower::Tower():
    m_eta(0.),
    m_phi(0.),
    m_x(0.),
    m_y(0.),
    m_energy(0.),
    m_calibratedEnergy(0.)
/*****************************************************************/
{
    // subdet 3 (EE)
    m_layerEnergies          .push_back(vector<float>(31));
    m_layerCalibratedEnergies.push_back(vector<float>(31));
    m_layerHits              .push_back(vector<int>  (31));
    // subdet 4 (HEF)
    m_layerEnergies          .push_back(vector<float>(13));
    m_layerCalibratedEnergies.push_back(vector<float>(13));
    m_layerHits              .push_back(vector<int>  (13));
    // subdet 5 (HEB)
    m_layerEnergies          .push_back(vector<float>(13));
    m_layerCalibratedEnergies.push_back(vector<float>(13));
    m_layerHits              .push_back(vector<int>  (13));

    // subdet
    for(unsigned s=0;s<m_layerEnergies.size();s++)
    {
        // layer
        for(unsigned l=0;l<m_layerEnergies[s].size();l++)
        {
            m_layerEnergies[s][l] = 0.;
            m_layerCalibratedEnergies[s][l] = 0.;
            m_layerHits[s][l] = 0;
        }
    }
}

/*****************************************************************/
Tower::Tower(const Tower& tower):
    m_eta(tower.eta()),
    m_phi(tower.phi()),
    m_x(tower.x()),
    m_y(tower.y()),
    m_energy(tower.energy()),
    m_calibratedEnergy(tower.calibratedEnergy())
/*****************************************************************/
{
    // subdet 3 (EE)
    m_layerEnergies          .push_back(vector<float>(31));
    m_layerCalibratedEnergies.push_back(vector<float>(31));
    m_layerHits              .push_back(vector<int>  (31));
    // subdet 4 (HEF)
    m_layerEnergies          .push_back(vector<float>(13));
    m_layerCalibratedEnergies.push_back(vector<float>(13));
    m_layerHits              .push_back(vector<int>  (13));
    // subdet 5 (HEB)
    m_layerEnergies          .push_back(vector<float>(13));
    m_layerCalibratedEnergies.push_back(vector<float>(13));
    m_layerHits              .push_back(vector<int>  (13));

    // subdet
    for(unsigned s=0;s<m_layerEnergies.size();s++)
    {
        // layer
        for(unsigned l=0;l<m_layerEnergies[s].size();l++)
        {
            m_layerEnergies[s][l] = tower.layerEnergy(l, s+3);
            m_layerCalibratedEnergies[s][l] = tower.layerCalibratedEnergy(l, s+3);
            m_layerHits[s][l] = tower.layerNHits(l, s+3);
        }
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
    //cout<<"Phi before = "<<m_phi<<"\n";
    //cout<<" + "<<hit.phi()<<"\n";
    double phidiff = TVector2::Phi_mpi_pi(hit.phi() - m_phi);
    double phiold = m_phi;

    m_eta *= m_energy;
    m_phi *= m_energy;
    m_x *= m_energy;
    m_y *= m_energy;

    m_energy += hit.energy();
    m_eta += hit.eta()*hit.energy();
    m_phi += (phiold + phidiff)*hit.energy();
    m_x += hit.x()*hit.energy();
    m_y += hit.y()*hit.energy();
    m_eta /= m_energy;
    m_phi /= m_energy;
    m_phi = TVector2::Phi_mpi_pi(m_phi);
    m_x /= m_energy;
    m_y /= m_energy;
    int layer = hit.layer();
    int subdet = hit.subdet();
    m_layerEnergies[subdet-3][layer] += hit.energy();
    m_layerHits[subdet-3][layer]++;

    //cout<<" Phi after = "<<m_phi<<"\n";
}


/*****************************************************************/
float Tower::layersEnergy(unsigned l1, unsigned l2, unsigned subdet) const
/*****************************************************************/
{

    double energy = 0.;
    for(unsigned l=l1; l<=l2; l++)
    {
        energy += (double)m_layerEnergies[subdet-3][l];
    }
    return (float)energy;
}

/*****************************************************************/
float Tower::layersCalibratedEnergy(unsigned l1, unsigned l2, unsigned subdet) const
/*****************************************************************/
{

    double energy = 0.;
    for(unsigned l=l1; l<=l2; l++)
    {
        energy += (double)m_layerCalibratedEnergies[subdet-3][l];
    }
    return (float)energy;
}

/*****************************************************************/
int Tower::nLayers(unsigned subdet) const
/*****************************************************************/
{

    int n = 0;
    for(unsigned l=0; l<m_layerEnergies[subdet-3].size(); l++)
    {
        if(m_layerEnergies[subdet-3][l]>0) n++;
    }
    return n;
}


/*****************************************************************/
void Tower::clear()
/*****************************************************************/
{
    m_eta = 0.;
    m_phi = 0.;
    m_x = 0.;
    m_y = 0.;
    m_energy = 0.;
    m_hits.clear();
    // subdet
    for(unsigned s=0;s<m_layerEnergies.size();s++)
    {
        // layer
        for(unsigned l=0;l<m_layerEnergies[s].size();l++)
        {
            m_layerEnergies[s][l] = 0.;
            m_layerCalibratedEnergies[s][l] = 0.;
            m_layerHits[s][l] = 0;
        }
    }
}



/*****************************************************************/
float Tower::longitudinalBarycenter(unsigned subdet) const
/*****************************************************************/
{
    double meanLayer = 0.;
    for(unsigned l=1;l<m_layerEnergies[subdet-3].size();l++)
    {
        meanLayer += (double)l * m_layerEnergies[subdet-3][l];
    }
    meanLayer /= m_energy;
    return meanLayer;

}
