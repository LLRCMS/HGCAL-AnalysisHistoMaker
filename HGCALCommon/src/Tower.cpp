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


using namespace AnHiMa;
using namespace std;

/*****************************************************************/
Tower::Tower():
    m_layerEnergies(31)
/*****************************************************************/
{
    for(unsigned l=0;l<m_layerEnergies.size();l++)
    {
        m_layerEnergies[l] = 0;
    }
}



/*****************************************************************/
Tower::~Tower()
/*****************************************************************/
{
}


/*****************************************************************/
void Tower::addHit(const SimHit& hit)
/*****************************************************************/
{
    m_hits.push_back(&hit);
    m_energy += hit.energy();
    int layer = hit.layer();
    m_layerEnergies[layer] += hit.energy();
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
