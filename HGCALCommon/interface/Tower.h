/**
 *  @file  Tower.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    26/10/2014
 *
 *  @internal
 *     Created :  26/10/2014
 * Last update :  26/10/2014 16:03:04
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef TOWER_H
#define TOWER_H

#include <vector>

#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"

namespace AnHiMa
{
    class Tower
    {
        public:
            Tower();
            ~Tower();

            void addHit(const SimHit& hit);

            void setEta(float eta) {m_eta = eta;}
            void setPhi(float phi) {m_phi = phi;}
            void setEnergy(float energy) {m_energy = energy;}

            float eta()    const {return m_eta;}
            float phi()    const {return m_phi;}
            float energy() const {return m_energy;}
            int nHits() const {return m_hits.size();}
            float layerEnergy(unsigned l) const {return m_layerEnergies[l];}
            float layersEnergy(unsigned l1, unsigned l2) const;
            float longitudinalBarycenter();

            void clear();


        private:
            float m_eta;
            float m_phi;
            float m_energy;

            std::vector<float> m_layerEnergies;
            std::vector<const SimHit*> m_hits;
    };
};

#endif
