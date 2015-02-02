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
#include <math.h>

#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"

namespace AnHiMa
{

    class Tower
    {
        public:
            Tower();
            Tower(const Tower& tower);
            ~Tower();

            void addHit(const SimHit& hit);

            void setEta(float eta) {m_eta = eta;}
            void setPhi(float phi) {m_phi = phi;}
            void setX(float x) {m_x = x;}
            void setY(float y) {m_y = y;}
            void setEnergy(float energy) {m_energy = energy;}
            void setCalibratedEnergy(float energy) {m_calibratedEnergy = energy;}
            void setLayerCalibratedEnergy(unsigned l, float energy, unsigned subdet=3) {m_layerCalibratedEnergies[subdet-3][l] = energy;}

            float eta()     const {return m_eta;}
            float phi()     const {return m_phi;}
            float x()     const {return m_x;}
            float y()     const {return m_y;}
            float energy()  const {return m_energy;}
            float calibratedEnergy()  const {return m_calibratedEnergy;}
            float et()      const {return m_energy/cosh(m_eta);}
            float calibratedEt()  const {return m_calibratedEnergy/cosh(m_eta);}
            int   nHits()   const {return m_hits.size();}
            int   nLayers(unsigned subdet=3) const;
            int   layerNHits(unsigned l, unsigned subdet=3)   const {return m_layerHits[subdet-3][l];}
            float layerEnergy(unsigned l, unsigned subdet=3) const {return m_layerEnergies[subdet-3][l];}
            float layersEnergy(unsigned l1, unsigned l2, unsigned subdet=3) const;
            float layerCalibratedEnergy(unsigned l, unsigned subdet=3) const {return m_layerCalibratedEnergies[subdet-3][l];}
            float layersCalibratedEnergy(unsigned l1, unsigned l2, unsigned subdet=3) const;
            float longitudinalBarycenter(unsigned subdet=3) const;
            const std::vector<const SimHit*>& hits() const {return m_hits;}

            void clear();


        private:
            float m_eta;
            float m_phi;
            float m_x;
            float m_y;
            float m_energy;
            float m_calibratedEnergy;

            std::vector< std::vector<float> > m_layerEnergies;
            std::vector< std::vector<float> > m_layerCalibratedEnergies;
            std::vector< std::vector<int> > m_layerHits;
            std::vector<const SimHit*> m_hits;
    };

    bool towerSort(const Tower* tw1, const Tower* tw2);
    bool towerSortE(const Tower* tw1, const Tower* tw2);

};

#endif
