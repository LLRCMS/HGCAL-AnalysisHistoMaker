/**
 *  @file  HGCALNavigator.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    19/10/2014
 *
 *  @internal
 *     Created :  19/10/2014
 * Last update :  19/10/2014 10:57:36
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef HGCALNAVIGATOR_H
#define HGCALNAVIGATOR_H


#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"


class HGCALNavigator
{
    public:
        HGCALNavigator();
        ~HGCALNavigator();

        bool initialize();

        bool valid(const DetId& id) const {return m_hgctopo->valid(id);}; 
        std::vector<DetId> north(const DetId& id){return m_hgctopo->north(id);};
        std::vector<DetId> south(const DetId& id){return m_hgctopo->south(id);};
        std::vector<DetId> east(const DetId& id) {return m_hgctopo->east(id);};
        std::vector<DetId> west(const DetId& id) {return m_hgctopo->west(id);};
        // FIXME: topo.up() and down() give no DetId... 
        //std::vector<DetId> up(const DetId& id)   {return m_hgctopo->up(id);};
        //std::vector<DetId> down(const DetId& id) {return m_hgctopo->down(id);};
        // temporary workaround to navigate up and down
        std::vector<HGCEEDetId> up(const HGCEEDetId& id, int nz=1);
        std::vector<HGCEEDetId> down(const HGCEEDetId& id, int nz=1);

        std::vector<HGCEEDetId> upProj(const HGCEEDetId& id, int nz=1);
        std::vector<HGCEEDetId> downProj(const HGCEEDetId& id, int nz=1);

        const HGCalGeometry* geometry() const {return m_hgcgeom;};
        const HGCalTopology* topology() const {return m_hgctopo;};
        const HGCalDDDConstants* dddConstants() const {return m_hgcdc;};


    private:
        DDCompactView* m_ddcv;
        HGCalDDDConstants* m_hgcdc;
        HGCalTopology* m_hgctopo;
        HGCalGeometry* m_hgcgeom;

};


#endif
