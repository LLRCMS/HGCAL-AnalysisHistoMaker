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

#include <string>
#include <vector>
#include <map>
#include <array>
#include <algorithm>


class HGCALNavigator
{
    public:
        HGCALNavigator();
        ~HGCALNavigator();

        bool initialize(const std::string& triggerLeftNavigationMapping, const std::string& triggerRightNavigationMapping, const std::string& triggerUpNavigationMapping, const std::string& triggerDownNavigationMapping);

        bool valid(const DetId& id, int subdet=3) const {return m_hgctopo.at(subdet)->valid(id);}; 
        std::vector<DetId> north(const DetId& id, int subdet=3)const {return m_hgctopo.at(subdet)->north(id);};
        std::vector<DetId> south(const DetId& id, int subdet=3)const {return m_hgctopo.at(subdet)->south(id);};
        std::vector<DetId> east(const DetId& id, int subdet=3) const {return m_hgctopo.at(subdet)->east(id);};
        std::vector<DetId> west(const DetId& id, int subdet=3) const {return m_hgctopo.at(subdet)->west(id);};
        // FIXME: topo.up() and down() give no DetId... 
        //std::vector<DetId> up(const DetId& id)   {return m_hgctopo->up(id);};
        //std::vector<DetId> down(const DetId& id) {return m_hgctopo->down(id);};
        // temporary workaround to navigate up and down
        std::vector<HGCEEDetId> up(const HGCEEDetId& id, int nz=1) const;
        std::vector<HGCEEDetId> down(const HGCEEDetId& id, int nz=1) const;

        std::vector<HGCEEDetId> upProj(const HGCEEDetId& id, int nz=1, double refEta=999., double refPhi=999.) const;
        std::vector<HGCEEDetId> downProj(const HGCEEDetId& id, int nz=1, double refEta=999., double refPhi=999.) const;

        // trigger navigation
        std::vector<HGCEEDetId> northTrigger(const HGCEEDetId& id) const;
        std::vector<HGCEEDetId> southTrigger(const HGCEEDetId& id) const;
        std::vector<HGCEEDetId> eastTrigger (const HGCEEDetId& id) const;
        std::vector<HGCEEDetId> westTrigger (const HGCEEDetId& id) const;

        const HGCalGeometry* geometry(int subdet=3) const {return m_hgcgeom.at(subdet);};
        const HGCalTopology* topology(int subdet=3) const {return m_hgctopo.at(subdet);};
        const HGCalDDDConstants* dddConstants(int subdet=3) const {return m_hgcdc.at(subdet);};


    private:
        bool loadTriggerLeftNavigation(const std::string& navigationMapping);
        bool loadTriggerRightNavigation(const std::string& navigationMapping);
        bool loadTriggerUpNavigation(const std::string& navigationMapping);
        bool loadTriggerDownNavigation(const std::string& navigationMapping);

        DDCompactView* m_ddcv;
        std::map<int, HGCalDDDConstants*> m_hgcdc;
        std::map<int, HGCalTopology*> m_hgctopo;
        std::map<int, HGCalGeometry*> m_hgcgeom;

        std::array< std::map<short, std::vector<std::pair<short,short>>>, 30> m_triggerLeftNavigationMapping;
        std::array< std::map<short, std::vector<std::pair<short,short>>>, 30> m_triggerRightNavigationMapping;
        std::array< std::map<short, std::vector<std::pair<short,short>>>, 30> m_triggerUpNavigationMapping;
        std::array< std::map<short, std::vector<std::pair<short,short>>>, 30> m_triggerDownNavigationMapping;

};


#endif
