/**
 *  @file  AnalysisLayerCalibration.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    20/02/2015
 *
 *  @internal
 *     Created :  20/02/2015
 * Last update :  20/02/2015 16:26:31
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#ifndef ANALYSISLAYERCALIBRATION_H
#define ANALYSISLAYERCALIBRATION_H

#include "AnHiMaHGCAL/Core/interface/IAnalysis.h"
#include "AnHiMaHGCAL/Core/interface/EventAware.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"

#include "AnHiMaHGCAL/HGCALCommon/interface/TriggerEGammaAlgorithm.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SuperCluster.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"

#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TH2F.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>




namespace AnHiMa
{

    class AnalysisLayerCalibration: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisLayerCalibration();
            ~AnalysisLayerCalibration();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos();
            void fillGenElectrons();
            

            TriggerEGammaAlgorithm m_egammaAlgo;

            std::vector<const GenParticle*> m_genElectrons;
            const GenParticle* m_genParticle;

            std::vector<Tower> m_electronCones;
            const Tower* m_matchedElectronCone;

    };
}


#endif
