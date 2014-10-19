/**
 *  @file  Event.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    05/29/2014
 *
 *  @internal
 *     Created :  05/29/2014
 * Last update :  05/29/2014 04:45:27 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#ifndef EVENT_H
#define EVENT_H

#include "AnHiMaHGCAL/Core/interface/IEvent.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"

#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <string>




namespace AnHiMa
{
    static const int MAXGEN=10000;
    static const int MAXHIT=10000000;
    static const int MAXHITREDUCED=10000;

    class Event: public IEvent
    {
        public:
            Event();
            ~Event();

            bool passSelection(int selection=0);
            void connectVariables(TChain* inputChain);
            void buildObjects();

            //-------------
            // Input variables
            //-------------
            int run;
            int lumi;
            int event;

            // gen particles
            int   ngen;
            int   gen_id    [MAXGEN];
            int   gen_status[MAXGEN];
            float gen_pt    [MAXGEN];
            float gen_eta   [MAXGEN];
            float gen_phi   [MAXGEN];
            float gen_en    [MAXGEN];

            // hits
            int   nhits;      
            int   hit_type [MAXHIT];
            int   hit_layer[MAXHIT];
            int   hit_sec  [MAXHIT];
            int   hit_bin  [MAXHIT];
            float hit_x    [MAXHIT];
            float hit_y    [MAXHIT];
            float hit_z    [MAXHIT];
            float hit_eta  [MAXHIT];
            float hit_phi  [MAXHIT];
            float hit_edep [MAXHIT];

            // harder scatter hits
            int   hard_nhits;      
            int   hard_hit_type [MAXHIT];
            int   hard_hit_layer[MAXHIT];
            int   hard_hit_sec  [MAXHIT];
            int   hard_hit_bin  [MAXHIT];
            float hard_hit_x    [MAXHIT];
            float hard_hit_y    [MAXHIT];
            float hard_hit_z    [MAXHIT];
            float hard_hit_eta  [MAXHIT];
            float hard_hit_phi  [MAXHIT];
            float hard_hit_edep [MAXHIT];


            const std::vector<SimHit>& simHits() const {return m_simHits;}
            const std::vector<SimHit>& hardSimHits() const {return m_hardSimHits;}
            double energySum() const {return m_energySum;}
            double hardEnergySum() const {return m_hardEnergySum;}

        private:
            std::vector<SimHit> m_simHits;
            std::vector<SimHit> m_hardSimHits;
            double m_energySum;
            double m_hardEnergySum;

    };
}

#endif
