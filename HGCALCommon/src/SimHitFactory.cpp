/**
 *  @file  SimHitFactory.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    22/10/2014
 *
 *  @internal
 *     Created :  22/10/2014
 * Last update :  22/10/2014 21:31:42
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */





#include "AnHiMaHGCAL/HGCALCommon/interface/SimHitFactory.h"

#include <iostream>

using namespace AnHiMa;
using namespace std;

/*****************************************************************/
bool AnHiMa::simHitSort(const SimHit* hit1, const SimHit* hit2)
/*****************************************************************/
{
    return hit1->energy() > hit2->energy();
};


/*****************************************************************/
void SimHitFactory::initialize(IEvent* event, TChain* chain)
/*****************************************************************/
{
    event->registerCallback((void*)this, SimHitFactory::callback);
    connectVariables(chain);

}

/*****************************************************************/
void SimHitFactory::initialize(IEvent* event, TChain* chain, HitType type)
/*****************************************************************/
{
    m_type = type;
    initialize(event, chain);
    m_data.push_back( vector<SimHit>() ); // HE
    m_data.push_back( vector<SimHit>() ); // FH
    m_data.push_back( vector<SimHit>() ); // BH
}



/*****************************************************************/
void SimHitFactory::connectVariables(TChain* chain)
/*****************************************************************/
{
    m_hit_n         = 0;
    m_hit_detid     = 0;
    m_hit_subdet    = 0;
    m_hit_cell      = 0;
    m_hit_sector    = 0;
    m_hit_subsector = 0;
    m_hit_layer     = 0;
    m_hit_zside     = 0;
    m_hit_energy    = 0;
    m_hit_eta       = 0;
    m_hit_phi       = 0;
    m_hit_x         = 0;
    m_hit_y         = 0;
    m_hit_z         = 0;

    switch(m_type)
    {
        case ALL:
            chain->SetBranchStatus("hit_n"        , true);
            chain->SetBranchStatus("hit_detid"    , true);
            chain->SetBranchStatus("hit_subdet"   , true);
            chain->SetBranchStatus("hit_cell"     , true);
            chain->SetBranchStatus("hit_sector"   , true);
            chain->SetBranchStatus("hit_subsector", true);
            chain->SetBranchStatus("hit_layer"    , true);
            chain->SetBranchStatus("hit_zside"    , true);
            chain->SetBranchStatus("hit_energy"   , true);
            chain->SetBranchStatus("hit_eta"      , true);
            chain->SetBranchStatus("hit_phi"      , true);
            chain->SetBranchStatus("hit_x"        , true);
            chain->SetBranchStatus("hit_y"        , true);
            chain->SetBranchStatus("hit_z"        , true);


            chain->SetBranchAddress("hit_n"        , &m_hit_n        );
            chain->SetBranchAddress("hit_detid"    , &m_hit_detid    );
            chain->SetBranchAddress("hit_subdet"   , &m_hit_subdet   );
            chain->SetBranchAddress("hit_cell"     , &m_hit_cell     );
            chain->SetBranchAddress("hit_sector"   , &m_hit_sector   );
            chain->SetBranchAddress("hit_subsector", &m_hit_subsector);
            chain->SetBranchAddress("hit_layer"    , &m_hit_layer    );
            chain->SetBranchAddress("hit_zside"    , &m_hit_zside    );
            chain->SetBranchAddress("hit_energy"   , &m_hit_energy   );
            chain->SetBranchAddress("hit_eta"      , &m_hit_eta      );
            chain->SetBranchAddress("hit_phi"      , &m_hit_phi      );
            chain->SetBranchAddress("hit_x"        , &m_hit_x        );
            chain->SetBranchAddress("hit_y"        , &m_hit_y        );
            chain->SetBranchAddress("hit_z"        , &m_hit_z        );
            break;
        case HARDINT:
            chain->SetBranchStatus("hard_hit_n"        , true);
            chain->SetBranchStatus("hard_hit_detid"    , true);
            chain->SetBranchStatus("hard_hit_subdet"   , true);
            chain->SetBranchStatus("hard_hit_cell"     , true);
            chain->SetBranchStatus("hard_hit_sector"   , true);
            chain->SetBranchStatus("hard_hit_subsector", true);
            chain->SetBranchStatus("hard_hit_layer"    , true);
            chain->SetBranchStatus("hard_hit_zside"    , true);
            chain->SetBranchStatus("hard_hit_energy"   , true);
            chain->SetBranchStatus("hard_hit_eta"      , true);
            chain->SetBranchStatus("hard_hit_phi"      , true);
            chain->SetBranchStatus("hard_hit_x"        , true);
            chain->SetBranchStatus("hard_hit_y"        , true);
            chain->SetBranchStatus("hard_hit_z"        , true);


            chain->SetBranchAddress("hard_hit_n"        , &m_hit_n        );
            chain->SetBranchAddress("hard_hit_detid"    , &m_hit_detid    );
            chain->SetBranchAddress("hard_hit_subdet"   , &m_hit_subdet   );
            chain->SetBranchAddress("hard_hit_cell"     , &m_hit_cell     );
            chain->SetBranchAddress("hard_hit_sector"   , &m_hit_sector   );
            chain->SetBranchAddress("hard_hit_subsector", &m_hit_subsector);
            chain->SetBranchAddress("hard_hit_layer"    , &m_hit_layer    );
            chain->SetBranchAddress("hard_hit_zside"    , &m_hit_zside    );
            chain->SetBranchAddress("hard_hit_energy"   , &m_hit_energy   );
            chain->SetBranchAddress("hard_hit_eta"      , &m_hit_eta      );
            chain->SetBranchAddress("hard_hit_phi"      , &m_hit_phi      );
            chain->SetBranchAddress("hard_hit_x"        , &m_hit_x        );
            chain->SetBranchAddress("hard_hit_y"        , &m_hit_y        );
            chain->SetBranchAddress("hard_hit_z"        , &m_hit_z        );
            break;
        default:
            break;
    }

}



/*****************************************************************/
void SimHitFactory::callback(void* object)
/*****************************************************************/
{
    SimHitFactory* myself = reinterpret_cast<SimHitFactory*>(object);
    myself->update();
}


/*****************************************************************/
void SimHitFactory::update()
/*****************************************************************/
{
    m_data[0].clear();
    m_data[1].clear();
    m_data[2].clear();
    
    for(int i=0;i<m_hit_n;i++)
    {
        SimHit hit;
        hit.setDetid     ( (unsigned)(*m_hit_detid)[i] );
        hit.setSubdet    ( (int)(*m_hit_subdet)[i] );
        hit.setCell      ( (int)(*m_hit_cell)[i] );
        hit.setSector    ( (int)(*m_hit_sector)[i] );
        hit.setSubsector ( (int)(*m_hit_subsector)[i] );
        hit.setLayer     ( (int)(*m_hit_layer)[i] );
        hit.setZside     ( (int)(*m_hit_zside)[i] );
        hit.setEnergy    ( (*m_hit_energy)[i] );
        hit.setEta       ( (*m_hit_eta)[i] );
        hit.setPhi       ( (*m_hit_phi)[i] );
        hit.setX         ( (*m_hit_x)[i] );
        hit.setY         ( (*m_hit_y)[i] );
        hit.setZ         ( (*m_hit_z)[i] );

        int iCol = hit.subdet() - 3;
        if(iCol>2 || iCol<0)
        {
            cerr<<"WARNING: Unknown subdetector '"<<hit.subdet()<<"'\n";
            continue;
        }
        m_data[iCol].push_back( hit );
    }
}
