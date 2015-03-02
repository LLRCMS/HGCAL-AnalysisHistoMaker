/**
 *  @file  TriggerHitFactory.cpp
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





#include "AnHiMaHGCAL/HGCALCommon/interface/TriggerHitFactory.h"

#include <utility>
#include <iostream>
#include <fstream>

#include "TVector2.h"

using namespace AnHiMa;
using namespace std;


/*****************************************************************/
void TriggerHitFactory::initialize(IEvent* event, TChain* chain)
/*****************************************************************/
{
    event->registerCallback((void*)this, TriggerHitFactory::callback);
    connectVariables(chain);

}

/*****************************************************************/
void TriggerHitFactory::initialize(IEvent* event, TChain* chain, const SimHitFactory& hitFactory, const string& triggerCellMapping)
/*****************************************************************/
{
    m_hitFactory = &hitFactory;
    initialize(event, chain);
    m_data.push_back( vector<SimHit>() ); // HE
    m_data.push_back( vector<SimHit>() ); // FH
    m_data.push_back( vector<SimHit>() ); // BH

    m_cellToTriggerCell.resize(30);

    // Read trigger cell mapping
    ifstream ifs(triggerCellMapping, std::ifstream::in);
    if(!ifs.is_open()) 
    {
        cerr<<"ERROR: TriggerHitFactory::initialize(): Cannot open trigger cell mapping '"<<triggerCellMapping<<"'\n";
        return;
    }
    short layer = 0;
    short cell = 0;
    short triggerCell = 0;
    short subsector = 0;
    while(!ifs.eof())
    {
        ifs >> layer >> cell >> triggerCell >> subsector;
        if(layer>=30 || layer<0) 
        {
            cerr<<"ERROR: TriggerHitFactory::initialize(): Found layer '"<<layer<<"' in trigger cell mapping. Discard it.\n";
            continue;
        }
        short newcell = (subsector==1 ? cell+1 : -(cell+1));
        auto ret = m_cellToTriggerCell[layer].insert( std::make_pair(newcell, triggerCell) );
        if(ret.second==false)
        {
            cerr<<"ERROR: TriggerHitFactory::initialize(): Duplicated entry ("<<layer<<","<<cell<<") in trigger cell mapping. Use first mapping found.\n";
        }
    }
    ifs.close();


}

/*****************************************************************/
void TriggerHitFactory::connectVariables(TChain* chain)
/*****************************************************************/
{
}


/*****************************************************************/
void TriggerHitFactory::callback(void* object)
/*****************************************************************/
{
    TriggerHitFactory* myself = reinterpret_cast<TriggerHitFactory*>(object);
    myself->update();
}

/*****************************************************************/
unsigned TriggerHitFactory::cellIndex(int zside, int layer, int sector, int cell) const
/*****************************************************************/
{
    int z = (zside==-1 ? 0 : 1);
    int l = layer-1;
    int s = sector-1;
    //int ss = (subsector==-1 ? 0 : 1);
    int c = cell;
    //return (unsigned)(c + ss*3600 + s*7200 + l*129600 + z*3888000);
    return (unsigned)(c + s*3600 + l*64800 + z*1944000);
}

/*****************************************************************/
void TriggerHitFactory::update()
/*****************************************************************/
{
    m_data[0].clear();
    m_data[1].clear();
    m_data[2].clear();
    
    for(int iCol=0;iCol<3;iCol++)
    {
        const vector<SimHit>& hits = m_hitFactory->data().at(iCol);
        if(iCol==0)
        {
            map<int, pair<short,SimHit> > m_triggerHits;
            for(const auto& hit : hits)
            {
                short layer = hit.layer();
                short cell = hit.cell()*hit.subsector();
                short triggerCell = m_cellToTriggerCell[layer].at(cell);
                int index = cellIndex(hit.zside(), hit.layer(), hit.sector(), triggerCell);
                SimHit triggerHit(hit);
                triggerHit.setCell(triggerCell);
                pair<short,SimHit> p(1,triggerHit);
                auto ret = m_triggerHits.insert( pair<int, pair<short,SimHit> >( index, p ));
                if(!ret.second) 
                {
                    ret.first->second.first += 1;
                    ret.first->second.second.setEnergy(hit.energy() + ret.first->second.second.energy());
                    ret.first->second.second.setEta(hit.eta() + ret.first->second.second.eta());
                    double phidiff = TVector2::Phi_mpi_pi(hit.phi() - ret.first->second.second.phi());
                    ret.first->second.second.setPhi(2*ret.first->second.second.phi() + phidiff);
                }
            }
            for(auto& id_n_hit : m_triggerHits)
            {
                id_n_hit.second.second.setEta(id_n_hit.second.second.eta()/(double)id_n_hit.second.first);
                id_n_hit.second.second.setPhi( TVector2::Phi_mpi_pi(id_n_hit.second.second.phi()/(double)id_n_hit.second.first) );
                m_data[iCol].push_back( id_n_hit.second.second );
            }
        }
        else
        {
            for(const auto& hit : hits)
            {
                m_data[iCol].push_back( hit );
            }
        }
    }
}
