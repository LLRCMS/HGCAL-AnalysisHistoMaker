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
#include <cassert>

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
void TriggerHitFactory::initialize(IEvent* event, TChain* chain, const SimHitFactory& hitFactory, const string& triggerCellMapping, const HGCALNavigator& navigator)
/*****************************************************************/
{
    m_hitFactory = &hitFactory;
    initialize(event, chain);
    m_data.push_back( vector<SimHit>() ); // HE
    m_data.push_back( vector<SimHit>() ); // FH
    m_data.push_back( vector<SimHit>() ); // BH


    array< map<short, vector<short>>, 30> triggerCellToCell; 

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
    
    for(; ifs>>layer>>cell>>triggerCell>>subsector; )
    {
        if(layer>=30 || layer<0) 
        {
            cerr<<"ERROR: TriggerHitFactory::initialize(): Found layer '"<<layer<<"' in trigger cell mapping. Discard it.\n";
            continue;
        }
        short newcell = (subsector==1 ? cell+1 : -(cell+1));
        auto ret = m_cellToTriggerCell[layer].insert( std::make_pair(newcell, triggerCell) );
        if(ret.second==false)
        {
            cerr<<"ERROR: TriggerHitFactory::initialize(): Duplicated entry ("<<layer<<","<<newcell<<") in trigger cell mapping. Use first mapping found.\n";
        }
        vector<short> vec(0);
        triggerCellToCell[layer].insert( std::make_pair(triggerCell, vec) );
        triggerCellToCell[layer].at(triggerCell).push_back(newcell);
    } 
    assert(ifs.eof());
    ifs.close();

    // check the number of cells included in trigger cells
    for(unsigned layer=1; layer<=30; layer++)
    {
        const auto& mapping = triggerCellToCell.at(layer-1);
        for(const auto& trigger_cells : mapping)
        {
            if(trigger_cells.second.size()>5) 
            {
                cout<<"WARNING: TriggerHitFactory::initialize(): layer "<<layer<<" trigger cell "<<trigger_cells.first<<" includes "<<trigger_cells.second.size()<<" cells\n";
            }
        }
    }

    // compute positions of trigger cells
    for(unsigned sec=1;sec<=18; sec++)
    {
        for(unsigned layer=1; layer<=30; layer++)
        {
            for(const auto& trigger_cells : triggerCellToCell[layer-1])
            {
                const short triggerCell = trigger_cells.first;
                const vector<short>& cells = trigger_cells.second;
                SimplePosition pos;
                pos.x   = 0.;
                pos.y   = 0.;
                for(const auto& c : cells)
                {
                    int subsec = (c>0 ? 1 : -1);
                    int cell = abs(c)-1;
                    HGCEEDetId detid(HGCEE, 1, layer, sec, subsec, cell);
                    const GlobalPoint point = navigator.geometry()->getPosition(detid);
                    //pos.eta += point.eta();
                    //pos.phi += (pos.phi + TVector2::Phi_mpi_pi(point.phi() - pos.phi));
                    pos.x   += point.x();
                    pos.y   += point.y();
                }
                //pos.eta /= (float)cells.size();
                //pos.phi  = TVector2::Phi_mpi_pi(pos.phi/(float)cells.size());
                pos.x   /= (float)cells.size();
                pos.y   /= (float)cells.size();
                m_triggerCellPositions[sec-1][layer-1].insert( std::make_pair(triggerCell, pos) );
            }
        }
    }

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
    
    // Loop on EE, F-HE, B-HE
    // Trigger cell grouping is only done in EE
    for(int iCol=0;iCol<3;iCol++)
    {
        const vector<SimHit>& hits = m_hitFactory->data().at(iCol);
        if(iCol==0) // EE
        {
            map<int, SimHit> m_triggerHits;
            for(const auto& hit : hits)
            {
                try
                {
                    short layer = hit.layer()-1;
                    short cell = (hit.cell()+1)*hit.subsector();
                    short triggerCell = m_cellToTriggerCell[layer].at(cell);
                    int index = cellIndex(hit.zside(), hit.layer(), hit.sector(), triggerCell);
                    SimHit triggerHit(hit);
                    triggerHit.setCell(triggerCell);
                    auto ret = m_triggerHits.insert( pair<int, SimHit>( index, triggerHit ) );
                    //cerr<<"Hit "<<layer<<","<<cell<<" -> "<<triggerCell<<"\n";
                    if(!ret.second) // if trigger cells already store, add energy
                    {
                        ret.first->second.setEnergy(hit.energy() + ret.first->second.energy()); // increase trigger cell energy
                    }
                }
                catch(std::out_of_range& e)
                {
                    cerr<<"ERROR: Cannot find the trigger cell for hit layer "<<hit.layer()<<" sector "<<hit.sector()<<" subsec "<<hit.subsector()<<" cell "<<hit.cell()<<"\n";
                    cerr<<"       This hit is discarded\n";
                }
            }
            // retrieve trigger cell position and fill in data buffer
            for(auto& id_hit : m_triggerHits) 
            {
                const SimHit& hit = id_hit.second;
                try
                {
                    const SimplePosition& pos = m_triggerCellPositions[hit.sector()-1][hit.layer()-1].at(hit.cell());
                    float x = (float)hit.zside() * pos.x;
                    float y = pos.y;
                    float z = hit.z();
                    GlobalPoint point(x, y, z);
                    float eta = point.eta();
                    float phi = point.phi();
                    id_hit.second.setEta(eta);
                    id_hit.second.setPhi(phi);
                    id_hit.second.setX(x);
                    id_hit.second.setY(y);
                }
                catch(std::out_of_range& e)
                {
                    cerr<<"ERROR: Cannot find the position of trigger hit layer "<<hit.layer()<<" sector "<<hit.sector()<<" cell "<<hit.cell()<<"\n";
                    cerr<<"       The hit position will be incorrect\n";
                }
                m_data[iCol].push_back( id_hit.second );
            }
        }
        // not EE. Just copy the cells (no grouping)
        else
        {
            for(const auto& hit : hits)
            {
                m_data[iCol].push_back( hit );
            }
        }
    }
}
