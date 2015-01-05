/**
 *  @file  ProjectedSimHitFactory.cpp
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





#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/ProjectedSimHitFactory.h"

#include <iostream>
#include <fstream>

using namespace AnHiMa;
using namespace std;


/*****************************************************************/
void ProjectedSimHitFactory::initialize(IEvent* event, TChain* chain)
/*****************************************************************/
{
    event->registerCallback((void*)this, ProjectedSimHitFactory::callback);
    connectVariables(chain);

}

/*****************************************************************/
void ProjectedSimHitFactory::initialize(IEvent* event, TChain* chain, const EventHGCAL* eventHGCAL, const string& projectionMappingFile, HitType type)
/*****************************************************************/
{
    m_type = type;
    m_event = eventHGCAL;
    fillProjectionMapping(projectionMappingFile);
    initialize(event, chain);
}



/*****************************************************************/
void ProjectedSimHitFactory::connectVariables(TChain* chain)
/*****************************************************************/
{
// dummy
}



/*****************************************************************/
void ProjectedSimHitFactory::callback(void* object)
/*****************************************************************/
{
    ProjectedSimHitFactory* myself = reinterpret_cast<ProjectedSimHitFactory*>(object);
    myself->update();
}


/*****************************************************************/
void ProjectedSimHitFactory::update()
/*****************************************************************/
{
    map<unsigned, pair<unsigned,SimHit*> > hitMap;
    m_data.clear();
    const vector<SimHit>* simhits;
    switch(m_type)
    {
        case ALL:
            simhits = &(m_event->simhits());
            break;
        case HARDINT:
            simhits = &(m_event->hardsimhits());
            break;
        default:
            simhits = &(m_event->simhits());
            break;
    }
        cerr<<"toto\n";
    for(const auto& hit : *simhits)
    {
        HGCEEDetId projId = projectCell( (HGCEEDetId)hit.detid() );
        cerr<<"hit "<<projId.rawId()<<"\n";
        if(projId==HGCEEDetId(0))
        {
            cerr<<"WARNING: Invalid projected cell\n";
        }
        auto itr = hitMap.find(projId.rawId());
        if(itr==hitMap.end())
        {
            double calibEnergy = m_calibrator.calibratedEnergy(hit.energy(), hit.eta(), hit.layer());
            const GlobalPoint& projPoint = m_event->hgcalNavigator().geometry()->getPosition(projId);
            SimHit hitProj;
            hitProj.setDetid     ( projId.rawId() );
            hitProj.setSubdet    ( projId.subdet() );
            hitProj.setCell      ( projId.cell() );
            hitProj.setSector    ( projId.sector() );
            hitProj.setSubsector ( projId.subsector() );
            hitProj.setLayer     ( 1 );
            hitProj.setZside     ( projId.zside() );
            hitProj.setEnergy    ( calibEnergy );
            hitProj.setEta       ( projPoint.eta() );
            hitProj.setPhi       ( projPoint.phi() );
            hitProj.setX         ( projPoint.x() );
            hitProj.setY         ( projPoint.y() );
            hitProj.setZ         ( projPoint.z() );

            m_data.push_back( hitProj );
            hitMap[projId.rawId()] = make_pair( m_data.size()-1, &(*(m_data.end()--)) );
        }
        else
        {
            double calibEnergy = m_calibrator.calibratedEnergy(hit.energy(), hit.eta(), hit.layer());
            itr->second.second->setEnergy(itr->second.second->energy() + calibEnergy);
        }

    }
}

/*****************************************************************/
HGCEEDetId ProjectedSimHitFactory::projectCell(const HGCEEDetId& cellid)
/*****************************************************************/
{
    int layer = cellid.layer();
    int cell  = cellid.cell();
    std::pair<unsigned,unsigned> layer_cell = make_pair(layer,cell);
    const auto& itr = m_projectionMapping.find(layer_cell);
    if(itr==m_projectionMapping.end())
    {
        cerr<<"WARNING: Cannot find projected cell for (layer,cell)=("<<layer<<","<<cell<<")\n";
        return DetId(0);
    }
    return HGCEEDetId(HGCEE, cellid.zside(), 1, cellid.sector(), cellid.subsector(), itr->second);
}


/*****************************************************************/
void ProjectedSimHitFactory::fillProjectionMapping(const string& inputFile)
/*****************************************************************/
{
    fstream input(inputFile,fstream::in);
    if(!input.is_open())
    {
        cerr<<"ERROR: Cannot open projection mapping file '"<<inputFile<<"'\n";
        return;
    }
    while(!input.eof())
    {
        int layer, cell, projcell;
        input >> layer >> cell >> projcell;
        m_projectionMapping[make_pair(layer,cell)] = projcell;
    }
    input.close();

}
