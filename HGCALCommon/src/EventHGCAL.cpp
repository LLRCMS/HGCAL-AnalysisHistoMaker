/**
 *  @file  EventHGCAL.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    05/29/2014
 *
 *  @internal
 *     Created :  05/29/2014
 * Last update :  05/29/2014 05:17:36 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include <iostream>
#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"
#include <TChain.h>
#include "TMath.h"


using namespace AnHiMa;
using namespace std;

/*****************************************************************/
EventHGCAL::EventHGCAL():IEvent(),
    m_sortedHits(5184000)
/*****************************************************************/
{

    m_hgcalNavigator.initialize();

}



/*****************************************************************/
EventHGCAL::~EventHGCAL()
/*****************************************************************/
{
}


/*****************************************************************/
void EventHGCAL::connectVariables(TChain* inputChain)
/*****************************************************************/
{
    IEvent::connectVariables(inputChain);

    // Branches to read
    inputChain->SetBranchStatus("run"       , true);
    inputChain->SetBranchStatus("lumi"      , true);
    inputChain->SetBranchStatus("event"     , true);
    inputChain->SetBranchStatus("npu"       , true);


    // Connect branches
    inputChain->SetBranchAddress("run"       , &m_run);
    inputChain->SetBranchAddress("lumi"      , &m_lumi);
    inputChain->SetBranchAddress("event"     , &m_event);
    inputChain->SetBranchAddress("npu"       , &m_npu);


    // keep this order of initialization
    m_simhitFactoryAll.initialize(this, inputChain, SimHitFactory::ALL);
    m_simhitFactoryHard.initialize(this, inputChain, SimHitFactory::HARDINT);
    registerCallback((void*)this, EventHGCAL::callback);
    //m_towerFactory.initialize(this, inputChain, this);
    m_genparticleFactory.initialize(this, inputChain);

}


/*****************************************************************/
bool EventHGCAL::passSelection(int selection)
/*****************************************************************/
{
    return true;
}


/*****************************************************************/
void EventHGCAL::update()
/*****************************************************************/
{
    IEvent::update();

}

/*****************************************************************/
void EventHGCAL::callback(void* object)
/*****************************************************************/
{
    EventHGCAL* myself = reinterpret_cast<EventHGCAL*>(object);
    myself->buildHitsIndex();
}


/*****************************************************************/
void EventHGCAL::buildHitsIndex()
/*****************************************************************/
{
    // reset simhit pointers
    for(unsigned i=0; i<m_sortedHits.size(); i++)
    {
        //cerr<<"Resetting "<<i<<"\n";
        m_sortedHits[i] = NULL;
    }
    // sort hits according to their index
    for(auto itr=simhits().cbegin(); itr!=simhits().end(); itr++)
    {
        const SimHit& hit = *itr;
        unsigned index = cellIndex(hit.zside(), hit.layer(), hit.sector(), hit.subsector(), hit.cell());
        //cerr<<"Setting "<<index<<"\n";
        //cerr<<" ("<<hit.zside()<<","<< hit.layer()<<","<< hit.sector()<<","<< hit.subsector()<<","<< hit.cell()<<")\n";
        m_sortedHits[index] = &hit;
    }


    m_energySortedHits.clear();
    m_energySortedHits.reserve(simhits().size());
    for(auto itr=simhits().cbegin(); itr!=simhits().end(); itr++)
    {
        m_energySortedHits.push_back( &(*itr) );
    }
    sort(m_energySortedHits.begin(), m_energySortedHits.end(), AnHiMa::simHitSort);
}


/*****************************************************************/
unsigned EventHGCAL::cellIndex(int zside, int layer, int sector, int subsector, int cell) const
/*****************************************************************/
{
    int z = (zside==-1 ? 0 : 1);
    int l = layer-1;
    int s = sector-1;
    int ss = (subsector==-1 ? 0 : 1);
    int c = cell;
    return (unsigned)(c + ss*2400 + s*4800 + l*86400 + z*2592000);
}


/*****************************************************************/
unsigned EventHGCAL::cellIndex(const HGCEEDetId& detid) const
/*****************************************************************/
{
    return cellIndex(detid.zside(), detid.layer(), detid.sector(), detid.subsector(), detid.cell());
}


/*****************************************************************/
const SimHit* EventHGCAL::simhit(int zside, int layer, int sector, int subsector, int cell) const
/*****************************************************************/
{
    unsigned index = cellIndex(zside, layer, sector, subsector, cell);
    return m_sortedHits[index];
}


/*****************************************************************/
const SimHit* EventHGCAL::simhit(const HGCEEDetId& detid) const
/*****************************************************************/
{
    unsigned index = cellIndex(detid);
    return m_sortedHits[index];
}

