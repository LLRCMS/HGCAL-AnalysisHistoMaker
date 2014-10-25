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


    m_simhitFactory.initialize(this, inputChain);

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
}


/*****************************************************************/
unsigned EventHGCAL::cellIndex(int zside, int layer, int sector, int subsector, int cell)
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
unsigned EventHGCAL::cellIndex(const HGCEEDetId& detid)
/*****************************************************************/
{
    return cellIndex(detid.zside(), detid.layer(), detid.sector(), detid.subsector(), detid.cell());
}

