/**
 *  @file  TowerFactory.cpp
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





#include "AnHiMaHGCAL/HGCALCommon/interface/TowerFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"

using namespace AnHiMa;
using namespace std;



/*****************************************************************/
void TowerFactory::initialize(IEvent* event, TChain* chain)
/*****************************************************************/
{
    event->registerCallback((void*)this, TowerFactory::callback);
    connectVariables(chain);

}

/*****************************************************************/
void TowerFactory::initialize(IEvent* event, TChain* chain, EventHGCAL* hgcalEvent)
/*****************************************************************/
{
    m_hgcalEvent = hgcalEvent;
    initialize(event, chain);
}



/*****************************************************************/
void TowerFactory::connectVariables(TChain* chain)
/*****************************************************************/
{


}



/*****************************************************************/
void TowerFactory::callback(void* object)
/*****************************************************************/
{
    TowerFactory* myself = reinterpret_cast<TowerFactory*>(object);
    myself->update();
}


/*****************************************************************/
void TowerFactory::update()
/*****************************************************************/
{
    m_data.clear();
    float mip = 0.000055;

    vector<HGCEEDetId> idLayers(31);
    int lay = 15;
    int iz = 1;
    int subsec = 0;
    int sec = 1;
    //for (int izz=0; izz<=1; izz++) 
    //{
    //    int iz = (2*izz-1);
    //    for (int subsec=0; subsec<=1; ++subsec) 
    //    {
    //        for (int sec=1; sec<=18; ++sec) 
    //        {
                for (int cell=0; cell<2174; ++cell) 
                {
                    const HGCEEDetId detid(HGCEE,iz,lay,sec,subsec,cell);
                    if (m_hgcalEvent->hgcalNavigator().valid(detid)) 
                    {
                        //cout<<" cell z="<<iz<<", subsec="<<subsec<<", sec="<<sec<<", cell="<<cell<<"\n";
                        for(int l=1;l<=30;l++)
                        {
                            if(l==15) continue;
                            vector<HGCEEDetId> ids = m_hgcalEvent->hgcalNavigator().upProj(detid, l-15);
                            if(ids.size()==0)
                            {
                                cout<<"[WARNING] Cannot find cell in layer "<<l<<"\n";
                                idLayers[l] = (HGCEEDetId)DetId(0);
                            }
                            else
                            {
                                idLayers[l] = ids[0];
                            }
                        }

                        // build tower from the list of cells
                        Tower tower;
                        double etaMean = 0.;
                        double phiMean = 0.;
                        for(unsigned l=0; l<idLayers.size(); l++)
                        {
                            HGCEEDetId id = idLayers[l];
                            if(id==DetId(0)) continue;
                            const SimHit* hit = m_hgcalEvent->simhit(id);
                            if(!hit || hit->energy()<mip) continue;
                            //cout<<"   Adding cell eta="<<hit->eta()<<",phi="<<hit->phi()<<",layer="<<hit->layer()<<", E="<<hit->energy()<<"\n";
                            etaMean += hit->eta()*hit->energy();
                            phiMean += hit->phi()*hit->energy();
                            tower.addHit(*hit);
                        }
                        if(tower.energy()>0.)
                        {
                            tower.setEta(etaMean/tower.energy());
                            tower.setPhi(phiMean/tower.energy());
                            //cout<<"   Filling tower eta="<<tower.eta()<<", phi="<<tower.phi()<<", E="<<tower.energy()<<"\n";
                            m_data.push_back(tower);
                        }
                    }
                }
    //        }
    //    }
    //}

}
