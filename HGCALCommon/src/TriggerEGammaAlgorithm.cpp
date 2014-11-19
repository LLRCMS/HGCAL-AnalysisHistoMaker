/**
 *  @file  TriggerEGammaAlgorithm.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    17/11/2014
 *
 *  @internal
 *     Created :  17/11/2014
 * Last update :  17/11/2014 18:51:17
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include <fstream>

#include "AnHiMaHGCAL/HGCALCommon/interface/TriggerEGammaAlgorithm.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"

#include "AnHiMaHGCAL/HGCALCommon/interface/TowerCalibrator.h"

using namespace std;
using namespace AnHiMa;


/*****************************************************************/
TriggerEGammaAlgorithm::TriggerEGammaAlgorithm()
/*****************************************************************/
{
}



/*****************************************************************/
TriggerEGammaAlgorithm::~TriggerEGammaAlgorithm()
/*****************************************************************/
{
}



/*****************************************************************/
void TriggerEGammaAlgorithm::initialize(const EventHGCAL& event, const std::string& pileupParamsFile)
/*****************************************************************/
{
    // initialize vectors
    for(int s=0;s<12;s++)
    {
        for(int z=0;z<2;z++)
        {
            for(int u=0;u<2;u++)
            {
                int triggerRegion = s + 12*u + 24*z;
                for(int l=1;l<=30;l++)
                {
                    m_regionSimHits[make_pair(triggerRegion, l)] = vector<pair<HGCEEDetId,float> >(0);
                    m_regionHitsAboveTh[make_pair(triggerRegion, l)] = vector<const SimHit*>(0);
                }
            }
        }
    }
    // fill vectors
    for (int izz=0; izz<=1; izz++) 
    {
        int iz = (2*izz-1);
        for (int subsec=0; subsec<=1; ++subsec) 
        {
            int ssec = (2*subsec-1);
            for (int sec=1; sec<=18; ++sec) 
            {
                int sector30deg = (2*(sec-1) + (subsec==0 ? 1 : 0))/3;
                for (int lay=1; lay<=30; ++lay) 
                {
                    for (int cell=0; cell<4000; ++cell) 
                    {
                        const HGCEEDetId id(HGCEE,iz,lay,sec,ssec,cell);
                        if (event.hgcalNavigator().valid(id)) 
                        {
                            double eta = event.hgcalNavigator().geometry()->getPosition(id).eta();
                            //int up = (cell>=1200 ? 1 : 0);
                            int up = (fabs(eta)<=1.8 ? 1 : 0);
                            int triggerRegion = sector30deg + 12*up + 24*izz;
                            m_regionSimHits[make_pair(triggerRegion, lay)].push_back(make_pair(id,eta));
                        }
                    }
                }
            }
        }
    }

    // fill pileup params
    ifstream stream(pileupParamsFile);
    if(!stream.is_open())
    {
        cout<<"ERROR: Cannot open pileup parameters file '"<<pileupParamsFile<<"'\n";
        return;
    }
    while(!stream.eof())
    {
        int up, layer, subeta;
        float a, b;
        stream >> up >> layer >> subeta >> a >> b;
        int eta = subeta + 6*up;
        m_pileupParams[make_pair(eta,layer)] = make_pair(a,b);
    }
    stream.close();
}




/*****************************************************************/
void TriggerEGammaAlgorithm::fillPileupEstimators(const EventHGCAL& event)
/*****************************************************************/
{
    for(auto& region_hits : m_regionHitsAboveTh)
    {
        region_hits.second.clear();
    }

    for(const auto& region_id : m_regionSimHits)
    {
        pair<int,int> region_layer = region_id.first;
        for(const auto& id_eta : region_id.second)
        {
            const HGCEEDetId& id = id_eta.first;
            const SimHit* hit = event.simhit(id);
            if(hit && hit->energy()>mip) m_regionHitsAboveTh[region_layer].push_back(hit);
        }
    }
}


/*****************************************************************/
void TriggerEGammaAlgorithm::seeding(const EventHGCAL& event, vector<Tower>& seeds)
/*****************************************************************/
{
    // find hits in layers 15 > 5 mip
    //vector<const SimHit*> selectedHits;
    //for(const auto& hit : event.simhits())
    //{
    //    if(hit.energy()>5.*mip && hit.layer()==15) selectedHits.push_back(&hit);
    //}
    TowerCalibrator towerCalibrator;
    set<const SimHit*> usedHits;

    // build towers
    for(auto hit : event.sortedSeedingHits())
    {
        int layer = hit->layer();
        double refEta = hit->eta();
        double refPhi = hit->phi();

        // find projective cells in layers 15-18
        vector<HGCEEDetId> idLayers(31);
        idLayers[layer] = HGCEEDetId(hit->detid());
        for(int l=layer+1;l<=18;l++)
        {
            vector<HGCEEDetId> ids;
            ids = event.hgcalNavigator().upProj(idLayers[l-1], 1, refEta, refPhi);
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

        // find the 4x4 cells
        vector<HGCEEDetId> towerCells;
        for(int l=layer;l<=18;l++)
        {
            // center cell
            towerCells.push_back(idLayers[l]);
            // first ring (3x3)
            DetId id1 = event.hgcalNavigator().topology()->offsetBy(idLayers[l], -1,  1);
            DetId id2 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  0,  1);
            DetId id3 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  1,  1);
            DetId id4 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  1,  0);
            DetId id5 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  1, -1);
            DetId id6 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  0, -1);
            DetId id7 = event.hgcalNavigator().topology()->offsetBy(idLayers[l], -1, -1);
            DetId id8 = event.hgcalNavigator().topology()->offsetBy(idLayers[l], -1,  0);
            // finish 4x4
            DetId id9  = event.hgcalNavigator().topology()->offsetBy(idLayers[l], -1,  2);
            DetId id10 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  0,  2);
            DetId id11 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  1,  2);
            DetId id12 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  2,  2);
            DetId id13 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  2,  1);
            DetId id14 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  2,  0);
            DetId id15 = event.hgcalNavigator().topology()->offsetBy(idLayers[l],  2, -1);
            // Add cells
            towerCells.push_back((HGCEEDetId)id1);
            towerCells.push_back((HGCEEDetId)id2);
            towerCells.push_back((HGCEEDetId)id3);
            towerCells.push_back((HGCEEDetId)id4);
            towerCells.push_back((HGCEEDetId)id5);
            towerCells.push_back((HGCEEDetId)id6);
            towerCells.push_back((HGCEEDetId)id7);
            towerCells.push_back((HGCEEDetId)id8);
            towerCells.push_back((HGCEEDetId)id9);
            towerCells.push_back((HGCEEDetId)id10);
            towerCells.push_back((HGCEEDetId)id11);
            towerCells.push_back((HGCEEDetId)id12);
            towerCells.push_back((HGCEEDetId)id13);
            towerCells.push_back((HGCEEDetId)id14);
            towerCells.push_back((HGCEEDetId)id15);
        }
        // build towers from hits
        Tower tower;
        for(const auto& id : towerCells)
        {
            if(id==DetId(0)) continue;
            const SimHit* hit = event.simhit(id);
            if(!hit || hit->energy()<=mip) continue;
            if(usedHits.find(hit)!=usedHits.end()) continue;
            tower.addHit(*hit);
            usedHits.insert(hit);
        }
        // calibrate energy
        towerCalibrator.calibrate(tower);
        // filter seed candidates
        if(tower.calibratedEt()<0.5) continue;
        int nhits = tower.layerNHits(15)+tower.layerNHits(16)+tower.layerNHits(17)+tower.layerNHits(18);
        if(nhits<19) continue;
        seeds.push_back(tower);
    }
}
