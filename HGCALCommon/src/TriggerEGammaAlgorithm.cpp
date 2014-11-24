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

#include "TEnv.h"

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
void TriggerEGammaAlgorithm::initialize(const EventHGCAL& event, TEnv& params)
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
    string pileupParamsFile = params.GetValue("PileupParams", "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/src/AnHiMaHGCAL/ElectronClusterThreshold/data/thresholdParameters.txt");
    ifstream stream(pileupParamsFile);
    if(!stream.is_open())
    {
        cerr<<"ERROR: Cannot open pileup parameters file '"<<pileupParamsFile<<"'\n";
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

    // fill cluster sizes
    string clusterSizeFile = params.GetValue("ClusterSize", "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/src/AnHiMaHGCAL/HGCALCommon/data/clusterSizes.txt");
    stream.open(clusterSizeFile);
    if(!stream.is_open())
    {
        cerr<<"ERROR: Cannot open cluster size file '"<<clusterSizeFile<<"'\n";
        return;
    }
    while(!stream.eof())
    {
        int layer;
        float size;
        stream >> layer >> size;
        m_clusterSizes[layer] = size;
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
float TriggerEGammaAlgorithm::pileupThreshold(float eta, int layer, int nhits)
/*****************************************************************/
{
    int etaIndex = 0;
    if     (fabs(eta)>=1.8 && fabs(eta)<2.)  etaIndex = 0;
    else if(fabs(eta)>=2.  && fabs(eta)<2.2) etaIndex = 1;
    else if(fabs(eta)>=2.2 && fabs(eta)<2.4) etaIndex = 2;
    else if(fabs(eta)>=2.4 && fabs(eta)<2.6) etaIndex = 3;
    else if(fabs(eta)>=2.6 && fabs(eta)<2.8) etaIndex = 4;
    else if(fabs(eta)>=2.8 && fabs(eta)<3.)  etaIndex = 5;
    else if(fabs(eta)>=1.5 && fabs(eta)<1.6) etaIndex = 6;
    else if(fabs(eta)>=1.6 && fabs(eta)<1.7) etaIndex = 7;
    else if(fabs(eta)>=1.7 && fabs(eta)<1.8) etaIndex = 8;
    const auto itr = m_pileupParams.find(make_pair(etaIndex,layer));
    if(itr==m_pileupParams.end())
    {
        cerr<<"ERROR: cannot find pileup parameters for eta="<<eta<<",layer="<<layer<<"\n";
        return 1.;
    }
    double a = itr->second.first;
    double b = itr->second.second;
    return a*(double)nhits + b;
}

/*****************************************************************/
int TriggerEGammaAlgorithm::triggerRegionIndex(float eta, int zside, int sector, int subsector)
/*****************************************************************/
{
    int iz = (zside==-1 ? 0 : 1);
    int subsec = (subsector==-1 ? 1 : 0);
    int sector30deg = (2*(sector-1) + subsec)/3;
    int up = (abs(eta)<=1.8 ? 1 : 0);
    int triggerRegion = sector30deg + 12*up + 24*iz;
    return triggerRegion;
}

/*****************************************************************/
int TriggerEGammaAlgorithm::triggerRegionHits(int triggerRegion, int layer)
/*****************************************************************/
{
    const auto itr = m_regionHitsAboveTh.find(make_pair(triggerRegion, layer));
    if(itr==m_regionHitsAboveTh.end())
    {
        cerr<<"ERROR: cannot find number of hits in region "<<triggerRegion<<" layer "<<layer<<"\n";
        return 0;
    }
    return itr->second.size();
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
                cerr<<"[WARNING] Cannot find cell in layer "<<l<<"\n";
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
    } // end seeding hits loop

}

/*****************************************************************/
void TriggerEGammaAlgorithm::clustering(const EventHGCAL& event, const vector<Tower>& seeds, vector<Tower>& clusters)
/*****************************************************************/
{
    //cerr<<"Calling clustering\n";
    // count number of hits in each trigger region, needed for pileup threshold
    fillPileupEstimators(event);

    // find hits above pileup threshold
    map<int, vector<const SimHit*> > hitsAboveThreshold;
    for(int l=1;l<=30;l++)
    {
        hitsAboveThreshold[l] = vector<const SimHit*>();
    }
    for(const auto& hit : event.simhits())
    {
        int triggerRegion = triggerRegionIndex(hit.eta(), hit.zside(), hit.sector(), hit.subsector());
        int nhits = triggerRegionHits(triggerRegion, hit.layer());
        float threshold = pileupThreshold(hit.eta(), hit.layer(), nhits);
        //float threshold = 5.;
        if(hit.energy()>=threshold*mip) hitsAboveThreshold[hit.layer()].push_back(&hit);
    }


    // to calibrate clusters
    TowerCalibrator towerCalibrator;
    // set of hits already included in clusters
    set<const SimHit*> usedHits;

    // sort seeds before looping (to favor highest ET seeds for hit aggregation)
    vector<const Tower*> sortedSeeds;
    for(const auto& seed : seeds) sortedSeeds.push_back(&seed);
    sort(sortedSeeds.begin(), sortedSeeds.end(), AnHiMa::towerSort);
    
    // Build clusters around each seed
    for(const auto seed : sortedSeeds)
    {
        //cerr<<"seed ET="<<seed->calibratedEt()<<", eta="<<seed->eta()<<", phi="<<seed->phi()<<"\n";
        const SimHit* seedHit = seed->hits()[0];

        // first build a 1x1 tower
        int layer0 = seedHit->layer();
        double eta0 = seedHit->eta();
        double phi0 = seedHit->phi();

        vector<HGCEEDetId> idLayers(31);
        idLayers[layer0] = HGCEEDetId(seedHit->detid());
        for(int l=layer0+1;l<=30;l++)
        {
            vector<HGCEEDetId> ids;
            ids = event.hgcalNavigator().upProj(idLayers[l-1], 1, eta0, phi0);
            if(ids.size()==0)
            {
                cerr<<"[WARNING] Cannot find cell in layer "<<l<<"\n";
                idLayers[l] = (HGCEEDetId)DetId(0);
            }
            else
            {
                idLayers[l] = ids[0];
            }
        }
        for(int l=layer0-1;l>0;l--)
        {
            vector<HGCEEDetId> ids;
            ids = event.hgcalNavigator().downProj(idLayers[l+1], 1, eta0, phi0);
            if(ids.size()==0)
            {
                cerr<<"[WARNING] Cannot find cell in layer "<<l<<"\n";
                idLayers[l] = (HGCEEDetId)DetId(0);
            }
            else
            {
                idLayers[l] = ids[0];
            }
        }

        // then aggregate hits in each layer
        Tower cluster;
        for(const auto& id : idLayers)
        {
            if(id==HGCEEDetId(0)) continue;

            const GlobalPoint& point0 = event.hgcalNavigator().geometry()->getPosition(id);
            float x0 = point0.x();
            float y0 = point0.y();
            float z0 = point0.z();
            int layer = id.layer();
            float clusterSize = m_clusterSizes[layer];
            for(const auto hit: hitsAboveThreshold[layer])
            {
                if(usedHits.find(hit)!=usedHits.end()) continue;
                if(z0*hit->z()<0.) continue;
                float x = hit->x();
                float y = hit->y();
                double dx = (x-x0);
                double dy = (y-y0);
                float dr = sqrt(dx*dx+dy*dy);
                if(dr>clusterSize) continue;
                cluster.addHit(*hit);
                usedHits.insert(hit);
                //cerr<<"    Adding hit eta="<<hit->eta()<<", phi="<<hit->phi()<<"\n";
            }
        }
        if(cluster.energy()>0.)
        {
            // calibrate energy
            towerCalibrator.calibrate(cluster);
            clusters.push_back(cluster);
            //cerr<<"cluster ET="<<cluster.calibratedEt()<<", eta="<<cluster.eta()<<", phi="<<cluster.phi()<<"\n";
        }
        //cerr<<"cluster ET="<<cluster.calibratedEt()<<", eta="<<cluster.eta()<<", phi="<<cluster.phi()<<"\n";
        //cerr<<" NHits="<<cluster.nHits()<<"\n";
        //cerr<<" Layer energies = ";
        //for(int l=1;l<=30;l++) cerr<<cluster.layerEnergy(l)<<",";
        //cerr<<"\n";
    } // end seeds loop
    
}


/*****************************************************************/
void TriggerEGammaAlgorithm::superClustering(const vector<Tower>& clusters, vector<SuperCluster>& superClusters)
/*****************************************************************/
{
    //cerr<<"In superClustering()\n";
    vector<SuperCluster*> sortedSuperClusters;
    // Find seed clusters
    // Seed clusters are clusters with max Et in a 0.05x0.15 eta x phi window
    for(const auto& cluster : clusters)
    {
        double eta0 = cluster.eta();
        double phi0 = cluster.phi();
        double et0 = cluster.et();
        bool isSeedCluster = true;
        // find neighbour clusters
        for(const auto& subcluster : clusters)
        {
            double eta = subcluster.eta();
            double phi = subcluster.phi();
            double et = subcluster.et();
            double deta = (eta-eta0);
            double dphi = TVector2::Phi_mpi_pi(phi-phi0);
            if(fabs(deta)<1.e-5 && fabs(dphi)<1.e-5) continue;
            if(fabs(deta)>0.05 || fabs(dphi)>0.15) continue;
            if(et>et0) 
            {
                isSeedCluster = false;
                break;
            }
        }
        if(isSeedCluster) 
        {
            superClusters.push_back( SuperCluster() );
            superClusters.back().addCluster(&cluster);
            //cerr<<"seed cluster eta="<<cluster.eta()<<", phi="<<cluster.phi()<<", Et="<<cluster.calibratedEt()<<"\n";
        }
    }
    for(auto& sc : superClusters)
    {
        //cerr<<"seed cluster eta="<<sc.eta()<<", phi="<<sc.phi()<<", Et="<<sc.et()<<"\n";
        sortedSuperClusters.push_back(&sc);
    }
    sort(sortedSuperClusters.begin(), sortedSuperClusters.end(), AnHiMa::superClusterSort);


    std::set<const Tower*> usedClusters;
    // Gather subclusters
    for(auto sc : sortedSuperClusters)
    {
        double eta0 = sc->eta();
        double phi0 = sc->phi();
        //cerr<<"seed cluster eta="<<eta0<<", phi="<<phi0<<", Et="<<sc->et()<<"\n";
        for(const auto& subcluster : clusters)
        {
            if(usedClusters.find(&subcluster)!=usedClusters.end()) continue;
            double eta = subcluster.eta();
            double phi = subcluster.phi();
            double deta = (eta-eta0);
            double dphi = TVector2::Phi_mpi_pi(phi-phi0);
            // discard seed cluster
            if(fabs(deta)<1.e-5 && fabs(dphi)<1.e-5) continue;
            // outside window
            if(fabs(deta)>0.05 || fabs(dphi)>0.15) continue;
            //
            //cerr<<"   Add subcluster eta="<<eta<<", phi="<<phi<<", Et="<<subcluster.calibratedEt()<<"\n";
            sc->addCluster(&subcluster);
            usedClusters.insert(&subcluster);
        }
    }
    //cerr<<"End superClustering()\n";
}
