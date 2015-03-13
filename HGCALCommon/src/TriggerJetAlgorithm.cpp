/**
 *  @file  TriggerJetAlgorithm.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    14/01/2015
 *
 *  @internal
 *     Created :  14/01/2015
 * Last update :  14/01/2015 13:59:18
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */





#include <fstream>
#include <list>

#include "TEnv.h"
#include "TFile.h"

#include "AnHiMaHGCAL/HGCALCommon/interface/TriggerJetAlgorithm.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"

#include "AnHiMaHGCAL/HGCALCommon/interface/TowerCalibrator.h"

using namespace std;
using namespace AnHiMa;


/*****************************************************************/
TriggerJetAlgorithm::TriggerJetAlgorithm():
    m_truthCorrection(0)
/*****************************************************************/
{
}



/*****************************************************************/
TriggerJetAlgorithm::~TriggerJetAlgorithm()
/*****************************************************************/
{
}



/*****************************************************************/
void TriggerJetAlgorithm::initialize(const EventHGCAL& event, TEnv& params)
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
                // ECAL
                for(int l=1;l<=30;l++)
                {
                    m_regionSimHitsECAL[make_pair(triggerRegion, l)] = vector<pair<HGCEEDetId,float> >(0);
                    m_regionHitsAboveThECAL[make_pair(triggerRegion, l)] = vector<const SimHit*>(0);
                }
                // HCAL
                for(int l=1;l<=12;l++)
                {
                    m_regionSimHitsHCAL[make_pair(triggerRegion, l)] = vector<pair<HGCHEDetId,float> >(0);
                    m_regionHitsAboveThHCAL[make_pair(triggerRegion, l)] = vector<const SimHit*>(0);
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
                // ECAL
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
                            m_regionSimHitsECAL[make_pair(triggerRegion, lay)].push_back(make_pair(id,eta));
                        }
                    }
                }
                // HCAL
                for (int lay=1; lay<=12; ++lay) 
                {
                    for (int cell=0; cell<4000; ++cell) 
                    {
                        const HGCHEDetId id(HGCHEF,iz,lay,sec,ssec,cell);
                        if (event.hgcalNavigator().valid(id,4)) 
                        {
                            double eta = event.hgcalNavigator().geometry(4)->getPosition(id).eta();
                            //int up = (cell>=1200 ? 1 : 0);
                            int up = (fabs(eta)<=1.8 ? 1 : 0);
                            int triggerRegion = sector30deg + 12*up + 24*izz;
                            m_regionSimHitsHCAL[make_pair(triggerRegion, lay)].push_back(make_pair(id,eta));
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
    cout<<"INFO: PileupParams = "<<pileupParamsFile<<"\n";
    while(!stream.eof())
    {
        int up, subdet, layer, subeta;
        float a, b;
        stream >> up >> subdet >> layer >> subeta >> a >> b;
        int eta = subeta + 6*up;
        if(subdet==3) m_pileupParamsECAL[make_pair(eta,layer)] = make_pair(a,b);
        else if(subdet==4) m_pileupParamsHCAL[make_pair(eta,layer)] = make_pair(a,b);
    }
    stream.close();

    // initialize energy corrections
    int nPUsubtraction = params.GetValue("JetPileupSubtraction.N", 1);
    for(int i=1;i<nPUsubtraction+1;i++)
    {
        stringstream fileName, name;
        name << "JetPileupSubtraction." << i <<".Name";
        fileName << "JetPileupSubtraction." << i <<".File";
        string pileupSubtractionName = params.GetValue(name.str().c_str(), "");
        string pileupSubtractionFileName = params.GetValue(fileName.str().c_str(), "");
        TFile* pileupSubtractionFile = TFile::Open(pileupSubtractionFileName.c_str());
        assert(pileupSubtractionFile);
        cout<<"INFO: "<<fileName.str()<<" ("<<pileupSubtractionName<<") = "<<pileupSubtractionFileName<<"\n";
        GBRForest* pileupSubtraction_up = (GBRForest*)pileupSubtractionFile->Get("EBCorrection");
        GBRForest* pileupSubtraction_down = (GBRForest*)pileupSubtractionFile->Get("EECorrection");
        assert(pileupSubtraction_up);
        assert(pileupSubtraction_down);
        m_pileupSubtractions[pileupSubtractionName] = make_pair(pileupSubtraction_up, pileupSubtraction_down);
        pileupSubtractionFile->Close();
    }

    int nThresholdCorrections = params.GetValue("JetThresholdCorrection.N", 1);
    for(int i=1;i<nThresholdCorrections+1;i++)
    {
        stringstream fileName, name;
        name << "JetThresholdCorrection." << i <<".Name";
        fileName << "JetThresholdCorrection." << i <<".File";
        string jetThresholdCorrectionName = params.GetValue(name.str().c_str(), "");
        string jetThresholdCorrectionFileName = params.GetValue(fileName.str().c_str(), "");
        TFile* jetThresholdCorrectionFile = TFile::Open(jetThresholdCorrectionFileName.c_str());
        assert(jetThresholdCorrectionFile);
        cout<<"INFO: "<<fileName.str()<<" ("<<jetThresholdCorrectionName<<") = "<<jetThresholdCorrectionFileName<<"\n";
        GBRForest* jetThresholdCorrection_up = (GBRForest*)jetThresholdCorrectionFile->Get("EBCorrection");
        GBRForest* jetThresholdCorrection_down = (GBRForest*)jetThresholdCorrectionFile->Get("EECorrection");
        assert(jetThresholdCorrection_up);
        assert(jetThresholdCorrection_down);
        m_thresholdCorrections[jetThresholdCorrectionName] = make_pair(jetThresholdCorrection_up, jetThresholdCorrection_down);
        jetThresholdCorrectionFile->Close();
    }
    string truthCorrectionFileName = params.GetValue("JetTruthCorrectionFile", "");
    TFile* truthCorrectionFile = TFile::Open(truthCorrectionFileName.c_str());
    assert(truthCorrectionFile);
    cout<<"INFO: JetTruthCorrectionFile = "<<truthCorrectionFileName<<"\n";
    //
    m_truthCorrection = (GBRForest*)truthCorrectionFile->Get("EECorrection");
    
    truthCorrectionFile->Close();

}



/*****************************************************************/
void TriggerJetAlgorithm::fillPileupEstimators(const EventHGCAL& event)
/*****************************************************************/
{
    /// ECAL
    for(auto& region_hits : m_regionHitsAboveThECAL)
    {
        region_hits.second.clear();
    }

    for(const auto& region_id : m_regionSimHitsECAL)
    {
        pair<int,int> region_layer = region_id.first;
        for(const auto& id_eta : region_id.second)
        {
            const HGCEEDetId& id = id_eta.first;
            const SimHit* hit = event.simhit(id);
            if(hit && hit->energy()>mip) m_regionHitsAboveThECAL[region_layer].push_back(hit);
        }
    }
    /// HCAL
    for(auto& region_hits : m_regionHitsAboveThHCAL)
    {
        region_hits.second.clear();
    }

    for(const auto& region_id : m_regionSimHitsHCAL)
    {
        pair<int,int> region_layer = region_id.first;
        for(const auto& id_eta : region_id.second)
        {
            const HGCHEDetId& id = id_eta.first;
            const SimHit* hit = event.simhitFH(id);
            if(hit && hit->energy()>mip_HF) m_regionHitsAboveThHCAL[region_layer].push_back(hit);
        }
    }
}




/*****************************************************************/
float TriggerJetAlgorithm::pileupThreshold(float eta, int layer, int nhits, int subdet)
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
    double a = 0.;
    double b = 1.;
    if(subdet==3)
    {
        const auto itr = m_pileupParamsECAL.find(make_pair(etaIndex,layer));
        if(itr==m_pileupParamsECAL.end())
        {
            cerr<<"ERROR: cannot find pileup parameters for eta="<<eta<<",layer="<<layer<<"\n";
            return 1.;
        }
        a = itr->second.first;
        b = itr->second.second;
    }
    else if(subdet==4)
    {
        const auto itr = m_pileupParamsHCAL.find(make_pair(etaIndex,layer));
        if(itr==m_pileupParamsHCAL.end())
        {
            cerr<<"ERROR: cannot find pileup parameters for eta="<<eta<<",layer="<<layer<<"\n";
            return 1.;
        }
        a = itr->second.first;
        b = itr->second.second;
    }
    return a*(double)nhits + b;
}

/*****************************************************************/
int TriggerJetAlgorithm::triggerRegionIndex(float eta, int zside, int sector, int subsector)
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
int TriggerJetAlgorithm::triggerRegionHits(int triggerRegion, int layer, int subdet)
/*****************************************************************/
{
    int nhits = 0;
    if(subdet==3)
    {
        const auto itr = m_regionHitsAboveThECAL.find(make_pair(triggerRegion, layer));
        if(itr==m_regionHitsAboveThECAL.end())
        {
            cerr<<"ERROR: cannot find number of hits in region "<<triggerRegion<<" layer "<<layer<<"\n";
            return 0;
        }
        nhits = itr->second.size();
    }
    else if(subdet==4)
    {
        const auto itr = m_regionHitsAboveThHCAL.find(make_pair(triggerRegion, layer));
        if(itr==m_regionHitsAboveThHCAL.end())
        {
            cerr<<"ERROR: cannot find number of hits in region "<<triggerRegion<<" layer "<<layer<<"\n";
            return 0;
        }
        nhits = itr->second.size();
    }
    return nhits;
}

/*****************************************************************/
void TriggerJetAlgorithm::pileupSubtraction(vector<Tower>& jets, std::string type)
/*****************************************************************/
{
    for(auto& jet : jets)
    {
        if(jet.nHits()==0) 
        {
            //cout<<"WARNING: trying to correct jet without hits\n";
            continue;
        }
        int triggerRegion = triggerRegionIndex(jet.eta(), jet.hits()[0]->zside(), jet.hits()[0]->sector(), jet.hits()[0]->subsector());
        int nhits = 0;
        nhits += triggerRegionHits(triggerRegion, 1, 3);
        nhits += triggerRegionHits(triggerRegion, 2, 3);
        nhits += triggerRegionHits(triggerRegion, 3, 3);
        nhits += triggerRegionHits(triggerRegion, 4, 3);
        nhits += triggerRegionHits(triggerRegion, 5, 3);
        double eta = fabs(jet.eta());
        double energy = jet.calibratedEnergy();
        float* puInputs = new float[2];
        puInputs[0] = eta;
        puInputs[1] = nhits;
        double esub = 0.;
        esub = (eta<=1.8 ? m_pileupSubtractions[type].first->GetResponse(puInputs) : m_pileupSubtractions[type].second->GetResponse(puInputs))*(double)jet.nHits();
        jet.setCalibratedEnergy(energy+esub);

        delete[] puInputs;
    }
}

/*****************************************************************/
void TriggerJetAlgorithm::thresholdCorrection(vector<Tower>& jets, std::string type)
/*****************************************************************/
{
    for(auto& jet : jets)
    {
        if(jet.nHits()==0) 
        {
            //cout<<"WARNING: trying to correct jet without hits\n";
            continue;
        }
        //int triggerRegion = triggerRegionIndex(jet.eta(), jet.hits()[0]->zside(), jet.hits()[0]->sector(), jet.hits()[0]->subsector());
        //int nhits = 0;
        //nhits += triggerRegionHits(triggerRegion, 1);
        //nhits += triggerRegionHits(triggerRegion, 2);
        //nhits += triggerRegionHits(triggerRegion, 3);
        //nhits += triggerRegionHits(triggerRegion, 4);
        //nhits += triggerRegionHits(triggerRegion, 5);
        double eta = fabs(jet.eta());
        double energy = jet.calibratedEnergy();
        float* thInputs = new float[2];
        thInputs[0] = energy;
        thInputs[1] = eta;
        double ecorr = 1.;
        ecorr = (eta<=1.8 ? m_thresholdCorrections[type].first->GetResponse(thInputs) : m_thresholdCorrections[type].second->GetResponse(thInputs));
        jet.setCalibratedEnergy(energy*ecorr);

        delete[] thInputs;
    }
}

/*****************************************************************/
void TriggerJetAlgorithm::truthCorrection(vector<Tower>& jets, const vector<Tower>& cores)
/*****************************************************************/
{
    if(jets.size()!=cores.size())
    {
        cout<<"ERROR: numbers of full jets and core jets are different!\n";
        return;
    }
    for(unsigned i=0; i<jets.size(); i++)
    {
        auto& jet = jets[i];
        const auto& core = cores[i];
        if(jet.nHits()==0) 
        {
            continue;
        }
        double eta = fabs(jet.eta());
        double energy = jet.calibratedEnergy();
        double et = jet.calibratedEt();
        double ooc = (core.calibratedEt()>0. ? et/core.calibratedEt(): 999.);
        float* thInputs = new float[3];
        thInputs[0] = energy;
        thInputs[1] = eta;
        thInputs[2] = ooc;
        double ecorr = 1.;
        for(unsigned i=0;i<3;i++)
        {
            ecorr = m_truthCorrection->GetResponse(thInputs);
            thInputs[0] = energy/ecorr;
        }
        jet.setCalibratedEnergy(energy/ecorr);

        delete[] thInputs;
    }
}



/*****************************************************************/
void TriggerJetAlgorithm::hcalSeeding(const EventHGCAL& event, const vector<SimHit>& hits, vector<Tower>& seeds, double clusterSize)
/*****************************************************************/
{
    fillPileupEstimators(event);
    TowerCalibrator towerCalibrator;

    // find hits above threshold in FH
    list<const SimHit*> selectedHits;
    for(const auto& hit : hits)
    {
        int triggerRegion = triggerRegionIndex(hit.eta(), hit.zside(), hit.sector(), hit.subsector());
        int nhits = triggerRegionHits(triggerRegion, hit.layer(), hit.subdet());
        float threshold = pileupThreshold(hit.eta(), hit.layer(), nhits, hit.subdet());
        //cout<<"hit eta="<<hit.eta()<<", layer="<<hit.layer()<<", nhits="<<nhits<<", th="<<threshold<<"\n";
        if(hit.energy()>threshold*mip_HF) selectedHits.push_back(&hit);
    }
    selectedHits.sort(AnHiMa::simHitSort);

    // Build clusters around selected hits, using only selected hits
    // If a hit is clustered, remove it from the list for the next clusters
    for(auto it1=selectedHits.begin(); it1!=std::prev(selectedHits.end()) && it1!=selectedHits.end(); it1++)
    {
        Tower cluster;
        const SimHit* hit1 = *it1;
        //cout<<"Main hit E="<<hit1->energy()<<", eta="<<hit1->eta()<<",phi="<<hit1->phi()<<"\n";
        cluster.addHit(*hit1);
        //
        GlobalPoint point1(hit1->x(), hit1->y(), hit1->z());
        GlobalPoint proj1( GlobalPoint::Polar(point1.theta(), point1.phi(), point1.mag()*354.77/point1.z()) );
        //cout<<"   Projected on eta="<<proj1.eta()<<",phi="<<proj1.phi()<<", z="<<proj1.z()<<"\n";
        for(auto it2=std::next(it1); it2!=selectedHits.end();)
        {
            const SimHit* hit2 = *it2;
            GlobalPoint point2(hit2->x(), hit2->y(), hit2->z());
            GlobalPoint proj2( GlobalPoint::Polar(point2.theta(), point2.phi(), point2.mag()*354.77/point2.z()) );
            double dx = proj1.x()-proj2.x();
            double dy = proj1.y()-proj2.y();
            double dr = sqrt(dx*dx+dy*dy);
            if(dr<clusterSize) // cluster hit
            {
                //cout<<"     Adding hit E="<<hit2->energy()<<", eta="<<hit2->eta()<<",phi="<<hit2->phi()<<"\n";
                cluster.addHit(*hit2);
                it2 = selectedHits.erase(it2);
            }
            else // move forward
            {
                it2++;
            }
        }
        // calibrate energy
        towerCalibrator.calibrate(cluster);
        //cout<<"Cluster energy="<<cluster.energy()<<", eta="<<cluster.eta()<<",phi="<<cluster.phi()<<"\n";
        // filter cluster candidates
        //if(cluster.calibratedEt()<1.) continue;
        // fill cluster
        seeds.push_back(cluster);
    }

}


/*****************************************************************/
void TriggerJetAlgorithm::superClustering(const EventHGCAL& event, const vector<Tower>& seeds, vector<SuperCluster>& superClusters, double coneSize)
/*****************************************************************/
{
    vector<SuperCluster*> sortedSuperClusters;
    // Find seed clusters
    // Seed clusters are clusters with max Et in a <coneSize> eta x phi cone
    for(const auto& cluster : seeds)
    {
        double eta0 = cluster.eta();
        double phi0 = cluster.phi();
        double et0 = cluster.et();
        bool isSeedCluster = true;
        // find neighbour clusters
        for(const auto& subcluster : seeds)
        {
            double eta = subcluster.eta();
            double phi = subcluster.phi();
            double et = subcluster.et();
            double deta = (eta-eta0);
            double dphi = TVector2::Phi_mpi_pi(phi-phi0);
            double dr = sqrt(deta*deta + dphi*dphi);
            if(fabs(deta)<1.e-5 && fabs(dphi)<1.e-5) continue;
            if(dr>coneSize) continue;
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
        }
    }
    for(auto& sc : superClusters) sortedSuperClusters.push_back(&sc);
    sort(sortedSuperClusters.begin(), sortedSuperClusters.end(), AnHiMa::superClusterSort);

    std::set<const Tower*> usedClusters;
    // Gather subclusters
    for(auto sc : sortedSuperClusters)
    {
        double eta0 = sc->eta();
        double phi0 = sc->phi();
        for(const auto& subcluster : seeds)
        {
            if(usedClusters.find(&subcluster)!=usedClusters.end()) continue;
            double eta = subcluster.eta();
            double phi = subcluster.phi();
            double deta = (eta-eta0);
            double dphi = TVector2::Phi_mpi_pi(phi-phi0);
            double dr = sqrt(deta*deta + dphi*dphi);
            // discard seed cluster
            if(fabs(deta)<1.e-5 && fabs(dphi)<1.e-5) continue;
            // outside window
            if(dr>coneSize) continue;
            //
            sc->addCluster(&subcluster);
            usedClusters.insert(&subcluster);
        }
    }
    //
}

/*****************************************************************/
void TriggerJetAlgorithm::coneClustering(const EventHGCAL& event, const vector<SimHit>& hits, const vector<SuperCluster>& superClusters, vector<Tower>& cones, double coneSize, bool applyThreshold)
/*****************************************************************/
{
    TowerCalibrator towerCalibrator;

    for(const auto& sc : superClusters)
    {
        double eta0 = sc.weightedEta();
        double phi0 = sc.weightedPhi();
        // Build cone around the ROI
        Tower cone;
        for(const auto& hit: hits)
        {
            float eta = hit.eta();
            float phi = hit.phi();
            double deta = (eta-eta0);
            double dphi =  TVector2::Phi_mpi_pi(phi-phi0);
            float dr = sqrt(deta*deta+dphi*dphi);
            if(dr>coneSize) continue;
            int triggerRegion = triggerRegionIndex(hit.eta(), hit.zside(), hit.sector(), hit.subsector());
            int nhits = triggerRegionHits(triggerRegion, hit.layer(), hit.subdet());
            float threshold = pileupThreshold(hit.eta(), hit.layer(), nhits, hit.subdet());
            double m = mip;
            switch(hit.subdet())
            {
                case 3:
                    m = mip;
                    break;
                case 4:
                    m = mip_HF;
                    break;
                default:
                    break;
            }
            // pass PU threshold
            if(!applyThreshold || hit.energy()>threshold*m) cone.addHit(hit);
        }
        //if(cone.nHits()==0) continue;
        towerCalibrator.calibrate(cone);
        cones.push_back(cone);
    }

}


/*****************************************************************/
vector<double> TriggerJetAlgorithm::computeJetProfile(const Tower& jet, double eta0, double phi0)
/*****************************************************************/
{
    vector<double> energies(40, 0.);
    TowerCalibrator calib;
    for(const auto& hit : jet.hits())
    {
        double deta = hit->eta() - eta0;
        double dphi = TVector2::Phi_mpi_pi(hit->phi() - phi0);
        double dr = sqrt(deta*deta + dphi*dphi);
        double energy = calib.calibratedEnergy(hit->energy(), hit->eta(), hit->layer(), hit->subdet());
        for(unsigned i=1;i<=40;i++)
        {
            double coneSize = (double)i*0.01; // bins of 0.01
            if(dr<coneSize)
            {
                for(unsigned j=i;j<=40;j++) // fill all larger cones
                {
                    energies[j-1] += energy;
                }
                break;
            }
        }
    }
    for(unsigned i=0;i<energies.size();i++)
    {
        energies[i] = energies[i]/cosh(eta0);
    }
    return energies;
}

