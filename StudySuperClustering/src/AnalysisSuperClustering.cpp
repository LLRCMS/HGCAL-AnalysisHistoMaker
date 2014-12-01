/**
 *  @file  AnalysisSuperClustering.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    19/10/2014
 *
 *  @internal
 *     Created :  19/10/2014
 * Last update :  19/10/2014 13:44:32
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#include <iostream>
#include <fstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/StudySuperClustering/interface/AnalysisSuperClustering.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisSuperClustering::AnalysisSuperClustering():IAnalysis(),
    m_genParticle(0),
    m_matchedSuperCluster(0),
    m_matchedIdealCluster(0)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisSuperClustering::~AnalysisSuperClustering()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisSuperClustering::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    // Read parameters
    //string pileupParamsFile = m_reader.params().GetValue("PileupParams", "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/src/AnHiMaHGCAL/ElectronClusterThreshold/data/thresholdParameters.txt");

    m_egammaAlgo.initialize(event(), m_reader.params());
    return true;
}


/*****************************************************************/
void AnalysisSuperClustering::execute()
/*****************************************************************/
{
    event().update();
    cerr<<"Event "<<event().event()<<"\n";
    if(!event().passSelection()) return;

    m_seeds.clear();
    m_clusters.clear();
    m_superClusters.clear();
    // seeding
    m_egammaAlgo.seeding(event(), m_seeds);
    // clustering
    m_egammaAlgo.clustering(event(), m_seeds, m_clusters);
    // cluster energy corrections
    m_egammaAlgo.clusterPileupSubtraction(m_clusters);
    m_egammaAlgo.clusterThresholdCorrection(m_clusters);
    // super clustering
    m_egammaAlgo.superClustering(m_clusters, m_superClusters);
    // supercluster energy corrections
    m_egammaAlgo.superClusterCorrection(m_superClusters);


    for(int i=0;i<2;i++)
    {
        m_genParticle = &event().genparticles()[i];
        if(fabs(m_genParticle->Eta())<2.9 && fabs(m_genParticle->Eta())>1.6) 
        {
            // ideal clustering
            m_idealClusters.clear();
            m_egammaAlgo.idealClustering(event(), m_idealClusters, m_genParticle->Eta(), m_genParticle->Phi());
            m_matchedIdealCluster = (m_idealClusters.size()>0 ? &m_idealClusters[0] : 0);
            // match supercluster
            matchSuperCluster();
            fillHistos();
        }
    }
}

/*****************************************************************/
void AnalysisSuperClustering::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;


    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    //

    if(!m_matchedSuperCluster) return;
    // seed cluster
    const Tower* seedCluster = m_matchedSuperCluster->cluster(0);
    map<int, pair<double,double> >  sizes = clusterSizes(*seedCluster);
    for(auto layer_sizes : sizes)
    {
        int layer = layer_sizes.first;
        double dx = layer_sizes.second.first;
        double dy = layer_sizes.second.second;
        double dr = sqrt(dx*dx+dy*dy);
        m_histos.Fill1BinHisto(100+hoffset, layer, dx, weight, sysNum);
        m_histos.Fill1BinHisto(150+hoffset, layer, dy, weight, sysNum);
        m_histos.Fill1BinHisto(200+hoffset, layer, dx, dy, weight, sysNum);
        m_histos.Fill1BinHisto(250+hoffset, layer, seedCluster->layerNHits(layer), weight, sysNum);
        m_histos.Fill1BinHisto(300+hoffset, layer, dr, seedCluster->layerNHits(layer), weight, sysNum);
        //cout<<"Layer "<<layer_sizes.first<<" : "<<layer_sizes.second.first<<","<<layer_sizes.second.second<<"\n";
    }
    if(m_matchedSuperCluster->nClusters()>1)
    {
        for(int c=1;c<m_matchedSuperCluster->nClusters();c++)
        {
            const Tower* subCluster = m_matchedSuperCluster->cluster(c);
            map<int, pair<double,double> >  subsizes = clusterSizes(*subCluster);
            for(auto layer_sizes : subsizes)
            {
                int layer = layer_sizes.first;
                double dx = layer_sizes.second.first;
                double dy = layer_sizes.second.second;
                //double dr = sqrt(dx*dx+dy*dy);
                m_histos.Fill1BinHisto(500+hoffset, layer, dx, weight, sysNum);
                m_histos.Fill1BinHisto(550+hoffset, layer, dy, weight, sysNum);
                m_histos.Fill1BinHisto(600+hoffset, layer, dx, dy, weight, sysNum);
                m_histos.Fill1BinHisto(650+hoffset, layer, subCluster->layerNHits(layer), weight, sysNum);
            }
            // position wrt seed
            double dx = subCluster->x()-seedCluster->x();
            double dy = subCluster->y()-seedCluster->y();
            double dxy = sqrt(dx*dx+dy*dy);
            double deta =  subCluster->eta()-seedCluster->eta();
            double dphi =  subCluster->phi()-seedCluster->phi();
            //double dep = sqrt(deta*deta+dphi*dphi);
            m_histos.Fill1BinHisto(1000+hoffset, fabs(m_genParticle->Eta()), fabs(dx), weight, sysNum);
            m_histos.Fill1BinHisto(1010+hoffset, fabs(m_genParticle->Eta()), fabs(dy), weight, sysNum);
            m_histos.Fill1BinHisto(1020+hoffset, fabs(m_genParticle->Eta()), fabs(dx), fabs(dy), weight, sysNum);
            m_histos.Fill1BinHisto(1030+hoffset, fabs(m_genParticle->Eta()), dxy, weight, sysNum);
            m_histos.Fill1BinHisto(1040+hoffset, fabs(m_genParticle->Eta()), fabs(deta), weight, sysNum);
            m_histos.Fill1BinHisto(1050+hoffset, fabs(m_genParticle->Eta()), fabs(dphi), weight, sysNum);
        }
    }

    // difference with ideal cluster
    if(m_matchedIdealCluster)
    {
        double eta0 = m_matchedIdealCluster->eta();
        double phi0 = m_matchedIdealCluster->phi();
        double x0 = m_matchedIdealCluster->x();
        double y0 = m_matchedIdealCluster->y();
        for(const auto hit : m_matchedIdealCluster->hits())
        {
            bool found = false;
            for(int c=0;c<m_matchedSuperCluster->nClusters();c++)
            {
                const Tower* cluster = m_matchedSuperCluster->cluster(c);
                for(const auto hit2 : cluster->hits())
                {
                    if(hit->detid()==hit2->detid())
                    {
                        found = true;
                        break;
                    }
                }
                if(found==true) break;
            }
            if(!found)
            {
                int layer = hit->layer();
                double eta = hit->eta();
                double phi = hit->phi();
                double x = hit->x();
                double y = hit->y();
                double energy = hit->energy();
                double deta = (eta - eta0);
                double dphi = TVector2::Phi_mpi_pi(phi - phi0);
                double dx = (x - x0);
                double dy = (y - y0);
                // bad super-cluster
                if((m_matchedSuperCluster->et()-m_genParticle->Pt())/m_genParticle->Pt()<-0.7)
                {
                    m_histos.Fill1BinHisto(2200+hoffset, layer, deta, dphi, weight*energy, sysNum);
                    m_histos.Fill1BinHisto(2300+hoffset, layer, dx, dy, weight*energy, sysNum);
                }
                // good super-cluster
                else
                {
                    m_histos.Fill1BinHisto(2000+hoffset, layer, deta, dphi, weight*energy, sysNum);
                    m_histos.Fill1BinHisto(2100+hoffset, layer, dx, dy, weight*energy, sysNum);
                }
            }
        }
        double seedEta = m_matchedSuperCluster->cluster(0)->eta();
        double seedPhi = m_matchedSuperCluster->cluster(0)->phi();
        double seedX   = m_matchedSuperCluster->cluster(0)->x();
        double seedY   = m_matchedSuperCluster->cluster(0)->y();
        double seedDeta = (seedEta - eta0);
        double seedDphi = TVector2::Phi_mpi_pi(seedPhi - phi0);
        double seedDx = (seedX - x0);
        double seedDy = (seedY - y0);
        // bad super-cluster
        if((m_matchedSuperCluster->et()-m_genParticle->Pt())/m_genParticle->Pt()<-0.7)
        {
            m_histos.FillHisto(2502+hoffset, seedDeta, seedDphi, weight, sysNum);
            m_histos.FillHisto(2503+hoffset, seedDx, seedDy, weight, sysNum);
        }
        // good super-cluster
        else
        {
            m_histos.FillHisto(2500+hoffset, seedDeta, seedDphi, weight, sysNum);
            m_histos.FillHisto(2501+hoffset, seedDx, seedDy, weight, sysNum);
        }
             
    }
}




/*****************************************************************/
void AnalysisSuperClustering::matchSuperCluster()
/*****************************************************************/
{
    m_matchedSuperCluster = 0;

    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxEt = 0.;
    for(const auto& sc : m_superClusters )
    {
        double deta = (sc.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(sc.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && sc.et()>maxEt)
        {
            m_matchedSuperCluster = &sc;
            maxEt = sc.et();
        }
    }
}


/*****************************************************************/
map<int, pair<double,double> > AnalysisSuperClustering::clusterSizes(const Tower& cluster)
/*****************************************************************/
{
    map<int, pair<double,double> > clusterPositions;
    map<int, pair<double,double> > clusterPositionSquares;
    map<int, pair<double,double> > clusterSizes;
    map<int, double> layerEnergies;
    for(int l=1;l<=30;l++)
    {
        layerEnergies[l] = 0.;
        clusterSizes[l] = make_pair(0.,0.);
        clusterPositions[l] = make_pair(0.,0.);
        clusterPositionSquares[l] = make_pair(0.,0.);
    }
    // find cluster hits in each layer
    map<int, vector<const SimHit*> > clusterHits;
    for(const auto& hit : cluster.hits())
    {
        int layer = hit->layer();
        if(clusterHits.find(layer)==clusterHits.end()) clusterHits[layer] = vector<const SimHit*>();
        clusterHits[layer].push_back(hit);
    }

    // compute energy weighted x and y (and x**2 and y**2) and RMSs
    for(auto layer_hits : clusterHits)
    {
        int layer = layer_hits.first;
        vector<const SimHit*> hits = layer_hits.second;
        for(auto hit : hits)
        {
            clusterPositions[layer].first  += hit->x()*hit->energy();
            clusterPositions[layer].second += hit->y()*hit->energy();
            clusterPositionSquares[layer].first  += hit->x()*hit->x()*hit->energy();
            clusterPositionSquares[layer].second += hit->y()*hit->y()*hit->energy();
            layerEnergies[layer] += hit->energy();
        }
    }
    for(int l=1;l<=30;l++)
    {
        double energy = layerEnergies[l];
        if(energy>0.)
        {
            clusterPositions[l].first /= energy;
            clusterPositions[l].second /= energy;
            clusterPositionSquares[l].first /= energy;
            clusterPositionSquares[l].second /= energy;
            //
            clusterSizes[l].first = sqrt(max(0.,clusterPositionSquares[l].first - clusterPositions[l].first*clusterPositions[l].first));
            clusterSizes[l].second = sqrt(max(0.,clusterPositionSquares[l].second - clusterPositions[l].second*clusterPositions[l].second));
        }
    }
    return clusterSizes;
}



