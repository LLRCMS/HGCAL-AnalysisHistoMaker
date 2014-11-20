/**
 *  @file  AnalysisClustering.cpp
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

#include "AnHiMaHGCAL/StudyClustering/interface/AnalysisClustering.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisClustering::AnalysisClustering():IAnalysis(),
    m_genParticle(0)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisClustering::~AnalysisClustering()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisClustering::initialize(const string& parameterFile)
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
void AnalysisClustering::execute()
/*****************************************************************/
{
    event().update();
    if(!event().passSelection()) return;

    m_seeds.clear();
    m_clusters.clear();
    m_egammaAlgo.seeding(event(), m_seeds);
    m_egammaAlgo.clustering(event(), m_seeds, m_clusters);


    fillHistos();

}

/*****************************************************************/
void AnalysisClustering::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;


    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    //
    m_histos.FillHisto(10+hoffset, m_seeds.size(), weight, sysNum);

    for(const auto& tower : m_seeds)
    {
        m_histos.FillHisto(11+hoffset, fabs(tower.eta()), weight, sysNum);
        m_histos.FillHisto(12+hoffset, tower.phi(), weight, sysNum);
        m_histos.FillHisto(13+hoffset, fabs(tower.eta()), tower.phi(), weight, sysNum);
        m_histos.FillHisto(14+hoffset, tower.calibratedEt(), weight, sysNum);

        //m_histos.Fill1BinHisto(100+hoffset, event().event()-1, fabs(tower.eta()), tower.phi(), weight, sysNum);
    }

    for(const auto& cluster : m_clusters)
    {
        m_histos.FillHisto(21+hoffset, fabs(cluster.eta()), weight, sysNum);
        m_histos.FillHisto(22+hoffset, cluster.phi(), weight, sysNum);
        m_histos.FillHisto(24+hoffset, cluster.calibratedEt(), weight, sysNum);
    }

    for(int i=0;i<2;i++)
    {
        matchCluster(i);
        m_histos.FillHisto(50+hoffset, m_matchedClusters.size(), weight, sysNum);
        for(const auto cluster : m_matchedClusters)
        {
            double deta = (cluster->eta()-m_genParticle->Eta());
            double dphi = TVector2::Phi_mpi_pi(cluster->phi()-m_genParticle->Phi());
            double response = (cluster->calibratedEt() - m_genParticle->Et())/m_genParticle->Et();
            m_histos.FillHisto(51+hoffset, deta, weight, sysNum);
            m_histos.FillHisto(52+hoffset, dphi, weight, sysNum);
            m_histos.FillHisto(53+hoffset, deta, dphi, weight, sysNum);
            m_histos.FillHisto(54+hoffset, response, weight, sysNum);

            // position with respect to main cluster
            if(cluster!=m_matchedClusters[0])
            {
                double detaSub = (cluster->eta()-m_matchedClusters[0]->eta());
                double dphiSub = TVector2::Phi_mpi_pi(cluster->phi()-m_matchedClusters[0]->phi());
                m_histos.FillHisto(61+hoffset, detaSub, weight, sysNum);
                m_histos.FillHisto(62+hoffset, dphiSub, weight, sysNum);
                m_histos.FillHisto(63+hoffset, detaSub, dphiSub, weight, sysNum);
            }
        }
    }

}



/*****************************************************************/
void AnalysisClustering::matchCluster(int i)
/*****************************************************************/
{
    m_matchedClusters.clear();

    m_genParticle = &event().genparticles()[i];
    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    for(const auto& cluster : m_clusters )
    {
        double deta = (cluster.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cluster.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.3)
        {
            m_matchedClusters.push_back(&cluster);
        }
    }
    sort(m_matchedClusters.begin(), m_matchedClusters.end(), AnHiMa::towerSort);
    //cout<<"Clusters Et:";
    //for(const auto cluster : m_matchedClusters)
    //{
    //    cout<<cluster->calibratedEt()<<", ";
    //}
    //cout<<"\n";
}
