/**
 *  @file  Cluster.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    07/02/2014
 *
 *  @internal
 *     Created :  07/02/2014
 * Last update :  07/02/2014 05:04:59 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */

#include <math.h>
#include <iostream>

#include "TVector2.h"

#include "AnHiMaHGCAL/HGCALCommon/interface/Cluster.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Event.h"

using namespace std;
using namespace AnHiMa;


/*****************************************************************/
Cluster::Cluster():TLorentzVector(),
    m_sigmaEtaEta(0.),
    m_sigmaPhiPhi(0.),
    m_sigmaIzIz(0.),
    m_maxLayer(-1)
/*****************************************************************/
{
}

/*****************************************************************/
Cluster::~Cluster()
/*****************************************************************/
{
}


/*****************************************************************/
void Cluster::build(Event& event, double eta0, double phi0)
/*****************************************************************/
{
    double etaWidth = 0.05; // half widths
    double phiWidth = 0.2;
    double clusterEnergy = 0.;
    double clusterEta = 0.;
    double clusterPhi = 0.;
    double clusterLayer = 0.;

    vector<double> etas;
    vector<double> phis;
    vector<double> layers;
    vector<double> energies;
    map<int,double> layerEnergies;
    map<int, vector<double> > layerCellEnergies;
    // compute energy and position of the cluster
    for(int i=0;i<event.nhits;i++)
    {
        if(event.hit_type[i]!=0) continue;
        double eta = event.hit_eta[i];
        double phi = event.hit_phi[i];
        //double z   = event.hit_z[i];
        int layer = event.hit_layer[i];
        double energy = event.hit_edep[i]; 
        double deta = eta-eta0;
        double dphi = TVector2::Phi_mpi_pi(phi-phi0);
        if(fabs(deta)<etaWidth && fabs(dphi)<phiWidth) // box 
        {
            clusterEnergy += energy;
            clusterEta += eta*energy;
            clusterPhi += phi*energy;
            clusterLayer += (double)layer*energy;
            etas.push_back(eta);
            phis.push_back(phi);
            layers.push_back(layer);
            energies.push_back(energy);
            if(layerEnergies.find(layer)==layerEnergies.end()) layerEnergies[layer] = 0.;
            if(layerCellEnergies.find(layer)==layerCellEnergies.end()) layerCellEnergies[layer] = vector<double>();
            layerEnergies[layer] += energy;
            layerCellEnergies[layer].push_back(energy);
        }
    }
    if(clusterEnergy==0.) 
    {
        SetPtEtaPhiM(0.,0.,0.,0.);
        return;
    }
    clusterEta /= clusterEnergy;
    clusterPhi /= clusterEnergy;
    clusterLayer /= clusterEnergy;
    double clusterEt = clusterEnergy/cosh(clusterEta);
    SetPtEtaPhiM(clusterEt, clusterEta, clusterPhi, 0.);
    m_layerEnergies = layerEnergies;

    // compute shape variables
    double clusterSigmaIzIz = 0.;
    for(unsigned i=0;i<energies.size();i++)
    {
        double layer = layers[i];
        double energy = energies[i];
        double dlayer = (double)layer-clusterLayer;
        clusterSigmaIzIz += dlayer*dlayer*energy;
    }
    clusterSigmaIzIz   /= clusterEnergy;
    m_sigmaIzIz = sqrt(clusterSigmaIzIz);

    computeSigmaEtaEta(event, clusterEta, clusterPhi);
    computeSigmaPhiPhi(event, clusterEta, clusterPhi);

    // 
    for(auto itr=layerEnergies.cbegin(); itr!=layerEnergies.end(); itr++)
    {
        int layer = itr->first;
        double layerEnergy = itr->second;
        vector<double> cellEnergies = layerCellEnergies[layer];
        sort(cellEnergies.begin(), cellEnergies.end());
        reverse(cellEnergies.begin(), cellEnergies.end());
        double energy4 = 0.;
        for(int i=0;i<min(4,(int)cellEnergies.size());i++)
        {
            energy4 += cellEnergies[i];
        }
        m_layerEnergyCompactness[layer] = energy4/layerEnergy;
        m_layerEnergyOverTot[layer] = layerEnergy/clusterEnergy;
    }


    // find maximum layer
    int maxLayer = -1;
    double maxEnergy = 0.;
    for(auto itr=layerEnergies.cbegin();itr!=layerEnergies.end();itr++)
    {
        if(itr->second>maxEnergy)
        {
            maxLayer = itr->first;
            maxEnergy = itr->second;
        }
    }
    m_maxLayer = maxLayer;
    
}


/*****************************************************************/
void Cluster::computeSigmaEtaEta(Event& event, double eta0, double phi0)
/*****************************************************************/
{
    double etaWidth = 0.05; // half widths
    double phiWidth = 0.02;
    double clusterEnergy = 0.;

    vector<double> etas;
    vector<double> layers;
    vector<double> energies;
    map<int,double> layerEnergies;
    // compute energy and position of the cluster
    for(int i=0;i<event.nhits;i++)
    {
        if(event.hit_type[i]!=0) continue;
        double eta = event.hit_eta[i];
        double phi = event.hit_phi[i];
        //double z   = event.hit_z[i];
        int layer = event.hit_layer[i];
        double energy = event.hit_edep[i]; 
        double deta = eta-eta0;
        double dphi = TVector2::Phi_mpi_pi(phi-phi0);
        if(fabs(deta)<etaWidth && fabs(dphi)<phiWidth) // box 
        {
            clusterEnergy += energy;
            etas.push_back(eta);
            layers.push_back(layer);
            energies.push_back(energy);
            if(layerEnergies.find(layer)==layerEnergies.end()) layerEnergies[layer] = 0.;
            layerEnergies[layer] += energy;
        }
    }
    if(clusterEnergy==0.) return;

    map<int,double> layerEnergiesLog;
    // compute shape variables
    double clusterSigmaEtaEta = 0.;
    map<int,double> layerSigmaEtaEta;
    map<int,double> layerSigmaEta2Log;
    map<int,double> layerSigmaEta3;
    map<int,double> layerSigmaEta4;
    for(unsigned i=0;i<energies.size();i++)
    {
        double eta = etas[i];
        double layer = layers[i];
        double energy = energies[i];
        double deta = eta-eta0;
        clusterSigmaEtaEta += deta*deta*energy;

        if(layerEnergiesLog.find(layer)==layerEnergiesLog.end()) layerEnergiesLog[layer] = 0.;
        layerEnergiesLog[layer] += log(1.+energy/layerEnergies[layer]);
        //
        if(layerSigmaEtaEta.find(layer)==layerSigmaEtaEta.end()) layerSigmaEtaEta[layer] = 0.;
        if(layerSigmaEta3.find(layer)==layerSigmaEta3.end()) layerSigmaEta3[layer] = 0.;
        if(layerSigmaEta4.find(layer)==layerSigmaEta4.end()) layerSigmaEta4[layer] = 0.;
        if(layerSigmaEta2Log.find(layer)==layerSigmaEta2Log.end()) layerSigmaEta2Log[layer] = 0.;
        layerSigmaEtaEta[layer] += deta*deta*energy;
        layerSigmaEta3[layer] += deta*deta*deta*energy;
        layerSigmaEta4[layer] += deta*deta*deta*deta*energy;
        layerSigmaEta2Log[layer] += deta*deta*log(1.+energy/layerEnergies[layer]);
    }
    clusterSigmaEtaEta /= clusterEnergy;
    m_sigmaEtaEta = sqrt(clusterSigmaEtaEta);

    for(auto itr=layerEnergies.cbegin();itr!=layerEnergies.end();itr++)
    {
        int layer = itr->first;
        double energy = itr->second;
        double energyLog = layerEnergiesLog[layer];
        layerSigmaEtaEta[layer] /= energy;
        layerSigmaEta3[layer] /= energy;
        layerSigmaEta4[layer] /= energy;
        layerSigmaEta2Log[layer] /= energyLog;
        layerSigmaEtaEta[layer] = sqrt(layerSigmaEtaEta[layer]);

    }
    m_layerSigmaEtaEta = layerSigmaEtaEta;
    m_layerSigmaEta3 = layerSigmaEta3;
    m_layerSigmaEta4 = layerSigmaEta4;
    m_layerSigmaEta2Log = layerSigmaEta2Log;
}

/*****************************************************************/
void Cluster::computeSigmaPhiPhi(Event& event, double eta0, double phi0)
/*****************************************************************/
{
    double etaWidth = 0.02; // half widths
    double phiWidth = 0.2;
    double clusterEnergy = 0.;

    vector<double> etas;
    vector<double> phis;
    vector<double> layers;
    vector<double> energies;
    map<int,double> layerEnergies;
    // compute energy and position of the cluster
    for(int i=0;i<event.nhits;i++)
    {
        if(event.hit_type[i]!=0) continue;
        double eta = event.hit_eta[i];
        double phi = event.hit_phi[i];
        //double z   = event.hit_z[i];
        int layer = event.hit_layer[i];
        double energy = event.hit_edep[i]; 
        double deta = eta-eta0;
        double dphi = TVector2::Phi_mpi_pi(phi-phi0);
        if(fabs(deta)<etaWidth && fabs(dphi)<phiWidth) // box 
        {
            clusterEnergy += energy;
            etas.push_back(eta);
            phis.push_back(phi);
            layers.push_back(layer);
            energies.push_back(energy);
            if(layerEnergies.find(layer)==layerEnergies.end()) layerEnergies[layer] = 0.;
            layerEnergies[layer] += energy;
        }
    }
    if(clusterEnergy==0.) return;

    // compute shape variables
    double clusterSigmaPhiPhi = 0.;
    //double clusterSigmaIzIz = 0.;
    map<int,double> layerSigmaPhiPhi;
    for(unsigned i=0;i<energies.size();i++)
    {
        double phi = phis[i];
        double layer = layers[i];
        double energy = energies[i];
        double dphi = TVector2::Phi_mpi_pi(phi-phi0);
        clusterSigmaPhiPhi += dphi*dphi*energy;
        if(layerSigmaPhiPhi.find(layer)==layerSigmaPhiPhi.end()) layerSigmaPhiPhi[layer] = 0.;
        layerSigmaPhiPhi[layer] += dphi*dphi*energy;
    }
    clusterSigmaPhiPhi /= clusterEnergy;
    m_sigmaPhiPhi = sqrt(clusterSigmaPhiPhi);

    for(auto itr=m_layerEnergies.cbegin();itr!=m_layerEnergies.end();itr++)
    {
        int layer = itr->first;
        double energy = itr->second;
        layerSigmaPhiPhi[layer] /= energy;
        layerSigmaPhiPhi[layer] = sqrt(layerSigmaPhiPhi[layer]);
    }
    m_layerSigmaPhiPhi = layerSigmaPhiPhi;
}


/*****************************************************************/
double Cluster::layerEnergy(int layer)
/*****************************************************************/
{
    double energy = 0.;
    auto itr = m_layerEnergies.find(layer);
    if(itr!=m_layerEnergies.end()) energy = itr->second;

    return energy;
}

/*****************************************************************/
double Cluster::layerSigmaEtaEta(int layer)
/*****************************************************************/
{
    double sigmaEtaEta = 0.;
    auto itr = m_layerSigmaEtaEta.find(layer);
    if(itr!=m_layerSigmaEtaEta.end()) sigmaEtaEta = itr->second;

    return sigmaEtaEta;
}

/*****************************************************************/
double Cluster::layerSigmaEta3(int layer)
/*****************************************************************/
{
    double sigmaEta3 = 0.;
    auto itr = m_layerSigmaEta3.find(layer);
    if(itr!=m_layerSigmaEta3.end()) sigmaEta3 = itr->second;

    return sigmaEta3;
}

/*****************************************************************/
double Cluster::layerSigmaEta4(int layer)
/*****************************************************************/
{
    double sigmaEta4 = 0.;
    auto itr = m_layerSigmaEta4.find(layer);
    if(itr!=m_layerSigmaEta4.end()) sigmaEta4 = itr->second;

    return sigmaEta4;
}

/*****************************************************************/
double Cluster::layerSigmaEta2Log(int layer)
/*****************************************************************/
{
    double sigmaEta2Log = 0.;
    auto itr = m_layerSigmaEta2Log.find(layer);
    if(itr!=m_layerSigmaEta2Log.end()) sigmaEta2Log = itr->second;

    return sigmaEta2Log;
}


/*****************************************************************/
double Cluster::layerSigmaPhiPhi(int layer)
/*****************************************************************/
{
    double sigmaPhiPhi = 0.;
    auto itr = m_layerSigmaPhiPhi.find(layer);
    if(itr!=m_layerSigmaPhiPhi.end()) sigmaPhiPhi = itr->second;

    return sigmaPhiPhi;
}

/*****************************************************************/
double Cluster::layerEnergyCompactness(int layer)
/*****************************************************************/
{
    double energyCompactness = 0.;
    auto itr = m_layerEnergyCompactness.find(layer);
    if(itr!=m_layerEnergyCompactness.end()) energyCompactness = itr->second;

    return energyCompactness;
}

/*****************************************************************/
double Cluster::layerEnergyOverTot(int layer)
/*****************************************************************/
{
    double energyOverTot = 0.;
    auto itr = m_layerEnergyOverTot.find(layer);
    if(itr!=m_layerEnergyOverTot.end()) energyOverTot = itr->second;

    return energyOverTot;
}
