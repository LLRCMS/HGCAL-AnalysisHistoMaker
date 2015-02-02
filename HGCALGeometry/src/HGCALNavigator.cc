/**
 *  @file  HGCALNavigator.cc
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    19/10/2014
 *
 *  @internal
 *     Created :  19/10/2014
 * Last update :  19/10/2014 11:22:02
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */


#include <vector>
#include <string>

#include "TVector2.h"

#include "AnHiMaHGCAL/HGCALGeometry/interface/HGCALNavigator.h"
#include "AnHiMaHGCAL/HGCALGeometry/interface/GeometryConfigurationHack.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/PluginManager/interface/standard.h"

#include "DetectorDescription/Base/interface/DDdebug.h"
#include "DetectorDescription/Parser/interface/DDLParser.h"
#include "DetectorDescription/Core/interface/DDRoot.h"
#include "DetectorDescription/Core/interface/DDMaterial.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Core/interface/DDSpecifics.h"
#include "DetectorDescription/Base/interface/DDRotationMatrix.h"
#include "DetectorDescription/Core/src/Material.h"
#include "DetectorDescription/Core/src/Solid.h"
#include "DetectorDescription/Core/src/LogicalPart.h"
#include "DetectorDescription/Core/src/Specific.h"

#include "Geometry/FCalGeometry/interface/HGCalGeometryLoader.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"



using namespace std;


/*****************************************************************/
HGCALNavigator::HGCALNavigator()
/*****************************************************************/
{
    // initialize CMSSW plugin manager
    edmplugin::PluginManager::configure(edmplugin::standard::config());
}



/*****************************************************************/
HGCALNavigator::~HGCALNavigator()
/*****************************************************************/
{
    delete m_ddcv;
    for(auto& itr : m_hgcdc) delete itr.second;
    for(auto& itr : m_hgctopo) delete itr.second;
    for(auto& itr : m_hgcgeom) delete itr.second;
}



/*****************************************************************/
bool HGCALNavigator::initialize()
/*****************************************************************/
{
    // xml geometry files defined in https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_X_SLHC/Geometry/CMSCommonData/python/cmsExtendedGeometry2023HGCalMuonXML_cfi.py 
    vector<string> geovec;
    geovec.push_back("Geometry/CMSCommonData/data/materials.xml");
    geovec.push_back("Geometry/CMSCommonData/data/rotations.xml");
    geovec.push_back("Geometry/CMSCommonData/data/extend/cmsextent.xml");
    geovec.push_back("Geometry/CMSCommonData/data/PhaseI/cms.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cmsMother.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cmsTracker.xml");
    geovec.push_back("Geometry/CMSCommonData/data/eta3/etaMax.xml");
    geovec.push_back("Geometry/CMSCommonData/data/PhaseII/caloBase.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cmsCalo.xml");
    geovec.push_back("Geometry/CMSCommonData/data/PhaseII/muonBase.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cmsMuon.xml");
    geovec.push_back("Geometry/CMSCommonData/data/mgnt.xml");
    geovec.push_back("Geometry/CMSCommonData/data/PhaseII/beampipe.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cmsBeam.xml");
    geovec.push_back("Geometry/CMSCommonData/data/muonMB.xml");
    geovec.push_back("Geometry/CMSCommonData/data/muonMagnet.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cavern.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdMaterials.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdCommon.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdCylinder.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwd.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdDisks.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdInnerDisk1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdInnerDisk2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdInnerDisk3.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdInnerDisk4.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdInnerDisk5.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdInnerDisk6.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdInnerDisk7.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdOuterDisk1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdOuterDisk2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdOuterDisk3.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdOuterDisk4.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdOuterDisk5.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdOuterDisk6.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdOuterDisk7.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdOuterDisk8.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdOuterDisk9.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdOuterDisk10.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdblade1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdblade2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdblade3.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdblade4.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdblade5.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdblade6.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdblade7.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdblade8.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdblade9.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/pixfwdblade10.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarmaterial.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarladder.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarladderfull0.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarladderfull1.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarladderfull2.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarladderfull3.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarlayer.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarlayer0.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarlayer1.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarlayer2.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/pixbarlayer3.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/BarrelEndcap5D/pixbar.xml"); 
    geovec.push_back("Geometry/TrackerCommonData/data/trackermaterial.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/BarrelEndcap5D/tracker.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackerpixbar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseI/trackerpixfwd.xml");
    geovec.push_back("Geometry/EcalCommonData/data/ectkcable.xml");
    geovec.push_back("Geometry/EcalCommonData/data/PhaseII/eregalgo.xml");
    geovec.push_back("Geometry/EcalCommonData/data/ebalgo.xml");
    geovec.push_back("Geometry/EcalCommonData/data/ebcon.xml");
    geovec.push_back("Geometry/EcalCommonData/data/ebrot.xml");
    geovec.push_back("Geometry/EcalCommonData/data/eecon.xml");
    geovec.push_back("Geometry/EcalCommonData/data/PhaseII/escon.xml");
    geovec.push_back("Geometry/EcalCommonData/data/PhaseII/esalgo.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalrotations.xml");
    geovec.push_back("Geometry/HcalCommonData/data/PhaseII/NoHE/hcalalgo.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalbarrelalgo.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalouteralgo.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalforwardalgo.xml");
    geovec.push_back("Geometry/HcalCommonData/data/PhaseII/NoHE/hcalSimNumbering.xml");
    geovec.push_back("Geometry/HcalCommonData/data/PhaseII/NoHE/hcalRecNumbering.xml");
    geovec.push_back("Geometry/HcalCommonData/data/average/hcalforwardmaterial.xml");
    geovec.push_back("Geometry/HGCalCommonData/data/v5/hgcal.xml");
    geovec.push_back("Geometry/HGCalCommonData/data/v5/hgcalEE.xml");
    geovec.push_back("Geometry/HGCalCommonData/data/v5/hgcalHEsil.xml");
    geovec.push_back("Geometry/HGCalCommonData/data/v5/hgcalHEsci.xml");
    geovec.push_back("Geometry/HGCalCommonData/data/v5/hgcalCons.xml");
    geovec.push_back("Geometry/MuonCommonData/data/v1/mbCommon.xml");
    geovec.push_back("Geometry/MuonCommonData/data/v1/mb1.xml");
    geovec.push_back("Geometry/MuonCommonData/data/v1/mb2.xml");
    geovec.push_back("Geometry/MuonCommonData/data/v1/mb3.xml");
    geovec.push_back("Geometry/MuonCommonData/data/v1/mb4.xml");
    geovec.push_back("Geometry/MuonCommonData/data/design/muonYoke.xml");
    geovec.push_back("Geometry/MuonCommonData/data/PhaseII/mf.xml");
    geovec.push_back("Geometry/MuonCommonData/data/PhaseII/rpcf.xml");
    geovec.push_back("Geometry/MuonCommonData/data/v2/gemf.xml");
    geovec.push_back("Geometry/MuonCommonData/data/v4/gem11.xml");
    geovec.push_back("Geometry/MuonCommonData/data/v6/gem21.xml");
    geovec.push_back("Geometry/MuonCommonData/data/v2/csc.xml");
    geovec.push_back("Geometry/MuonCommonData/data/PhaseII/mfshield.xml");
    geovec.push_back("Geometry/MuonCommonData/data/PhaseII/me0.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/forward.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/v2/forwardshield.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/brmrotations.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/brm.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/totemMaterials.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/totemRotations.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/totemt1.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/totemt2.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/ionpump.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/castor.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/zdcmaterials.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/lumimaterials.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/zdcrotations.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/lumirotations.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/zdc.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/zdclumi.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/cmszdc.xml");
    geovec.push_back("Geometry/MuonCommonData/data/PhaseII/muonNumbering.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/PhaseII/Pixel10D/trackerStructureTopology.xml");
    geovec.push_back("Geometry/TrackerSimData/data/PhaseII/Pixel10D/trackersens.xml");
    geovec.push_back("Geometry/TrackerRecoData/data/PhaseII/Pixel10D/trackerRecoMaterial.xml");
    geovec.push_back("Geometry/EcalSimData/data/PhaseII/ecalsens.xml");
    geovec.push_back("Geometry/HcalCommonData/data/PhaseII/NoHE/hcalsenspmf.xml");
    geovec.push_back("Geometry/HcalSimData/data/hf.xml");
    geovec.push_back("Geometry/HcalSimData/data/hfpmt.xml");
    geovec.push_back("Geometry/HcalSimData/data/hffibrebundle.xml");
    geovec.push_back("Geometry/HcalSimData/data/CaloUtil.xml");
    geovec.push_back("Geometry/HGCalSimData/data/hgcsens.xml");
    geovec.push_back("Geometry/HGCalSimData/data/hgccons.xml");
    geovec.push_back("Geometry/HGCalSimData/data/hgcProdCuts.xml");
    geovec.push_back("Geometry/MuonSimData/data/PhaseII/muonSens.xml");
    geovec.push_back("Geometry/DTGeometryBuilder/data/dtSpecsFilter.xml");
    geovec.push_back("Geometry/CSCGeometryBuilder/data/cscSpecsFilter.xml");
    geovec.push_back("Geometry/CSCGeometryBuilder/data/cscSpecs.xml");
    geovec.push_back("Geometry/RPCGeometryBuilder/data/PhaseII/RPCSpecs.xml");
    geovec.push_back("Geometry/GEMGeometryBuilder/data/v5/GEMSpecs.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/brmsens.xml");
    geovec.push_back("Geometry/ForwardSimData/data/castorsens.xml");
    geovec.push_back("Geometry/ForwardSimData/data/zdcsens.xml");
    geovec.push_back("Geometry/HcalSimData/data/HcalProdCuts.xml");
    geovec.push_back("Geometry/EcalSimData/data/EcalProdCuts.xml");
    geovec.push_back("Geometry/TrackerSimData/data/PhaseII/Pixel10D/trackerProdCuts.xml");
    geovec.push_back("Geometry/TrackerSimData/data/trackerProdCutsBEAM.xml");
    geovec.push_back("Geometry/MuonSimData/data/PhaseII/muonProdCuts.xml");
    geovec.push_back("Geometry/ForwardSimData/data/CastorProdCuts.xml");
    geovec.push_back("Geometry/ForwardSimData/data/zdcProdCuts.xml");
    geovec.push_back("Geometry/ForwardSimData/data/ForwardShieldProdCuts.xml");
    geovec.push_back("Geometry/CMSCommonData/data/FieldParameters.xml");
    edm::ParameterSet parset; 
    parset.addParameter< vector<string> >("geomXMLFiles", geovec);
    parset.addParameter< string >("rootNodeName", "cms:OCMS");


    // Copy-paste from https://github.com/PFCal-dev/cmssw/blob/CMSSW_6_2_X_SLHC/GeometryReaders/XMLIdealGeometryESSource/src/XMLIdealGeometryESSource.cc 
    cout<<"INFO: Parsing geometry files: ";
    bool userNS = false;
    string rootNodeName(parset.getParameter<string>("rootNodeName"));
    GeometryConfigurationHack geoConfig(parset);

    DDName ddName(rootNodeName);
    DDLogicalPart rootNode(ddName);
    DDRootDef::instance().set(rootNode);
    m_ddcv = new DDCompactView(rootNode);
    DDLParser parser(*m_ddcv); 
    parser.getDDLSAX2FileHandler()->setUserNS(userNS);
    int result = parser.parse(geoConfig);
    if (result != 0) 
    {
        cout<<"Failed!\n";
        return false;
    }
    if( !rootNode.isValid() ){
        cout<<"No valid node named "<<rootNodeName<<"\n";
        return false;
    }
    m_ddcv->lockdown();
    cout<<" OK\n";
    string nameEE("HGCalEESensitive");
    string nameHEF("HGCalHESiliconSensitive");
    cout<<"INFO: Building HGCAL topology: ";
    // Copy-paste from https://github.com/PFCal-dev/cmssw/blob/CMSSW_6_2_X_SLHC/Geometry/HGCalCommonData/plugins/HGCalNumberingInitialization.cc 
    m_hgcdc[3] = new HGCalDDDConstants(*m_ddcv, nameEE);
    m_hgcdc[4] = new HGCalDDDConstants(*m_ddcv, nameHEF);
    if(!m_hgcdc[3] || !m_hgcdc[4])
    {
        cout<<" Failed to build HGCalDDDConstants\n";
        return false;
    }
    // Copy-paste from https://github.com/PFCal-dev/cmssw/blob/CMSSW_6_2_X_SLHC/Geometry/CaloEventSetup/plugins/HGCalTopologyBuilder.cc 
    m_hgctopo[3] = new HGCalTopology(*m_hgcdc[3], HGCEE, false);
    m_hgctopo[4] = new HGCalTopology(*m_hgcdc[4], HGCHEF, false);
    if(!m_hgctopo[3] || !m_hgctopo[4])
    {
        cout<<" Failed!\n";
        return false;
    }
    cout<<"OK\n";

    // Copy-paste from https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_X_SLHC/Geometry/FCalGeometry/plugins/HGCalGeometryESProducer.cc 
    cout<<"INFO: Building HGCAL geometry: ";
    HGCalGeometryLoader builder; 
    m_hgcgeom[3] = builder.build(*m_hgctopo[3]);
    m_hgcgeom[4] = builder.build(*m_hgctopo[4]);
    if(!m_hgcgeom[3] || !m_hgcgeom[4])
    {
        cout<<" Failed!\n";
        return false;
    }
    cout<<"OK\n";

    return true;

}


// temporary workaround to navigate up and down
/*****************************************************************/
std::vector<HGCEEDetId> HGCALNavigator::up(const HGCEEDetId& id, int nz) const
/*****************************************************************/
{
    std::vector<HGCEEDetId> vCells;
    int layer = id.layer()+nz;
    if(layer<=0 || layer>30) return vCells;

    HGCalDDDConstants* hgcdc = m_hgcdc.at(id.subdet());
    HGCalTopology* hgctopo = m_hgctopo.at(id.subdet());

    std::pair<float,float> xy = hgcdc->locateCell(id.cell(), id.layer(), id.subsector(), true);
    std::pair<int,int>     kcell = hgcdc->assignCell(xy.first, xy.second, layer, id.subsector(), true);
    int cell = kcell.second;
    if(!hgctopo->valid(HGCEEDetId(id.subdet(), id.zside(), layer, id.sector(), id.subsector(), cell))) return vCells;

    vCells.push_back( HGCEEDetId(id.subdet(), id.zside(), layer, id.sector(), id.subsector(), cell) );
    return vCells;
}

/*****************************************************************/
std::vector<HGCEEDetId> HGCALNavigator::down(const HGCEEDetId& id, int nz) const
/*****************************************************************/
{
    std::vector<HGCEEDetId> vCells;
    int layer = id.layer()-nz;
    if(layer<=0 || layer>30) return vCells;

    HGCalDDDConstants* hgcdc = m_hgcdc.at(id.subdet());
    HGCalTopology* hgctopo = m_hgctopo.at(id.subdet());

    std::pair<float,float> xy = hgcdc->locateCell(id.cell(), id.layer(), id.subsector(), true);
    std::pair<int,int>     kcell = hgcdc->assignCell(xy.first, xy.second, layer, id.subsector(), true);
    int cell = kcell.second;
    if(!hgctopo->valid(HGCEEDetId(id.subdet(), id.zside(), layer, id.sector(), id.subsector(), cell))) return vCells;

    vCells.push_back( HGCEEDetId(id.subdet(), id.zside(), layer, id.sector(), id.subsector(), cell) );
    return vCells;
}

/*****************************************************************/
std::vector<HGCEEDetId> HGCALNavigator::upProj(const HGCEEDetId& id, int nz, double refEta, double refPhi) const
/*****************************************************************/
{
    std::vector<HGCEEDetId> vCellsProj;
    int layer = id.layer()+nz;
    if(layer<=0 || layer>30) return vCellsProj;

    HGCalTopology* hgctopo = m_hgctopo.at(id.subdet());
    HGCalGeometry* hgcgeom = m_hgcgeom.at(id.subdet());

    std::vector<HGCEEDetId> vCells = up(id, nz);
    if(vCells.size()==0) // try simple navigation
    {
        if( !hgctopo->valid(HGCEEDetId(id.subdet(), id.zside(), layer, id.sector(), id.subsector(), id.cell())) ) return vCellsProj; // Simple navigation doesn't work neither
        vCells.push_back( HGCEEDetId(id.subdet(), id.zside(), layer, id.sector(), id.subsector(), id.cell()) );   
    }

    // Try to find the cell in the up layer with the smallest DeltaR, among a 3x9 grid
    std::vector<DetId> ids;
    // center
    ids.push_back( vCells[0]); // 0
    // first ring
    ids.push_back( hgctopo->offsetBy(ids[0], -1,  1) ); // 1
    ids.push_back( hgctopo->offsetBy(ids[0],  0,  1) ); // 2
    ids.push_back( hgctopo->offsetBy(ids[0],  1,  1) ); // 3
    ids.push_back( hgctopo->offsetBy(ids[0],  1,  0) ); // 4
    ids.push_back( hgctopo->offsetBy(ids[0],  1, -1) ); // 5
    ids.push_back( hgctopo->offsetBy(ids[0],  0, -1) ); // 6
    ids.push_back( hgctopo->offsetBy(ids[0], -1, -1) ); // 7
    ids.push_back( hgctopo->offsetBy(ids[0], -1,  0) ); // 8
    // second ring
    //ids.push_back( hgctopo->offsetBy(ids[0], -2,  2) ); // 9
    ids.push_back( hgctopo->offsetBy(ids[0], -1,  2) ); // 10
    ids.push_back( hgctopo->offsetBy(ids[0],  0,  2) ); // 11
    ids.push_back( hgctopo->offsetBy(ids[0],  1,  2) ); // 12
    //ids.push_back( hgctopo->offsetBy(ids[0],  2,  2) ); // 13
    //ids.push_back( hgctopo->offsetBy(ids[0],  2,  1) ); // 14
    //ids.push_back( hgctopo->offsetBy(ids[0],  2,  0) ); // 15
    //ids.push_back( hgctopo->offsetBy(ids[0],  2, -1) ); // 16
    //ids.push_back( hgctopo->offsetBy(ids[0],  2, -2) ); // 17
    ids.push_back( hgctopo->offsetBy(ids[0],  1, -2) ); // 18
    ids.push_back( hgctopo->offsetBy(ids[0],  0, -2) ); // 19
    ids.push_back( hgctopo->offsetBy(ids[0], -1, -2) ); // 20
    //ids.push_back( hgctopo->offsetBy(ids[0], -2, -2) ); // 21
    //ids.push_back( hgctopo->offsetBy(ids[0], -2, -1) ); // 22
    //ids.push_back( hgctopo->offsetBy(ids[0], -2,  0) ); // 23
    //ids.push_back( hgctopo->offsetBy(ids[0], -2,  1) ); // 24
    // third ring
    //ids.push_back( hgctopo->offsetBy(ids[0], -3,  3) ); // 25
    //ids.push_back( hgctopo->offsetBy(ids[0], -2,  3) ); // 26
    ids.push_back( hgctopo->offsetBy(ids[0], -1,  3) ); // 27
    ids.push_back( hgctopo->offsetBy(ids[0],  0,  3) ); // 28
    ids.push_back( hgctopo->offsetBy(ids[0],  1,  3) ); // 29
    //ids.push_back( hgctopo->offsetBy(ids[0],  2,  3) ); // 30
    //ids.push_back( hgctopo->offsetBy(ids[0],  3,  3) ); // 31
    //ids.push_back( hgctopo->offsetBy(ids[0],  3,  2) ); // 32
    //ids.push_back( hgctopo->offsetBy(ids[0],  3,  1) ); // 33
    //ids.push_back( hgctopo->offsetBy(ids[0],  3,  0) ); // 34
    //ids.push_back( hgctopo->offsetBy(ids[0],  3, -1) ); // 35
    //ids.push_back( hgctopo->offsetBy(ids[0],  3, -2) ); // 36
    //ids.push_back( hgctopo->offsetBy(ids[0],  3, -3) ); // 37
    //ids.push_back( hgctopo->offsetBy(ids[0],  2, -3) ); // 38
    ids.push_back( hgctopo->offsetBy(ids[0],  1, -3) ); // 39
    ids.push_back( hgctopo->offsetBy(ids[0],  0, -3) ); // 40
    ids.push_back( hgctopo->offsetBy(ids[0], -1, -3) ); // 41
    //ids.push_back( hgctopo->offsetBy(ids[0], -2, -3) ); // 42
    //ids.push_back( hgctopo->offsetBy(ids[0], -3, -3) ); // 43
    //ids.push_back( hgctopo->offsetBy(ids[0], -3, -2) ); // 44
    //ids.push_back( hgctopo->offsetBy(ids[0], -3, -1) ); // 45
    //ids.push_back( hgctopo->offsetBy(ids[0], -3,  0) ); // 46
    //ids.push_back( hgctopo->offsetBy(ids[0], -3,  1) ); // 47
    //ids.push_back( hgctopo->offsetBy(ids[0], -3,  2) ); // 48
    // fourth ring
    ids.push_back( hgctopo->offsetBy(ids[0], -1,  4) );
    ids.push_back( hgctopo->offsetBy(ids[0],  0,  4) );
    ids.push_back( hgctopo->offsetBy(ids[0],  1,  4) );
    ids.push_back( hgctopo->offsetBy(ids[0],  1, -4) );
    ids.push_back( hgctopo->offsetBy(ids[0],  0, -4) );
    ids.push_back( hgctopo->offsetBy(ids[0], -1, -4) );


    if(refEta==999. && refPhi==999.)
    {
        const GlobalPoint pos ( std::move( hgcgeom->getPosition(id ) ) ); 
        refEta = (double)pos.eta();
        refPhi = (double)pos.phi();
    }
    std::vector<GlobalPoint> poss;
    std::vector<double> drs;
    for(unsigned i=0;i<ids.size();i++)
    {
        if(ids[0]==DetId(0)) 
        {
            drs.push_back(99999.);
            continue;
        }
        poss.push_back( std::move( hgcgeom->getPosition(ids[i]) ) );
        double deta = (double)poss[i].eta()-refEta;
        double dphi = TVector2::Phi_mpi_pi((double)poss[i].phi()-refPhi);
        double dr = sqrt(deta*deta + dphi*dphi);
        drs.push_back(dr);
    }


    double mindr = 99999.;
    unsigned minIndex = 0;
    for(unsigned i=0;i<drs.size();i++)
    {
        //cout<<"Cell #"<<i<<": DR="<<drs[i]<<"\n";
        if(drs[i]<mindr)
        {
            mindr = drs[i];
            minIndex = i;
        }
    }
    //cout<<"Min DR="<<mindr<<" ("<<minIndex<<")\n";
    vCellsProj.push_back(ids[minIndex]);

    return vCellsProj;
}

/*****************************************************************/
std::vector<HGCEEDetId> HGCALNavigator::downProj(const HGCEEDetId& id, int nz, double refEta, double refPhi) const
/*****************************************************************/
{
    return upProj(id, -nz, refEta, refPhi);
}
