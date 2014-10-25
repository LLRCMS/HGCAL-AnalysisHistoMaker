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
    delete m_hgcdc;
    delete m_hgctopo;
}



/*****************************************************************/
bool HGCALNavigator::initialize()
/*****************************************************************/
{
    // xml geometry files defined in https://github.com/PFCal-dev/cmssw/blob/CMSSW_6_2_X_SLHC/Geometry/CMSCommonData/python/cmsExtendedGeometry2023HGCalMuonXML_cfi.py 
    vector<string> geovec;
    geovec.push_back("Geometry/CMSCommonData/data/materials.xml");
    geovec.push_back("Geometry/CMSCommonData/data/rotations.xml");
    geovec.push_back("Geometry/CMSCommonData/data/extend/cmsextent.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cms.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cmsMother.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cmsTracker.xml");
    geovec.push_back("Geometry/CMSCommonData/data/PhaseII/caloBase.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cmsCalo.xml");
    geovec.push_back("Geometry/CMSCommonData/data/muonBase.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cmsMuon.xml");
    geovec.push_back("Geometry/CMSCommonData/data/mgnt.xml");
    geovec.push_back("Geometry/CMSCommonData/data/beampipe.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cmsBeam.xml");
    geovec.push_back("Geometry/CMSCommonData/data/muonMB.xml");
    geovec.push_back("Geometry/CMSCommonData/data/muonMagnet.xml");
    geovec.push_back("Geometry/CMSCommonData/data/cavern.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdMaterials.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdCommon.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdPlaq.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdPlaq1x2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdPlaq1x5.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdPlaq2x3.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdPlaq2x4.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdPlaq2x5.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdPanelBase.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdPanel.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdBlade.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdNipple.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdDisk.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwdCylinder.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixfwd.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixbarmaterial.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixbarladder.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixbarladderfull.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixbarladderhalf.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixbarlayer.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixbarlayer0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixbarlayer1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixbarlayer2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/pixbar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibtidcommonmaterial.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibmaterial.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibmodpar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibmodule0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibmodule0a.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibmodule0b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibmodule2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstringpar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring0ll.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring0lr.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring0ul.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring0ur.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring1ll.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring1lr.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring1ul.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring1ur.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring2ll.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring2lr.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring2ul.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring2ur.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring3ll.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring3lr.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring3ul.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring3ur.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibstring3.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tiblayerpar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tiblayer0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tiblayer1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tiblayer2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tiblayer3.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tib.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidmaterial.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidmodpar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidmodule0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidmodule0r.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidmodule0l.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidmodule1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidmodule1r.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidmodule1l.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidmodule2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidringpar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidring0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidring0f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidring0b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidring1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidring1f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidring1b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidring2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tid.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidf.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tidb.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibtidservices.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibtidservicesf.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tibtidservicesb.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobmaterial.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobmodpar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobmodule0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobmodule2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobmodule4.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrodpar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod0c.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod0l.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod0h.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod1l.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod1h.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod2c.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod2l.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod2h.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod3l.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod3h.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod3.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod4c.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod4l.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod4h.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod4.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod5l.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod5h.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tobrod5.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tob.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmaterial.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodpar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule0r.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule0s.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule1r.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule1s.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule3.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule4.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule4r.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule4s.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule5.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecmodule6.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetpar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring1.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring2.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring3.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring4.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring5.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring6.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring0f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring1f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring2f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring3f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring4f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring5f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring6f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring0b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring1b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring2b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring3b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring4b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring5b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecring6b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetalf.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetalb.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetal0.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetal0f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetal0b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetal3.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetal3f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetal3b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetal6f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetal6b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetal8f.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecpetal8b.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecwheel.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecwheela.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecwheelb.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecwheelc.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecwheeld.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecwheel6.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecservices.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tecbackplate.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tec.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackermaterial.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/tracker.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackerpixbar.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackerpixfwd.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackertibtidservices.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackertib.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackertid.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackertob.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackertec.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackerbulkhead.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackerother.xml");
    geovec.push_back("Geometry/EcalCommonData/data/PhaseII/eregalgo.xml");
    geovec.push_back("Geometry/EcalCommonData/data/ebalgo.xml");
    geovec.push_back("Geometry/EcalCommonData/data/ebcon.xml");
    geovec.push_back("Geometry/EcalCommonData/data/ebrot.xml");
    geovec.push_back("Geometry/EcalCommonData/data/eecon.xml");
    geovec.push_back("Geometry/EcalCommonData/data/ectkcable.xml");
    geovec.push_back("Geometry/EcalCommonData/data/PhaseII/esalgo.xml");
    geovec.push_back("Geometry/EcalCommonData/data/PhaseII/escon.xml");
    geovec.push_back("Geometry/CMSCommonData/data/eta3/etaMax.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalrotations.xml");
    geovec.push_back("Geometry/HcalCommonData/data/PhaseII/NoHE/hcalalgo.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalcablealgo.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalbarrelalgo.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalouteralgo.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalforwardalgo.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalforwardfibre.xml");
    geovec.push_back("Geometry/HcalCommonData/data/hcalforwardmaterial.xml");
    geovec.push_back("Geometry/HGCalCommonData/data/v5/hgcal.xml");
    geovec.push_back("Geometry/HGCalCommonData/data/v5/hgcalEE.xml");
    geovec.push_back("Geometry/HGCalCommonData/data/v5/hgcalHEsil.xml");
    geovec.push_back("Geometry/HGCalCommonData/data/v5/hgcalHEsci.xml");
    geovec.push_back("Geometry/HGCalCommonData/data/v5/hgcalCons.xml");
    geovec.push_back("Geometry/MuonCommonData/data/mbCommon.xml");
    geovec.push_back("Geometry/MuonCommonData/data/mb1.xml");
    geovec.push_back("Geometry/MuonCommonData/data/mb2.xml");
    geovec.push_back("Geometry/MuonCommonData/data/mb3.xml");
    geovec.push_back("Geometry/MuonCommonData/data/mb4.xml");
    geovec.push_back("Geometry/MuonCommonData/data/muonYoke.xml");
    geovec.push_back("Geometry/MuonCommonData/data/mf.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/forward.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/forwardshield.xml");
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
    geovec.push_back("Geometry/ForwardCommonData/data/plt.xml");   
    geovec.push_back("Geometry/ForwardCommonData/data/zdc.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/zdclumi.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/cmszdc.xml");
    geovec.push_back("Geometry/FP420CommonData/data/fp420.xml");
    geovec.push_back("Geometry/FP420CommonData/data/zzzrectangle.xml");
    geovec.push_back("Geometry/FP420CommonData/data/materialsfp420.xml");
    geovec.push_back("Geometry/FP420CommonData/data/FP420Rot.xml");
    geovec.push_back("Geometry/FP420CommonData/data/cmsfp420.xml");
    geovec.push_back("Geometry/MuonCommonData/data/muonNumbering.xml");
    geovec.push_back("Geometry/TrackerCommonData/data/trackerStructureTopology.xml");
    geovec.push_back("Geometry/TrackerSimData/data/trackersens.xml");
    geovec.push_back("Geometry/TrackerRecoData/data/trackerRecoMaterial.xml");
    geovec.push_back("Geometry/EcalSimData/data/PhaseII/ecalsens.xml");
    geovec.push_back("Geometry/HcalCommonData/data/PhaseII/NoHE/hcalsens.xml");
    geovec.push_back("Geometry/HcalCommonData/data/PhaseII/NoHE/hcalSimNumbering.xml");
    geovec.push_back("Geometry/HcalCommonData/data/PhaseII/NoHE/hcalRecNumbering.xml");
    geovec.push_back("Geometry/HcalSimData/data/CaloUtil.xml");
    geovec.push_back("Geometry/HcalSimData/data/hf.xml"); 
    geovec.push_back("Geometry/HcalSimData/data/hffibre.xml"); 
    geovec.push_back("Geometry/HGCalSimData/data/hgcsens.xml"); 
    geovec.push_back("Geometry/HGCalSimData/data/hgccons.xml"); 
    geovec.push_back("Geometry/MuonSimData/data/muonSens.xml");
    geovec.push_back("Geometry/DTGeometryBuilder/data/dtSpecsFilter.xml");
    geovec.push_back("Geometry/CSCGeometryBuilder/data/cscSpecsFilter.xml");
    geovec.push_back("Geometry/CSCGeometryBuilder/data/cscSpecs.xml");
    geovec.push_back("Geometry/RPCGeometryBuilder/data/RPCSpecs.xml");
    geovec.push_back("Geometry/ForwardCommonData/data/brmsens.xml");
    geovec.push_back("Geometry/ForwardSimData/data/totemsensT1.xml");
    geovec.push_back("Geometry/ForwardSimData/data/totemsensT2.xml");
    geovec.push_back("Geometry/ForwardSimData/data/castorsens.xml");
    geovec.push_back("Geometry/ForwardSimData/data/pltsens.xml");
    geovec.push_back("Geometry/ForwardSimData/data/zdcsens.xml");
    geovec.push_back("Geometry/FP420SimData/data/fp420sens.xml");
    geovec.push_back("Geometry/HcalSimData/data/HcalProdCuts.xml");
    geovec.push_back("Geometry/EcalSimData/data/EcalProdCuts.xml");
    geovec.push_back("Geometry/HGCalSimData/data/hgcProdCuts.xml"); 
    geovec.push_back("Geometry/TrackerSimData/data/trackerProdCuts.xml");
    geovec.push_back("Geometry/TrackerSimData/data/trackerProdCutsBEAM.xml");
    geovec.push_back("Geometry/MuonSimData/data/muonProdCuts.xml");
    geovec.push_back("Geometry/ForwardSimData/data/CastorProdCuts.xml");
    geovec.push_back("Geometry/ForwardSimData/data/zdcProdCuts.xml");
    geovec.push_back("Geometry/ForwardSimData/data/ForwardShieldProdCuts.xml");
    geovec.push_back("Geometry/FP420SimData/data/FP420ProdCuts.xml");
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
    string name("HGCalEESensitive");
    cout<<"INFO: Building HGCAL topology: ";
    // Copy-paste from https://github.com/PFCal-dev/cmssw/blob/CMSSW_6_2_X_SLHC/Geometry/HGCalCommonData/plugins/HGCalNumberingInitialization.cc 
    m_hgcdc = new HGCalDDDConstants(*m_ddcv, name);
    if(!m_hgcdc)
    {
        cout<<" Failed to build HGCalDDDConstants\n";
        return false;
    }
    // Copy-paste from https://github.com/PFCal-dev/cmssw/blob/CMSSW_6_2_X_SLHC/Geometry/CaloEventSetup/plugins/HGCalTopologyBuilder.cc 
    m_hgctopo = new HGCalTopology(*m_hgcdc, HGCEE, false);
    if(!m_hgctopo)
    {
        cout<<" Failed!\n";
        return false;
    }
    cout<<"OK\n";
    return true;

}


// temporary workaround to navigate up and down
/*****************************************************************/
std::vector<HGCEEDetId> HGCALNavigator::up(const HGCEEDetId& id)
/*****************************************************************/
{
    std::vector<HGCEEDetId> vCells;
    int layer = id.layer()+1;
    if(layer<=0 || layer>30) return vCells;

    std::pair<float,float> xy = m_hgcdc->locateCell(id.cell(), id.layer(), id.subsector(), true);
    std::pair<int,int>     kcell = m_hgcdc->assignCell(xy.first, xy.second, layer, id.subsector(), true);
    int cell = kcell.second;
    vCells.push_back( HGCEEDetId(id.subdet(), id.zside(), layer, id.sector(), id.subsector(), cell) );
    return vCells;
}

/*****************************************************************/
std::vector<HGCEEDetId> HGCALNavigator::down(const HGCEEDetId& id)
/*****************************************************************/
{
    std::vector<HGCEEDetId> vCells;
    int layer = id.layer()-1;
    if(layer<=0 || layer>30) return vCells;

    std::pair<float,float> xy = m_hgcdc->locateCell(id.cell(), id.layer(), id.subsector(), true);
    std::pair<int,int>     kcell = m_hgcdc->assignCell(xy.first, xy.second, layer, id.subsector(), true);
    int cell = kcell.second;
    vCells.push_back( HGCEEDetId(id.subdet(), id.zside(), layer, id.sector(), id.subsector(), cell) );
    return vCells;
}
