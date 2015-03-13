from AnhimaBatchLauncher import AnhimaBatchLauncher
import glob

batch = AnhimaBatchLauncher()
batch.name = "Identification_Zee_140PU"
batch.exe = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/bin/slc6_amd64_gcc472/identification.exe"
batch.baseDir = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/"
batch.inputFiles.extend(glob.glob("/data_CMS/cms/sauvan/HGCAL/Zee_140PU/6_2_0_SLHC20_15.03.10/*/Zee_140PU_*.root"))
batch.tree = "HGC"
batch.outputDirectory = "/data_CMS/cms/sauvan/HGCAL/histos/TriggerIdentification/Zee_140PU_HalfLayers/"
batch.outputFile = "identification_zee_140PU.root"
batch.histoParameters = "../histos.par"
batch.histoTag = "Histos"
batch.nFilesPerJob = 3

batch.scram_arch = "slc6_amd64_gcc472"
batch.cmsswDir   = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/"
batch.queue      = "cms"

## e/g algo
batch.additionalParameters["PileupParams"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/thresholdParameters.txt"
batch.additionalParameters["ClusterSize"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/clusterSizes.txt"
batch.additionalParameters["PileupSubtractionFile"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_EG_halfLayers_PUsubtraction_eta_nhits_results.root"
batch.additionalParameters["ThresholdCorrectionFile"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_EG_halfLayers_thresholdCorr_E_eta_results.root"
batch.additionalParameters["SuperClusterCorrectionFile"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_EG_halfLayers_superClusterCorr_eta_phi_ncl_results.root"
batch.additionalParameters["UseHalfLayers"] = "true"

batch.additionalParameters["BDT.N"] = "4"
batch.additionalParameters["BDT.1.Name"] = "BDTG"
batch.additionalParameters["BDT.1.File"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/ele_pu_firstLayerMaxLayerEratio.xml"
batch.additionalParameters["BDT.2.Name"] = "BDTG_MinBias_HalfLayers"
batch.additionalParameters["BDT.2.File"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/ele_minbias_HalfLayers_firstLayerMaxLayerEratio.xml"
batch.additionalParameters["BDT.3.Name"] = "BDTG_QCD"
batch.additionalParameters["BDT.3.File"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/ele_qcd_firstLayerMaxLayerEratio.xml"
batch.additionalParameters["BDT.4.Name"] = "BDTG_QCD_HalfLayers"
batch.additionalParameters["BDT.4.File"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/ele_qcd_HalfLayers_firstLayerMaxLayerEratio.xml"

## jet algo
batch.additionalParameters["JetPileupSubtraction.N"] = "2"
batch.additionalParameters["JetPileupSubtraction.1.Name"] = "0.2"
batch.additionalParameters["JetPileupSubtraction.1.File"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_Jets_PUsubtraction_eta_nhits_results.root"
batch.additionalParameters["JetPileupSubtraction.2.Name"] = "0.1"
batch.additionalParameters["JetPileupSubtraction.2.File"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_Jets_Core_PUsubtraction_eta_nhits_results.root"
#
batch.additionalParameters["JetThresholdCorrection.N"] = "2"
batch.additionalParameters["JetThresholdCorrection.1.Name"] = "0.2"
batch.additionalParameters["JetThresholdCorrection.1.File"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_Jets_thresholdCorr_E_eta_results.root"
batch.additionalParameters["JetThresholdCorrection.2.Name"] = "0.1"
batch.additionalParameters["JetThresholdCorrection.2.File"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_Jets_Core_thresholdCorr_E_eta_results.root"
#
batch.additionalParameters["JetTruthCorrectionFile"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_Jets_truthCorr_ET_eta_ooc_results.root"


### analysis
batch.additionalParameters["Sample"] = "signal_zee"
batch.additionalParameters["EtCut"] = "10."

batch.additionalParameters["ProjectionMapping"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/projectionMapping.txt"



