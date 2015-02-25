from AnhimaBatchLauncher import AnhimaBatchLauncher
import glob

batch = AnhimaBatchLauncher()
batch.name = "Electron_140PU"
batch.exe = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/bin/slc6_amd64_gcc472/calibration.exe"
batch.baseDir = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/"
#batch.inputFiles.extend(glob.glob("/data_CMS/cms/sauvan/HGCAL/ElectronGun_140PU/CMSSW_6_2_0_SLHC20_14.11.27/Electron_140PU_*.root"))
batch.inputFiles.extend(glob.glob("/data_CMS/cms/sauvan/HGCAL/ElectronGun_140PU/CMSSW_6_2_0_SLHC20_pt20to50_14.11.29/Electron_pt20to50_140PU_0_*.root"))
batch.tree = "HGC"
batch.outputDirectory = "/data_CMS/cms/sauvan/HGCAL/histos/TriggerCalibration/Electron_140PU/"
batch.outputFile = "electron_140PU.root"
batch.histoParameters = "../histos.par"
batch.histoTag = "Histos"
batch.nFilesPerJob = 2

batch.scram_arch = "slc6_amd64_gcc472"
batch.cmsswDir   = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/"
batch.queue      = "cms"


## e/g algo
batch.additionalParameters["PileupParams"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/thresholdParameters.txt"
batch.additionalParameters["ClusterSize"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/clusterSizes.txt"
batch.additionalParameters["PileupSubtractionFile"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_EG_PUsubtraction_eta_nhits_results.root"
batch.additionalParameters["ThresholdCorrectionFile"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_EG_thresholdCorr_E_eta_results.root"
batch.additionalParameters["SuperClusterCorrectionFile"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_EG_superClusterCorr_eta_phi_ncl_results.root"
batch.additionalParameters["UseHalfLayers"] = "false"

batch.additionalParameters["ProjectionMapping"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/projectionMapping.txt"

