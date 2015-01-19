from AnhimaBatchLauncher import AnhimaBatchLauncher
import glob

batch = AnhimaBatchLauncher()
batch.name = "MinBias_140PU"
batch.exe = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/bin/slc6_amd64_gcc472/jets.exe"
batch.baseDir = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/"
batch.inputFiles.extend(glob.glob("/data_CMS/cms/sauvan/HGCAL/MinBias_140PU/6_2_0_SLHC20_15.01.17/*/MinBias_*.root"))
batch.tree = "HGC"
batch.outputDirectory = "/data_CMS/cms/sauvan/HGCAL/histos/StudyJets/MinBias_140PU/"
batch.outputFile = "MinBias_140PU.root"
batch.histoParameters = "../histos.par"
batch.histoTag = "Histos"
batch.nFilesPerJob = 1

batch.scram_arch = "slc6_amd64_gcc472"
batch.cmsswDir   = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/"
batch.queue      = "cms"

batch.additionalParameters["PileupParams"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/thresholdParameters.txt"
batch.additionalParameters["ClusterSize"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/clusterSizes.txt"
batch.additionalParameters["PileupSubtractionFile"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_PUsubtraction_eta_nhits_results.root"
batch.additionalParameters["ThresholdCorrectionFile"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_thresholdCorr_E_eta_results.root"
batch.additionalParameters["SuperClusterCorrectionFile"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/GBR_Trigger_HGCAL_superClusterCorr_eta_phi_ncl_results.root"
batch.additionalParameters["ProjectionMapping"] = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/projectionMapping.txt"



