from AnhimaBatchLauncher import AnhimaBatchLauncher
import glob

batch = AnhimaBatchLauncher()
batch.name = "MinBias_140PU"
batch.exe = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/bin/slc6_amd64_gcc472/pileup.exe"
batch.baseDir = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/"
batch.inputFiles.extend(glob.glob("/data_CMS/cms/sauvan/HGCAL/MinBias_140PU/6_2_0_SLHC19_14.10.28/MinBias_140PU_*.root"))
batch.tree = "HGC"
batch.outputDirectory = "/data_CMS/cms/sauvan/HGCAL/histos/StudyPileup/MinBias_140PU/"
batch.outputFile = "minbias_140PU.root"
batch.histoParameters = "../histos.par"
batch.histoTag = "Histos"
batch.nFilesPerJob = 1

batch.scram_arch = "slc6_amd64_gcc472"
batch.cmsswDir   = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/"
batch.queue      = "cms"



