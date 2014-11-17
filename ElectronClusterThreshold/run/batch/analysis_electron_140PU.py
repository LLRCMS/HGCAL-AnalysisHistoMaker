from AnhimaBatchLauncher import AnhimaBatchLauncher
import glob

batch = AnhimaBatchLauncher()
batch.name = "Electron_140PU"
batch.exe = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/bin/slc6_amd64_gcc472/seeding.exe"
batch.baseDir = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/"
batch.inputFiles.extend(glob.glob("/data_CMS/cms/sauvan/HGCAL/ElectronGun_140PU/CMSSW_6_2_0_SLHC19_14.10.28/Electron_140PU_*.root"))
batch.tree = "HGC"
batch.outputDirectory = "/data_CMS/cms/sauvan/HGCAL/histos/StudySeeding/Electron_140PU/"
batch.outputFile = "electron_140PU.root"
batch.histoParameters = "../histos.par"
batch.histoTag = "Histos"
batch.nFilesPerJob = 5

batch.scram_arch = "slc6_amd64_gcc472"
batch.cmsswDir   = "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/"
batch.queue      = "cms"

batch.additionalParameters["Sample"] = "electronPU"



