
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Fall15_GluGlu_HToInvisible_M600_13TeV_powheg_pythia8'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'makingBacon_MC_25ns_MINIAOD.py'
config.JobType.outputFiles = ['Output.root']

config.Data.inputDataset = '/GluGlu_HToInvisible_M600_13TeV_powheg_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 20000
config.Data.outLFNDirBase = '/store/group/cmst3/group/monojet/production/08/Fall15_GluGlu_HToInvisible_M600_13TeV_powheg_pythia8/'
config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3'

config.Site.storageSite = 'T2_CH_CERN'
