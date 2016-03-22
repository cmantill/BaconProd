
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'XX-LABEL-XX'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'XX-CONFIG-XX'
config.JobType.outputFiles = ['Output.root']

config.Data.inputDataset = 'XX-DATASET-XX'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 20000
config.Data.outLFNDirBase = '/store/group/cmst3/group/monojet/production/singleb/XX-LABEL-XX/'
config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3'

config.Site.storageSite = 'T2_CH_CERN'
