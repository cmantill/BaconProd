from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'XX-LABEL-XX'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'XX-CONFIG-XX'
config.JobType.outputFiles = ['Output.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 5000
#config.Data.lumiMask = 'notFinishedLumis.json'
#config.Data.lumiMask = 'lumimask.json'
config.Data.inputDataset = 'XX-DATASET-XX'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 35000
#config.Data.unitsPerJob = 10000
#config.Data.unitsPerJob = 5000
config.Data.outLFNDirBase = '/store/group/lpcbacon/15/XX-LABEL-XX/'
config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3'
config.Site.storageSite = 'T3_US_FNALLPC'
