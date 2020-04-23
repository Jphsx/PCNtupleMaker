from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()


#pyCfgParams = ['outputFile=Run2018D.root']

config.General.requestName = 'DPG_SingleMuon2018D'
config.General.workArea = '../crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.outputFiles      = ['Run2018D.root']


config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/home/t3-ku/janguian/CMSSW_10_6_11_patch1/src/MaterialStudy/PCNtupleMaker/test/conversions_cfg.py'

config.Data.inputDataset = '/SingleMuon/Run2018D-PromptReco-v2/AOD'
#config.Data.useParent = True
config.Data.inputDBS = 'global'

config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

#config.Data.splitting = 'EventAwareLumiBased'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 50
config.Data.outLFNDirBase = '/store/user/jsingera/DPG/PC/Run2018'
config.Data.publication = True
config.Data.outputDatasetTag = 'DPG_SingleMuon2018D'

config.Site.storageSite = 'T2_US_Nebraska'
