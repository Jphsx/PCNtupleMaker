from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'TnP_SingleMuon2017C_zmumu'
config.General.workArea = '../crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True


config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/home/t3-ku/janguian/CMSSW_10_2_5/src/MuonAnalysis/TagAndProbe/test/susy/zmumu/KUTnP_Data.py'

config.Data.inputDataset = '/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD'
config.Data.useParent = True
config.Data.inputDBS = 'global'
#config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 50000
config.Data.outLFNDirBase = '/store/user/jsingera/TnP'
config.Data.publication = True
config.Data.outputDatasetTag = 'TnP_SingleMuon2017C_zmumu'

config.Site.storageSite = 'T2_US_Nebraska'
