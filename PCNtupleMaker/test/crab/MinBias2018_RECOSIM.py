from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()


#pyCfgParams = ['outputFile=Run2018D.root']

config.General.requestName = 'DPG_MinBias2018'
config.General.workArea = '../crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

#config.JobType.outputFiles      = ['Run2018C.root']


config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/home/t3-ku/janguian/CMSSW_10_6_11_patch1/src/MaterialStudy/PCNtupleMaker/test/conversionsMC_RECOSIM_cfg.py'

config.Data.inputDataset = '/minBias2018/sroychow-tkmaterial_y2018mcRECO-f2bc6d20807befae26d85e37f259503a/USER'
#config.Data.useParent = True
config.Data.inputDBS = 'phys03'

#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

#config.Data.lumiMask = ' Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt


#config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 50000
#config.Data.unitsPerJob = 50
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 5



config.Data.outLFNDirBase = '/store/user/jsingera/DPG/PC/MinBias2018'
config.Data.publication = True
config.Data.outputDatasetTag = 'DPG_MinBias2018'

config.Site.storageSite = 'T2_US_Nebraska'
