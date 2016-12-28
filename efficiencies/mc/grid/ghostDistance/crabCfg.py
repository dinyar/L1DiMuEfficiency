from WMCore.Configuration import Configuration

config = Configuration()

requestName = 'L1Ntuple-GhostDistance'
dataset     = '/SingleMuon/Run2016B-MuTau-PromptReco-v2/RAW-RECO'
# splitting   = 'EventAwareLumiBased'
splitting   = 'LumiBased'
output      = '/store/user/dinyar/cancel_out_studies/'

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = requestName
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = '../../l1NtupleRECO_RAW2DIGI.py'
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
config.Data.splitting = splitting
# config.Data.useParent = True
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = output
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt'
config.Data.runRange = '275125'
config.Data.ignoreLocality = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

