from WMCore.Configuration import Configuration

config = Configuration()

requestName = 'L1Ntuple-Run2016G-280017-280385-ZeroBias-aggressiveTuning'
dataset     = '/ZeroBias/Run2016G-23Sep2016-v1/AOD'
# splitting   = 'EventAwareLumiBased'
splitting   = 'LumiBased'
output      = '/store/user/dinyar/cancel_out_studies/'

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = requestName
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = '../../l1NtupleData_RAW2DIGI.py'
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
config.Data.splitting = splitting
config.Data.useParent = True
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = output
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T_MuonPhys.txt'
config.Data.runRange = '280017-280385' # 15. Juni to 25. June for E run. -- 4. September to 9. September for G run.
config.Data.ignoreLocality = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

