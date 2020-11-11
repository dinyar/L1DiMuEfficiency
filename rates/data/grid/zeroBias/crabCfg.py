from WMCore.Configuration import Configuration

config = Configuration()

requestName = 'L1Ntuple-Run2018D-320838-321479-ZeroBias-emtfTracks'
dataset     = '/ZeroBias/Run2018D-PromptReco-v2/AOD'
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
config.Data.unitsPerJob = 15
config.Data.outLFNDirBase = output
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-321777_13TeV_PromptReco_Collisions18_JSON_MuonPhys.txt'
# config.Data.runRange = '280017-280385' # 4. September to 9. September for G run.
# config.Data.runRange = '279993-280017' # 3. September to 4. September (?)
# config.Data.runRange = '279588-279991' # 25. August to 3. September
# config.Data.runRange = '279588-279887' # 25. August to 3. September
config.Data.runRange = '320838-321479' # 4. August to 20. August
config.Data.ignoreLocality = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T2_CH_CERN','T2_FR_*','T2_IT_*','T2_DE_*']

