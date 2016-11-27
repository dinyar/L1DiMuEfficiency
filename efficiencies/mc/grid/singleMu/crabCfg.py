from WMCore.Configuration import Configuration

config = Configuration()

requestName = 'L1Ntuple-SingleMuMinusFlatPt3To70_EtaPhiRestricted-tuning_v1'
dataset     = '/SingleMuMinusFlatPt3To70_EtaPhiRestricted/RunIISpring16DR80-NoPURAW_NZS_withHLT_80X_mcRun2_asymptotic_v14-v1/GEN-SIM-RAW'
splitting   = 'FileBased'
output      = '/store/user/dinyar/cancel_out_studies/'

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = requestName
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = '../../l1NtupleMC_RAW2DIGI_GEN.py'
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
config.Data.splitting = splitting
config.Data.useParent = False
config.Data.unitsPerJob = 6
config.Data.outLFNDirBase = output
config.Data.ignoreLocality = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

