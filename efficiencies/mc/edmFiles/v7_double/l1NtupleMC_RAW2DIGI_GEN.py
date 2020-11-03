# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --conditions 80X_mcRun2_asymptotic_v14 -s L1REPACK:FullMC,RAW2DIGI --mc --no_output --no_exec -n 1000 --era=Run2_2016 --filein root://xrootd-cms.infn.it//store/mc/RunIISpring16DR80/SingleMuMinusFlatPt3To70_EtaPhiRestricted/GEN-SIM-RAW/NoPURAW_NZS_withHLT_80X_mcRun2_asymptotic_v14-v1/50000/001ABF88-BD5C-E611-A4F4-FA163EC5F1D1.root --processName reL1T --python_filename l1NtupleMC_RAW2DIGI_GEN.py --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMUGEN_MC
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('reL1T',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1EmulatorRepack_FullMC_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/mc/RunIISpring16DR80/SingleMuMinusFlatPt3To70_EtaPhiRestricted/GEN-SIM-RAW/NoPURAW_NZS_withHLT_80X_mcRun2_asymptotic_v14-v1/50000/001ABF88-BD5C-E611-A4F4-FA163EC5F1D1.root'),
    fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/mc/RunIISpring16DR80/JPsiToMuMu_Pt20to120_EtaPhiRestricted-pythia8-gun/GEN-SIM-RAW/NoPURAW_NZS_withHLT_80X_mcRun2_asymptotic_v14-v1/40000/00235558-C25C-E611-BA58-20CF3027A5A7.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('--conditions nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_v14', '')

# Path and EndPath definitions
process.L1RePack_step = cms.Path(process.SimL1Emulator)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.L1RePack_step,process.raw2digi_step,process.endjob_step)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleRAWEMUGEN_MC 

#call to customisation function L1NtupleRAWEMUGEN_MC imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleRAWEMUGEN_MC(process)

# End of customisation functions

