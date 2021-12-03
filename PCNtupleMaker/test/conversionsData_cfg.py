import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

#from Configuration.Eras.Era_Run2_2017_cff import Run2_2017
#from Configuration.Eras.Modifier_run2_nanoAOD_92X_cff import run2_nanoAOD_92X
#process = cms.Process('NANO',Run2_2017,run2_nanoAOD_92X)
process = cms.Process("miniflatntuple")
options = VarParsing.VarParsing('analysis')
#options.inputFiles = 
options.outputFile = "defaultout.root"
options.maxEvents = 100
options.parseArguments()

process.load("FWCore.MessageLogger.MessageLogger_cfi")

"""
service = MessageLogger {
      vstring destinations =  {"debug.txt"}
      PSet debug.txt = { 
         string threshold = "DEBUG" 
         PSet DEBUG = { int32 limit = -1}    
      }
      vstring debugModules = { "stuff"} 
}
"""


process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['phase1_2017_realistic']

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.maxEvents = cms.untracked.PSet( input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = cms.untracked.vstring(options.inputFiles)


#from MaterialStudy.PCNtupleMaker.genparticles_cff import *
#process.load("MaterialStudy.PCNtupleMaker.genparticles_cff")
#from MaterialStudy.PCNtupleMaker.simparticles_cff import *
#process.load("MaterialStudy.PCNtupleMaker.simparticles_cff")



from MaterialStudy.PCNtupleMaker.globals_cff import *
process.load("MaterialStudy.PCNtupleMaker.globals_cff")

from MaterialStudy.PCNtupleMaker.vertices_cff import *
process.load("MaterialStudy.PCNtupleMaker.vertices_cff")

from MaterialStudy.PCNtupleMaker.convs_cff import *
process.load("MaterialStudy.PCNtupleMaker.convs_cff")



#process.load("PhysicsTools.NanoAOD.nano_cff")

#from KUSoftMVA.MuonAnalysis.genparticles_cff import *
#process.load("KUSoftMVA.MuonAnalysis.genparticles_cff")

#from KUSoftMVA.MuonAnalysis.muons_cff import *
#process.load("KUSoftMVA.MuonAnalysis.muons_cff")

#from KUSoftMVA.MuonAnalysis.qparts_cff import *
#process.load("KUSoftMVA.MuonAnalysis.qparts_cff")

#process.Path = cms.Path(genParticleSequence+ genParticleTables + muonSequence + muonTables + muonMC +qTables + qMC )
#process.Path = cms.Path( qTables )
#process.Path = cms.Path( globalTables + vertexTables + convTables )
#process.Path = cms.Path(genParticleSequence + genParticleTable + vertexTables * testseq)
#process.Path = cms.Path( genParticleTable + vertexTables + testseq)
#process.Path = cms.Path(genParticleTable + vertexTables * testseq)
#process.Path = cms.Path(simParticleTables * vertexTables * convTables)
#process.Path = cms.Path( vertexTables * convTables)
process.Path = cms.Path( vertexTables * testseq )
#process.Path = cms.Path(process.nanoSequenceMC)
#for data:
#process.nanoPath = cms.Path(process.nanoSequence)
#process.GlobalTag.globaltag = autoCond['run2_data']

process.out = cms.OutputModule("NanoAODOutputModule", 
#process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = cms.untracked.vstring(
        'drop *',
        "keep nanoaodFlatTable_*Table_*_*",     # event data
       # "keep edmTriggerResults_*_*_*",  # event data
       # "keep String_*_genModel_*",  # generator model data
       # "keep nanoaodMergeableCounterTable_*Table_*_*", # accumulated per/run or per/lumi data
       # "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    ),
  #  outputCommands = process.NanoAODEDMEventContent.outputCommands,
   #compressionLevel = cms.untracked.int32(9),
    #compressionAlgorithm = cms.untracked.string("LZMA"),

)
#process.out1 = cms.OutputModule("NanoAODOutputModule",
#    fileName = cms.untracked.string('lzma.root'),
 #   outputCommands = process.NanoAODEDMEventContent.outputCommands,
 #   compressionLevel = cms.untracked.int32(9),
 #   compressionAlgorithm = cms.untracked.string("LZMA"),

#)
#process.end = cms.EndPath(process.out+process.out1)  '
process.end = cms.EndPath(process.out)
