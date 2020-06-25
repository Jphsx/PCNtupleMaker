import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *



##################### User floats producers, selectors ##########################


##################### Tables for final output and docs ##########################
simTrackTable = cms.EDProducer("SimpleSimTrackFlatTableProducer",
	src = cms.InputTag("g4SimHits"),
	cut = cms.string(""),
	name= cms.string("SimTrk"),
	doc= cms.string("Sim track collection"),
	singleton = cms.bool(False),
	extension = cms.bool(False),
	variables = cms.PSet(
		pt = Var("momentum.Pt()", float, precision=8, doc="pt of sim track"),
		phi = Var("momentum.Phi()", float, precision=8, doc="phi of sim track"),
		eta = Var("momentum.Eta()", float, precision=8, doc="eta of sim track"),
		pdgId = Var("type()",int,doc="PDG id"),
		charge = Var("charge()",float,doc="charge of sim track"),
		trackId = Var("trackId()",int,doc="unique integer id of sim track"),
		simvtx_Idx = Var("vertIndex()",int,doc="index of origin sim vertex")
	)

)
simVertexTable = cms.EDProducer("SimpleSimVertexFlatTableProducer",
	src = cms.InputTag("g4SimHits"),
	cut = cms.string(""),
	name = cms.string("SimVtx"),
	doc=cms.string("Sim Vertex collection"),
	singleton = cms.bool(False),
	extension = cms.bool(False),
	variables = cms.PSet(
		x = Var("position.x()", float, precision=8, doc="x position of sim vertex"),
		y = Var("position.y()", float, precision=8, doc="y position of sim vertex"),
		z = Var("position.z()", float, precision=8, doc="z postiion of sim vertex"),
		tof = Var("position.t()", float, precision=8, doc="time of flight of sim vertex"),
		simtrk_parent_Id = Var("parentIndex()", int, doc="tiD of the parent sim track WARNING: NOT A VECTOR INDEX"),
		processType = Var("processType()", int, doc="type of process associated with vtx")
	)
)

simParticleTables = cms.Sequence(simTrackTable +simVertexTable)
