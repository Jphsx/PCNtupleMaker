import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *



##################### User floats producers, selectors ##########################


##################### Tables for final output and docs ##########################
vertexTable = cms.EDProducer("VertexTableProducer",
#    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    goodPvCut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"), 
   # svSrc = cms.InputTag("slimmedSecondaryVertices"),
   # svCut = cms.string(""),
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(3),
    pvName = cms.string("PV"),
   # svName = cms.string("SV"),
   # svDoc  = cms.string("secondary vertices from IVF algorithm"),
)



pvCandidateTable =  cms.EDProducer("SimpleVertexTableProducer",
   
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(position().Z() ) <= 24 && position().Rho() <= 2"),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("PV"),
    doc = cms.string("all Primary Vertices after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), 
    variables = cms.PSet(
        X   = Var("position().X()", float, doc = "secondary vertex X position, in cm",precision=10),
        Y   = Var("position().Y()", float, doc = "secondary vertex Y position, in cm",precision=10),
        Z   = Var("position().Z()", float, doc = "secondary vertex Z position, in cm",precision=14),
        ndof   = Var("ndof()", float, doc = "number of degrees of freedom",precision=8),
        normalizedChi2   = Var("normalizedChi2()", float, doc = "reduced chi2, i.e. chi/ndof",precision=8),
	chi2 = Var("chi2()", float, doc ="chi squared"),
	cov_00 = Var("covariance(0,0)", float, doc = " pv covrariance element xx "),
	cov_01 = Var("covariance(0,1)", float, doc = " pv covrariance element xy "),
	cov_02 = Var("covariance(0,2)", float, doc = " pv covrariance element xz "),
	cov_11 = Var("covariance(1,1)", float, doc = " pv covrariance element yy "),
	cov_12 = Var("covariance(1,2)", float, doc = " pv covrariance element yz "),
	cov_22 = Var("covariance(2,2)", float, doc = " pv covrariance element zz "),
		

    ),
)
#pvCandidateTable.variables.pt.precision=10
#pvCandidateTable.variables.phi.precision=12


#before cross linking
#vertexSequence = cms.Sequence()
#after cross linkining
#vertexTables = cms.Sequence( vertexTable+svCandidateTable)


#i may need to do some crosslinking
vertexTables = cms.Sequence( pvCandidateTable )
