import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Modifier_run2_muon_2016_cff import run2_muon_2016
from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
from Configuration.Eras.Modifier_run2_nanoAOD_94X2016_cff import run2_nanoAOD_94X2016
from Configuration.Eras.Modifier_run2_nanoAOD_94XMiniAODv1_cff import run2_nanoAOD_94XMiniAODv1
from Configuration.Eras.Modifier_run2_nanoAOD_94XMiniAODv2_cff import run2_nanoAOD_94XMiniAODv2
from Configuration.Eras.Modifier_run2_nanoAOD_102Xv1_cff import run2_nanoAOD_102Xv1
from PhysicsTools.NanoAOD.common_cff import *



#TODO build a preselection selector
#convs = cms.EDFilter("ObjectSelector",
#	src = cms.InputTag("allConversions"),
#	cut = cms.string(""),
#	
#	)


convTable = cms.EDProducer("SimpleConversionTableProducer",
    src = cms.InputTag("allConversions"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Conv"),
    doc  = cms.string("all conversions after basic selection"), #(" + convs.cut.value()+")"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons 
    
    variables = cms.PSet( 
	isConverted = Var("isConverted()",bool,doc="Bool flagging objects having track size >0"),
	nTracks = Var("nTracks()", int, doc="Number of tracks= 0,1,2"),
	#MVAout = Var("MVAout()", float, doc="the value  of the TMVA output"),#UNUSED
	pairInvariantMass = Var("pairInvariantMass()", float, doc="if nTracks=2 returns the pair invariant mass. Original tracks are used here"),
	pairCotThetaSeparation = Var("pairCotThetaSeparation()",float, doc="Delta cot(Theta) where Theta is the angle in the (y,z) plane between the two tracks. Original tracks are used"),
	pairMomentum_Px = Var("pairMomentum().X()", float, doc="Conversion tracks momentum from the tracks inner momentum, x-component"),
	pairMomentum_Py = Var("pairMomentum().Y()", float, doc="Conversion tracks momentum from the tracks inner momentum, y-component"),
	pairMomentum_Pz = Var("pairMomentum().Z()", float, doc="Conversion tracks momentum from the tracks inner momentum, z-component"),
	refittedPair4Momentum_Px = Var("refittedPair4Momentum().Px()", float, doc="Conversion track pair 4-momentum from the tracks refitted with vertex constraint, Px"),
	refittedPair4Momentum_Py = Var("refittedPair4Momentum().Py()", float, doc="Conversion track pair 4-momentum from the tracks refitted with vertex constraint, Py"),
	refittedPair4Momentum_Pz = Var("refittedPair4Momentum().Pz()", float, doc="Conversion track pair 4-momentum from the tracks refitted with vertex constraint, Pz"),
	refittedPair4Momentum_E = Var("refittedPair4Momentum().E()", float, doc="Conversion track pair 4-momentum from the tracks refitted with vertex constraint, E"),
	refittedPair4Momentum_M = Var("refittedPair4Momentum().M()", float, doc="Conversion track pair 4-momentum from the tracks refitted with vertex constraint, M"),
	refittedPairMomentum_Px = Var("refittedPairMomentum().X()", float, doc="Conversion tracks momentum from the tracks refitted with vertex constraint, Px"),
     	refittedPairMomentum_Py = Var("refittedPairMomentum().Y()", float, doc="Conversion tracks momentum from the tracks refitted with vertex constraint, Py"),
	refittedPairMomentum_Pz = Var("refittedPairMomentum().Z()", float, doc="Conversion tracks momentum from the tracks refitted with vertex constraint, Pz"),
	EoverP = Var("EoverP()",float,doc="Super Cluster energy divided by track pair momentum if Standard seeding method. If a pointer to two (or more clusters) is stored in the conversion, this method returns the energy sum of clusters divided by the track pair momentum. Track innermost momentum is used here",precision=10),
	EoverPrefittedTracks = Var("EoverPrefittedTracks", float, doc="Track momentum refitted with vertex constraint is used"),
        distOfMinimumApproach = Var("distOfMinimumApproach()", float, doc="Dist of minimum approach between tracks"),
	dPhiTracksAtVtx = Var("dPhiTracksAtVtx()",float, doc="deltaPhi tracks at innermost point"),
#	dPhiTracksAtEcal = Var("dPhiTracksAtEcal()",float,doc="deltaPhi tracks at ECAL"), #this one breaks the code? null maybe?
#	dEtaTracksAtEcal = Var("dEtaTracksAtEcal()",float,doc="deltaEta tracks at ECAL"), #also breaks
	dxy = Var("dxy()",float,doc="transverse impact parameter, computed with respect to given beamspot (0,0,0) and refitted pair momentum"),
	dz = Var("dz()", float, doc="longitudinal impact parameter, computed with respect to given beamspot (0,0,0) and refitted pair momentum"),
	lxy = Var("lxy()", float, doc="transverse decay length, computed with respect to given beamspot (0,0,0) and refitted pair momentum"),
	lz = Var("lz()", float, doc="longitudinal decay length, computed with respect to given beamspot (0,0,0) and refitted pair momentum"),
	zOfPrimaryVertexFromTracks = Var("zOfPrimaryVertexFromTracks()",float,doc="z position of intersection with beamspot in rz plane (possible tilt of beamspot is neglected)"),
	#ecalImpactPosition_X_Tk0 = Var("ecalImpactPosition().at(0).X()",float,doc="positions of the track extrapolation at the ECAL front face, track index 0"), #breaks code
	tracksSigned_d0_Tk0 = Var("tracksSigned_d0.at(0)",float,doc=" signed transverse impact parameter for each track, track index 0"),		
	tracksInnerPosition_X_Tk0 = Var("tracksInnerPosition().at(0).X()",float,doc="the X position of the innermost hit of track index 0"),
	tracksInnerPosition_Y_Tk0 = Var("tracksInnerPosition().at(0).Y()",float,doc="the Y position of the innermost hit of track index 0"),
	tracksInnerPosition_Z_Tk0 = Var("tracksInnerPosition().at(0).Z()",float,doc="the Z position of the innermost hit of track index 0"),
	tracksPout_Px_Tk0 = Var("tracksPout().at(0).X()",float,doc="track Px measured at the outermost hit, track index 0"),
	tracksPout_Py_Tk0 = Var("tracksPout().at(0).Y()",float,doc="track Py measured at the outermost hit, track index 0"),
	tracksPout_Pz_Tk0 = Var("tracksPout().at(0).Z()",float,doc="track Pz measured at the outermost hit, track index 0"),
	tracksPin_Px_Tk0 = Var("tracksPin().at(0).X()",float,doc="track Px measured at the innermost hit, track index 0"),
	tracksPin_Py_Tk0 = Var("tracksPin().at(0).Y()",float,doc="track Py measured at the innermost hit, track index 0"),
	tracksPin_Pz_Tk0 = Var("tracksPin().at(0).Z()",float,doc="track Pz measured at the innermost hit, track index 0"),
	nHitsBeforeVtx_Tk0 = Var("nHitsBeforeVtx().at(0)",int,doc="number of hits before the vertex along each track trajector, track index 0"),
 	dlClosestHitToVtx_Tk0 = Var("dlClosestHitToVtx().at(0).value()",float,doc="signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 0"),	
	dlClosestHitToVtx_err_Tk0 =  Var("dlClosestHitToVtx().at(0).error()",float,doc=" Error of signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 0"),
	dlClosestHitToVtx_sig_Tk0 =  Var("dlClosestHitToVtx().at(0).significance()",float,doc="Significance of signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 0"),
	

	tracksSigned_d0_Tk1 = Var("tracksSigned_d0.at(1)",float,doc=" signed transverse impact parameter for each track, track index 1"),		
	tracksInnerPosition_X_Tk1 = Var("tracksInnerPosition().at(1).X()",float,doc="the X position of the innermost hit of track index 1"),
	tracksInnerPosition_Y_Tk1 = Var("tracksInnerPosition().at(1).Y()",float,doc="the Y position of the innermost hit of track index 1"),
	tracksInnerPosition_Z_Tk1 = Var("tracksInnerPosition().at(1).Z()",float,doc="the Z position of the innermost hit of track index 1"),
	tracksPout_Px_Tk1 = Var("tracksPout().at(1).X()",float,doc="track Px measured at the outermost hit, track index 1"),
	tracksPout_Py_Tk1 = Var("tracksPout().at(1).Y()",float,doc="track Py measured at the outermost hit, track index 1"),
	tracksPout_Pz_Tk1 = Var("tracksPout().at(1).Z()",float,doc="track Pz measured at the outermost hit, track index 1"),
	tracksPin_Px_Tk1 = Var("tracksPin().at(1).X()",float,doc="track Px measured at the innermost hit, track index 1"),
	tracksPin_Py_Tk1 = Var("tracksPin().at(1).Y()",float,doc="track Py measured at the innermost hit, track index 1"),
	tracksPin_Pz_Tk1 = Var("tracksPin().at(1).Z()",float,doc="track Pz measured at the innermost hit, track index 1"),
	nHitsBeforeVtx_Tk1 = Var("nHitsBeforeVtx().at(1)",int,doc="number of hits before the vertex along each track trajector, track index 1"),
 	dlClosestHitToVtx_Tk1 = Var("dlClosestHitToVtx().at(1).value()",float,doc="signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 1"),	
	dlClosestHitToVtx_err_Tk1 =  Var("dlClosestHitToVtx().at(1).error()",float,doc=" Error of signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 1"),
	dlClosestHitToVtx_sig_Tk1 =  Var("dlClosestHitToVtx().at(1).significance()",float,doc="Significance of signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 1"),

	nSharedHits = Var("nSharedHits()",int, doc="number of shared hits between the two tracks"),	
	algo = Var("algo()",int, doc="Conversion Track algorithm/provenance undefined=0, ecalSeeded=1, trackerOnly=2, mixed=3, pflow=4, algoSize=5"),
	#quality = Var("quality()",int, doc="generalTracksOnly=0 arbitratedEcalSeeded=1, arbitratedMerged=2, arbitratedMergedEcalGeneral=3, gsfTracksOpenOnly=4, highPurity=8, highEfficiency=9, ecalMatched1Track=10, ecalMatched2Track=11"), #not worth getting to work

	vtx_X = Var("conversionVertex().position().X()",float,doc="x component of the reco conversion vertex"),
	vtx_Y = Var("conversionVertex().position().Y()",float,doc="y component of the reco conversion vertex"),
	vtx_Z = Var("conversionVertex().position().Z()",float,doc="z component of the reco conversion vertex"),
	vtx_cov_00 = Var("conversionVertex().covariance(0,0)",float,doc="vertex covariance element xx"),
	vtx_cov_01 = Var("conversionVertex().covariance(0,1)",float,doc="vertex covariance element xy"),
	vtx_cov_02 = Var("conversionVertex().covariance(0,2)",float,doc="vertex covariance element xz"),
	vtx_cov_11 = Var("conversionVertex().covariance(1,1)",float,doc="vertex covariance element yy"),
	vtx_cov_12 = Var("conversionVertex().covariance(1,2)",float,doc="vertex covariance element yz"),
	vtx_cov_22 = Var("conversionVertex().covariance(2,2)",float,doc="vertex covariance element zz"),
	vtx_chi2 = Var("conversionVertex().chi2()",float,doc="chi-squared"),
	vtx_normalizedChi2 = Var("conversionVertex().normalizedChi2()",float,doc="chi-squared divided by n.d.o.f"),
	vtx_ndof = Var("conversionVertex().ndof()",float,doc="Number of degrees of freedom, tracks may contribute to the vertex with fractional weights. The ndof is then equal to the sum of the track weights. see e.g. CMS NOTE-2006/032, CMS NOTE-2004/002"),
	
	Tk0_pt = Var("tracks().at(0).pt()", float, doc="transverse momentum of track index 0"),		
	Tk0_eta = Var("tracks().at(0).eta()", float, doc="pseudorapidity of track index 0"),
	Tk0_phi = Var("tracks().at(0).phi()", float, doc="azimuth of track index 0"),
        Tk0_charge = Var("tracks().at(0).charge()", int, doc="charge of track index 0"),
        Tk0_chi2 = Var("tracks().at(0).chi2()", float, doc="chi-squared of the fit, track index 0"),
	Tk0_normalizedChi2 = Var("tracks().at(0).normalizedChi2()", float, doc="chi-squared divided by n.d.o.f. (or chi-squared * 1e6 if n.d.o.f. is zero), track index 0"),
	Tk0_ndof = Var("tracks().at(0).ndof()", float, doc="number of degrees of freedom of the fit, track index 0"),
	Tk0_found = Var("tracks().at(0).found()", int, doc="number of valid hits on track, track index 0"),
	Tk0_lost = Var("tracks().at(0).lost()", int, doc="numberof lost hits (invalid) on track, track index 0"),
        Tk0_quality = Var("tracks().at(0).qualityMask()", int, doc="track index 0 quality: undefQuality = -1, loose = 0, tight = 1, highPurity = 2, confirmed = 3 ( means found by more than one iteration), goodIterative = 4 ( meaningless ), looseSetWithPV = 5, highPuritySetWithPV = 6, discarded = 7 ( because a better track found. kept in the collection for reference....) "),
	Tk0_algo = Var("tracks().at(0).algo()", int, doc="track index 0 algorithm: see https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackReco/interface/TrackBase.h for a list of algorithms"),		


	Tk1_pt = Var("tracks().at(1).pt()", float, doc="transverse momentum of track index 1"),		
	Tk1_eta = Var("tracks().at(1).eta()", float, doc="pseudorapidity of track index 1"),
	Tk1_phi = Var("tracks().at(1).phi()", float, doc="azimuth of track index 1"),
        Tk1_charge = Var("tracks().at(1).charge()", int, doc="charge of track index 1"),
        Tk1_chi2 = Var("tracks().at(1).chi2()", float, doc="chi-squared of the fit, track index 1"),
	Tk1_normalizedChi2 = Var("tracks().at(1).normalizedChi2()", float, doc="chi-squared divided by n.d.o.f. (or chi-squared * 1e6 if n.d.o.f. is zero), track index 1"),
	Tk1_ndof = Var("tracks().at(1).ndof()", float, doc="number of degrees of freedom of the fit, track index 1"),
	Tk1_found = Var("tracks().at(1).found()", int, doc="number of valid hits on track, track index 1"),
	Tk1_lost = Var("tracks().at(1).lost()", int, doc="numberof lost hits (invalid) on track, track index 1"),
        Tk1_quality = Var("tracks().at(1).qualityMask()", int, doc="track index 1 quality: undefQuality = -1, loose = 0, tight = 1, highPurity = 2, confirmed = 3 ( means found by more than one iteration), goodIterative = 4 ( meaningless ), looseSetWithPV = 5, highPuritySetWithPV = 6, discarded = 7 ( because a better track found. kept in the collection for reference....) "),
	Tk1_algo = Var("tracks().at(1).algo()", int, doc="track index 1 algorithm: see https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackReco/interface/TrackBase.h for a list of algorithms"),	
#	convVtxIdx = Var("vtxmatch:convVtxIdx",int, doc="index of nearest SimVertex"),

	),
	externalVariables = cms.PSet(
		convVtxIdx = ExtVar( cms.InputTag("vtxmatch:convVtxIdx"),int,doc="index of nearest sim vertex"),
		vtxdl = ExtVar( cms.InputTag("vtxmatch:vtxdl"),float,doc="distance between reco conversion and nearest sim vertex"),
	),

)

# vtx matching #########################
vtxmatch = cms.EDProducer("convVtxSimVtxMatcher",
	src= cms.InputTag("allConversions"),
	simsrc = cms.InputTag("g4SimHits"),
)
#this table doesnt work
vtxidxtable = cms.EDProducer("SimpleIntTableProducer",
	src= cms.InputTag("vtxmatch:convVtxIdx"),
	cut= cms.string(""),
	name = cms.string("TEST"),
	doc = cms.string("doc"),
	singleton = cms.bool(False),
	extension = cms.bool(False),
	variables = cms.PSet(
		convVtxIdx = Var("convVtxIdx",int,doc="test"),
				
	)
)

############################

c0tks = cms.EDProducer("convTrkProducer",
	src = cms.InputTag("allConversions"),
	trkCandIdx = cms.int32(0),
)
c1tks = cms.EDProducer("convTrkProducer",
	src = cms.InputTag("allConversions"),
	trkCandIdx = cms.int32(1),
)

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *
c0TrackCandidates = cms.EDProducer("ChargedCandidateProducer",
    src = cms.InputTag("c0tks:convTrks"),
    particleType = cms.string('e+')
)

Tk0Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("c0TrackCandidates"),
    #src = cms.InputTag("GenParticles"),
    cut = cms.string(""), #we should not filter after pruning
    name= cms.string("Tk0"),
    doc = cms.string("conversion track idx 0 particles "),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the taus
    variables = cms.PSet(
         pt  = Var("pt",  float, precision=8),
    )
)


#TODO MC MATCHING
convMCMatchForTable = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = cms.InputTag("c0TrackCandidates"),                         # final reco collection
    matched     = cms.InputTag("finalGenParticles"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),               # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.3),              # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),    # False = just match input in order; True = pick lowest deltaR pair first
)

#convMCMatchForTable = cms.EDProducer("MCTruthDeltaRMatcher",
#	src = cms.InputTag("c0tks"),
#	distMin = cms.double(0.15),
#	matchPDGId = cms.vint32(),
#	matched = cms.InputTag("finalGenParticles"),
#)

c0MCTable = cms.EDProducer("CandMCMatchTableProducer",
    src     = cms.InputTag("c0TrackCandidates"),
    mcMap   = cms.InputTag("convMCMatchForTable"),
    objName = cms.string("Tk0"),
    objType = cms.string("Electron"),
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 qparts"),
)



#convTables = cms.Sequence( convs + convTable)
convTables = cms.Sequence( vtxmatch *  convTable)
#convTablesMC = cms.Sequence( convTable + convMCMatchForTable + qMCTable #)
#testseq = cms.Sequence(convTable + c0tks +  c0TrackCandidates+Tk0Table + convMCMatchForTable  +c0MCTable) 
#testseq = cms.Sequence(c0tks + c0TrackCandidates + Tk0Table + convTable + convMCMatchForTable + c0MCTable)
testseq = cms.Sequence(convTable)

