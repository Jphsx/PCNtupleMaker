import FWCore.ParameterSet.Config as cms

process = cms.Process("miniflatntuple")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)
process.CandVars = cms.PSet(
    charge = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('electric charge'),
        expr = cms.string('charge'),
        mcOnly = cms.bool(False),
        precision = cms.int32(-1),
        type = cms.string('int')
    ),
    eta = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('eta'),
        expr = cms.string('eta'),
        mcOnly = cms.bool(False),
        precision = cms.int32(12),
        type = cms.string('float')
    ),
    mass = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('mass'),
        expr = cms.string('mass'),
        mcOnly = cms.bool(False),
        precision = cms.int32(10),
        type = cms.string('float')
    ),
    pdgId = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('PDG code assigned by the event reconstruction (not by MC truth)'),
        expr = cms.string('pdgId'),
        mcOnly = cms.bool(False),
        precision = cms.int32(-1),
        type = cms.string('int')
    ),
    phi = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('phi'),
        expr = cms.string('phi'),
        mcOnly = cms.bool(False),
        precision = cms.int32(12),
        type = cms.string('float')
    ),
    pt = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('pt'),
        expr = cms.string('pt'),
        mcOnly = cms.bool(False),
        precision = cms.int32(-1),
        type = cms.string('float')
    )
)

process.HFRecalParameterBlock = cms.PSet(
    HFdepthOneParameterA = cms.vdouble(
        0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
        0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
        0.058939, 0.125497
    ),
    HFdepthOneParameterB = cms.vdouble(
        -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
        2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
        0.000425, 0.000209
    ),
    HFdepthTwoParameterA = cms.vdouble(
        0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
        0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
        0.051579, 0.086593
    ),
    HFdepthTwoParameterB = cms.vdouble(
        -2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
        1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
        0.000157, -3e-06
    )
)

process.P3Vars = cms.PSet(
    eta = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('eta'),
        expr = cms.string('eta'),
        mcOnly = cms.bool(False),
        precision = cms.int32(12),
        type = cms.string('float')
    ),
    phi = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('phi'),
        expr = cms.string('phi'),
        mcOnly = cms.bool(False),
        precision = cms.int32(12),
        type = cms.string('float')
    ),
    pt = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('pt'),
        expr = cms.string('pt'),
        mcOnly = cms.bool(False),
        precision = cms.int32(-1),
        type = cms.string('float')
    )
)

process.P4Vars = cms.PSet(
    eta = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('eta'),
        expr = cms.string('eta'),
        mcOnly = cms.bool(False),
        precision = cms.int32(12),
        type = cms.string('float')
    ),
    mass = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('mass'),
        expr = cms.string('mass'),
        mcOnly = cms.bool(False),
        precision = cms.int32(10),
        type = cms.string('float')
    ),
    phi = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('phi'),
        expr = cms.string('phi'),
        mcOnly = cms.bool(False),
        precision = cms.int32(12),
        type = cms.string('float')
    ),
    pt = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('pt'),
        expr = cms.string('pt'),
        mcOnly = cms.bool(False),
        precision = cms.int32(-1),
        type = cms.string('float')
    )
)

process.PTVars = cms.PSet(
    phi = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('phi'),
        expr = cms.string('phi'),
        mcOnly = cms.bool(False),
        precision = cms.int32(12),
        type = cms.string('float')
    ),
    pt = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('pt'),
        expr = cms.string('pt'),
        mcOnly = cms.bool(False),
        precision = cms.int32(-1),
        type = cms.string('float')
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.convTable = cms.EDProducer("SimpleConversionTableProducer",
    cut = cms.string(''),
    doc = cms.string('all conversions after basic selection'),
    extension = cms.bool(False),
    name = cms.string('Conv'),
    singleton = cms.bool(False),
    src = cms.InputTag("allConversions"),
    variables = cms.PSet(
        EoverP = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Super Cluster energy divided by track pair momentum if Standard seeding method. If a pointer to two (or more clusters) is stored in the conversion, this method returns the energy sum of clusters divided by the track pair momentum. Track innermost momentum is used here'),
            expr = cms.string('EoverP()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        EoverPrefittedTracks = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Track momentum refitted with vertex constraint is used'),
            expr = cms.string('EoverPrefittedTracks'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk0_algo = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track index 0 algorithm: see https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackReco/interface/TrackBase.h for a list of algorithms'),
            expr = cms.string('tracks().at(0).algo()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        Tk0_charge = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('charge of track index 0'),
            expr = cms.string('tracks().at(0).charge()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        Tk0_chi2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('chi-squared of the fit, track index 0'),
            expr = cms.string('tracks().at(0).chi2()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk0_eta = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('pseudorapidity of track index 0'),
            expr = cms.string('tracks().at(0).eta()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk0_found = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('number of valid hits on track, track index 0'),
            expr = cms.string('tracks().at(0).found()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        Tk0_lost = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('numberof lost hits (invalid) on track, track index 0'),
            expr = cms.string('tracks().at(0).lost()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        Tk0_ndof = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('number of degrees of freedom of the fit, track index 0'),
            expr = cms.string('tracks().at(0).ndof()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk0_normalizedChi2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('chi-squared divided by n.d.o.f. (or chi-squared * 1e6 if n.d.o.f. is zero), track index 0'),
            expr = cms.string('tracks().at(0).normalizedChi2()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk0_phi = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('azimuth of track index 0'),
            expr = cms.string('tracks().at(0).phi()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk0_pt = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('transverse momentum of track index 0'),
            expr = cms.string('tracks().at(0).pt()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk0_quality = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track index 0 quality: undefQuality = -1, loose = 0, tight = 1, highPurity = 2, confirmed = 3 ( means found by more than one iteration), goodIterative = 4 ( meaningless ), looseSetWithPV = 5, highPuritySetWithPV = 6, discarded = 7 ( because a better track found. kept in the collection for reference....) '),
            expr = cms.string('tracks().at(0).qualityMask()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        Tk1_algo = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track index 1 algorithm: see https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackReco/interface/TrackBase.h for a list of algorithms'),
            expr = cms.string('tracks().at(1).algo()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        Tk1_charge = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('charge of track index 1'),
            expr = cms.string('tracks().at(1).charge()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        Tk1_chi2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('chi-squared of the fit, track index 1'),
            expr = cms.string('tracks().at(1).chi2()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk1_eta = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('pseudorapidity of track index 1'),
            expr = cms.string('tracks().at(1).eta()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk1_found = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('number of valid hits on track, track index 1'),
            expr = cms.string('tracks().at(1).found()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        Tk1_lost = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('numberof lost hits (invalid) on track, track index 1'),
            expr = cms.string('tracks().at(1).lost()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        Tk1_ndof = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('number of degrees of freedom of the fit, track index 1'),
            expr = cms.string('tracks().at(1).ndof()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk1_normalizedChi2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('chi-squared divided by n.d.o.f. (or chi-squared * 1e6 if n.d.o.f. is zero), track index 1'),
            expr = cms.string('tracks().at(1).normalizedChi2()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk1_phi = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('azimuth of track index 1'),
            expr = cms.string('tracks().at(1).phi()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk1_pt = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('transverse momentum of track index 1'),
            expr = cms.string('tracks().at(1).pt()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        Tk1_quality = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track index 1 quality: undefQuality = -1, loose = 0, tight = 1, highPurity = 2, confirmed = 3 ( means found by more than one iteration), goodIterative = 4 ( meaningless ), looseSetWithPV = 5, highPuritySetWithPV = 6, discarded = 7 ( because a better track found. kept in the collection for reference....) '),
            expr = cms.string('tracks().at(1).qualityMask()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        algo = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion Track algorithm/provenance undefined=0, ecalSeeded=1, trackerOnly=2, mixed=3, pflow=4, algoSize=5'),
            expr = cms.string('algo()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        dPhiTracksAtVtx = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('deltaPhi tracks at innermost point'),
            expr = cms.string('dPhiTracksAtVtx()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        distOfMinimumApproach = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Dist of minimum approach between tracks'),
            expr = cms.string('distOfMinimumApproach()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        dlClosestHitToVtx_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 0'),
            expr = cms.string('dlClosestHitToVtx().at(0).value()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        dlClosestHitToVtx_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 1'),
            expr = cms.string('dlClosestHitToVtx().at(1).value()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        dlClosestHitToVtx_err_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string(' Error of signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 0'),
            expr = cms.string('dlClosestHitToVtx().at(0).error()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        dlClosestHitToVtx_err_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string(' Error of signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 1'),
            expr = cms.string('dlClosestHitToVtx().at(1).error()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        dlClosestHitToVtx_sig_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Significance of signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 0'),
            expr = cms.string('dlClosestHitToVtx().at(0).significance()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        dlClosestHitToVtx_sig_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Significance of signed decay length with uncertainty from nearest hit on track to the conversion vtx positions, track index 1'),
            expr = cms.string('dlClosestHitToVtx().at(1).significance()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        dxy = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('transverse impact parameter, computed with respect to given beamspot (0,0,0) and refitted pair momentum'),
            expr = cms.string('dxy()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        dz = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('longitudinal impact parameter, computed with respect to given beamspot (0,0,0) and refitted pair momentum'),
            expr = cms.string('dz()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        isConverted = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Bool flagging objects having track size >0'),
            expr = cms.string('isConverted()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('bool')
        ),
        lxy = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('transverse decay length, computed with respect to given beamspot (0,0,0) and refitted pair momentum'),
            expr = cms.string('lxy()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        lz = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('longitudinal decay length, computed with respect to given beamspot (0,0,0) and refitted pair momentum'),
            expr = cms.string('lz()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        nHitsBeforeVtx_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('number of hits before the vertex along each track trajector, track index 0'),
            expr = cms.string('nHitsBeforeVtx().at(0)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nHitsBeforeVtx_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('number of hits before the vertex along each track trajector, track index 1'),
            expr = cms.string('nHitsBeforeVtx().at(1)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nSharedHits = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('number of shared hits between the two tracks'),
            expr = cms.string('nSharedHits()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        nTracks = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Number of tracks= 0,1,2'),
            expr = cms.string('nTracks()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('int')
        ),
        pairCotThetaSeparation = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Delta cot(Theta) where Theta is the angle in the (y,z) plane between the two tracks. Original tracks are used'),
            expr = cms.string('pairCotThetaSeparation()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        pairInvariantMass = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('if nTracks=2 returns the pair invariant mass. Original tracks are used here'),
            expr = cms.string('pairInvariantMass()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        pairMomentum_Px = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion tracks momentum from the tracks inner momentum, x-component'),
            expr = cms.string('pairMomentum().X()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        pairMomentum_Py = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion tracks momentum from the tracks inner momentum, y-component'),
            expr = cms.string('pairMomentum().Y()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        pairMomentum_Pz = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion tracks momentum from the tracks inner momentum, z-component'),
            expr = cms.string('pairMomentum().Z()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        refittedPair4Momentum_E = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion track pair 4-momentum from the tracks refitted with vertex constraint, E'),
            expr = cms.string('refittedPair4Momentum().E()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        refittedPair4Momentum_M = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion track pair 4-momentum from the tracks refitted with vertex constraint, M'),
            expr = cms.string('refittedPair4Momentum().M()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        refittedPair4Momentum_Px = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion track pair 4-momentum from the tracks refitted with vertex constraint, Px'),
            expr = cms.string('refittedPair4Momentum().Px()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        refittedPair4Momentum_Py = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion track pair 4-momentum from the tracks refitted with vertex constraint, Py'),
            expr = cms.string('refittedPair4Momentum().Py()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        refittedPair4Momentum_Pz = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion track pair 4-momentum from the tracks refitted with vertex constraint, Pz'),
            expr = cms.string('refittedPair4Momentum().Pz()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        refittedPairMomentum_Px = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion tracks momentum from the tracks refitted with vertex constraint, Px'),
            expr = cms.string('refittedPairMomentum().X()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        refittedPairMomentum_Py = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion tracks momentum from the tracks refitted with vertex constraint, Py'),
            expr = cms.string('refittedPairMomentum().Y()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        refittedPairMomentum_Pz = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Conversion tracks momentum from the tracks refitted with vertex constraint, Pz'),
            expr = cms.string('refittedPairMomentum().Z()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksInnerPosition_X_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('the X position of the innermost hit of track index 0'),
            expr = cms.string('tracksInnerPosition().at(0).X()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksInnerPosition_X_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('the X position of the innermost hit of track index 1'),
            expr = cms.string('tracksInnerPosition().at(1).X()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksInnerPosition_Y_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('the Y position of the innermost hit of track index 0'),
            expr = cms.string('tracksInnerPosition().at(0).Y()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksInnerPosition_Y_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('the Y position of the innermost hit of track index 1'),
            expr = cms.string('tracksInnerPosition().at(1).Y()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksInnerPosition_Z_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('the Z position of the innermost hit of track index 0'),
            expr = cms.string('tracksInnerPosition().at(0).Z()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksInnerPosition_Z_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('the Z position of the innermost hit of track index 1'),
            expr = cms.string('tracksInnerPosition().at(1).Z()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPin_Px_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Px measured at the innermost hit, track index 0'),
            expr = cms.string('tracksPin().at(0).X()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPin_Px_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Px measured at the innermost hit, track index 1'),
            expr = cms.string('tracksPin().at(1).X()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPin_Py_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Py measured at the innermost hit, track index 0'),
            expr = cms.string('tracksPin().at(0).Y()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPin_Py_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Py measured at the innermost hit, track index 1'),
            expr = cms.string('tracksPin().at(1).Y()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPin_Pz_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Pz measured at the innermost hit, track index 0'),
            expr = cms.string('tracksPin().at(0).Z()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPin_Pz_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Pz measured at the innermost hit, track index 1'),
            expr = cms.string('tracksPin().at(1).Z()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPout_Px_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Px measured at the outermost hit, track index 0'),
            expr = cms.string('tracksPout().at(0).X()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPout_Px_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Px measured at the outermost hit, track index 1'),
            expr = cms.string('tracksPout().at(1).X()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPout_Py_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Py measured at the outermost hit, track index 0'),
            expr = cms.string('tracksPout().at(0).Y()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPout_Py_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Py measured at the outermost hit, track index 1'),
            expr = cms.string('tracksPout().at(1).Y()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPout_Pz_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Pz measured at the outermost hit, track index 0'),
            expr = cms.string('tracksPout().at(0).Z()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksPout_Pz_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('track Pz measured at the outermost hit, track index 1'),
            expr = cms.string('tracksPout().at(1).Z()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksSigned_d0_Tk0 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string(' signed transverse impact parameter for each track, track index 0'),
            expr = cms.string('tracksSigned_d0.at(0)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        tracksSigned_d0_Tk1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string(' signed transverse impact parameter for each track, track index 1'),
            expr = cms.string('tracksSigned_d0.at(1)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_X = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('x component of the reco conversion vertex'),
            expr = cms.string('conversionVertex().position().X()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_Y = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('y component of the reco conversion vertex'),
            expr = cms.string('conversionVertex().position().Y()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_Z = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('z component of the reco conversion vertex'),
            expr = cms.string('conversionVertex().position().Z()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_chi2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('chi-squared'),
            expr = cms.string('conversionVertex().chi2()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_cov_00 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('vertex covariance element xx'),
            expr = cms.string('conversionVertex().covariance(0,0)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_cov_01 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('vertex covariance element xy'),
            expr = cms.string('conversionVertex().covariance(0,1)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_cov_02 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('vertex covariance element xz'),
            expr = cms.string('conversionVertex().covariance(0,2)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_cov_11 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('vertex covariance element yy'),
            expr = cms.string('conversionVertex().covariance(1,1)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_cov_12 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('vertex covariance element yz'),
            expr = cms.string('conversionVertex().covariance(1,2)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_cov_22 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('vertex covariance element zz'),
            expr = cms.string('conversionVertex().covariance(2,2)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_ndof = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Number of degrees of freedom, tracks may contribute to the vertex with fractional weights. The ndof is then equal to the sum of the track weights. see e.g. CMS NOTE-2006/032, CMS NOTE-2004/002'),
            expr = cms.string('conversionVertex().ndof()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        vtx_normalizedChi2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('chi-squared divided by n.d.o.f'),
            expr = cms.string('conversionVertex().normalizedChi2()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        zOfPrimaryVertexFromTracks = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('z position of intersection with beamspot in rz plane (possible tilt of beamspot is neglected)'),
            expr = cms.string('zOfPrimaryVertexFromTracks()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        )
    )
)


process.genTable = cms.EDProducer("SimpleGenEventFlatTableProducer",
    cut = cms.string(''),
    doc = cms.string('Generator information'),
    extension = cms.bool(False),
    name = cms.string('Generator'),
    singleton = cms.bool(True),
    src = cms.InputTag("generator"),
    variables = cms.PSet(
        binvar = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('MC generation binning value'),
            expr = cms.string('?hasBinningValues()?binningValues()[0]:-1'),
            mcOnly = cms.bool(False),
            precision = cms.int32(14),
            type = cms.string('float')
        ),
        id1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('id of first parton'),
            expr = cms.string('?hasPDF?pdf().id.first:-1'),
            mcOnly = cms.bool(False),
            precision = cms.int32(6),
            type = cms.string('int')
        ),
        id2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('id of second parton'),
            expr = cms.string('?hasPDF?pdf().id.second:-1'),
            mcOnly = cms.bool(False),
            precision = cms.int32(6),
            type = cms.string('int')
        ),
        scalePDF = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('Q2 scale for PDF'),
            expr = cms.string('?hasPDF?pdf().scalePDF:-1'),
            mcOnly = cms.bool(False),
            precision = cms.int32(14),
            type = cms.string('float')
        ),
        weight = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('MC generator weight'),
            expr = cms.string('weight()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(14),
            type = cms.string('float')
        ),
        x1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('x1 fraction of proton momentum carried by the first parton'),
            expr = cms.string('?hasPDF?pdf().x.first:-1'),
            mcOnly = cms.bool(False),
            precision = cms.int32(14),
            type = cms.string('float')
        ),
        x2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('x2 fraction of proton momentum carried by the second parton'),
            expr = cms.string('?hasPDF?pdf().x.second:-1'),
            mcOnly = cms.bool(False),
            precision = cms.int32(14),
            type = cms.string('float')
        ),
        xpdf1 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('x*pdf(x) for the first parton'),
            expr = cms.string('?hasPDF?pdf().xPDF.first:-1'),
            mcOnly = cms.bool(False),
            precision = cms.int32(14),
            type = cms.string('float')
        ),
        xpdf2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('x*pdf(x) for the second parton'),
            expr = cms.string('?hasPDF?pdf().xPDF.second:-1'),
            mcOnly = cms.bool(False),
            precision = cms.int32(14),
            type = cms.string('float')
        )
    )
)


process.puTable = cms.EDProducer("NPUTablesProducer",
    pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    src = cms.InputTag("slimmedAddPileupInfo"),
    zbins = cms.vdouble(
        0.0, 1.7, 2.6, 3.0, 3.5, 
        4.2, 5.2, 6.0, 7.5, 9.0, 
        12.0
    )
)


process.pvCandidateTable = cms.EDProducer("SimpleVertexTableProducer",
    cut = cms.string('!isFake && ndof > 4 && abs(position().Z() ) <= 24 && position().Rho() <= 2'),
    doc = cms.string('all Primary Vertices after basic selection'),
    extension = cms.bool(False),
    name = cms.string('PV'),
    singleton = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices"),
    variables = cms.PSet(
        X = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('secondary vertex X position, in cm'),
            expr = cms.string('position().X()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        Y = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('secondary vertex Y position, in cm'),
            expr = cms.string('position().Y()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(10),
            type = cms.string('float')
        ),
        Z = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('secondary vertex Z position, in cm'),
            expr = cms.string('position().Z()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(14),
            type = cms.string('float')
        ),
        chi2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('chi squared'),
            expr = cms.string('chi2()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        cov_00 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string(' pv covrariance element xx '),
            expr = cms.string('covariance(0,0)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        cov_01 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string(' pv covrariance element xy '),
            expr = cms.string('covariance(0,1)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        cov_02 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string(' pv covrariance element xz '),
            expr = cms.string('covariance(0,2)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        cov_11 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string(' pv covrariance element yy '),
            expr = cms.string('covariance(1,1)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        cov_12 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string(' pv covrariance element yz '),
            expr = cms.string('covariance(1,2)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        cov_22 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string(' pv covrariance element zz '),
            expr = cms.string('covariance(2,2)'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            type = cms.string('float')
        ),
        ndof = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('number of degrees of freedom'),
            expr = cms.string('ndof()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(8),
            type = cms.string('float')
        ),
        normalizedChi2 = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('reduced chi2, i.e. chi/ndof'),
            expr = cms.string('normalizedChi2()'),
            mcOnly = cms.bool(False),
            precision = cms.int32(8),
            type = cms.string('float')
        )
    )
)


process.randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")


process.rhoTable = cms.EDProducer("GlobalVariablesTableProducer",
    variables = cms.PSet(
        fixedGridRhoFastjetAll = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('rho from all PF Candidates, used e.g. for JECs'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            src = cms.InputTag("fixedGridRhoFastjetAll"),
            type = cms.string('double')
        ),
        fixedGridRhoFastjetCentral = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('rho from all PF Candidates for central region, used e.g. for JECs'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            src = cms.InputTag("fixedGridRhoFastjetCentral"),
            type = cms.string('double')
        ),
        fixedGridRhoFastjetCentralCalo = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('rho from calo towers with |eta| < 2.5, used e.g. egamma PFCluster isolation'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            src = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
            type = cms.string('double')
        ),
        fixedGridRhoFastjetCentralChargedPileUp = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('rho from charged PF Candidates for central region, used e.g. for JECs'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            src = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
            type = cms.string('double')
        ),
        fixedGridRhoFastjetCentralNeutral = cms.PSet(
            compression = cms.string('none'),
            doc = cms.string('rho from neutral PF Candidates with |eta| < 2.5, used e.g. for rho corrections of some lepton isolations'),
            mcOnly = cms.bool(False),
            precision = cms.int32(-1),
            src = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
            type = cms.string('double')
        )
    )
)


process.vertexTable = cms.EDProducer("VertexTableProducer",
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(3),
    goodPvCut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
    pvName = cms.string('PV'),
    pvSrc = cms.InputTag("offlinePrimaryVertices")
)


process.out = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string('defaultout_numEvent100.root'),
    outputCommands = cms.untracked.vstring(
        'drop *', 
        'keep nanoaodFlatTable_*Table_*_*'
    )
)


process.DQMStore = cms.Service("DQMStore",
    LSbasedMode = cms.untracked.bool(False),
    collateHistograms = cms.untracked.bool(False),
    enableMultiThread = cms.untracked.bool(False),
    forceResetOnBeginLumi = cms.untracked.bool(False),
    referenceFileName = cms.untracked.string(''),
    saveByLumi = cms.untracked.bool(False),
    verbose = cms.untracked.int32(0),
    verboseQT = cms.untracked.int32(0)
)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring(
        'FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'
    ),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(10000)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring(
        'warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'
    ),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    CTPPSFastRecHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(1357987)
    ),
    LHCTransport = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(87654321)
    ),
    MuonSimHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(987346)
    ),
    VtxSmeared = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(98765432)
    ),
    ecalPreshowerRecHit = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(6541321)
    ),
    ecalRecHit = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(654321)
    ),
    externalLHEProducer = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(234567)
    ),
    famosPileUp = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    fastSimProducer = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(13579)
    ),
    fastTrackerRecHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(24680)
    ),
    g4SimHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(11)
    ),
    generator = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(123456789)
    ),
    hbhereco = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(541321)
    ),
    hfreco = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(541321)
    ),
    hiSignal = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(123456789)
    ),
    hiSignalG4SimHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(11)
    ),
    hiSignalLHCTransport = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(88776655)
    ),
    horeco = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(541321)
    ),
    l1ParamMuons = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(6453209)
    ),
    mix = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(12345)
    ),
    mixData = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(12345)
    ),
    mixGenPU = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    mixRecoTracks = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    mixSimCaloHits = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    paramMuons = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(54525)
    ),
    saveFileName = cms.untracked.string(''),
    simBeamSpotFilter = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(87654321)
    ),
    simMuonCSCDigis = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(11223344)
    ),
    simMuonDTDigis = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(1234567)
    ),
    simMuonRPCDigis = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(1234567)
    ),
    simSiStripDigiSimLink = cms.PSet(
        engineName = cms.untracked.string('MixMaxRng'),
        initialSeed = cms.untracked.uint32(1234567)
    )
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring(
        'HCAL', 
        'ZDC', 
        'CASTOR', 
        'EcalBarrel', 
        'EcalEndcap', 
        'EcalPreshower', 
        'TOWER'
    )
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer",
    appendToDataLabel = cms.string('')
)


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.MuonNumberingInitialization = cms.ESProducer("MuonNumberingInitialization")


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    compatibiltyWith11 = cms.untracked.bool(True),
    useDDD = cms.untracked.bool(False)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer")


process.XMLFromDBSource = cms.ESProducer("XMLIdealGeometryESProducer",
    label = cms.string('Extended'),
    rootDDName = cms.string('cms:OCMS')
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    MergePosition = cms.untracked.bool(False),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiPixelQualityFromDbRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        )
    ),
    siPixelQualityLabel = cms.string('')
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGainRcd')
        ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )
    ),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiStripDetVOffRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )
    ),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('106X_mc2017_realistic_v7'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.HcalTimeSlewEP = cms.ESSource("HcalTimeSlewEP",
    appendToDataLabel = cms.string('HBHE'),
    timeSlewParametersM2 = cms.VPSet(
        cms.PSet(
            slope = cms.double(-3.178648),
            tmax = cms.double(16.0),
            tzero = cms.double(23.960177)
        ), 
        cms.PSet(
            slope = cms.double(-1.5610227),
            tmax = cms.double(10.0),
            tzero = cms.double(11.977461)
        ), 
        cms.PSet(
            slope = cms.double(-1.075824),
            tmax = cms.double(6.25),
            tzero = cms.double(9.109694)
        )
    ),
    timeSlewParametersM3 = cms.VPSet(
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(15.5),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-3.2),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(32.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        )
    )
)


process.HepPDTESSource = cms.ESSource("HepPDTESSource",
    pdtFileName = cms.FileInPath('SimGeneral/HepPDTESSource/data/pythiaparticle.tbl')
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HBRecalibration = cms.bool(False),
    HBmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHB.txt'),
    HBreCalibCutoff = cms.double(20.0),
    HERecalibration = cms.bool(False),
    HEmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHE.txt'),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalParameterBlock = cms.PSet(
        HFdepthOneParameterA = cms.vdouble(
            0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
            0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
            0.058939, 0.125497
        ),
        HFdepthOneParameterB = cms.vdouble(
            -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
            2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
            0.000425, 0.000209
        ),
        HFdepthTwoParameterA = cms.vdouble(
            0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
            0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
            0.051579, 0.086593
        ),
        HFdepthTwoParameterB = cms.vdouble(
            -2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
            1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
            0.000157, -3e-06
        )
    ),
    HFRecalibration = cms.bool(False),
    SiPMCharacteristics = cms.VPSet(
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(36000)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(2500)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(0)
        )
    ),
    hb = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.19),
        gainWidth = cms.vdouble(0.0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.285),
        pedestalWidth = cms.double(0.809),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.49, 1.8, 7.2, 37.9),
        qieSlope = cms.vdouble(0.912, 0.917, 0.922, 0.923),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(8)
    ),
    hbUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(150),
            intlumiToNeutrons = cms.double(367000000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(-5)
        ),
        recoShape = cms.int32(206),
        zsThreshold = cms.int32(16)
    ),
    he = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.23),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.163),
        pedestalWidth = cms.double(0.9698),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.38, 2.0, 7.6, 39.6),
        qieSlope = cms.vdouble(0.912, 0.916, 0.92, 0.922),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(9)
    ),
    heUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(75),
            intlumiToNeutrons = cms.double(29200000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(5)
        ),
        recoShape = cms.int32(206),
        zsThreshold = cms.int32(16)
    ),
    hf = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(9.354),
        pedestalWidth = cms.double(2.516),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(-0.87, 1.4, 7.8, -29.6),
        qieSlope = cms.vdouble(0.359, 0.358, 0.36, 0.367),
        qieType = cms.int32(0),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    hfUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(13.33),
        pedestalWidth = cms.double(3.33),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(0.0697, -0.7405, 12.38, -671.9),
        qieSlope = cms.vdouble(0.297, 0.298, 0.298, 0.313),
        qieType = cms.int32(1),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    ho = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.006, 0.0087),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(201),
        pedestal = cms.double(12.06),
        pedestalWidth = cms.double(0.6285),
        photoelectronsToAnalog = cms.double(4.0),
        qieOffset = cms.vdouble(-0.44, 1.4, 7.1, 38.5),
        qieSlope = cms.vdouble(0.907, 0.915, 0.92, 0.921),
        qieType = cms.int32(0),
        recoShape = cms.int32(201),
        zsThreshold = cms.int32(24)
    ),
    iLumi = cms.double(-1.0),
    killHE = cms.bool(False),
    testHEPlan1 = cms.bool(False),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths'),
    useHBUpgrade = cms.bool(False),
    useHEUpgrade = cms.bool(False),
    useHFUpgrade = cms.bool(False),
    useHOUpgrade = cms.bool(True),
    useIeta18depth1 = cms.bool(True),
    useLayer0Weight = cms.bool(False)
)


process.prefer("es_hardcode")

process.convTables = cms.Sequence(process.convTable)


process.globalTablesMC = cms.Sequence(process.puTable+process.genTable)


process.vertexTables = cms.Sequence(process.pvCandidateTable)


process.globalTables = cms.Sequence(process.rhoTable)


process.Path = cms.Path(process.vertexTables+process.convTables)


process.end = cms.EndPath(process.out)


