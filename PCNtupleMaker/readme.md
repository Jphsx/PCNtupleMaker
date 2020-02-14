




To setup: uses CMSSW 10_6_8

grab the release
cmsrel 10_6_8

clone the repository into CMSSW_10_6_8/src/

setup environment
cmsenv

build the Ntuplizer
scram b

if the files to be run over are not local, set up grid
voms-proxy-init --voms cms

authenticate grid

To run interactively use configuration scripts in the test directory
cmsRun DYMC2018_cfg.py 


To run over different files, modify or create a new cfg script and change the input source
	fileNames = cms.untracked.vstring(
		'file:/pathtolocalfile/file1.root',
		'file:/pathtolocalfile/file2.root',
	),

If the files are not local, here is syntax for 2017C single muon data file
	fileNames = cms.untracked.vstring(
		'/store/data/Run2017C/SingleMuon/MINIAOD/17Nov2017-v1/70000/FE39A6C4-7FDA-E711-8C6B-02163E014496.root',
	        '/store/data/Run2017C/SingleMuon/MINIAOD/17Nov2017-v1/70000/FCDAA686-73DA-E711-BBE6-02163E013747.root',
	),		



To add things to the tree 

	need to modify NtupleMakerPhotonConversions.h and .cc files

	to add a new collection e.g. ecal photons 
		find the proper collection name in the datatype being used -- edmDumpEventContent someaodfile.root
		add in the new collection via "tokens" which "consume" some edm input tag
		the tokens should be defined  as  class members in the header

		to access the collection objects in the token use a edm::Handle this is a wrapper for the collection object, access the object datamembers via ->
		place the token collections in the handle via iEvent.getByToken(thetoken, thehandle)


	to add variables to the output tree (from some collection )
		declare the variable in the header e.g. std::vector<double> photonpt;

		create the associated branch in the .cc file  e.g. outputTree->Branch("photonpt", "std::vector<double>", &photonpt);
	
		fill the variable somewhere in the analyze function

	after modifying C based files remember to scram b in some upper-directory to build or scram b clean to clean up if necessary

		
