#cmsRun conversions_cfg.py print inputFiles_load=./dasquery/bottomup/SingleMuon2017/testaodlist.list outputFile=./OutputFiles/SingleMuon2017.root maxEvents=20000



#cmsRun conversionsMC_cfg.py print inputFiles_load=./dasquery/bottomup/SingleMuon2017/testaodmclist.list outputFile=./OutputFiles/Jpsimc2017.root maxEvents=5000


#cmsRun conversionsMC_cfg.py print inputFiles_load=./dasquery/bottomup/SingleMuon2017/testaodmclist.list outputFile=./OutputFiles/Jpsimc2017.root maxEvents=100


cmsRun conversionsMC_RECOSIM_cfg.py print inputFiles_load=./InputFiles/recosimlist.list outputFile=./OutputFiles/recotest.root maxEvents=100
