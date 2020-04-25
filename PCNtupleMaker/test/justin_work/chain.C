


void chain(){
std::string PATH = "/home/t3-ku/janguian/jsingera/DPG/PC/Run2018C/SingleMuon/DPG_SingleMuon2018C/200425_171633/0000/"; 
std::vector<std::string> names{
"defaultout_numEvent100_1.root" }; /*,    "defaultout_numEvent100_28.root ", "defaultout_numEvent100_46.root",   "defaultout_numEvent100_64.root",   "defaultout_numEvent100_82.root", 
"defaultout_numEvent100_10.root",   "defaultout_numEvent100_280.root", "defaultout_numEvent100_460.root",  "defaultout_numEvent100_640.root",  "defaultout_numEvent100_820.root",
"defaultout_numEvent100_100.root",  "defaultout_numEvent100_281.root", "defaultout_numEvent100_461.root" , "defaultout_numEvent100_641.root",  "defaultout_numEvent100_821.root",
"defaultout_numEvent100_101.root",  "defaultout_numEvent100_282.root",  "defaultout_numEvent100_462.root",  "defaultout_numEvent100_642.root",  "defaultout_numEvent100_822.root",
"defaultout_numEvent100_102.root",  "defaultout_numEvent100_283.root",  "defaultout_numEvent100_463.root",  "defaultout_numEvent100_643.root",  "defaultout_numEvent100_823.root",
"defaultout_numEvent100_103.root",  "defaultout_numEvent100_284.root",  "defaultout_numEvent100_464.root",  "defaultout_numEvent100_644.root",  "defaultout_numEvent100_824.root",
"defaultout_numEvent100_104.root",  "defaultout_numEvent100_285.root",  "defaultout_numEvent100_465.root",  "defaultout_numEvent100_645.root",  "defaultout_numEvent100_825.root",
"defaultout_numEvent100_105.root",  "defaultout_numEvent100_286.root",  "defaultout_numEvent100_466.root",  "defaultout_numEvent100_646.root",  "defaultout_numEvent100_826.root"
};*/

TChain* c = new TChain("Events");
for(int i=0; i<names.size(); i++){
	c->AddFile( (PATH+names.at(i)).c_str() );

}

c->Process("convsel.C");

}
