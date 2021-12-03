


void runproof(){

//defaultout_numEvent100_64.root

TProof::Open("workers=16");
//TProof::Open("workers=1");
TChain *chain = new TChain("Events", " ");

//
//chain->Add("./defaultout_numEvent100_64.root");
//chain->Add("/home/t3-ku/janguian/jsingera/DPG/PC/MinBias2018/minBias2018/DPG_MinBias2018/200625_195959/0000/defaultout_numEvent100_1*.root");
//chain->Add("/home/t3-ku/janguian/jsingera/DPG/PC/MinBias2018/minBias2018/DPG_MinBias2018/200625_195959/0000/defaultout_numEvent100_2*.root");
//chain->Add("/home/t3-ku/janguian/jsingera/DPG/PC/MinBias2018/minBias2018/DPG_MinBias2018/200625_195959/0000/defaultout_numEvent100_3*.root");
//chain->Add("/home/t3-ku/janguian/jsingera/DPG/PC/MinBias2018/minBias2018/DPG_MinBias2018/200625_195959/0000/defaultout_numEvent100_4*.root");
//chain->Add("/home/t3-ku/janguian/jsingera/DPG/PC/MinBias2018/minBias2018/DPG_MinBias2018/200625_195959/0000/defaultout_numEvent100_5*.root");
//chain->Add("/home/t3-ku/janguian/jsingera/DPG/PC/MinBias2018/minBias2018/DPG_MinBias2018/200625_195959/0000/defaultout_numEvent100_6*.root");
//chain->Add("/home/t3-ku/janguian/jsingera/DPG/PC/MinBias2018/minBias2018/DPG_MinBias2018/200625_195959/0000/defaultout_numEvent100_7*.root");
//chain->Add("/home/t3-ku/janguian/jsingera/DPG/PC/MinBias2018/minBias2018/DPG_MinBias2018/200625_195959/0000/defaultout_numEvent100_8*.root");
chain->Add("/home/t3-ku/janguian/jsingera/DPG/PC/MinBias2018/minBias2018/DPG_MinBias2018/200625_195959/0000/defaultout_numEvent100_9*.root");

chain->SetProof();

chain->Process("PCNtuple.C+");

//TProofLog *pl = TProof::Mgr("")->GetSessionLogs();
//pl->Save("*", "log.txt");
//TProofLog *pl = TProofMgrLite::GetSessionLogs();
//pl->Save("*","log.txt");

}
