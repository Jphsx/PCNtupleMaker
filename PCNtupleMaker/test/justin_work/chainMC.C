
#include "mc1.C"
#include "mc2.C"
#include "mc3.C"

void chainMC(){
mc1* d1 = new mc1();
std::vector<std::string> names = d1->files;
std::string path = d1->PATH;


mc2* d2 = new mc2();
std::vector<std::string> names2 = d2->files;
std::string path2 = d2->PATH;

mc3* d3 = new mc3();
std::vector<std::string> names3 = d3->files;
std::string path3 = d3->PATH;
//std::cout<<d1->PATH<<std::endl;
//for(int i=0; i<5; i++){
//	std::cout<<names[i]<<std::endl;
//}

int MAXSIZE1= names.size();
int MAXSIZE2= names2.size();
int MAXSIZE3= names3.size();
//MAXSIZE = 1;
//MAXSIZE1 = 1;
//MAXSIZE2 = 1;
//MAXSIZE3 = 1;
TChain* cmc = new TChain("Events");
for(int i=0; i<MAXSIZE1; i++){
	cmc->AddFile( (path+names.at(i)).c_str() );

}
for(int i=0; i<MAXSIZE2; i++){
	cmc->AddFile( (path2+names2.at(i)).c_str() );
}
for(int i=0; i<MAXSIZE3; i++){
	cmc->AddFile( (path3+names3.at(i)).c_str() );
}


cmc->Process("convselMC_noMatch.C");

}
