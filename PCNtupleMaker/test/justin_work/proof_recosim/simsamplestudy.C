#define PCNtuple_cxx
// The class definition in PCNtuple.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("PCNtuple.C")
// root> T->Process("PCNtuple.C","some options")
// root> T->Process("PCNtuple.C+")
//


#include "PCNtuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <iostream>
#include <iomanip>
#include <TLorentzVector.h>
void PCNtuple::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void PCNtuple::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   nPC = new TH1D("nconv","Number of Reconstructed Conversions",21,-0.5,21.5);
   fOutput->Add(nPC);
   vtxdl = new TH1D("vtxdl"," ",50,-0.5,4.5);
   fOutput->Add(vtxdl);
   xy10 = new TH2D("xy10","",200,-10,10,200,-10,10);
   xy25 = new TH2D("xy25","",500,-25,25,500,-25,25);
   r10 = new TH1D("r10","",100,0,10.0);
   r25 = new TH1D("r25","",250,0,25.0);
   fOutput->Add(xy10);
   fOutput->Add(xy25);
   fOutput->Add(r10);
   fOutput->Add(r25);


  //pure sim
  xy10s = new TH2D("xy10s","xy of sim pc",200,-10,10,200,-10,10);
  xy25s = new TH2D("xy25s","xy of sim pc",500,-25,25,500,-25,25);
  r10s = new TH1D("r10s"," r of sim pc",100,0,10);
  r25s = new TH1D("r25s"," r of sim pc", 250,0,25);
  fOutput->Add(xy10s);
  fOutput->Add(xy25s);
  fOutput->Add(r10s);
  fOutput->Add(r25s);



   xy10sim = new TH2D("xy10sim","xy of true sim pc matched to reco pc",200,-10,10,200,-10,10);
   xy10sim_cut = new TH2D("xysim10cut","xy of true sim pc matched to reco pc with cut",200,-10,10,200,-10,10);
   xy25sim = new TH2D("xysim25","xy of true sim pc matched to reco pc",500,-25,25,500,-25,25);
   xy25sim_cut = new TH2D("xysim25cut","xy of true sim pc matched to reco pc with cut",500,-25,25,500,-25,25);
   r10sim = new TH1D("r10sim","r of true sim pc matched to reco pc",100,0,10);
   r10sim_cut = new TH1D("r10simcut","r of true sim pc matched to reco pc with cut",100,0,10);
   r25sim = new TH1D("r25sim", "r of true sim pc matched to reco pc",250,0,25);   
   r25sim_cut = new TH1D("r25simcut","r of true sim pc matched to reco pc with cut",250,0,25);
   
   fOutput->Add(xy10sim);
   fOutput->Add(xy10sim_cut);
   fOutput->Add(xy25sim);
   fOutput->Add(xy25sim_cut);
   fOutput->Add(r10sim);
   fOutput->Add(r10sim_cut);
   fOutput->Add(r25sim);
   fOutput->Add(r25sim_cut);
  
//   numsimphot = new TH1D("numsimphot","num of sim track photons",1,0,1);
//   numsimconv = new TH1D("numsimconv","num of sim conversions",1,0,1);
   nump14     = new TH1D("nump14","num of sim vertex process14",1,0,1);
   numchild   = new TH1D("numchild","num children at sim conv",5,-0.5,4.5);

 //  fOutput->Add(numsimphot);
//   fOutput->Add(numsimconv);
   fOutput->Add(nump14);
   fOutput->Add(numchild);

	numprimarypart = new TH1D( "numprimarypart"," ",101,-0.5,200.5);
	numprimaryphot = new TH1D( "numprimaryphot"," ",101,-0.5,200.5);
	numprimaryq = new TH1D("numprimaryq","",101,-0.5,200.5);
	numsecondaryvtx = new TH1D( "numsecondaryvtx"," ",16,-0.5,15.5);
	secondarypt=new TH1D( "secondarypt"," ",100,0,4);
        secondaryeta=new TH1D( "secondaryeta"," ",100,-10,10);
    fOutput->Add(numprimarypart);
    fOutput->Add(numprimaryphot);
    fOutput->Add(numsecondaryvtx);
    fOutput->Add(secondarypt);
    fOutput->Add(secondaryeta);
    fOUtput->Add(numprimaryq);

        nsimphot_raw = new TH1D("nsimphotraw"," ",201,-0.5,200.5);	
	nsimconv_raw = new TH1D("nsimconvraw"," ",201,-0.5,200.5);
	nsimconv_cut = new TH1D("nsimconvcut"," ",201,-0.5,200.5);
        simphotpt_raw = new TH1D("simphotptraw"," ",40,0,12);
//	simphotpt_cut = new TH1D("simphotptcut"," ",40,0,12);
	simconvpt_raw = new TH1D("simconvptraw"," ",40,0,12);
	simconvpt_cut = new TH1D("simconvptcut"," ",40,0,12);
	simconvR_raw = new TH1D("simconvRraw"," ",300,0,30);
	simconvR_cut = new TH1D("simconvRcut"," ",300,0,30);
	simconvptR_cut= new TH2D("simconvptR"," ",20,0,20,5,0,30);

	numsimconv_cut=new TH1D("numsimconvcut","",1,0,1);
       numsimconv_raw=new TH1D("numsimconvraw","",1,0,1);
	numsimphot_raw =new TH1D("numsimphotraw","",1,0,1);

	fOutput->Add(nsimphot_raw);
	fOutput->Add(nsimconv_raw);
	fOutput->Add(nsimconv_cut);
	fOutput->Add(simphotpt_raw);
	fOutput->Add(simconvpt_raw);
	fOutput->Add(simconvpt_cut);
	fOutput->Add(simconvR_raw);
	fOutput->Add(simconvR_cut);
	fOutput->Add(simconvptR_cut);
	fOutput->Add(numsimconv_cut);
	fOutput->Add(numsimconv_raw);
	fOutput->Add(numsimphot_raw);






   TString option = GetOption();

}

Bool_t PCNtuple::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);

   nPC->Fill(*nConv);
//auto & 

   //std::cout<<"test"<<std::endl;
 
   for(int i=0; i<Conv_vtxdl.GetSize(); i++){
	vtxdl->Fill(Conv_vtxdl[i]);	
   }
  double x,y,r;
   for(int i=0; i< *nConv; i++){
	x = Conv_vtx_X[i];
	y = Conv_vtx_Y[i];
	r = sqrt(x*x + y*y);
        xy10->Fill(x,y);
	xy25->Fill(x,y);
	r10->Fill(r);
	r25->Fill(r);
	
   }  

 std::cout<< std::setprecision(12);
  if(ctr<20){

	std::cout<<"Event "<<*event<<std::endl;
	//find all process type 14, and find all children of that vertex
//	double TOTAL_E=0;
//	TLorentzVector sum;
	for(int i =0; i< *nSimVtx; i++){
		//if( SimVtx_processType[i] == 14 ){
		if( SimVtx_processType[i] == 0){
			std::cout<<"simvtx "<<i<<" ptype "<<SimVtx_processType[i]<< " (x,y,z) "<<SimVtx_x[i]<<" "<<SimVtx_y[i]<<" "<<SimVtx_z[i]<<" ptid "<< SimVtx_simtrk_parent_tid[i];
			//find all sim tracks that map to this vertex
			if( SimVtx_simtrk_parent_tid[i] == -1){
				std::cout<< " ppdg -1"<<std::endl;
			}
			else{
				for(int j=0; j<*nSimTrk; j++){
					if( SimVtx_simtrk_parent_tid[i] == SimTrk_trackId[j]){
						std::cout<<" ppdg "<<SimTrk_pdgId[j]<<std::endl;
					}
				}

			}


			//about to find children -- 0 total E
			TLorentzVector sum;
			TLorentzVector p;
			//assume 0 mass :(
			for(int j=0; j< *nSimTrk; j++){
				if( SimTrk_simvtx_Idx[j] == i ){
					std::cout<<"child of vtx "<<i<<" pdg: "<<SimTrk_pdgId[j]<<" tid: "<<SimTrk_trackId[j]<<std::endl;
					p.SetPtEtaPhiM(SimTrk_pt[j],SimTrk_eta[j],SimTrk_phi[j],0);
					sum = sum + p; 
				}
			}	
			std::cout<<"total approx E: "<<sum.E();	
			std::cout<<std::endl;
		}
	}


}

 //cout number of particles and photons in primary vtx0
 int nphotons=0;
 int nparticles=0;
 int nq=0;
 int nsec=0;
 int pvtx, ptid;
  for(int i=0; i<*nSimTrk; i++){
	if(SimTrk_simvtx_Idx[i] == 0){//this is the primary int
		nparticles++;
		if( SimTrk_pdgId[i] == 22){
			nphotons++;
		}
		if( SimTrk_pdgId[i] != 22){
			nq++;
		}
	}
  }
  numprimarypart->Fill(nparticles);
  numprimaryphot->Fill(nphotons);
 numprimaryq->Fill(nq);
 //loop for plots
 for( int i=1; i<*nSimVtx; i++){
	if( SimVtx_processType[i] ==0 && SimVtx_simtrk_parent_tid[i] == -1){
		nsec++;
		//find the child photons of this vtx
		for(int j=0; j<*nSimTrk; j++){
			if( SimTrk_simvtx_Idx[j] == i){
				secondarypt->Fill(SimTrk_pt[j] );
				secondaryeta->Fill( SimTrk_eta[j] );
			}
		}	

	}	
 }
 numsecondaryvtx->Fill(nsec);


std::vector<bool> vtxmask(*nSimVtx);
std::vector<int> vtxc1(*nSimVtx);//index of sim track child 1
std::vector<int> vtxc2(*nSimVtx);// index of sim track child 1
std::vector<int> childholder;
std::vector<double> costpc(*nSimVtx);
std::vector<double> ptpc(*nSimVtx);
TLorentzVector gamma;
//test code create simvertex mask for 2 child p14s
int numchild=0;
ptid=-1;
for(int i=0; i<*nSimVtx; i++){
  vtxmask[i]= false;
  
  vtxc1[i] = -1;
  vtxc2[i] = -1;
  costpc[i] = -2;
  ptpc[i] = 0;
  if( SimVtx_processType[i] != 14){
        continue;
  }
  else{
        numchild=0;
	ptid = SimVtx_simtrk_parent_tid[i]; 
        for(int j=0; j<*nSimTrk; j++){
		//compute cos theta for acceptance cut
		if( SimTrk_trackId[j] == ptid ){
			gamma.SetPtEtaPhiM(SimTrk_pt[j],SimTrk_eta[j],SimTrk_phi[j],0.0);
			costpc[i] = gamma.CosTheta();						           ptpc[i] = gamma.Pt();
		}	

                if( SimTrk_simvtx_Idx[j] == i){
                        numchild++;
                        childholder.push_back(j);
                }
                if(numchild == 2){
                        vtxmask[i] = true;
                        vtxc1[i] = childholder[0];
                        vtxc2[i] = childholder[1];
                        childholder.clear();
                        break;
                }

        }
        childholder.clear();

  }

}

//analyze with masks 
//count all phot
double abzvtx;
double abcost;
int nphotons_raw=0;
int nphotons_cut=0;
int numconv_raw =0;
int numconv_cut =0;
for(int i=0; i<*nSimTrk; i++){
	if( SimTrk_pdgId[i] == 22){
		numsimphot_raw->Fill(0);	
		nphotons_raw++;
		simphotpt_raw->Fill(SimTrk_pt[i]);
	}

}
nsimphot_raw->Fill(nphotons);

double R;
for(int i=0; i<*nSimVtx; i++){

	if( vtxmask[i]){
		numconv_raw++;
		numsimconv_raw->Fill(0);

		simconvpt_raw->Fill( ptpc[i] );
		R = sqrt( SimVtx_x[i]*SimVtx_x[i] + SimVtx_y[i]*SimVtx_y[i] );
		simconvR_raw->Fill( R);	 

		abzvtx = abs( SimVtx_z[i] );
        	abcost = abs( costpc[i]);

		if( abcost < 0.85 && abzvtx < 25.0){
			numsimconv_cut->Fill(0);
			numconv_cut++;
			simconvpt_cut->Fill(ptpc[i] );
			simconvR_cut->Fill(R);
			simconvptR_cut->Fill(ptpc[i],R);

		}
		

	}
	
	
}
nsimconv_raw->Fill(numconv_raw);
nsimconv_cut->Fill(numconv_cut);
  
  
ctr++;   

   return kTRUE;
}

void PCNtuple::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void PCNtuple::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   TFile* output = new TFile("out.root","RECREATE");
   output->WriteTObject(nPC);
   output->WriteTObject(vtxdl);
   output->WriteTObject(xy10);
   output->WriteTObject(xy25);
   output->WriteTObject(r10);
   output->WriteTObject(r25);
  output->WriteTObject(xy10s);
  output->WriteTObject(xy25s);
  output->WriteTObject(r10s);
  output->WriteTObject(r25s);
 // output->WriteTObject(numsimphot);
 // output->WriteTObject(numsimconv);
  output->WriteTObject(nump14);
  output->WriteTObject(numchild);
 output->WriteTObject(numprimarypart);
 output->WriteTObject(numprimaryphot);
 output->WriteTObject(numsecondaryvtx);
 output->WriteTObject(secondarypt);
 output->WriteTObject(secondaryeta);

 output->WriteTObject(numsimconv_raw);
 output->WriteTObject(numsimconv_cut);
 output->WriteTObject(nsimconv_raw);
 output->WriteTObject(nsimphot_raw);
 output->WriteTObject(simphotpt_raw);
 output->WriteTObject(simconvpt_raw);
 output->WriteTObject(simconvpt_cut);
 output->WriteTObject(simconvR_raw);
 output->WriteTObject(simconvR_cut);
 output->WriteTObject(simconvptR_cut);
 output->WriteTObject(numsimphot_raw);
   output->Write(); 
   output->Close();



}
