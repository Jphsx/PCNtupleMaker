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
   //do we need to add these ints?
//   fOutput->Add(n_sim_pc);
//   fOutput->Add(n_sim_pc_rec);
//   fOutput->Add(n_sim_pc_rec_cut);
   



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
  

//if(ctr<20){ 
//	std::cout<<"test event "<<*event<<std::endl;   
//sim stuff with no cut
   //identify true sim photon  conversions
   //count the number of simtrk from a vtx, if 2 see if they are e,p pair
   //look at process type also
   int ivtx,jvtx;
  // bool ep,em;
   int p_tid;
   double sx,sy,sz, sr;
   //since some vtx seem identical...?
  // double ivtxposX, ivtxposY, ivtxposZ;
   std::vector<int> vtxidxs;
   for(int i=0; i< *nSimTrk; i++){
	ivtx = SimTrk_simvtx_Idx[i];
	vtxidxs.push_back(i);
	for(int j=i+1; j< *nSimTrk; j++){
		jvtx = SimTrk_simvtx_Idx[j];
		if(ivtx == jvtx){
			vtxidxs.push_back(j);
		}
	
		
	}
	//look through vtx idx, are there 2 children?
	//do triplets happen?
	//is the mother track a photon
	if( vtxidxs.size() == 2 || vtxidxs.size() == 3){
//		std::cout<<"invtxids ";
		//get parent tid of sim vertex at ivtx	
		p_tid = SimVtx_simtrk_parent_tid[ivtx];
//		std::cout<<"p_tid "<<p_tid<<" ";
		for(int k=0; k<*nSimTrk; k++){
			if( SimTrk_trackId[k] == p_tid ){
				//std::cout<<"found match parent";	
				if( SimTrk_pdgId[k] == 22 ){
					//print out this conversion
//					std::cout<<"found conversion event:"<< *event <<std::endl;
//					std::cout<<"number of vtx children: "<<vtxidxs.size()<<" ";
//					std::cout<<"ptid "<<p_tid<<" ";
//					for(int l=0; l<vtxidxs.size(); l++){
//						std::cout<< SimTrk_pdgId[ vtxidxs[l]]<<" ";
//					}
//					std::cout<<std::endl;	

					//in here we have found a true sim conversion
					n_sim_pc++;
					//do analysis to sim
					//ivtx is the location of conv
				       sx= SimVtx_x[ivtx];
				       sy= SimVtx_y[ivtx];
				       sz= SimVtx_z[ivtx];
				       sr= sqrt(sx*sx + sy*sy);
					xy10s->Fill(sx,sy);
					xy25s->Fill(sx,sy);
					r10s->Fill(sr);
					r25s->Fill(sr);	
				}
			}
		}		
		
	}
	vtxidxs.clear();
   }// end looking for true sim
//endctr
//}
ctr++;
 // for(int i=0; i< *nSimVtx; i++){
	
//   }   
   

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

   output->Write(); 
   output->Close();



}
