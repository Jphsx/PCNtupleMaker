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
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>
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
 
  thetacheck= new TH1D("thetacheck","",20,-3.14,3.14);
  fOutput->Add(thetacheck);
/*

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
  
   numsimphot = new TH1D("numsimphot","num of sim track photons",1,0,1);
   numsimconv = new TH1D("numsimconv","num of sim conversions",1,0,1);
   nump14     = new TH1D("nump14","num of sim vertex process14",1,0,1);
   numchild   = new TH1D("numchild","num children at sim conv",5,-0.5,4.5);

   fOutput->Add(numsimphot);
   fOutput->Add(numsimconv);
   fOutput->Add(nump14);
   fOutput->Add(numchild);

   p14 = new TH1D("maskp14","p14 from mask",1,0,1);
   fOutput->Add(p14);
   npc = new TH1D("npc","number of reco conv",1,0,1);
   fOutput->Add(npc);

   recomatchp14 = new TH1D("recomatchp14", "reco conv matched to sim conv",1,0,1);
   recomatchany = new TH1D("recomatchany", "reco conv matched to sim vtx",1,0,1);
   recomatchp14alltrk = new TH1D("recomatchp14alltrk","reco conv matched to sim con and children",1,0,1);
   recomatchp14onetrk = new TH1D("recomatchp14onetrk","reco conv matched to sim conv and 1 child",1,0,1);
   recomatchp14dup = new TH1D("recomatchp14dup","number of recos matched single sim conv",7,-0.5,6.5);
   recomatchbg = new TH1D("recomatchbg", "reco conv matched to !sim conv", 1,0,1);
   fOutput->Add(recomatchp14);
   fOutput->Add(recomatchany);
   fOutput->Add(recomatchp14alltrk);
   fOutput->Add(recomatchp14onetrk);
   fOutput->Add(recomatchp14dup);


   p14_pt = new TH1D("maskp14pt","p14 from mask",21,-.5,20.5);
   fOutput->Add(p14_pt);
   npc_pt = new TH1D("npcpt","number of reco conv",21,-.5,20.5);
   fOutput->Add(npc_pt);

   recomatchp14_pt = new TH1D("recomatchp14pt", "reco conv matched to sim conv",21,-0.5,20.5);
   recomatchany_pt = new TH1D("recomatchanypt", "reco conv matched to sim vtx",21,-0.5,20.5);
   recomatchp14alltrk_pt = new TH1D("recomatchp14alltrkpt","reco conv matched to sim con and children",21,-0.5,20.5);
   recomatchp14onetrk_pt = new TH1D("recomatchp14onetrkpt","reco conv matched to sim conv and 1 child",21,-0.5,20.5);
     recomatchbg_pt = new TH1D("recomatchbgpt", "reco conv matched to !sim conv", 21,-0.5,20.5);
   fOutput->Add(recomatchp14_pt);
   fOutput->Add(recomatchany_pt);
   fOutput->Add(recomatchp14alltrk_pt);
   fOutput->Add(recomatchp14onetrk_pt);
  

  
   p14_r = new TH1D("maskp14r","p14 from mask",6,0,30);
   fOutput->Add(p14_r);
   npc_r = new TH1D("npcr","number of reco conv",6,0,30);
   fOutput->Add(npc_r);


   recomatchp14_r = new TH1D("recomatchp14r", "reco conv matched to sim conv",6,0,30);
   recomatchany_r = new TH1D("recomatchanyr", "reco conv matched to sim vtx",6,0,30);
   recomatchp14alltrk_r = new TH1D("recomatchp14alltrkr","reco conv matched to sim con and children",6,0,30);
   recomatchp14onetrk_r = new TH1D("recomatchp14onetrkr","reco conv matched to sim conv and 1 child",6,0,30);
     recomatchbg_r = new TH1D("recomatchbgr", "reco conv matched to !sim conv", 6,0,20.5);
   fOutput->Add(recomatchp14_r);
   fOutput->Add(recomatchany_r);
   fOutput->Add(recomatchp14alltrk_r);
   fOutput->Add(recomatchp14onetrk_r);



   p14_ptr = new TH2D("maskp14ptr","p14 from mask",21,-.5,20.5,6,0,30);
   fOutput->Add(p14_ptr);
   npc_ptr = new TH2D("npcptr","number of reco conv",21,-.5,20.5,6,0,30);
   fOutput->Add(npc_ptr);


   recomatchp14_ptr = new TH2D("recomatchp14ptr", "reco conv matched to sim conv",21,-.5,20.5,6,0,30);
   recomatchany_ptr = new TH2D("recomatchanyptr", "reco conv matched to sim vtx",21,-.5,20.5,6,0,30);
   recomatchp14alltrk_ptr = new TH2D("recomatchp14alltrkptr","reco conv matched to sim con and children",21,-.5,20.5,6,0,30);
   recomatchp14onetrk_ptr = new TH2D("recomatchp14onetrkptr","reco conv matched to sim conv and 1 child",21,-.5,20.5,6,0,30);
     recomatchbg_ptr = new TH2D("recomatchbgptr", "reco conv matched to !sim conv",21,-.5,20.5, 6,0,20.5);
   fOutput->Add(recomatchp14_ptr);
   fOutput->Add(recomatchany_ptr);
   fOutput->Add(recomatchp14alltrk_ptr);
   fOutput->Add(recomatchp14onetrk_ptr);
*/


   totnpc_raw = new TH1D("totnpcraw","",1,0,1);
   totnpc_cut = new TH1D("totnpccut","",1,0,1);
   xypc_cut = new TH2D("xypccut","",500,-25,25,500,-25,25);
   rpc_cut = new TH1D("rpccut","",250,0,25.0);
   
   ptdenom = new TH1D( "ptdenom","",20,0,20);
   ptnum = new TH1D( "ptnum","",20,0,20);
   rdenom = new TH1D( "rdenom","",30,0,30);
   rnum = new TH1D( "rnum","",30,0,30);
   ptdenom->Sumw2(true);
   ptnum->Sumw2(true);
   rdenom->Sumw2(true);
   rnum->Sumw2(true);
   fOutput->Add(ptdenom);
   fOutput->Add(ptnum);
   fOutput->Add(rdenom);
   fOutput->Add(rnum);
   fOutput->Add(totnpc_raw);
   fOutput->Add(totnpc_cut);
   fOutput->Add(xypc_cut);
   fOutput->Add(rpc_cut);
  
   xypc_raw = new TH2D("xypcraw","",500,-25,25,500,-25,25);
   rpc_raw = new TH1D("rpccut","",250,0,25.0);
   fOutput->Add(xypc_raw);
   fOutput->Add(rpc_raw);


   ptrnum = new TH2D("ptrnum","",20,0,20,30,0,30);
   ptrdenom = new TH2D("ptrdenom","",20,0,20,30,0,30);
   ptrnum->Sumw2(true);
   ptrdenom->Sumw2(true);
   fOutput->Add(ptrdenom);
   fOutput->Add(ptrnum);
  
   totsig_cut = new TH1D("totsigcut","",1,0,1);
   totbg_cut= new TH1D("totbgcut","",1,0,1);
   totunmatch_cut= new TH1D("totunmatchcut","",1,0,1);
   totsig_raw= new TH1D("totsigraw","",1,0,1);
   totbg_raw= new TH1D("totbgraw","",1,0,1);
   totunmatch_raw= new TH1D("totunmatchraw","",1,0,1);

   fOutput->Add(totsig_cut);
   fOutput->Add(totbg_cut);
   fOutput->Add(totunmatch_cut);
   fOutput->Add(totsig_raw);
   fOutput->Add(totbg_raw);
   fOutput->Add(totunmatch_raw);




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
/*
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
*/

/////MASKING//////////////////////
std::vector<bool> vtxmask(*nSimVtx);
std::vector<TLorentzVector*> photmask(*nSimVtx);
std::vector<int> vtxc1(*nSimVtx);//index of sim track child 1
std::vector<int> vtxc2(*nSimVtx);// index of sim track child 1
std::vector<int> childholder;
//test code create simvertex mask for 2 child p14s
int numchild=0;
int ptid=-1;
for(int i=0; i<*nSimVtx; i++){
  vtxmask[i]= false;
  vtxc1[i] = -1;
  vtxc2[i] = -1;
  photmask[i] = new TLorentzVector();
  if( SimVtx_processType[i] != 14){
	continue;
  }
  else{
	numchild=0;
	ptid = SimVtx_simtrk_parent_tid[i];
	for(int j=0; j<*nSimTrk; j++){
		if( SimTrk_trackId[j] == ptid ){
			photmask[i]->SetPtEtaPhiM(SimTrk_pt[j],SimTrk_eta[j],SimTrk_phi[j],0.0);
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
///END MASK//////


	const double RERRCUT = 0.25;
	const double COSTCUT = 0.85;
	const double ZCUT = 25.0;
	const double FITPROBCUT = 0.010;
    const double MASSCUT = 0.15;
double vxx,vxy,vyy,vzz; //variances
double varsum_r, varsum_phi; //intermediate calculation variables
double rerr,phierr,zerr;//errors
double sphi,cphi;
double fitprob;

// We now have various "centers" to compare to for radial coordinates.

	//beam pipe displacement (in cm) from Anna's DPF2019 talk
    const double x0bpdata =  0.171;
    const double y0bpdata = -0.176;
    const double x0bpmc = 0.0;
    const double y0bpmc = 0.0;

    //BPIX center displacement (in cm) (Anna's 16-Dec-2019 talk page 3)
    const double x0data =  0.086;
    const double y0data = -0.102;
    const double x0mc = 0.0;
    const double y0mc = 0.0;

    //Pixel support displacement (in cm) (Anna's December talk page 4)
    const double x0psdata = -0.080;
    const double y0psdata = -0.318;
    const double x0psmc = 0.0;
    const double y0psmc = 0.0;


double x0,y0,x0bp,y0bp,x0ps,y0ps;

       x0bp = x0bpmc;
       y0bp = y0bpmc;
       x0 = x0mc;
       y0 = y0mc;
       x0ps = x0psmc;
       y0ps = y0psmc;
double x,y,z,r,theta,phi,rnominal,rps,rho,phip;

//make num for eff here
for(int i=0; i<*nConv; i++){
//        vcuts.push_back(false);
		x = Conv_vtx_X[i];
		y = Conv_vtx_Y[i];
		z = Conv_vtx_Z[i];
        double px = Conv_refittedPairMomentum_Px[i];
        double py = Conv_refittedPairMomentum_Py[i];
        double pz = Conv_refittedPairMomentum_Pz[i];
        double E = sqrt(px*px + py*py + pz*pz);
        double pt = sqrt(px*px + py*py); 
	TVector3 temp(px,py,pz);
		r = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
		phi = atan2(y-y0, x-x0);
		//theta = atan2(pt,pz);
		theta = temp.Theta();
		thetacheck->Fill(theta);

		rho  =  sqrt( (x-x0bp)*(x-x0bp) + (y-y0bp)*(y-y0bp)) ;
  	        phip =  atan2(y-y0bp, x-x0bp);

		rps = sqrt( (x-x0ps)*(x-x0ps) + (y-y0ps)*(y-y0ps) );

                rnominal = sqrt( x*x + y*y );

		vxx = Conv_vtx_cov_00[i];
		vxy = Conv_vtx_cov_01[i];
		vyy = Conv_vtx_cov_11[i];
		vzz = Conv_vtx_cov_22[i];

		cphi = cos(phi);
		sphi = sin(phi);

		// This is the correct one
		varsum_r   = cphi*cphi*vxx + 2.0*sphi*cphi*vxy + sphi*sphi*vyy;
		varsum_phi = sphi*sphi*vxx - 2.0*sphi*cphi*vxy + cphi*cphi*vyy;
		rerr = sqrt(varsum_r);
		phierr = sqrt(varsum_phi)/r;
		zerr = sqrt(vzz);

	 	fitprob = TMath::Prob(Conv_vtx_chi2[i], 3);


		int nhits0 = Conv_nHitsBeforeVtx_Tk0[i]; 
		int nhits1 = Conv_nHitsBeforeVtx_Tk1[i];
		bool nohits = false;
		if(nhits0 == 0 && nhits1==0){
			nohits = true;
		}


		//plots
		totnpc_raw->Fill(0);
		 if( vtxmask[ Conv_convVtxIdx[i] ] && Conv_vtxdl[i] <2 ){
                                        totsig_raw->Fill(0);    
                                }
                                if( !vtxmask[ Conv_convVtxIdx[i] ] && Conv_vtxdl[i] <2){
                                        totbg_raw->Fill(0);
                                }
                                if( !vtxmask[ Conv_convVtxIdx[i] ] && Conv_vtxdl[i] >=2){
                                        totunmatch_raw->Fill(0);
                                }





		if( (rerr < RERRCUT) && (abs(cos(theta)) < COSTCUT) && (abs(z) < ZCUT) && (fitprob > FITPROBCUT) && (nohits) ){
			//we have passed cuts
				totnpc_cut->Fill(0);
				xypc_cut->Fill(x,y);
				rpc_cut->Fill(r);

				if( vtxmask[ Conv_convVtxIdx[i] ] && Conv_vtxdl[i] <2 ){
					totsig_cut->Fill(0);					
					ptnum->Fill(pt);
					rnum->Fill(r);
					ptrnum->Fill(pt,r);	
				}
				if( !vtxmask[ Conv_convVtxIdx[i] ] && Conv_vtxdl[i] <2){
					totbg_cut->Fill(0);
				}
				if( !vtxmask[ Conv_convVtxIdx[i] ] && Conv_vtxdl[i] >=2){
					totunmatch_cut->Fill(0);
				}
		}

}
//construct denoms
//
double rs;
for(int i=0; i<*nSimVtx; i++){

	//after cuts denom eff
	if( vtxmask[i] ){
		if( (abs( photmask[i]->CosTheta() )< COSTCUT) && (abs( SimVtx_x[i])  < ZCUT) ){
			ptdenom->Fill( photmask[i]->Pt() );
			rs = sqrt(SimVtx_x[i]*SimVtx_x[i] + SimVtx_y[i]*SimVtx_y[i]) ;
			rdenom->Fill( rs);
			ptrdenom->Fill(photmask[i]->Pt(),rs);

		}

	}
}

/////////////////////begin eff tests
//loop over reco convs
//   TTreeReaderArray<Int_t> Conv_Tk0_Idx = {fReader, "Conv_Tk0_Idx"};
//   TTreeReaderArray<Int_t> Conv_Tk1_Idx = {fReader, "Conv_Tk1_Idx"};
//   TTreeReaderArray<Int_t> Conv_convVtxIdx = {fReader, "Conv_convVtxIdx"};
/*
double maxdl = 2.0;
int vidx = -1;
bool t1=false;
bool t0=false;
double pt;//r is already defined
for(int i=0; i<*nConv; i++){
	vidx = Conv_convVtxIdx[i];
   	npc->Fill(0);
	pt = sqrt( Conv_refittedPairMomentum_Px[i]*Conv_refittedPairMomentum_Px[i] + Conv_refittedPairMomentum_Py[i]*Conv_refittedPairMomentum_Py[i]);
        r = sqrt( Conv_vtx_X[i]*Conv_vtx_X[i] + Conv_vtx_Y[i]*Conv_vtx_Y[i] );
	npc_pt->Fill(pt);
  	npc_r->Fill(r);
	npc_ptr->Fill(pt,r);
	if( (Conv_Tk1_Idx[i] == vtxc1[vidx]) || (Conv_Tk1_Idx[i] == vtxc2[vidx])){
		t1=true;
	}
	if( (Conv_Tk0_Idx[i] == vtxc1[vidx]) || (Conv_Tk0_Idx[i] == vtxc2[vidx])){
		t0=true;
	}

	if( vtxmask[vidx] && Conv_vtxdl[i] <= maxdl ){
		recomatchp14->Fill(0);
		recomatchp14_pt->Fill(pt);
		recomatchp14_r->Fill(r);
		recomatchp14_ptr->Fill(pt,r);
		//do tracks match
		if( t1 && t0){
			recomatchp14alltrk->Fill(0);
			recomatchp14alltrk_pt->Fill(pt);
			recomatchp14alltrk_r->Fill(r);
			recomatchp14alltrk_ptr->Fill(pt,r);
		}
		if( t1 || t0){
			recomatchp14onetrk->Fill(0);
			recomatchp14onetrk_pt->Fill(pt);
			recomatchp14onetrk_r->Fill(r);
			recomatchp14onetrk_ptr->Fill(pt,r);
		}


		

	}	
	if( Conv_vtxdl[i] <= maxdl ){
		recomatchany->Fill(0);	
		recomatchany_pt->Fill(pt);
		recomatchany_r->Fill(r);
		recomatchany_ptr->Fill(pt,r);
	}
	if( !vtxmask[vidx] && Conv_vtxdl[i] <= maxdl ){
		recomatchbg->Fill(0);
		recomatchbg_pt->Fill(pt);
		recomatchbg_r->Fill(r);
		recomatchbg_ptr->Fill(pt,r);
	}



 
}
*/
/*
//find dups
std::vector<int> vtxdups(*nSimVtx);
for(int i=0; i<vtxdups.size();i++){
	vtxdups.at(i)=0;
}
for(int i=0; i<*nConv; i++){
	vidx = Conv_convVtxIdx[i];
	if( vtxmask[vidx] && (Conv_vtxdl[i] <= maxdl) ){
		vtxdups[vidx]= vtxdups[vidx]+1;
	}

}
TLorentzVector part1;
TLorentzVector part2;
TLorentzVector part3;
double me = 0.511 * 10e-3; //electron mass units gev
for(int i=0; i<*nSimVtx; i++){
	if(vtxmask[i]){
		int j1 = vtxc1[i];
		int j2 = vtxc2[i];
		r = sqrt(SimVtx_x[i]*SimVtx_x[i] + SimVtx_y[i]*SimVtx_y[i]);
		part1.SetPtEtaPhiM(SimTrk_pt[j1],SimTrk_eta[j1],SimTrk_phi[j1],me);
		part2.SetPtEtaPhiM(SimTrk_pt[j2],SimTrk_eta[j2],SimTrk_phi[j2],me);
		part3 = part1+part2;
	
		recomatchp14dup->Fill(vtxdups[i]);
		p14->Fill(0);
		p14_pt->Fill(part3.Pt());
		p14_r->Fill(r);
		p14_ptr->Fill(part3.Pt(),r);
	}
}

*/







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
/*   output->WriteTObject(nPC);
   output->WriteTObject(vtxdl);
   output->WriteTObject(xy10);
   output->WriteTObject(xy25);
   output->WriteTObject(r10);
   output->WriteTObject(r25);
  output->WriteTObject(xy10s);
  output->WriteTObject(xy25s);
  output->WriteTObject(r10s);
  output->WriteTObject(r25s);
  output->WriteTObject(numsimphot);
  output->WriteTObject(numsimconv);
  output->WriteTObject(nump14);
  output->WriteTObject(numchild);
   output->WriteTObject(recomatchp14);
   output->WriteTObject(recomatchany);
   output->WriteTObject(recomatchp14alltrk);
   output->WriteTObject(recomatchp14onetrk);
   output->WriteTObject(recomatchp14dup);
   output->WriteTObject(recomatchbg);
  output->WriteTObject(p14);
  output->WriteTObject(npc);
 
    output->WriteTObject(recomatchp14_pt);
   output->WriteTObject(recomatchany_pt);
   output->WriteTObject(recomatchp14alltrk_pt);
   output->WriteTObject(recomatchp14onetrk_pt);
   output->WriteTObject(recomatchbg_pt);
  output->WriteTObject(p14_pt);
  output->WriteTObject(npc_pt);

output->WriteTObject(recomatchp14_r);
   output->WriteTObject(recomatchany_r);
   output->WriteTObject(recomatchp14alltrk_r);
   output->WriteTObject(recomatchp14onetrk_r);
   output->WriteTObject(recomatchbg_r);
  output->WriteTObject(p14_r);
  output->WriteTObject(npc_r);

output->WriteTObject(recomatchp14_ptr);
   output->WriteTObject(recomatchany_ptr);
   output->WriteTObject(recomatchp14alltrk_ptr);
   output->WriteTObject(recomatchp14onetrk_ptr);
   output->WriteTObject(recomatchbg_ptr);
  output->WriteTObject(p14_ptr);
  output->WriteTObject(npc_ptr);
*/
   output->WriteTObject(totnpc_raw);
   output->WriteTObject(totnpc_cut);
   output->WriteTObject(xypc_cut);
   output->WriteTObject(rpc_cut);
   output->WriteTObject(xypc_raw);
   output->WriteTObject(rpc_raw);
   output->WriteTObject(ptdenom);
   output->WriteTObject(ptnum);
   output->WriteTObject(rdenom);
   output->WriteTObject(rnum);
   output->WriteTObject(ptrnum);
   output->WriteTObject(ptrdenom);  

   output->WriteTObject(totsig_cut);
   output->WriteTObject(totbg_cut);
   output->WriteTObject(totunmatch_cut);
   output->WriteTObject(totsig_raw);
   output->WriteTObject(totbg_raw);
   output->WriteTObject(totunmatch_raw);


   output->WriteTObject(thetacheck);
   output->Write(); 
   output->Close();



}
