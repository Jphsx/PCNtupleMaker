#define convselMC_cxx
// The class definition in convselMC.h has been generated automatically
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
// root> T->Process("convselMC.C")
// root> T->Process("convselMC.C","some options")
// root> T->Process("convselMC.C+")
//


#include "convselMC.h"
#include <TH2.h>
#include <TStyle.h>

void convselMC::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void convselMC::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t convselMC::Process(Long64_t entry)
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

   
   for(int i=0; i< *nGenPart; i++){
	if(GenPart_pdgId[i] == 22 ){
	 	int ndaughter = GenPart_nDaughter[i];
/* blockstart 	if(ndaughter >0 ){
			std::cout<<" pdg 22 with ndaughter "<<ndaughter<< " ; "<< GenPart_pdgId[ GenPart_daughter0Idx[i] ];
		
		
		if(ndaughter > 1 ){
			std::cout<<" "<<GenPart_pdgId[ GenPart_daughter1Idx[i] ];
		}
		
			std::cout<<std::endl;
		}
 blockend*/
		if(ndaughter == 2){
			std::cout<<"Found Gen conversion in event: "<<*event<<std::endl;
			std::cout<<"Event contains nConv: "<<*nConv<<std::endl;
			int d1pdg= GenPart_pdgId[ GenPart_daughter0Idx[i] ];
			int d2pdg= GenPart_pdgId[ GenPart_daughter1Idx[i] ];
			int d1idx = GenPart_daughter0Idx[i];
			int d2idx = GenPart_daughter1Idx[i];
			if( abs(d1pdg) == 11 ){
				std::cout<<"Converted Photon and Daughter Vtx: pdg (x,y,z) (pt,eta,phi)"<<std::endl;
				std::cout<<" 22 "<<GenPart_vtx_X[i]<<" "<<GenPart_vtx_Y[i]<<" "<<GenPart_vtx_Z[i]<<" "<<GenPart_pt[i]<<" "<<GenPart_eta[i]<<" "<<GenPart_phi[i]<<std::endl;
				std::cout<<d1pdg<<" "<<GenPart_vtx_X[d1idx]<<" "<<GenPart_vtx_Y[d1idx]<<" "<<GenPart_vtx_Z[d1idx]<<" "<<GenPart_pt[d1idx]<<" "<<GenPart_eta[d1idx]<<" "<<GenPart_phi[d1idx]<<std::endl; 
				std::cout<<d2pdg<<" "<<GenPart_vtx_X[d2idx]<<" "<<GenPart_vtx_Y[d2idx]<<" "<<GenPart_vtx_Z[d2idx]<<" "<<GenPart_pt[d2idx]<<" "<<GenPart_eta[d2idx]<<" "<<GenPart_phi[d2idx]<<std::endl;
		
				std::cout<<"Scanning for nearby conversion vertices"<<std::endl;
				double gx,gy,gz,gpt;
				gx = GenPart_vtx_X[i];
				gy = GenPart_vtx_Y[i];
				gz = GenPart_vtx_Z[i];
				gpt = GenPart_pt[i];
				double vx,vy,vz,cpx,cpy,cpt,refcpx,refcpy,refcpt;
				double d;
				for(int j=0; j<*nConv; j++){
					vx= Conv_vtx_X[j];
					vy= Conv_vtx_Y[j];
					vz= Conv_vtx_Z[j];
					cpx = Conv_pairMomentum_Px[j];
					cpy = Conv_pairMomentum_Py[j];
					cpt = sqrt( cpx*cpx + cpy*cpy);
					refcpx = Conv_refittedPairMomentum_Px[j];
					refcpy = Conv_refittedPairMomentum_Py[j];
					refcpt = sqrt(refcpx*refcpx + refcpy*refcpy);						

					d = sqrt( (vx-gx)*(vx-gx) + (vy-gy)*(vy-gy) + (vz-gz)*(vz-gz) );
					//if(vx-gx <0.2 && vy-gy<0.2 &&vz-gz<0.2){{
			//		if(d<5){
			//			std::cout<<"Found nearby conversion Vertex: "<<vx<<" "<<vy<<" "<<vz<<std::endl;
			//		}
					std::cout<<"Reconstructed conversion vertex (x,y,z) (conv_pt, refitted_conv_pt): "<<vx<<" "<<vy<<" "<<vz<<" "<<cpt<<" "<<refcpt<<std::endl;
					std::cout<<"T0 pt,eta,phi "<<Conv_Tk0_pt[j]<<" "<<Conv_Tk0_eta[j]<<" "<<Conv_Tk0_phi[j] <<std::endl;
					std::cout<<"T1 pt,eta,phi "<<Conv_Tk1_pt[j]<<" "<<Conv_Tk1_eta[j]<<" "<<Conv_Tk1_phi[j]<<std::endl<<std::endl;
											
	
				}	
				

			} 


			
		}
		
	
	}


 	/*if( abs(GenPart_pdgId[i]) == 11){
		int motherpdg = GenPart_pdgId[ GenPart_genPartIdxMother[i] ];
		if(motherpdg == 22){
			std::cout<<" pdg "<<GenPart_pdgId[i]<<" with mother 22 "<<std::endl;;
		}	
		if(GenPart_isConvertedPhoton[i]){
			std::cout<<" pdg "<<GenPart_pdgId[i]<<" is converted photon"<<std::endl;
		}
	} */

   }//end genpart scan


/*	std::cout<<"begin new scan, match any conv vtx to gen vtx"<<std::endl<<std::endl;
	double vx,vy,vz;
	double gx,gy,gz;
	double d;
	for(int i=0; i<*nConv; i++){
		vx= Conv_vtx_X[i];
                vy= Conv_vtx_Y[i];
                vz= Conv_vtx_Z[i];
		std::cout<<"conv vtx: "<<vx<<" "<<vy<<" "<<vz<<std::endl;
		for(int j=0; j< *nGenPart; j++){
			 gx = GenPart_vtx_X[j];
                         gy = GenPart_vtx_Y[j];
                         gz = GenPart_vtx_Z[j];
			d = sqrt( (vx-gx)*(vx-gx) + (vy-gy)*(vy-gy) + (vz-gz)*(vz-gz) );
			if(d<1){
			 std::cout<<"Found " <<gx<<" "<<gy<<" "<<gz<<std::endl;
			}

		}		
				


	}
*/
if(*nConv > 0){
std::cout<<"Event: "<<*event<<" with nConv: "<<*nConv<<std::endl;
}

double pt1,eta1,phi1,pt2,eta2,phi2;
double gpt,geta,gphi, gpdg;
double dr,dp;
//Begin scan with DR matching for conv tracks
for(int i=0; i<*nConv; i++){
                pt1= Conv_Tk0_pt[i];
                eta1= Conv_Tk0_eta[i];
                phi1= Conv_Tk0_phi[i];
		pt2= Conv_Tk1_pt[i];
                eta2= Conv_Tk1_eta[i];
                phi2= Conv_Tk1_phi[i];
		std::cout<<"Conversion "<<i<<std::endl;
		std::cout<<"c vtx (x,y,z) : "<<Conv_vtx_X[i]<<" "<<Conv_vtx_Y[i]<<" "<<Conv_vtx_Z[i]<<std::endl;
                std::cout<<"Tk0 (pt,eta,phi) : "<<pt1<<" "<<eta1<<" "<<phi1<<std::endl;
		std::cout<<"Tk1 (pt,eta,phi) : "<<pt2<<" "<<eta2<<" "<<phi2<<std::endl;
		std::cout<<"Matches for Tk0 "<<std::endl;
                for(int j=0; j< *nGenPart; j++){
                         gpt = GenPart_pt[j];
                         geta = GenPart_eta[j];
                         gphi = GenPart_phi[j];
			 gpdg = GenPart_pdgId[j];
                        dr = dR(geta,gphi,eta1,phi1);
			dp = abs(gpt - pt1);
                        if((dr<0.1 && dp<0.2 && gpdg != 21) || (abs(gpdg)==11 && dr<2 && dp<0.6  ) ){
                         std::cout<<"Found gen Track (pdg, dr, pt, eta, phi) "<<gpdg<<" | "<<dr<<" "<<dp<<" | "<<gpt<<" "<<geta<<" "<<gphi<<std::endl;
			 int midx = GenPart_genPartIdxMother[j];

			 std::cout<<"Ancestry (Idx, pdg): "<<midx<<" "<<GenPart_pdgId[midx]<<" | "<<GenPart_vtx_X[midx]<<" "<<GenPart_vtx_Y[midx]<<" "<<GenPart_vtx_Z[midx]<<std::endl;
                        }

                }
		std::cout<<std::endl;
		std::cout<<"Matches for Tk1 "<<std::endl;
                for(int j=0; j< *nGenPart; j++){
                         gpt = GenPart_pt[j];
                         geta = GenPart_eta[j];
                         gphi = GenPart_phi[j];
                         gpdg = GenPart_pdgId[j];
                        dr = dR(geta,gphi,eta2,phi2);
			dp = abs(gpt - pt2);
                        if((dr<0.2 && dp<0.3 && gpdg != 21) || (abs(gpdg)==11 && dr<2 && dp<0.6)    ){
                         std::cout<<"Found gen Track (pdg, dr, pt, eta, phi) " <<gpdg<<" | "<<dr<<" "<<dp<<" | "<<gpt<<" "<<geta<<" "<<gphi<<std::endl;
			 std::cout<<"match candidate vtx (x,y,z)" <<GenPart_vtx_X[j]<<" "<<GenPart_vtx_Y[j]<<" "<<GenPart_vtx_Z[j]<<std::endl;
			 int midx = GenPart_genPartIdxMother[j];
			 std::cout<<"Ancestry (Idx, pdg), (x,y,z): "<<midx<<" "<<GenPart_pdgId[midx]<<" | "<<GenPart_vtx_X[midx]<<" "<<GenPart_vtx_Y[midx]<<" "<<GenPart_vtx_Z[midx]<<std::endl; 
                        }

                }        
		std::cout<<std::endl;       
                                

}





   return kTRUE;
}

void convselMC::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void convselMC::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
