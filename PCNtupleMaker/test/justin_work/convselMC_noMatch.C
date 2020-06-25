#define convselMC_noMatch_cxx
// The class definition in convselMC_noMatch.h has been generated automatically
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
// root> T->Process("convselMC_noMatch.C")
// root> T->Process("convselMC_noMatch.C","some options")
// root> T->Process("convselMC_noMatch.C+")
//


#include "convselMC_noMatch.h"
#include <TH2.h>
#include <TStyle.h>

void convselMC_noMatch::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void convselMC_noMatch::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t convselMC_noMatch::Process(Long64_t entry)
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
//test plot
/*
  for(int i=0; i< *nConv; i++){
    xyhist_25.Fill(Conv_vtx_X[i],Conv_vtx_Y[i]);
    xyhist_10.Fill(Conv_vtx_X[i],Conv_vtx_Y[i]);
  }
*/	
double x,y,z,r,phi,theta;
double vxx, vxy, vyy, vzz;
double varsum_r,varsum_phi, cphi,sphi;
double rerr,phierr,zerr;	
double fitprob;

 for(int i=0; i< *nConv; i++){
	        x = Conv_vtx_X[i];
                y = Conv_vtx_Y[i];
                z = Conv_vtx_Z[i];

                r = sqrt( x*x + y*y );
                phi = atan2(y, x);
                theta = atan2( r,z);

		thetatest.Fill(theta);
		phitest.Fill(phi);

               // FillTH2(ind_xyHist, -x, y);
               // FillTH2(ind_rphiHist, r, phi);
               // FillTH2(ind_rzHist, z, r);


                vxx = Conv_vtx_cov_00[i];
                vxy = Conv_vtx_cov_01[i];
                vyy = Conv_vtx_cov_11[i];
                vzz = Conv_vtx_cov_22[i];
                cphi = cos(phi);
                sphi = sin(phi);

                // This is the correct one
                varsum_r   = cphi*cphi*vxx + 2.0*sphi*cphi*vxy + sphi*sphi*vyy;
                varsum_phi = sphi*sphi*vxx - 2.0*sphi*cphi*vxy + cphi*cphi*vyy;



//              phierr = sqrt(varsum_phi); 1/r factor?
//              rerr = r*sqrt(varsum_r); no r factor?
                rerr = sqrt(varsum_r);
                phierr = sqrt(varsum_phi)/r;

                zerr = sqrt(vzz);

 //               FillTH1(ind_rerrHist, rerr);
 //               FillTH1(ind_phierrHist, phierr);
 //               FillTH1(ind_zerrHist, zerr);


                fitprob = TMath::Prob(Conv_vtx_chi2[i], Conv_vtx_ndof[i]);

		if(x<25 && y<25){
			//grahams cuts
			if(rerr < 0.25 && abs(z)<25 && cos(theta)<0.85 && fitprob>0.01){ 
				xyhist_25.Fill(x,y);
			}
		}
		if(x<10 && y<10){
			if(rerr < 0.25 && abs(z)<25 && cos(theta)<0.85 && fitprob>0.01){
                                xyhist_10.Fill(x,y);
                        }

		}
		if(r<25){
                        if(rerr < 0.25 && abs(z)<25 && cos(theta)<0.85 && fitprob>0.01){
                                rhist_25.Fill(r);
				if(*nPV == 10){
					rhist_25_PV10.Fill(r);
				}
				if(*nPV == 25){
					rhist_25_PV25.Fill(r);
				}
				if(*nPV == 40){
					rhist_25_PV40.Fill(r);
				}

                        }

                }


		
		int nhits0 =Conv_nHitsBeforeVtx_Tk0[i];
		int nhits1 =Conv_nHitsBeforeVtx_Tk1[i];
		//begin ghost elim cuts
		if(x<25 && y<25){
                        //grahams cuts
                        if(rerr < 0.25 && abs(z)<25 && cos(theta)<0.85 && fitprob>0.01 && nhits0==0 &&nhits1==0){
                                xyhist_25_gst.Fill(x,y);
                        }
                }
                if(x<10 && y<10){
                        if(rerr < 0.25 && abs(z)<25 && cos(theta)<0.85 && fitprob>0.01 && nhits0==0 &&nhits1==0){
                                xyhist_10_gst.Fill(x,y);
                        }

                }

		if(r<25){
                        if(rerr < 0.25 && abs(z)<25 && cos(theta)<0.85 && fitprob>0.01 && nhits0==0 &&nhits1==0){
                                rhist_25_gst.Fill(r);
                                if(*nPV == 10){
                                        rhist_25_PV10_gst.Fill(r);
                                }
                                if(*nPV == 25){
                                        rhist_25_PV25_gst.Fill(r);
                                }
                                if(*nPV == 40){
                                        rhist_25_PV40_gst.Fill(r);
                                }

                        }

                }

		double sigt0=Conv_dlClosestHitToVtx_sig_Tk0[i];
                double sigt1=Conv_dlClosestHitToVtx_sig_Tk1[i];
                if(x<25 && y<25){
                        //grahams cuts
                        if(rerr < 0.25 && abs(z)<25 && cos(theta)<0.85 && fitprob>0.01 && nhits0==0 &&nhits1==0 && abs(sigt0)<1. && abs(sigt1)<1.){
                                xyhist_25_q1.Fill(x,y);
                        }
                }
                if(x<10 && y<10){
                        if(rerr < 0.25 && abs(z)<25 && cos(theta)<0.85 && fitprob>0.01 && nhits0==0 &&nhits1==0 && abs(sigt0)<1. && abs(sigt1)<1.){
                                xyhist_10_q1.Fill(x,y);
                        }

                }

                if(r<25){
                        if(rerr < 0.25 && abs(z)<25 && cos(theta)<0.85 && fitprob>0.01 && nhits0==0 &&nhits1==0 && abs(sigt0)<1. && abs(sigt1)<1.){
                                rhist_25_q1.Fill(r);
                                if(*nPV == 10){
                                        rhist_25_PV10_q1.Fill(r);
                                }
                                if(*nPV == 25){
                                        rhist_25_PV25_q1.Fill(r);
                                }
                                if(*nPV == 40){
                                        rhist_25_PV40_q1.Fill(r);
                                }

                        }

                }	
		


 }






   return kTRUE;
}

void convselMC_noMatch::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void convselMC_noMatch::Terminate()
{

  TFile *MyFile = new TFile("plotmcnew.root","RECREATE");
  MyFile->WriteObject(&xyhist_25, xyhist_25.GetName());
  MyFile->WriteObject(&xyhist_10, xyhist_10.GetName());
  MyFile->WriteObject(&rhist_25, rhist_25.GetName());
  MyFile->WriteObject(&rhist_25_PV10, rhist_25_PV10.GetName());
  MyFile->WriteObject(&rhist_25_PV25, rhist_25_PV25.GetName());
  MyFile->WriteObject(&rhist_25_PV40, rhist_25_PV40.GetName());

 MyFile->WriteObject(&xyhist_25_gst, xyhist_25_gst.GetName());
  MyFile->WriteObject(&xyhist_10_gst, xyhist_10_gst.GetName());
  MyFile->WriteObject(&rhist_25_gst, rhist_25_gst.GetName());
  MyFile->WriteObject(&rhist_25_PV10_gst, rhist_25_PV10_gst.GetName());
  MyFile->WriteObject(&rhist_25_PV25_gst, rhist_25_PV25_gst.GetName());
  MyFile->WriteObject(&rhist_25_PV40_gst, rhist_25_PV40_gst.GetName());

 MyFile->WriteObject(&xyhist_25_q1, xyhist_25_q1.GetName());
   MyFile->WriteObject(&xyhist_10_q1, xyhist_10_q1.GetName());
  MyFile->WriteObject(&rhist_25_q1, rhist_25_q1.GetName());
  MyFile->WriteObject(&rhist_25_PV10_q1, rhist_25_PV10_q1.GetName());
  MyFile->WriteObject(&rhist_25_PV25_q1, rhist_25_PV25_q1.GetName());
  MyFile->WriteObject(&rhist_25_PV40_q1, rhist_25_PV40_q1.GetName());



  MyFile->Write(); 
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  //TCanvas* c1 = new TCanvas();
 // xyhist_25.Draw("COLZ");
  
 // TCanvas* c2 = new TCanvas();
 // xyhist_10.Draw("COLZ");
 // thetatest.Draw();
 // TCanvas* c3 = new TCanvas();
 //  phitest.Draw();
}
