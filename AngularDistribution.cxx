//g++ -O2 -std=c++0x AngularDistribution.cxx -L$PROGDIR/SharcAnalysis -lSharcAnalysis -I$PROGDIR/SharcAnalysis -L$PROGDIR/TigressAnalysis -lTigressAnalysis -I$PROGDIR/TigressAnalysis  `grsi-config --cflags --all-libs --root` -o angdist.out

#include <Rtypes.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TChain.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TBox.h>


#include <string>
#include <sstream>
#include <ctype.h>
#include <map>
#include <algorithm>

#include "TSharcAnalysis.h"
#include "TTigressAnalysis.h"

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//
//  ____AngularDistribution____
//  
//  Uses SharcAnalysis and TigressAnalysis libraries to create
//  angular distributions from RedwoodMats histogram file by
//  making projections of an ExcVsTheta matrix and extracting
//  the counts within a specified excitation energy window
//
// ___Input Description___
// A standardized input file is required to create an angular 
// distribution. This describes all energy and angular gates 
// to be used, including background shapes. The input can also
// be specified and included using the CNTSPEC option, which is
// useful for elastic channels when the peaks overlap
//
// A sample input file is shown below :-
/*
INFILE: Results_RedwoodMats.root
OUTFILE: Results_1229.root					
CSFILE: AngDist_1229.txt							// angular distribution output file

FRAME:     	Cm
REACTION:		dp
NORMALIZATION: 7.93e-4 0.07e-4

EXC: 				1229

GAM: 		 		415
GAMPEAK: 		400 425
GAMBGLO:	 	380 400 
GAMBGHI:	 	425 445

EXCBINSZ: 	40.0
GAMBINSZ:   1.0
THETABINSZ: 4.0

//bgspec: theta[mid] bglo1 bglo2 bghi1 bghi2 const slope
BGSPEC:			94   			0			0			0			0		 			// REJECT
BGSPEC:			98   			0			0			0			0		 			// REJECT
BGSPEC:			102 			0			0			0			0		 			// REJECT
BGSPEC:			106 			0			0			0			0		 			// REJECT
*/
//
// The full list of accepted input specifications can be found 
// in the structure below entitled ' SteffenOptions '
//
// An efficiency correction is automatically carried out if available
// as is a branching ratio correction, see TTigressAnalysis docs.
//
// ___Notes___
/*
 Further work : -
 
 Known bugs : -
*/
//
// Made By Steffen Cruz, 2015-2016
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

TList *list, *acclist, *cntlist;
const unsigned long npos = std::string::npos;

struct SteffenOptions {

	std::string infile;
	std::string outfile;	
	std::string csfile;
		
	std::string cntfile;		
	std::string cnthist;		
	
	std::string badstripsfile;
	
	std::string reaction;	
	std::string frame;
	
	Double_t Exc;
	Double_t ExcLo;
	Double_t ExcHi;
	Double_t ExcSig;

	Double_t Gam;
	Double_t GamPeak[2];
	Double_t GamBgLo[2];
	Double_t GamBgHi[2];
	
	Double_t ExcBinSz;
	Double_t GamBinSz;
	Double_t ThetaBinSz;
	Double_t Norm;	
	Double_t NormErr;	
	
	Double_t MaxRelErr;	// discard cross sections with > max rel err

	std::vector<double> cnttheta; // set of theta values which have usr defined counts
			
	std::vector<double> bgtheta; // index of theta in this vector is map key
	std::map<int,std::vector<double> > bgspec;
	
	Bool_t fit;
};

SteffenOptions *info;
void ClearOptions(SteffenOptions *opt){
	opt->infile        = "";
	opt->outfile       = "";
	opt->csfile 			 = "";				
	opt->cntfile			 = "";	
	opt->cnthist 			 = "";								
	opt->badstripsfile = "";	
	opt->frame         = "";
	opt->reaction      = "";
	opt->Exc           = 0.;
	opt->ExcLo         = 0.;
	opt->ExcHi         = 0.;
	opt->Gam           = 0.;
	opt->ExcSig			   = 180.;
	opt->ExcBinSz      = 0.;
	opt->GamBinSz      = 0.;
	opt->ThetaBinSz    = 0.;	
	opt->Norm					 = 1.;
	opt->NormErr 		 	 = 0.;	
	opt->GamPeak[0]    = 0.;
	opt->GamPeak[1]    = 0.;
	opt->GamBgLo[0]    = 0.;
	opt->GamBgLo[1]    = 0.;
	opt->GamBgHi[0]    = 0.;
	opt->GamBgHi[1]    = 0.;
	opt->MaxRelErr     = 0.;
	
	opt->bgtheta.clear();		
	std::vector<double>   emptyvec;
	opt->bgspec[0]     = emptyvec;
	opt->fit						= false;
	opt->cnttheta.clear();	
	return;
}
void PrintOptions(SteffenOptions *opt){

	TCanvas *optcanvas = new TCanvas("Options","Options",580,600);
	TPaveText *pt = new TPaveText(0.1,0.1,0.9,0.9);
	pt->SetTextAlign(11);
	pt->SetTextFont(82);
  pt->SetShadowColor(0);
	
	std::string msg;
	msg.assign("\n\n Input Options :-");
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> infile        = ' %s ' ",opt->infile.c_str()));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());
	msg.assign(Form("\n\t-> outfile       = ' %s ' ",opt->outfile.c_str()));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> csfile        = ' %s ' ",opt->csfile.c_str()));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	
	msg.assign(Form("\n\t-> reaction      = ' %s ' ",opt->reaction.c_str()));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> frame         = ' %s ' ",opt->frame.c_str()));	
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> Exc Energy    = %.2f ",opt->Exc));
	if(opt->ExcLo){
	  pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	  msg.assign(Form("\n\t-> ExcLo Energy  = %.2f ",opt->ExcLo));
	}
	if(opt->ExcHi){
    pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
    msg.assign(Form("\n\t-> ExcHi Energy  = %.2f ",opt->ExcHi));	
	}
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> Exc Sigma     = %.2f ",opt->ExcSig));	
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());		
	msg.assign(Form("\n\t-> Normalization = %.2e ",opt->Norm));
	if(opt->NormErr)
    msg.append(Form("+/- %.2e ",opt->NormErr));
  pt->AddText(msg.c_str()); printf("%s",msg.c_str());		
		
	if(opt->Gam>0){
		msg.assign(Form("\n\t-> Gam Energy    = %.2f ",opt->Gam));
		pt->AddText(msg.c_str()); printf("%s",msg.c_str());			
		msg.assign(Form("\n\t    * Gam Peak Rng  = [ %.1f - %.1f ] *",opt->GamPeak[0],opt->GamPeak[1]));
		pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
		msg.assign(Form("\n\t    * Gam BgLo Rng  = [ %.1f - %.1f ] *",opt->GamBgLo[0],opt->GamBgLo[1]));
		pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
		msg.assign(Form("\n\t    * Gam BgHi Rng  = [ %.1f - %.1f ] *",opt->GamBgHi[0],opt->GamBgHi[1]));
		pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	}
	msg.assign(Form("\n\t-> Exc Bin Sz    = %.1f ",opt->ExcBinSz));	
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> Gam Bin Sz    = %.1f ",opt->GamBinSz));	
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());		
	msg.assign(Form("\n\t-> Theta Bin Sz  = %.1f ",opt->ThetaBinSz));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	if(opt->MaxRelErr){
    msg.assign(Form("\n\t-> Max Rel. Err. = %.2f [ %.2f %% ]",opt->MaxRelErr,opt->MaxRelErr*100));
    pt->AddText(msg.c_str()); printf("%s",msg.c_str());		
  }
  
	if(opt->fit){
		msg.assign(Form("\n\t-> Fit           = %s ",opt->fit?"TRUE":"FALSE"));	
		pt->AddText(msg.c_str()); printf("%s",msg.c_str());				
	}
	
	if(opt->cntfile.length() && opt->cnthist.length()){
		msg.assign(Form("\n\t-> cntfile       = ' %s ' ",opt->cntfile.c_str()));
		pt->AddText(msg.c_str()); printf("%s",msg.c_str());		
		msg.assign(Form("\n\t-> cnthist       = ' %s ' ",opt->cnthist.c_str()));
		pt->AddText(msg.c_str()); printf("%s",msg.c_str());			
		if(opt->cnttheta.size()){
			pt->AddText(" ");
			msg.assign("\n\n\t  * Count Specification For Theta Projections *\n");
			pt->AddText(msg.c_str()); printf("%s",msg.c_str());			
			msg.assign("\n theta     counts +/- counts_err");
			pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	
			for(int i=0; i<(int)opt->cnttheta.size(); i++){
				msg.assign(Form("\n  %5.1f   - will be taken from %s -",opt->cnttheta.at(i),opt->cntfile.c_str()));
				pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
			}			
		}
	}
	
	if(opt->bgtheta.size()){
		pt->AddText(" ");
		msg.assign("\n\n\t  * Background Specification For Theta Projections *\n");
		pt->AddText(msg.c_str()); printf("%s",msg.c_str());			
		msg.assign("\n\t theta    x1lo    x1hi    x2lo   x2hi   const  linear  quad ");
		pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	
		for(int i=0; i<(int)opt->bgtheta.size(); i++){
			msg.assign("\n\t");						
			for(int j=0; j<(int)opt->bgspec[i].size(); j++)
				msg.append(Form("  %5.1f",opt->bgspec[i].at(j)));
				pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
		}		
	}


	
	printf("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - \n\n");
	pt->AddText("");	
	pt->Draw();
	
	list->Add(optcanvas);
}
SteffenOptions *SetOptions(std::string);

Int_t GetIndex(Double_t theta, std::vector<double> v){

	for(int i=0; i<(int)v.size(); i++)
		if(v.at(i)==theta)
			return i;
	
	return -1;		
}
TH1D *MakeScaleHist(Double_t binsz, Double_t scale, Double_t err, const char *name);

TH1D *MakeCountsHist(TH2F *, TList *);

Bool_t ExtractCounts(TH1D *, Double_t, UInt_t, UInt_t, Double_t&, Double_t&);
Double_t PolBg(Double_t *x, Double_t *par);

TH1D *AngularDistribution(std::string);
Bool_t InitVars(std::string);


TH1D *AngularDistribution(std::string OptionsFile){

	if(!InitVars(OptionsFile))
		return 0;	
						
	// acceptance    ////////////////////////////////////////////////////////////
	TH1D *hacc = (TH1D*)acclist->FindObject(Form("SharcCorrection%s",info->frame.c_str()));	
	// make a rebinned clone
	TH1D *hcorr = TSharcAnalysis::RebinAcceptance(hacc,info->ThetaBinSz);
	acclist->Add(hcorr);
	// sine correction ensures solid angles integrate to 4 pi
	TH1D *hsin = TSharcAnalysis::MakeSin(info->ThetaBinSz);
	acclist->Add(hsin);	
		
		
	// hexctheta ////////////////////////////////////////////////////////////
	TH2F *hexctheta = (TH2F*)list->FindObject(Form("ExcTheta%s",info->frame.c_str()));
	// rebin theta & energy, and zoom into energy range
	hexctheta->RebinX(info->ThetaBinSz);
	hexctheta->RebinY(info->ExcBinSz/hexctheta->GetYaxis()->GetBinWidth(0));	
	hexctheta->GetYaxis()->SetRangeUser(info->Exc-1000.,info->Exc+1200.);
		
		
	// counts  ////////////////////////////////////////////////////////////			
	TList *list2 = new TList;  // fill list with projections
  TH1D *hcounts = MakeCountsHist(hexctheta,list2); // get counts versus theta
  hcounts->Scale(1/info->ThetaBinSz); // average the counts instead of rebinning	  
  list->Add(hcounts);	

	if(info->Gam){
	  TCanvas *ceff = TTigressAnalysis::SetEfficiencyCurve();
	  list->Add(ceff);
		Double_t gameff = TTigressAnalysis::Efficiency(info->Gam);		
		Double_t gamerr = TTigressAnalysis::EfficiencyError(info->Gam);			
		printf("\n\n\n\n TIGRESS efficiency for %.1f keV gamma gate = %.2f%% +/- %.2f%% [rel. error = %.2f]\n",info->Gam,gameff*100,gamerr*100,gamerr/gameff);
		TH1D *hgeff = MakeScaleHist(info->ThetaBinSz,1/gameff,gamerr/(gameff*gameff),"TigEfficiency");
	  list->Add(hgeff);
	  hcounts->Scale(1./gameff); // scale without uncertainty and apply after SFRESCO fit.
	//	hcounts->Multiply(hgeff); // we don't want additional uncertainty yet

		Double_t br = TTigressAnalysis::BranchingRatio(info->Exc,info->Gam);		
    if(br>0.01){
      printf("\n Branching ratio for this transition = %.4f%% -> correction = %.4f\n\n",br*100.0,1/br);
      hcounts->Scale(1/br); // scale without uncertainty and apply after SFRESCO fit.
    } else
      printf("\n\t Error :  No branching ratio correction can be applied..\n");
    
	}

	TH1D *hcorrcounts = (TH1D*)hcounts->Clone(Form("CorrectedCountsVsTheta%s",info->frame.c_str()));	
	hcorrcounts->SetTitle(Form("Corrected Counts Vs Theta %s; Theta Cm [deg]; Counts Divided By Fractional Coverage",info->frame.c_str()));
  list->Add(hcorrcounts);
	hcorrcounts->Multiply(hcorr); // save each intermediate step in the process


	// sigma   ////////////////////////////////////////////////////////////
	TH1D *hsigma = (TH1D*)hcounts->Clone(Form("CrossSectionVsTheta%s",info->frame.c_str()));
	hsigma->SetTitle(Form("Differential Cross Section; Theta %s [deg]; Cross section [%s]",info->frame.c_str(),info->Norm>0?"mb/sr":"arb."));	
  list->Add(hsigma);	
  
	hsigma->Multiply(hcorr); // acceptance correction applied to counts
	hsigma->Divide(hsin);  	// sine correction
	TH1D *hnorm = MakeScaleHist(info->ThetaBinSz,info->Norm,info->NormErr,"Normalization");
	list->Add(hnorm);
	hsigma->Scale(info->Norm); // we don't want uncertainty in our result yet
//	hsigma->Multiply(hnorm);
	
  TGraphErrors *gsigma = new TGraphErrors();
	gsigma->SetName(Form("CrossSectionVsTheta%s_Graph",info->frame.c_str()));
	gsigma->SetTitle(Form("Differential Cross Section; Theta %s [deg]; Cross section [%s]",info->frame.c_str(),info->Norm>0?"mb/sr":"arb."));	
  list->Add(gsigma);
 
  std::ofstream cfile(info->csfile.c_str());
 	Double_t theta, sigval, cntval, cnterr, corval, corerr, relerr; 
  for(int i=1; i<hsigma->GetNbinsX(); i++){
  	 
  	 theta = hsigma->GetBinCenter(i);
     sigval = hsigma->GetBinContent(i);
     
     cntval = hcounts->GetBinContent(i);
     cnterr = hcounts->GetBinError(i);
     
     corval = hcorr->GetBinContent(i);
     corerr = hcorr->GetBinError(i);
     relerr = cnterr/cntval+corerr/corval;
          
     if(sigval==0 || relerr>info->MaxRelErr){
        if(sigval){
          printf("\n %i.\t theta = %3.1f \t sigma = %5.3e +/- %5.3e",i,theta,sigval,sigval*relerr);   
          printf("\n ---- relative error for theta = %3.1f is too large, %.3f%% ---\n",theta,relerr*100.0);   
				}
				hsigma->SetBinContent(i,0.0);
     		hsigma->SetBinError(i,0.0);     		
        continue;
		 } 			
     printf("\n %i.\t theta = %3.1f \t sigma = %5.3e +/- %5.3e",i,theta,sigval,sigval*relerr);   
        
     hsigma->SetBinError(i,sigval*relerr);
     gsigma->SetPoint(gsigma->GetN(),hsigma->GetBinCenter(i),sigval);
     gsigma->SetPointError(gsigma->GetN()-1,hsigma->GetBinWidth(i)/2,sigval*relerr);

     if(cfile.is_open()) cfile << theta << "\t" << sigval << "\t" << sigval*relerr << "\n"; 
  }
  
  if(cfile.is_open()) {
		printf("\n\n\t-> Cross Section file ' %s ' has been created.",info->csfile.c_str());
		cfile.close();		
	}
		
	TFile *file = new TFile(info->outfile.c_str(),"RECREATE");

  printf("\n\t-> Output file ' %s ' has been created.\n\n",file->GetName());

	file->mkdir("Results");
	file->cd("Results");	
	list->Write();

	file->cd("");
	file->mkdir("Acceptance");
	file->cd("Acceptance");
	acclist->Write();
	
	file->cd("");
	file->mkdir("Projections");
	file->cd("Projections");
	list2->Write();

	if(info->cntfile.length()){
		file->cd("");
		file->mkdir(info->cntfile.c_str());
		file->cd(info->cntfile.c_str());
		cntlist->Write();
	}
	
  file->Close();
  
  new TCanvas;
  hsigma->Draw();
  
	return hsigma;
}


Bool_t InitVars(std::string OptionsFile){
	
	list = new TList;
	
	info = SetOptions(OptionsFile);
	PrintOptions(info);		

	gStyle->SetOptFit(0112);
						
	TReaction *r;
	if(info->reaction.find("dp")!=npos){
		r = new TReaction("sr95","d","p","sr96",510.9,0,true);
	} else if(info->reaction.find("pp")!=npos){
		r = new TReaction("sr95","p","p","sr95",510.9,0,true); 
	} else if(info->reaction.find("dd")!=npos){
		r = new TReaction("sr95","d","d","sr95",510.9,0,true); 
	} else{
		printf("\nReaction type '%s' not recognised. Fail. \n\n",info->reaction.c_str());
		return 0;
	}	
	printf("Reaction: %s\n\n",r->GetNameFull());
	printf("\n\n Frame ' %s ' Will Be Calculated\n",info->frame.c_str());
	
	// get acceptance correction with dead strips
	TSharcAnalysis::SetTarget(0.0,0.0,0.0,4.5,"cd2",0.5,0.5,0.5);
	const char *stripsfile = "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis/BadStrips.txt";
	
	acclist = TSharcAnalysis::GetAcceptanceList(r,stripsfile,30); // high resolution
  
	// get exc theta mat
	TFile *file = new TFile(info->infile.c_str(),"READ"); // NO BEAM	
	
	// set default excitation energy range if not specified
	if(!info->ExcLo)
	  info->ExcLo = info->Exc-400.0;	  
	if(!info->ExcHi)
	  info->ExcHi = info->Exc+400.0;		
		
	if(!info->MaxRelErr)
	  info->MaxRelErr = 1.0;
	  	
	TH2F *hexctheta;
	if(info->Gam>0.0) {
		TTigressAnalysis::LoadHistos(info->infile.c_str(),info->reaction.c_str());
		if(!TTigressAnalysis::InitLevelsGammas(100,false))
			return false;
		
		TTigressAnalysis::SetGamBinSz(info->GamBinSz);
		TTigressAnalysis::SetExcBinSz(info->ExcBinSz);		
		list->Add(TTigressAnalysis::AnalyzeGammas(info->GamPeak[0],info->GamPeak[1],info->GamBgLo[0],info->GamBgLo[1],info->GamBgHi[0],info->GamBgHi[1],info->ExcLo,info->ExcHi));
		hexctheta = TTigressAnalysis::ExcThetaGated(info->GamPeak[0],info->GamPeak[1],info->GamBgLo[0],info->GamBgLo[1],info->GamBgHi[0],info->GamBgHi[1]); //?		
	} else { 
		printf("\n\n No Gamma Gate will be used.\n");
		hexctheta = (TH2F*)file->FindObjectAny(Form("ExcTheta%s_%s",info->frame.c_str(),info->reaction.c_str()));	
		if(!hexctheta){
			printf("\n\t Error :  Could not find ExcTheta mat!\n");
			return false;
			}
	}	
	hexctheta->SetNameTitle(Form("ExcTheta%s",info->frame.c_str()),Form("Excitation Energy Versus Theta %s",info->frame.c_str()));	
	list->Add(hexctheta);
	
	if(info->cntfile.length()){ // if required locate the counts file 
		TFile *cntfile = new TFile(info->cntfile.c_str(),"READ");
		if(!cntfile){
			printf("\n\tCouldn't locate counts file ' %s '.\n",info->cntfile.c_str());
			return false;
		}
		TH1D *hcntfile = (TH1D*) cntfile->Get(info->cnthist.c_str());
		if(!hcntfile){
			printf("\n\tCouldn't locate counts hist ' %s ' inside file ' %s '.\n",info->cnthist.c_str(),info->cntfile.c_str());
			return false;
		}		
		list->Add(hcntfile);
		printf("\n Counts histogram ' %s ' has been read from file ' %s ' and added to results list.\n",info->cnthist.c_str(),info->cntfile.c_str());
	
		cntlist = new TList;
		TIter next(cntfile->GetListOfKeys());
		TKey *key;
		while ((key = (TKey*)next())) {
			TClass *cl = gROOT->GetClass(key->GetClassName());
			if (cl->InheritsFrom("TH1D")){
				TH1D *obj = (TH1D*)key->ReadObj();
				cntlist->Add(obj);
			}
		}
	}
	
  return true;
}


TH1D *MakeCountsHist(TH2F *h2, TList *list2){
	
	TH1D *htmpx = h2->ProjectionX(), *htmpy = h2->ProjectionY("ExcitationEnergy");
	list2->Add(htmpy);
	Int_t nprows = 4, npcols = 5, npads = nprows*npcols;
	Int_t ncanvas = ceil(180.0/info->ThetaBinSz/((double)npads));
	TCanvas *c[ncanvas];
	for(int i=0; i<ncanvas; i++){
		c[i] = new TCanvas(Form("Projections_Panel%i",i),Form("Theta Projections Canvas: Panel %i",i),1100,700);
		list2->Add(c[i]);
	}
	
	TH1D *hcounts,*h1[180];	
	hcounts = new TH1D(Form("CountsVsTheta%s",info->frame.c_str()),"",h2->GetNbinsX(),0,180);
	hcounts->SetTitle(Form("CountsVsTheta%s exc=%.1f; Theta %s [deg]; Counts",info->frame.c_str(),info->Exc,info->frame.c_str()));

	Double_t theta, tot_counts, bg_counts, counts, err;		
	Double_t binwid = htmpy->GetBinWidth(0);
	TAxis *ax = htmpy->GetXaxis();
  UInt_t peakbin = ax->FindBin(info->Exc);
	UInt_t minbin  = ax->FindBin(info->ExcLo);
	UInt_t maxbin  = ax->FindBin(info->ExcHi);
	UInt_t n = 0;
	printf("\nUsing Exc.Rng. = %.1f-%.1f keV -> Min Bin: %i\t Max Bin: %i\t Peak Bin: %i\n",info->ExcLo,info->ExcHi,minbin,maxbin,peakbin);
			
	printf("\n Projecting theta bins :-");	
	for(int i=htmpx->FindFirstBinAbove(); i<=htmpx->FindLastBinAbove(); i++){

		theta = htmpx->GetBinCenter(i);
		
	// read from cntfile if one is present and if this theta has been set to manual
	
		if(info->cntfile.length() && info->cnthist.length() && GetIndex(theta,info->cnttheta)>=0){ 
			TH1D *hcntfile = (TH1D*) list->FindObject(info->cnthist.c_str());
			UInt_t cntbin = hcntfile->GetXaxis()->FindBin(theta);
			counts = hcntfile->GetBinContent(cntbin);
			err = hcntfile->GetBinError(cntbin);
			
			printf("\n\n %i.\t THETA = %3.1f ..",i,theta);
			printf("\n\t *** USER SPECIFIED COUNTS :   peak counts = %5.1f+/-%4.1f\n",counts,err); 		
			hcounts->SetBinContent(i,counts);
			hcounts->SetBinError(i,err);	
			continue;		
		}
		
    h1[n] = h2->ProjectionY(Form("ExcProj_Theta%s%.1f",info->frame.c_str(),theta),i,i);
		h1[n]->SetTitle(h1[n]->GetName());		
		tot_counts = h1[n]->Integral();
		printf("\n\n %i.\t THETA = %3.1f \ttot_counts = %4.1f",i,theta,tot_counts); 
		if(tot_counts==0)
			continue;
	 // no counts in a projection either be due to zero acceptance or zero statistics		
	 // we should take all non-zero slices as this means that there is detector coverage			
			
		Bool_t success = ExtractCounts(h1[n],theta,minbin,maxbin,counts,err);
		if(success){
			printf("\n\t --->  bg_counts = %4.0f .. peak counts = %5.1f+/-%4.1f\n",tot_counts-counts,counts,err); 		
			hcounts->SetBinContent(i,counts);
			hcounts->SetBinError(i,err);
		} else 
			printf("\n\t***\t Counts were not extracted.. Skipping !! ***\n");

		list2->Add(h1[n]);						
		n++;
	}

	Int_t cnum = -1, pad;
	for(int j=0; j<(int)n; j++){
		
		pad = j%npads;		
		if(pad==0){
			cnum++;
			c[cnum]->Clear();
			c[cnum]->Divide(npcols,nprows);
		}			
		
		c[cnum]->cd(pad+1);
		h1[j]->Draw();

		TF1 *f1 = h1[j]->GetFunction("bg_pol");
		if(!f1) continue;		
		
		TBox *bp = new TBox(f1->GetParameter(3),h1[j]->GetMinimum()*1.05,f1->GetParameter(4),h1[j]->GetMaximum()*1.05);
		bp->SetFillColor(3);
		bp->SetFillStyle(3003); 
		bp->Draw();		
	}
	
	for(int i = cnum+1; i<ncanvas; i++){
		printf("\n\t\t^^^^ Removing canvas %i",i);
		list2->Remove(c[i]);
		c[i]->Close();
	}
			
	return hcounts;
}


Bool_t ExtractCounts(TH1D *h, Double_t theta, UInt_t binlo, UInt_t binhi, Double_t& counts, Double_t& err){
	
	counts = h->Integral(binlo,binhi);
	Double_t xmin = h->GetBinCenter(binlo);
	Double_t xmax = h->GetBinCenter(binhi);	
	Double_t binwid = h->GetBinWidth(0);
		
	TF1 *fbg = new TF1("bg_pol",PolBg,xmin-5*binwid,xmax+5*binwid,5); // fitted background
	fbg->SetParNames("const","slope","quad","x1hi","x2lo");
	fbg->SetParameters(0,0,0,0);
	
	// set background manually
	Double_t thta, x1lo, x1hi, x2lo, x2hi, bg_const=0.0, bg_slope=0.0, bg_quad=0.0;	
	Int_t indx = GetIndex(theta,info->bgtheta);
	if(indx>-1){
		x1lo = info->bgspec[indx].at(1);
		x1hi = info->bgspec[indx].at(2);
		x2lo = info->bgspec[indx].at(3);
		x2hi = info->bgspec[indx].at(4);
		
		if(!x1lo && !x1hi && !x2lo && !x2hi) // user reject point
			return false;
			
		if((x1lo < x1hi) && (x2lo < x2hi)){ // use specified ranges to fit bg
			fbg->SetRange(x1lo,x2hi);
			xmin = x1hi;
			xmax = x2lo;
			fbg->SetLineColor(kMagenta);
			printf("\n\t User Specified Fit Range x1=(%.0f)-(%.0f) & x2=(%.0f)-(%.0f)",x1lo,x1hi,x2lo,x2hi); 
			
			// re-evaluate counts using updated integral range
			counts = h->Integral(h->GetXaxis()->FindBin(xmin),h->GetXaxis()->FindBin(xmax));	
		}
	} 
	
	if(info->Gam>0){ // if a gamma gate has been applied don't subtract background
    fbg->FixParameter(0,0.0);
    fbg->FixParameter(1,0.0);
    fbg->FixParameter(2,0.0);
    fbg->SetLineColor(kOrange);			
  }	
  else if(indx>-1){ // use specified parameters to fix bg
			Int_t polorder = info->bgspec[indx].size()-6;
		if(polorder>=0){
      if(polorder==2){
        bg_quad  = info->bgspec[indx].at(7);	
        if(bg_quad) fbg->FixParameter(2,bg_quad);
      } else fbg->FixParameter(2,0.0);
      if(polorder==1){
        bg_slope = info->bgspec[indx].at(6);	
        if(bg_slope) fbg->FixParameter(1,bg_slope);
      } else fbg->FixParameter(1,0.0);
      fbg->SetLineColor(kOrange);		
      bg_const = info->bgspec[indx].at(5);			
      if(bg_const) fbg->FixParameter(0,bg_const);
      printf("\n\t*** User specified background *** const = %.1f   slope = %.1f   quad = %.1f",bg_const,bg_slope,bg_quad); 											
    } else {
      printf("\n\t Background will be fitted using a pol1 function in specified range");
      fbg->FixParameter(2,0.0);
    }
  } else { // default to linear bg if no spec
    fbg->FixParameter(2,0.0);							
		printf("\n\t Default background will be used for pol1. Using range x1=(%.0f)-(%.0f) & x2=(%.0f)-(%.0f)",xmin-5*binwid,xmin,xmax,xmax+5*binwid); 
	}

	if(info->fit){ // if fit is selected, do a pol1+gaus fit
	
		printf("\n\t  FULL FIT   * Spectrum will be fitted with a pol1 + gaus function *");
		const char *gform = "[0]/(sqrt(2*3.14159)*[2])*exp(-pow(x-[1],2.0)/(2*[2]*[2]))";
		TF1 *peakfit = new TF1("peak_fit",Form("%s+[3]+[4]*x",gform),xmin,xmax);
		peakfit->SetNpx(1000);
		
		peakfit->SetParNames("peak_content","peak_mean","peak_sigma","bg_const","bg_linear");
		peakfit->SetParameters(h->Integral(binlo,binhi)*binwid,info->Exc,info->ExcSig,0,0);
		peakfit->SetParLimits(0,0,h->Integral(binlo,binhi)*binwid); // content
		peakfit->SetParLimits(1,xmin,xmax); // mean position
		peakfit->SetParLimits(2,0.25*info->ExcSig,4.0*info->ExcSig);	// sigma
		
		peakfit->SetLineColor(kRed);
		h->Fit(peakfit,"QR","SAME");
		
		Double_t chi2 = peakfit->GetChisquare()/peakfit->GetNDF();
		printf("\n\t  FULL FIT   * Chi2/NDF = %.3f ",chi2);
		if(chi2 > 10.0){
			printf("!! CHI2 TOO LARGE [>%.3f]",10.0);
			return false;
		} 
		
		TF1 *gaus = new TF1("gaus",gform,info->Exc-1000,info->Exc+1000);
		gaus->SetParameters(peakfit->GetParameter(0),peakfit->GetParameter(1),peakfit->GetParameter(2));
		
		gaus->SetLineColor(kGreen);
		h->GetListOfFunctions()->Add(gaus);
		
		counts = peakfit->GetParameter(0)/binwid;
		err = peakfit->GetParError(0)/binwid;		
		return true; // RETURNS AT THIS POINT
	} 

	fbg->FixParameter(3,xmin);
	fbg->FixParameter(4,xmax);			
	// determine linear background by fitting a pol1 to edges of gaussian and subtracting the area under the line

	TFitResultPtr result = h->Fit(fbg,"RQ","SAME");
	
	Int_t p = result;
	if(p) // make sure that it's still added to the list of functions
		h->GetListOfFunctions()->Add(fbg);


	bg_const = fbg->GetParameter(0);
	bg_slope = fbg->GetParameter(1); 
	bg_quad  = fbg->GetParameter(2); 	

	Double_t bg_counts = fbg->Integral(xmin,xmax)/binwid;
	Double_t bg_err = fbg->IntegralError(xmin,xmax)/binwid;

	if(bg_counts<0 || bg_counts+bg_err>counts){ // this means fail
		printf("\n\t\t Warning :  bg_counts = %.0f+/-%.0f\tcounts = %.0f",bg_counts,bg_err,counts);
		err = 0;
		return false;
	} else if(counts<0)
		counts = 0;
		
	printf("\n\t Fit result :  bg_const = %.3f  bg_slope = %.3f  bg_quad = %.3f\n\t->  bg_counts = %.1f +/- %.1f",bg_const,bg_slope,bg_quad,bg_counts,bg_err);
	Double_t err_lo = sqrt(counts-bg_counts-bg_err);
	Double_t err_hi = sqrt(counts-bg_counts+bg_err); 	
	// assume that peak is zero at xmax and xmin 			
	counts-=bg_counts;
	
	// error = sqrt( fit err ^ 2	+	statistical error ^2 )
	Double_t err2 = pow(bg_err,2.0) + bg_counts + counts;		
	err = sqrt(err2);
	//err = 0.5*(err_lo+err_hi);
		
	TF1 *fbg2 = new TF1("bg_pol2","pol2",xmin,xmax);
	fbg2->SetLineColor(3);		
	fbg2->SetParameters(bg_const,bg_slope,bg_quad);
	fbg2->SetRange(xmin,xmax);		
	h->GetListOfFunctions()->Add(fbg2);
		
	return true;
}


TH1D *MakeScaleHist(Double_t binsz, Double_t scale, Double_t err, const char *name){

  Int_t nbins = (Int_t)180.0/binsz;
  TH1D *hscale = new TH1D(name,name,nbins,0,180);
  hscale->SetTitle(Form("Scale Histogram for %s; Theta %s [deg]; Scale",name,info->frame.c_str()));
  for(int i=1; i<=nbins; i++){
    hscale->SetBinContent(i,scale);
    hscale->SetBinError(i,err);
  }
  return hscale;
}    


Double_t PolBg(Double_t *x, Double_t *par){

	Double_t val =  par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
			
	if(x[0]<par[3] || x[0]>par[4]) // only evaluate outer region for bg estimation
		return val;
	else{
		TF1::RejectPoint();
		return val;
	}	

}


void trim(std::string * line, const std::string trimChars = " \f\n\r\t\v") {
//Removes the the string "trimCars" from  the string 'line'
   if (line->length() == 0)
      return;
   std::size_t found = line->find_first_not_of(trimChars);
   if (found != npos)
      *line = line->substr(found, line->length());
   found = line->find_last_not_of(trimChars);
   if (found != npos)
      *line = line->substr(0, found + 1);
   return;
}


SteffenOptions *SetOptions(std::string OptionsFile) {

  SteffenOptions *info = new SteffenOptions;
  ClearOptions(info);

	if(OptionsFile.length()==0)
		return info;

	std::ifstream optfile(OptionsFile.c_str());
	if(!optfile) {
		printf("could not open file.\n");
		return info;
	}

	std::string line;
	int linenumber = 0;
	bool brace_open = false;

  while (std::getline(optfile, line)) {
		linenumber++;
		trim(&line);
		int comment = line.find("//");
		if ((const unsigned long)comment != npos) 
			 line = line.substr(0, comment);
		
		if (!line.length())
			 continue;

		int colon = line.find(":");
    if((const unsigned long)colon == npos)
      continue;

		std::string type = line.substr(0, colon);
		line = line.substr(colon + 1, line.length());
		trim(&line);
		std::istringstream ss(line);
		int j = 0;
		while (type[j]) {
			char c = *(type.c_str() + j);
			c = toupper(c);
			type[j++] = c;
		}

		if(type.compare("INFILE")==0) {
			 info->infile = line;
		} else if(type.compare("OUTFILE")==0) {
			 info->outfile = line;
		} else if(type.compare("CSFILE")==0) {
			 info->csfile = line; 
		} else if(type.compare("CNTFILE")==0) {
			 info->cntfile = line; 
		} else if(type.compare("CNTHIST")==0) {
			 info->cnthist = line;  			 
		} else if(type.compare("BADSTRIPSFILE")==0) {
			 info->badstripsfile = line;
		} else if(type.compare("FRAME")==0) {
			 info->frame = line;		 
		} else if(type.compare("REACTION")==0) {
			 info->reaction = line;
		} else if(type.compare("EXC")==0) {
			double tempd; ss>>tempd;
			 info->Exc = tempd;
		} else if(type.compare("EXCLO")==0) {
			double tempd; ss>>tempd;
			info->ExcLo = tempd;
    } else if(type.compare("EXCHI")==0) {
			double tempd; ss>>tempd;
			info->ExcHi = tempd;							 
		} else if(type.compare("EXCSIG")==0) {
			double tempd; ss>>tempd;
			 info->ExcSig = tempd;			 
		} else if(type.compare("GAM")==0) {
			 double tempd; ss>>tempd;
			 info->Gam = tempd;
		} else if(type.compare("GAMPEAK")==0) {
			 double value;
			 std::vector<double> tempvec;
			 while (ss >> value) {   tempvec.push_back(value); }
			 if(tempvec.size()==2) {
					info->GamPeak[0] = tempvec.at(0);
					info->GamPeak[1] = tempvec.at(1);
			 }               
		} else if(type.compare("GAMBGLO")==0) {
			 double value;
			 std::vector<double> tempvec;
			 while (ss >> value) {   tempvec.push_back(value); }
			 if(tempvec.size()==2) {
					info->GamBgLo[0] = tempvec.at(0);
					info->GamBgLo[1] = tempvec.at(1);
			 }
		} else if(type.compare("GAMBGHI")==0) {
			 double value;
			 std::vector<double> tempvec;
			 while (ss >> value) {   tempvec.push_back(value); }
			 if(tempvec.size()==2) {
					info->GamBgHi[0] = tempvec.at(0);
					info->GamBgHi[1] = tempvec.at(1);
			 }
    } else if(type.compare("MAXRELERR")==0) {
			double tempd; ss>>tempd;
			info->MaxRelErr = tempd;				 
		} else if(type.compare("EXCBINSZ")==0) {
			 double tempd; ss>>tempd;
			 info->ExcBinSz = tempd;
		} else if(type.compare("GAMBINSZ")==0) {
			 double tempd; ss>>tempd;
			 info->GamBinSz = tempd;			 
		} else if(type.compare("THETABINSZ")==0) {
			 double tempd; ss>>tempd;
			 info->ThetaBinSz = tempd;
		} else if(type.compare("NORMALIZATION")==0) {
			 double value;		
			 std::vector<double> tempvec;
			 while (ss >> value) {   tempvec.push_back(value); }		
			 if(tempvec.size()>0)
         info->Norm = tempvec.at(0);	
			 if(tempvec.size()>1)
         info->NormErr = tempvec.at(1);	         	
		} else if(type.compare("FIT")==0) {
			 if(strcmp(line.c_str(),"TRUE")==0)
			 	info->fit = true; 			 
		} else if(type.compare("BGSPEC")==0) {
			 double value;
			 std::vector<double> tempvec;
			 while (ss >> value) {   tempvec.push_back(value); }
			 if(tempvec.size()>=5) {
					info->bgtheta.push_back(tempvec.at(0)); 
					info->bgspec[info->bgtheta.size()-1] = tempvec;
			 } 	
		}		else if(type.compare("CNTSPEC")==0) {
			 double value;
			 while (ss >> value) {   info->cnttheta.push_back(value); }
		}				
  }
  
  if(info->cnttheta.size() && (!info->cntfile.length() || !info->cnthist.length())){
  	printf("\n\t Error : Counts cannot be manually specified without a CNTFILE to read from.\n");
		ClearOptions(info);  	
	}
	return info;
}
                                                   

int main(int argc, char **argv){

	if(argc<2){
		printf("\n\n\tAn Options File Must Also Be Included!!\n\n");
		return 1;
	}
	
	AngularDistribution(argv[1]);
	
	return 0;
}


