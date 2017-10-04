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

#include "GRootCommands.h"

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
	std::string tigeffdata;
	std::string nndcdata;
	
  Int_t A         ;
  Double_t BeamE  ;
	Double_t TargetX;
	Double_t TargetY;
	Double_t TargetZ;	
	Double_t TargetThick;	
	
	std::string reaction;	
	std::string frame;
	
	Double_t Exc;
	Double_t ExcLo;
	Double_t ExcHi;
	Double_t ExcViewRng;
	Double_t ExcSig;

	Double_t Gam;
	Double_t GamPeak[2];
	Double_t GamBgLo[2];
	Double_t GamBgHi[2];
	Double_t TigAbsEff[3];
	
	Double_t ExcBinSz;
	Double_t GamBinSz;
	Double_t ThetaBinSz;
	Double_t Norm;	
	Double_t NormErr;	
	
	Double_t MaxRelErr;	// discard cross sections with > max rel err

	std::vector<double> cnttheta; // set of theta values which have usr defined counts
			
	std::vector<double> bgtheta; // index of theta in this vector is map key
	std::map<int,std::vector<double> > bgspec;
	
	std::vector<double> FitPeakMean;
	double FitPeakSigma;
	
	Bool_t fit;
	Bool_t expbg;
};
void ClearOptions(SteffenOptions *opt){
  opt->infile        = "";
	opt->outfile       = "";
	opt->csfile 			 = "";				
	opt->cntfile			 = "";	
	opt->cnthist 			 = "";								
	opt->badstripsfile = "";	
	opt->tigeffdata    = "";	
	opt->nndcdata      = "";	

  opt->A             = 0;
  opt->BeamE         = 0; 
	opt->TargetX       = 0.0;
	opt->TargetY       = 0.0;
	opt->TargetZ       = 0.0;
	opt->TargetThick   = 0.0;	
	opt->frame         = "";
	opt->reaction      = "";			
	opt->Exc           = 0.;
	opt->ExcLo         = 0.;
	opt->ExcHi         = 0.;
	opt->ExcViewRng    = 0.;	
	opt->Gam           = 0.;
	opt->ExcSig			   = 0.; // used to be 180
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
	opt->TigAbsEff[0]  = -1.;
	opt->TigAbsEff[1]  = -1.;
	opt->TigAbsEff[2]  = -1.;		
	opt->MaxRelErr     = 0.;
	
	opt->bgtheta.clear();		
	std::vector<double>   emptyvec;
	opt->bgspec[0]      = emptyvec;
	opt->fit						= false;
	opt->expbg					= false;
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
	msg.assign(Form("\n\t-> badstripsfile = ' %s ' ",opt->badstripsfile.c_str()));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());		

	msg.assign(Form("\n\t-> A             =  %i  ",opt->A));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> BeamE         =  %.1f ",opt->BeamE));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> TargetX       =  %.1f ",opt->TargetX));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> TargetY       =  %.1f ",opt->TargetY));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> TargetZ       =  %.1f ",opt->TargetZ));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());					
	msg.assign(Form("\n\t-> TargetThick   =  %.2f ",opt->TargetThick));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> reaction      = ' %s ' ",opt->reaction.c_str()));
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());			
	msg.assign(Form("\n\t-> frame         = ' %s ' ",opt->frame.c_str()));	
	pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	msg.assign(Form("\n\t-> Exc Energy    = %.2f ",opt->Exc));
	if(opt->ExcLo){
	  pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	  msg.assign(Form("\n\t-> ExcLo       = %.2f ",opt->ExcLo));
	}
	if(opt->ExcHi){
    pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
    msg.assign(Form("\n\t-> ExcHi       = %.2f ",opt->ExcHi));	
	}
	if(opt->ExcViewRng){
	  pt->AddText(msg.c_str()); printf("%s",msg.c_str());	
	  msg.assign(Form("\n\t-> ExcViewRng  = %.2f ",opt->ExcViewRng));
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
	  msg.assign(Form("\n\t-> Tig Eff. Data = ' %s ' ",opt->tigeffdata.c_str()));		
	  pt->AddText(msg.c_str()); printf("%s",msg.c_str());			
	  if(opt->TigAbsEff[0]){
      msg.assign(Form("\n\t-> Tig Abs. Eff. = %.3f+/-%.3f %% @ %.3f keV ",
        opt->TigAbsEff[1],opt->TigAbsEff[2],opt->TigAbsEff[0]));		
      pt->AddText(msg.c_str()); printf("%s",msg.c_str());			  
	  }
	  msg.assign(Form("\n\t-> NNDC Data     = ' %s ' ",opt->nndcdata.c_str()));		
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
		if(opt->FitPeakMean.size()==3){
      msg.assign(Form("\n\t-> FitPeakMean   = %.1f\t%.1f\t%.1f ",opt->FitPeakMean[0],opt->FitPeakMean[1],opt->FitPeakMean[2]));	
      pt->AddText(msg.c_str()); printf("%s",msg.c_str());				
      msg.assign(Form("\n\t-> FitPeakSigma  = %.1f ",opt->FitPeakSigma));	
      pt->AddText(msg.c_str()); printf("%s",msg.c_str());					
		}
	}
	if(opt->expbg){
		msg.assign(Form("\n\t-> Exp. Bg.        = %s ",opt->expbg?"TRUE":"FALSE"));	
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
		msg.assign("\n\t theta    x1lo    x1hi    x2lo   x2hi   const  slope  quad ");
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

SteffenOptions *info;

struct TriplePeakResults {

  double content[3];
  double height[3];  
  double mean[3];
  double sigma[3]; 
  double cnterr[3];
  //double bg[3]; //? 
  double sum;
  double sumerr;
  
  void Print(UInt_t indx){
    printf("\n\n TRIPLE PEAK RESULT : ");
    printf("\n\tcontent[%i] = %.1f +/- %.1f",indx,content[indx],cnterr[indx]);
    printf("\n\tmean[%i] = %.1f",indx,mean[indx]);
    printf("\n\tsigma[%i] = %.1f\n\n",indx,sigma[indx]);
  }
};
TriplePeakResults GetTripleFitResults(TH1 *hist) {
  TriplePeakResults tpr;
  TF1 *f = (TF1*)hist->GetListOfFunctions()->FindObject("triple_gaus_fit");
  if(!f) return tpr;
  double xlow,xhigh;
  f->GetRange(xlow,xhigh);
  
  TAxis *xax = hist->GetXaxis();
  tpr.sum = hist->Integral(xax->FindBin(xlow),xax->FindBin(xhigh));
  tpr.sumerr = 1/sqrt(tpr.sum);
  
  TF1 g("tempgaus","gaus",xlow,xhigh);
  g.SetParameters(f->GetParameter(0),f->GetParameter(3),f->GetParameter(6));
  tpr.content[0] = g.Integral(xlow,xhigh)/hist->GetBinWidth(0);  
  tpr.height[0]  = g.GetParameter(0);
  tpr.mean[0]  = g.GetParameter(1);
  tpr.sigma[0]  = g.GetParameter(2);
  tpr.cnterr[0] = sqrt(hist->Integral(xax->FindBin(tpr.mean[0]-3*tpr.sigma[0]),xax->FindBin(tpr.mean[0]+3*tpr.sigma[0])));  
  //tpr.cnterr[0] =  tpr.content[0]*sqrt(pow(tpr.cnterr[0]/tpr.content[0],2)+)
  
  g.SetParameters(f->GetParameter(1),f->GetParameter(4),f->GetParameter(6));
  tpr.content[1] = g.Integral(xlow,xhigh)/hist->GetBinWidth(0);
  tpr.height[1]  = g.GetParameter(0);
  tpr.mean[1]  = g.GetParameter(1);
  tpr.sigma[1]  = g.GetParameter(2);
  tpr.cnterr[1] = sqrt(hist->Integral(xax->FindBin(tpr.mean[1]-3*tpr.sigma[1]),xax->FindBin(tpr.mean[1]+3*tpr.sigma[1])));  
 
  g.SetParameters(f->GetParameter(2),f->GetParameter(5),f->GetParameter(6));
  tpr.content[2] = g.Integral(xlow,xhigh)/hist->GetBinWidth(0);
  tpr.height[2]  = g.GetParameter(0);
  tpr.mean[2]  = g.GetParameter(1);
  tpr.sigma[2]  = g.GetParameter(2);
  tpr.cnterr[2] = sqrt(hist->Integral(xax->FindBin(tpr.mean[2]-3*tpr.sigma[2]),xax->FindBin(tpr.mean[2]+3*tpr.sigma[2]))); 
 
  return tpr;
}

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
Double_t ExpBg(Double_t *x, Double_t *par);

Bool_t ScanEnergyWindow(Double_t exc, Int_t ex_min, Int_t ex_max, Int_t ex_step, Double_t ex_lo=0);
TH1D *AngularDistribution(std::string);
Bool_t InitVars(std::string);


Bool_t ScanEnergyWindow(Double_t exc, Int_t ex_min, Int_t ex_max, Int_t ex_step, Double_t ex_lo){

  Double_t ex_hi;
  for(int i=ex_min; i<=ex_max; i+=ex_step){
  
    std::string template_name = Form("Input_%.0f.txt",exc);
    std::ifstream file_tmp(template_name);
    if(!file_tmp.is_open()){
      printf("\n\t Error :  Couldn't find input file ' %s '\n\n",template_name.c_str());
      return false;
    }  
  
    ex_hi = (double)i;
    
    std::string input_name = Form("Input_%.0f_ExcHi.txt",exc);
    std::ofstream options_file(input_name);
  
    std::string line;
    while (std::getline(file_tmp, line)) {
    
      if(line.find("OUTFILE:")!=npos) // set root file
        line.assign(Form("OUTFILE: Results_%.0f_ExcHi%.0f.root",exc,ex_hi));
      if(line.find("CSFILE:")!=npos)  // set angdist file
        line.assign(Form("CSFILE: AngDist_%.0f_ExcHi%.0f.txt",exc,ex_hi));		  
      if(line.find("EXCHI:")!=npos)   // set ExcHi
        line.assign(Form("EXCHI: %.0f",ex_hi));			  	
      if(line.find("EXCLO:")!=npos && ex_lo)
        line.assign(Form("EXCLO: %.0f",ex_lo));			  	
      
      options_file << line << "\n";
    }

    options_file.close();
    AngularDistribution(input_name.c_str()); // extract angular distribution
  }

  return true;
}


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
	hexctheta->GetYaxis()->SetRangeUser(info->Exc-info->ExcViewRng,info->Exc+info->ExcViewRng);
		
		
	// counts  ////////////////////////////////////////////////////////////			
	TList *list2 = new TList;  // fill list with projections
  TH1D *hcountsraw = MakeCountsHist(hexctheta,list2); // get counts versus theta
  list->Add(hcountsraw);
  
  TH1D *hcounts = (TH1D*)hcountsraw->Clone("EffBrCorrectedCountsVsThetaCm");
  // very important step! If I group together five degrees i'll hav five times more cross
  // section, however this should not produce an angular distribution with 5xSF
  hcounts->Scale(1/info->ThetaBinSz); // average the counts instead of rebinning	  
  list->Add(hcounts);	

	if(info->Gam){
	  
	  TCanvas *ceff = TTigressAnalysis::SetEfficiencyCurve(info->tigeffdata.c_str(),info->TigAbsEff[0],info->TigAbsEff[1],info->TigAbsEff[2]);
	  list->Add(ceff);

    Double_t gameff, gamerr; 
		gameff = TTigressAnalysis::Efficiency(info->Gam);		
	  gamerr = TTigressAnalysis::EfficiencyError(info->Gam);			

		printf("\n\n\n\n TIGRESS efficiency for %.1f keV gamma gate = %.2f%% +/- %.2f%% [rel. error = %.2f]\n",info->Gam,gameff*100,gamerr*100,gamerr/gameff);
		TH1D *hgeff = MakeScaleHist(info->ThetaBinSz,1/gameff,gamerr/(gameff*gameff),"TigEfficiency");
	  list->Add(hgeff);
	  
	  hcounts->Scale(1./gameff); // scale without uncertainty and apply after SFRESCO fit.
	//	hcounts->Multiply(hgeff); // we don't want additional uncertainty yet

    // Load levels and gammas
    TTigressAnalysis::SetVerbose(false);
    TTigressAnalysis::LoadLevelsGammas(info->nndcdata.c_str(),100);
    
		Double_t br = TTigressAnalysis::BranchingRatio(info->Exc,info->Gam);		
    if(br>0.01){
      printf("\n Branching ratio for this transition = %.4f%% -> correction = %.4f\n\n",br*100.0,1/br);
      hcounts->Scale(1/br); // scale without uncertainty and apply after SFRESCO fit.
    } else
      printf("\n\t Error :  No branching ratio correction can be applied..\n");
    
	}

	TH1D *hcorrcounts = (TH1D*)hcounts->Clone(Form("CorrectedCountsVsTheta%s",info->frame.c_str()));	
	hcorrcounts->SetTitle(Form("Corrected Counts Vs #theta; #theta_{%s} [#circ]; Counts Divided By Fractional Coverage",info->frame.c_str()));
  list->Add(hcorrcounts);
	hcorrcounts->Multiply(hcorr); // save each intermediate step in the process
  
  Int_t a2 = info->A;
  if(info->reaction.find("dp")!=npos)
    a2+=1;
  else if(info->reaction.find("dp")!=npos)
    a2-=1;    
  std::string title = Form("^{%i}Sr(%c,%c)^{%i}Sr, E_{exc}=%.0f keV ",info->A,info->reaction[0],
    info->reaction[1],a2,info->Exc);
  if(info->Gam)
    title+=Form("#gamma=%.0f keV",info->Gam);
    
	// sigma   ////////////////////////////////////////////////////////////
	TH1D *hsigma = (TH1D*)hcounts->Clone(Form("CrossSectionVsTheta%s",info->frame.c_str()));
	hsigma->SetTitle(Form("%s; #theta_{%s} [#circ]; #frac{d#sigma}{d#Omega} [%s]",
	    title.c_str(),info->frame.c_str(),info->Norm>0?"mb/sr":"arb."));
  list->Add(hsigma);	
  
	hsigma->Multiply(hcorr); // acceptance correction applied to counts
	hsigma->Divide(hsin);  	// sine correction
	TH1D *hnorm = MakeScaleHist(info->ThetaBinSz,info->Norm,info->NormErr,"Normalization");
	list->Add(hnorm);
	hsigma->Scale(info->Norm); // we don't want uncertainty in our result yet
//	hsigma->Multiply(hnorm);
	
  TGraphErrors *gsigma = new TGraphErrors();
	gsigma->SetName(Form("CrossSectionVsTheta%s_Graph",info->frame.c_str()));
	gsigma->SetTitle(Form("%s; #theta_{%s} [#circ]; #frac{d#sigma}{d#Omega} [%s]",
	    title.c_str(),info->frame.c_str(),info->Norm>0?"mb/sr":"arb."));
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

	//gStyle->SetOptFit(0112);
						
	TReaction *r;
	// excitation energy is included 
	if(info->reaction.find("dp")!=npos){
		r = new TReaction(Form("sr%i",info->A),"d","p",Form("sr%i",info->A+1),info->BeamE,info->Exc*1e-3,true);
	} else if(info->reaction.find("pp")!=npos){
		r = new TReaction(Form("sr%i",info->A),"p","p",Form("sr%i",info->A),info->BeamE,info->Exc*1e-3,true); 
	} else if(info->reaction.find("dd")!=npos){
		r = new TReaction(Form("sr%i",info->A),"d","d",Form("sr%i",info->A),info->BeamE,info->Exc*1e-3,true); 
	} else if(info->reaction.find("dt")!=npos){
		r = new TReaction(Form("sr%i",info->A),"d","t",Form("sr%i",info->A-1),info->BeamE,info->Exc*1e-3,true);
	} else{
		printf("\nReaction type '%s' not recognised. Fail. \n\n",info->reaction.c_str());
		return 0;
	}	
	printf("Reaction: %s\n\n",r->GetNameFull());
	printf("\n\n Frame ' %s ' Will Be Calculated\n",info->frame.c_str());
	
	// get acceptance correction with dead strips
	TSharcAnalysis::SetTarget(info->TargetX,info->TargetY,info->TargetZ,info->TargetThick,"cd2",0.5,0.5,0.5);
	//const char *stripsfile = "/Users/steffencruz/Desktop/Steffen/Work/PhD/TRIUMF/CodesAndTools/SharcAnalysis/BadStrips.txt";
	
	acclist = TSharcAnalysis::GetAcceptanceList(r,info->badstripsfile.c_str(),30); // high resolution
  
	// get exc theta mat
	TFile *file = new TFile(info->infile.c_str(),"READ"); // NO BEAM	
	
	// set default excitation energy range if not specified
	if(!info->ExcLo)
	  info->ExcLo = info->Exc-400.0;	  
	if(!info->ExcHi)
	  info->ExcHi = info->Exc+400.0;		
		
	if(!info->ExcViewRng)
	  info->ExcViewRng = 2000.0;	 
	  		
	if(!info->MaxRelErr)
	  info->MaxRelErr = 1.0;
	  	
	TH2F *hexctheta;
	if(info->Gam>0.0) {
		TTigressAnalysis::LoadHistos(info->infile.c_str(),info->reaction.c_str());

//		if(!TTigressAnalysis::InitLevelsGammas(100,false))
//			return false;
		
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
	hcounts->SetTitle(Form("CountsVsTheta%s exc=%.1f; #theta_{%s} [#circ]; Counts",info->frame.c_str(),info->Exc,info->frame.c_str()));

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
			
			printf("\n\n %i.\t THETA = %3.1f [bin=%u]..",i,theta,cntbin);
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

		TF1 *f1;
    if(info->expbg)
  		 f1 = h1[j]->GetFunction("bg_exp");
  	else 	 
  		 f1 = h1[j]->GetFunction("bg_pol");
  	
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
		
	TF1 *fbg;
	if(info->expbg){
	  fbg = new TF1("bg_exp",ExpBg,xmin-5*binwid,xmax+5*binwid,5); // fitted background
	  fbg->SetParNames("bg_const","exp_slope","exp_mean","x1hi","x2lo");
    fbg->SetParameters(0,0,0,0);
  }else {
	  fbg = new TF1("bg_pol",PolBg,xmin-5*binwid,xmax+5*binwid,5); // fitted background
	  fbg->SetParNames("bg_const","bg_slope","bg_quad","x1hi","x2lo");  
    fbg->SetParameters(0,0,0,0);
	}
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
      if(info->expbg){
        printf("\n\t Background will be fitted using an exp function in specified range");
      } else {
        printf("\n\t Background will be fitted using a pol1 function in specified range");
        fbg->FixParameter(2,0.0);      
      }
    }
  } else { // default to linear bg if no spec
    fbg->FixParameter(2,0.0);							
		printf("\n\t Default background will be used for pol1. Using range x1=(%.0f)-(%.0f) & x2=(%.0f)-(%.0f)",xmin-5*binwid,xmin,xmax,xmax+5*binwid); 
	}

	if(info->fit){ // if fit is selected, do a pol1+gaus fit
	
    if(info->FitPeakMean.size()==3){
      std::vector<double> mean = info->FitPeakMean;
      
      UInt_t indx=0;
      // Grab appropriate peak from triple fit
      for(int i=0; i<(int)mean.size(); i++)
        if(mean.at(i)==info->Exc)
          indx = i;
     /*
     xmax = 2500.0;
     xmin = -2000.0;
     Double_t expgr = 0.0002;
      if(theta<=22.0){ // uqqq specific settings
        h->GetXaxis()->SetRangeUser(-1000,1200);
        xmin = -2000;        
        xmax = 1200.0;
        if(theta==22)
          
        expgr = 0.002;
      } else if(theta==26)
        expgr = 0.001;
       */ 
     
     TF1 *peakfit;     
     if(!info->Gam){
       Double_t expgr = 0.0002;     
       // constrain exponential background slope to be the same in each section
        if(theta<=22.0){ // uqqq specific settings
          h->GetXaxis()->SetRangeUser(-1000,1200);
          expgr = 0.002;
        } else if(theta==26)
          expgr = 0.001;    
          
        peakfit = DoTripleGausFit(h,xmin,xmax,mean[0],mean[1],mean[2],info->FitPeakSigma,expgr);         
      } else {
        h->GetXaxis()->SetRangeUser(info->ExcLo,info->ExcHi);
        for(int i=1; i<=h->GetNbinsX(); i++)
          h->SetBinError(i,sqrt(h->GetBinContent(i)));     
        // last argument is the paximum centroid shift from specified position
        peakfit = DoTripleGausFitNoBg(h,xmin,xmax,mean[0],mean[1],mean[2],info->FitPeakSigma,100.0);          
      }
      
      TriplePeakResults res = GetTripleFitResults(h);
      res.Print(indx);
      // choose output counts from 0, 352 or 680 peak.
      counts = res.content[indx];
      err = res.cnterr[indx];	      
      return true;
    }

		printf("\n\t  FULL FIT   * Spectrum will be fitted with a pol1 + gaus function *");
		const char *gform = "[0]/(sqrt(2*3.14159)*[2])*exp(-pow(x-[1],2.0)/(2*[2]*[2]))";
		TF1 *peakfit = new TF1("peak_fit",Form("%s+[3]+[4]*x",gform),info->Exc-1000,info->Exc+2000);//xmin,xmax);
		peakfit->SetLineWidth(1);
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
		
		TF1 *gaus = new TF1("gaus",gform,info->Exc-2*info->ExcSig,info->Exc+2*info->ExcSig);
		gaus->SetParameters(peakfit->GetParameter(0),peakfit->GetParameter(1),peakfit->GetParameter(2));

		gaus->SetLineWidth(1);		
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

	Double_t bg_counts = fbg->Integral(xmin,xmax)/binwid;
	Double_t bg_err = fbg->IntegralError(xmin,xmax)/binwid;

	if(bg_counts<0 || bg_counts+bg_err>counts){ // this means fail
		printf("\n\t\t Warning :  bg_counts = %.0f+/-%.0f\tcounts = %.0f",bg_counts,bg_err,counts);
		err = 0;
		return false;
	} else if(counts<0)
		counts = 0;

	Double_t err_lo = sqrt(counts-bg_counts-bg_err);
	Double_t err_hi = sqrt(counts-bg_counts+bg_err); 	
	// assume that peak is zero at xmax and xmin 			
	counts-=bg_counts;
	
	// error = sqrt( fit err ^ 2	+	statistical error ^2 )
	Double_t err2 = pow(bg_err,2.0) + bg_counts + counts;		
	err = sqrt(err2);
	//err = 0.5*(err_lo+err_hi);
		
	bg_const = fbg->GetParameter(0);
	bg_slope = fbg->GetParameter(1); 
	bg_quad  = fbg->GetParameter(2); 		
	printf("\n\t Fit result :  %s = %.3e  %s = %.3e  %s = %.3e\n\t->  bg_counts = %.1f +/- %.1f",
	  fbg->GetParName(0),bg_const,fbg->GetParName(1),bg_slope,fbg->GetParName(2),bg_quad,bg_counts,bg_err);
			
	TF1 *fbg2;
	if(info->expbg)
	  fbg2 = new TF1("bg_exp2","[0]+TMath::Exp([1]*(x-[2]))",xmin,xmax);
	else 
	  fbg2 = new TF1("bg_pol2","pol2",xmin,xmax);

  fbg2->SetParameters(bg_const,bg_slope,bg_quad);	  
	fbg2->SetLineColor(3);			
	fbg2->SetRange(xmin,xmax);		
	h->GetListOfFunctions()->Add(fbg2);
		
	return true;
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

Double_t ExpBg(Double_t *x, Double_t *par){

	Double_t val =  par[0] + TMath::Exp(par[1]*(x[0]-par[2]));
			
	if(x[0]<par[3] || x[0]>par[4]) // only evaluate outer region for bg estimation
		return val;
	else{
		TF1::RejectPoint();
		return val;
	}	

}

TH1D *MakeScaleHist(Double_t binsz, Double_t scale, Double_t err, const char *name){

  Int_t nbins = (Int_t)180.0/binsz;
  TH1D *hscale = new TH1D(name,name,nbins,0,180);
  hscale->SetTitle(Form("Scale Histogram for %s; #theta_{%s} [#circ]; Scale",name,info->frame.c_str()));
  for(int i=1; i<=nbins; i++){
    hscale->SetBinContent(i,scale);
    hscale->SetBinError(i,err);
  }
  return hscale;
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
		} else if(type.compare("TIGEFFDATA")==0) {
			 info->tigeffdata = line;			 
		} else if(type.compare("NNDCDATA")==0) {
			 info->nndcdata = line;			 			 
		} else if(type.compare("A")==0) {
			int tempd; ss>>tempd;
			 info->A = tempd;			 
		} else if(type.compare("BEAME")==0) {
			double tempd; ss>>tempd;
			 info->BeamE = tempd;
		} else if(type.compare("TARGPOS")==0) {
			 double value;
			 std::vector<double> tempvec;
			 while (ss >> value) {   tempvec.push_back(value); }
			 if(tempvec.size()==3) {
					info->TargetX = tempvec.at(0);
					info->TargetY = tempvec.at(1);
					info->TargetZ = tempvec.at(2);					
			 } 
		} else if(type.compare("TARGTHICK")==0) {
			double tempd; ss>>tempd;
			 info->TargetThick = tempd;			  	 
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
		} else if(type.compare("EXCVIEWRNG")==0) {
			double tempd; ss>>tempd;
			info->ExcViewRng = tempd;								 
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
	  } else if(type.compare("FITPEAKMEAN")==0) { // 3 NUMBERS
			 double value;
			 std::vector<double> tempvec;
			 while (ss >> value) {   tempvec.push_back(value); }
			 if(tempvec.size()==3) 
					info->FitPeakMean = tempvec;	  		 	
	  } else if(type.compare("FITPEAKSIGMA")==0) { // 1 NUMBER
			 double value; ss >> value;
			 info->FitPeakSigma = value; 
		} else if(type.compare("BGSPEC")==0) {
			 double value;
			 std::vector<double> tempvec;
			 while (ss >> value) {   tempvec.push_back(value); }
			 if(tempvec.size()>=5) {
					info->bgtheta.push_back(tempvec.at(0)); 
					info->bgspec[info->bgtheta.size()-1] = tempvec;
			 } 	
		}		else if(type.compare("COUNTSPEC")==0) {
			 double value;
			 while (ss >> value) {   info->cnttheta.push_back(value); }				
		} else if(type.compare("EXPBG")==0) {
			 if(strcmp(line.c_str(),"TRUE")==0)
			 	info->expbg = true; 	
    } else if(type.compare("TIGABSEFF")==0) {
			 double value;
			 std::vector<double> tempvec;
			 while (ss >> value) {   tempvec.push_back(value); }
			 if(tempvec.size()==3) {
					info->TigAbsEff[0] = tempvec.at(0);
					info->TigAbsEff[1] = tempvec.at(1);
					info->TigAbsEff[2] = tempvec.at(2);					
			 } 
    }
  }
  
  if(info->cnttheta.size() && (!info->cntfile.length() || !info->cnthist.length())){
  	printf("\n\t Error : Counts cannot be manually specified without a CNTFILE to read from.\n");
		ClearOptions(info);  	
	}
	
	if(info->FitPeakMean.size()==3 && info->FitPeakSigma==0){
  	printf("\n\t Error : Number of fit peak means must be 3 AND SIGMA MUST BE SPECIFIED.\n");
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



