#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;

void CMSSub_AddV2()
{
    TH1::SetDefaultSumw2();

    double BW2D = (9.9/33)*((2-1.0/16)*3.1416/31);
    double errfactor = sqrt(1.0);
    
    TFile* file[12];
    file[11] = TFile::Open("../MB0_10/ppCorr_pt033.root");

    file[0] = TFile::Open("../MB10_20/ppCorr_pt033.root");
    file[1] = TFile::Open("../MB20_30/ppCorr_pt033.root");
    file[2] = TFile::Open("../MB30_40/ppCorr_pt033.root");
    file[3] = TFile::Open("../MB40_60/ppCorr_pt033.root");
    file[4] = TFile::Open("../MB60_85/ppCorr_pt033.root");
    file[5] = TFile::Open("../MB85_95/ppCorr_pt033.root");
    file[6] = TFile::Open("../MB95_110/ppCorr_pt033.root");
    file[7] = TFile::Open("../MB110_above/ppCorr_pt033.root");
    //file[8] = TFile::Open("../HM120_130/ppCorr_pt033.root");
    //file[9] = TFile::Open("../HM130_150/ppCorr_pt033.root");
    //file[10] = TFile::Open("../HM150_above/ppCorr_pt033.root");
    //file[4] = TFile::Open("ppCorr_ref_60_90_dz10_v30.root");
    //file[4] = TFile::Open("ppCorr_ref_60_90_dz10_GplusPP_A.root");
    //file[5] = TFile::Open("ppCorr_ref_90_100_dz10_v30.root");
    //file[6] = TFile::Open("ppCorr_ref_100_110_dz10_v30.root");
    //file[7] = TFile::Open("ppCorr_ref_110_120_dz10_v30.root");
    //file[8] = TFile::Open("ppCorr_ref_120_130_dz10_v30.root");
    //file[9] = TFile::Open("ppCorr_ref_130_150_dz10_v30.root");
    //file[10] = TFile::Open("ppCorr_ref_150_above_dz10_v30.root");
    
    TCanvas* pm = new TCanvas("pm","pm",1200,800);
    TCanvas* c = new TCanvas("c","c",1200,800);
    TCanvas* c1 = new TCanvas("c1","c1",1200,800);
    TCanvas* c2 = new TCanvas("c2","c2",1200,800);
    pm->Divide(4,3);
    c1->Divide(4,3);
    c2->Divide(4,3);
    //high
    
    TH2D* signal[11];
    TH2D* background[11];
    TH1D* mult[11];
    TH1D* ntrk[11];
    TH1D* hPt[11];
    
    signal[0] = (TH2D*)file[0]->Get("pp_MB10_GplusPP/signal");
    background[0] = (TH2D*)file[0]->Get("pp_MB10_GplusPP/background");
    mult[0] = (TH1D*)file[0]->Get("pp_MB10_GplusPP/mult");
    ntrk[0] = (TH1D*)file[0]->Get("pp_MB10_GplusPP/mult_good");
    hPt[0] = (TH1D*)file[0]->Get("pp_MB10_GplusPP/pT");
    signal[1] = (TH2D*)file[1]->Get("pp_MB20_GplusPP/signal");
    background[1] = (TH2D*)file[1]->Get("pp_MB20_GplusPP/background");
    mult[1] = (TH1D*)file[1]->Get("pp_MB20_GplusPP/mult");
    ntrk[1] = (TH1D*)file[1]->Get("pp_MB20_GplusPP/mult_good");
    hPt[1] = (TH1D*)file[1]->Get("pp_MB20_GplusPP/pT");
    signal[2] = (TH2D*)file[2]->Get("pp_MB30_GplusPP/signal");
    background[2] = (TH2D*)file[2]->Get("pp_MB30_GplusPP/background");
    mult[2] = (TH1D*)file[2]->Get("pp_MB30_GplusPP/mult");
    ntrk[2] = (TH1D*)file[2]->Get("pp_MB30_GplusPP/mult_good");
    hPt[2] = (TH1D*)file[2]->Get("pp_MB30_GplusPP/pT");
    signal[3] = (TH2D*)file[3]->Get("pp_MB40_GplusPP/signal");
    background[3] = (TH2D*)file[3]->Get("pp_MB40_GplusPP/background");
    mult[3] = (TH1D*)file[3]->Get("pp_MB40_GplusPP/mult");
    ntrk[3] = (TH1D*)file[3]->Get("pp_MB40_GplusPP/mult_good");
    hPt[3] = (TH1D*)file[3]->Get("pp_MB40_GplusPP/pT");
    signal[4] = (TH2D*)file[4]->Get("pp_MB60_GplusPP/signal");
    background[4] = (TH2D*)file[4]->Get("pp_MB60_GplusPP/background");
    mult[4] = (TH1D*)file[4]->Get("pp_MB60_GplusPP/mult");
    ntrk[4] = (TH1D*)file[4]->Get("pp_MB60_GplusPP/mult_good");
    hPt[4] = (TH1D*)file[4]->Get("pp_MB60_GplusPP/pT");
    signal[5] = (TH2D*)file[5]->Get("pp_MB85_GplusPP/signal");
    background[5] = (TH2D*)file[5]->Get("pp_MB85_GplusPP/background");
    mult[5] = (TH1D*)file[5]->Get("pp_MB85_GplusPP/mult");
    ntrk[5] = (TH1D*)file[5]->Get("pp_MB85_GplusPP/mult_good");
    hPt[5] = (TH1D*)file[5]->Get("pp_MB85_GplusPP/pT");
    signal[6] = (TH2D*)file[6]->Get("pp_MB95_GplusPP/signal");
    background[6] = (TH2D*)file[6]->Get("pp_MB95_GplusPP/background");
    mult[6] = (TH1D*)file[6]->Get("pp_MB95_GplusPP/mult");
    ntrk[6] = (TH1D*)file[6]->Get("pp_MB95_GplusPP/mult_good");
    hPt[6] = (TH1D*)file[6]->Get("pp_MB95_GplusPP/pT");
    signal[7] = (TH2D*)file[7]->Get("pp_MB110_GplusPP/signal");
    background[7] = (TH2D*)file[7]->Get("pp_MB110_GplusPP/background");
    mult[7] = (TH1D*)file[7]->Get("pp_MB110_GplusPP/mult");
    ntrk[7] = (TH1D*)file[7]->Get("pp_MB110_GplusPP/mult_good");
    hPt[7] = (TH1D*)file[7]->Get("pp_MB110_GplusPP/pT");
    /*signal[7] = (TH2D*)file[7]->Get("signal");
    background[7] = (TH2D*)file[7]->Get("background");
    mult[7] = (TH1D*)file[7]->Get("mult");
    ntrk[7] = (TH1D*)file[7]->Get("mult_good");
    hPt[7] = (TH1D*)file[7]->Get("pT");*/
    /*signal[8] = (TH2D*)file[8]->Get("pp_HM120_GplusPP/signal");
    background[8] = (TH2D*)file[8]->Get("pp_HM120_GplusPP/background");
    mult[8] = (TH1D*)file[8]->Get("pp_HM120_GplusPP/mult");
    ntrk[8] = (TH1D*)file[8]->Get("pp_HM120_GplusPP/mult_good");
    hPt[8] = (TH1D*)file[8]->Get("pp_HM120_GplusPP/pT");
    signal[9] = (TH2D*)file[9]->Get("pp_HM130_GplusPP/signal");
    background[9] = (TH2D*)file[9]->Get("pp_HM130_GplusPP/background");
    mult[9] = (TH1D*)file[9]->Get("pp_HM130_GplusPP/mult");
    ntrk[9] = (TH1D*)file[9]->Get("pp_HM130_GplusPP/mult_good");
    hPt[9] = (TH1D*)file[9]->Get("pp_HM130_GplusPP/pT");*/
    /*signal[10] = (TH2D*)file[10]->Get("pp_HM150_GplusPP/signal");
    background[10] = (TH2D*)file[10]->Get("pp_HM150_GplusPP/background");
    mult[10] = (TH1D*)file[10]->Get("pp_HM150_GplusPP/mult");
    ntrk[10] = (TH1D*)file[10]->Get("pp_HM150_GplusPP/mult_good");*/

    
    double Nassoc[11];
    double Nassocfit[11];
    double nEvent[11];
    double Bz[11];
    double srye[11];
    double sry[11];
    double lrye[11];
    double lry[11];
    double subye[11];
    double suby[11];
    double V2[11];
    double V2e[11];
    double V3[11];
    double V3e[11];
    double V1[11];
    double V1e[11];
    double N[11];
   TH1D* lrs[11];
    TH1D* lrs1[11];
    TH1D* lrb[11];
    TH1D* lrb1[11];
    TH1D* srsks[11];
    TH1D* alrs[11];
    //TH1D* alrs1[11];
    double pTmean[11];
    double pTmeanerr[11];
    
    for(unsigned int i=0;i<8;i++)
    {
        mult[i]->GetXaxis()->SetRangeUser(2,300);
        Nassoc[i] = mult[i]->GetMean(1)-0.5;
        N[i] = ntrk[i]->GetMean(1);
        pTmean[i] = hPt[i]->GetMean(1);
        pTmeanerr[i] = hPt[i]->GetMeanError(1);
        
        TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.6,2.2);
	fit->SetLineColor(2);
        fit->SetParameters(1,1,1);
        
        nEvent[i] = mult[i]->Integral(3,10000);
        Bz[i] = background[i]->GetBinContent(background[i]->FindBin(0,0));
  c->cd();      
        //sr Y
        (TH1D*)srsks[i] = signal[i]->ProjectionY(Form("srsks%d",i),14,20);
        TH1D* srbks = background[i]->ProjectionY("srbks",14,20);
        srsks[i]->Divide(srbks);
        srsks[i]->Scale(Bz[i]/nEvent[i]/BW2D);
        
        srsks[i]->Fit("fit","R");
        srsks[i]->Fit("fit","R");
        srsks[i]->Fit("fit","R");
        srsks[i]->Fit("fit","R");
        srsks[i]->Fit("fit","R");
        
        double srmks = fit->GetMinimum(0.6,2.2);
        double srmksx = fit->GetMinimumX(0.6,2.2);
        TF1* mfsrks = new TF1("mfsrks","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mfsrks->SetParameter(0,-srmks);
	TH1D* srs = srsks[i]->Clone();
        srs->Add(mfsrks);
	c1->cd(i+1);
	srsks[i]->Draw(); c->cd();
        sry[i] = srs->IntegralAndError(srs->FindBin(0.0),srs->FindBin(srmksx),srye[i],"width");
        double bin0yield = srs->GetBinContent(srs->FindBin(0.0))*0.19635;
        sry[i] = sry[i]*2 - bin0yield;
        
        //lr Y
        TF1* fit1 = new TF1("fit1","[0]*x^2+[1]*x+[2]",0,2.0);
        fit1->SetLineColor(2);
        fit1->SetParameters(1,1,1);

	(TH1D*)alrs[i] = signal[i]->ProjectionY(Form("alrs_%d",i),1,10);
        alrs1 = signal[i]->ProjectionY("alrs1",24,33);
        alrb = background[i]->ProjectionY("alrb",1,10);
        alrb1 = background[i]->ProjectionY("alrb1",24,33);
        alrs[i]->Add(alrs1);
        alrb->Add(alrb1);
        alrs[i]->Divide(alrb);
        alrs[i]->Scale(Bz[i]/nEvent[i]/BW2D);
        
        alrs[i]->Fit("fit1","R");
        alrs[i]->Fit("fit1","R");
        alrs[i]->Fit("fit1","R");
        alrs[i]->Fit("fit1","R");
        alrs[i]->Fit("fit1","R");

        double lrmks = fit1->GetMinimum(0,2.0);
	double lrmksx = fit1->GetMinimumX(0,2.0);
        TF1* mflrks = new TF1("mflrks","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mflrks->SetParameter(0,-lrmks);
	TH1D* alr = alrs[i]->Clone();
        alr->Add(mflrks); 
	c2->cd(i+1);
	alrs[i]->Draw(); c->cd();
        //lry[i] = alr->IntegralAndError(2,14,lrye[i],"width");
        lry[i] = alr->IntegralAndError(alr->FindBin(0.0),alr->FindBin(lrmksx),lrye[i],"width");
        double bin0yield = alr->GetBinContent(alr->FindBin(0.0))*0.19635;
        lry[i] = lry[i]*2 - bin0yield;
        

        suby[i] = sry[i] - lry[i];
        subye[i] = sqrt(srye[i]*srye[i]+lrye[i]*lrye[i]);
        
        //V2
        //TF1* fit2 = new TF1("fit2","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        TF1* fit2 = new TF1("fit2","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
	fit2->SetLineColor(2);
        fit2->SetParNames("N","V1","V2","V3","V4");
        fit2->SetParameters(10,1,1,1,1);
        
        pm->cd(i+1);

        (TH1D*)lrs[i] = signal[i]->ProjectionY(Form("lrs%d",i),1,10);
        (TH1D*)lrs1[i] = signal[i]->ProjectionY(Form("lrs1%d",i),24,33);
        (TH1D*)lrb[i] = background[i]->ProjectionY(Form("lrb%d",i),1,10);
        (TH1D*)lrb1[i] = background[i]->ProjectionY(Form("lrb1%d",i),24,33);
        lrs[i]->Add(lrs1[i]);
        lrb[i]->Add(lrb1[i]);
        lrs[i]->Divide(lrb[i]);
        lrs[i]->Scale(Bz[i]/nEvent[i]/BW2D);
        lrs[i]->Fit("fit2","R");
        lrs[i]->Fit("fit2","R");
        lrs[i]->Fit("fit2","R");
        lrs[i]->Fit("fit2","R");
        
        Nassocfit[i] = fit2->GetParameter(0);

        //TF1* f1 = new TF1("f1",Form("%f*2.0*%f*cos(2.0*x)",Nassocfit[i],0.0025),-0.5*TMath::Pi(),1.5*TMath::Pi());
        TF1* f1 = new TF1("f1",Form("%f*2.0*%f*cos(2.0*x)",Nassocfit[i],N[i]*3e-05),-0.5*TMath::Pi(),1.5*TMath::Pi());
        lrs[i]->Add(f1);
        
        lrs[i]->Fit("fit2","R");
        lrs[i]->Fit("fit2","R");
        lrs[i]->Fit("fit2","R");
        lrs[i]->Fit("fit2","R");
        
        V2[i] = fit2->GetParameter(2);
        V2e[i] = fit2->GetParError(2)*errfactor;
        V3[i] = fit2->GetParameter(3);
        V3e[i] = fit2->GetParError(3)*errfactor;
        V1[i] = fit2->GetParameter(1);
        V1e[i] = fit2->GetParError(1)*errfactor;
        Nassocfit[i] = fit2->GetParameter(0);
    }
 c->cd();   
    //lowN
    TH2D* signal_low = file[0]->Get("pp_MB10_GplusPP/signal");
    TH2D* background_low = file[0]->Get("pp_MB10_GplusPP/background");
    TH1D* mult_low = file[0]->Get("pp_MB10_GplusPP/mult");
    TH1D* ntrk_low = file[0]->Get("pp_MB10_GplusPP/mult_good");
    N_low = ntrk_low->GetMean(1);

    mult_low->GetXaxis()->SetRangeUser(2,300);
    double Nassoc_low = mult_low->GetMean(1)-0.5;
    
    TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.6,2.2);
    fit->SetParameters(1,1,1);
    
    double nEvent_low = mult_low->Integral(3,10000);
    double Bz_low = background_low->GetBinContent(background_low->FindBin(0,0));
    
    //sr Y
    TH1D* srsks_low = signal_low->ProjectionY("srsks",14,20);
    TH1D* srbks = background_low->ProjectionY("srbks",14,20);
    srsks_low->Divide(srbks);
    srsks_low->Scale(Bz_low/nEvent_low/BW2D);
    
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    
    double srmks = fit->GetMinimum(0.6,2.2);
    double srmksx = fit->GetMinimumX(0.6,2.2);
    TF1* mfsrks = new TF1("mfsrks","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    mfsrks->SetParameter(0,-srmks);
    srsks_low->Add(mfsrks);
    double srye_low;
    double sry_low = srsks_low->IntegralAndError(srsks_low->FindBin(0.0),srsks_low->FindBin(srmksx),srye_low,"width");
    double bin0yield = srsks_low->GetBinContent(srsks_low->FindBin(0.0))*0.19635;
    sry_low = sry_low*2 - bin0yield;
    
    //lr Y
    TF1* fit1 = new TF1("fit1","[0]*x^2+[1]*x+[2]",0,2.0);
        fit1->SetLineColor(2);
        fit1->SetParameters(1,1,1);

    alrs_low = signal_low->ProjectionY("alrs",1,10);
    alrs1 = signal_low->ProjectionY("alrs1",24,33);
    alrb = background_low->ProjectionY("alrb",1,10);
    alrb1 = background_low->ProjectionY("alrb1",24,33);
    alrs_low->Add(alrs1);
    alrb->Add(alrb1);
    alrs_low->Divide(alrb);
    alrs_low->Scale(Bz_low/nEvent_low/BW2D);
    
    alrs_low->Fit("fit1","R");
    alrs_low->Fit("fit1","R");
    alrs_low->Fit("fit1","R");
    alrs_low->Fit("fit1","R");
    alrs_low->Fit("fit1","R");
    
    double lrmks = fit1->GetMinimum(0,2.0);
    double lrmksx = fit1->GetMinimumX(0,2.0);
    TF1* mflrks = new TF1("mflrks","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    mflrks->SetParameter(0,-lrmks);
    alrs_low->Add(mflrks);
    double lrye_low;
    double lry_low = alrs_low->IntegralAndError(2,14,lrye_low,"width");   
    double lry_low = alrs_low->IntegralAndError(alrs_low->FindBin(0.0),alrs_low->FindBin(lrmksx),lrye_low,"width");
    double bin0yield = alrs_low->GetBinContent(alrs_low->FindBin(0.0))*0.19635;
    lry_low = lry_low*2 - bin0yield;


    suby_low = sry_low - lry_low;
    subye_low = sqrt(srye_low*srye_low+lrye_low*lrye_low);

    //TF1* fit2 = new TF1("fit2","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    TF1* fit2 = new TF1("fit2","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit2->SetParNames("N","V1","V2","V3","V4");
    fit2->SetParameters(10,1,1,1,1);
    
    lrs_low = signal_low->ProjectionY("lrs_low",1,10);
    lrs1_low = signal_low->ProjectionY("lrs1_low",24,33);
    lrb_low = background_low->ProjectionY("lrb_low",1,10);
    lrb1_low = background_low->ProjectionY("lrb1_low",24,33);
    lrs_low->Add(lrs1_low);
    lrb_low->Add(lrb1_low);
    lrs_low->Divide(lrb_low);
    lrs_low->Scale(Bz_low/nEvent_low/BW2D);
    lrs_low->Fit("fit2","R");
    lrs_low->Fit("fit2","R");
    lrs_low->Fit("fit2","R");
    lrs_low->Fit("fit2","R");
    
    double V2_low = fit2->GetParameter(2)+N_low*3e-05;
    double V2e_low = fit2->GetParError(2)*errfactor;
    
    double v2_low = sqrt(V2_low);
    double v2e_low = sqrt(V2_low)*(V2e_low/V2_low)/2;
    
    double V3_low = fit2->GetParameter(3);
    double V3e_low = fit2->GetParError(3)*errfactor;

    double V1_low = fit2->GetParameter(1);
    double V1e_low = fit2->GetParError(1)*errfactor;

    double v3_low = 0;
    double v3e_low = 0;
    
    double Nassoc_lowfit = fit2->GetParameter(0);
    
    double V2sub[11];
    double V2sube[11];
    double v2sub[11];
    double v2sube[11];
    double v2[11];
    double v2e[11];

    double V3sub[11];
    double V3sube[11];
    double v3sub[11];
    double v3sube[11];
    double v3[11];
    double v3e[11];
    
    double V1sub[11];
    double V1sube[11];
    
    for(unsigned int i=0;i<8;i++)
    {
        V1sub[i] = V1[i] - V1_low*Nassoc_lowfit/Nassocfit[i]*suby[i]/suby_low;
        V1sube[i] = sqrt(V1e[i]*V1e[i] + V1e_low*Nassoc_lowfit/Nassocfit[i]*suby[i]/suby_low*V1e_low*Nassoc_lowfit/Nassocfit[i]*suby[i]/suby_low);
        
        V2sub[i] = V2[i] - V2_low*Nassoc_lowfit/Nassocfit[i]*suby[i]/suby_low;
        V2sube[i] = sqrt(V2e[i]*V2e[i] + V2e_low*Nassoc_lowfit/Nassocfit[i]*suby[i]/suby_low*V2e_low*Nassoc_lowfit/Nassocfit[i]*suby[i]/suby_low);
        
        v2sub[i] = sqrt(V2sub[i]);
        v2sube[i] = sqrt(V2sub[i])*(V2sube[i]/V2sub[i])/2;
        if(V2sub[i]<0.0000005){
                v2sub[i]=-0.001;
                v2sube[i]=0;
                //continue;
        }

        v2[i] = sqrt(V2[i]);
        v2e[i] = sqrt(V2[i])*(V2e[i]/V2[i])/2;

        V3sub[i] = V3[i] - V3_low*Nassoc_lowfit/Nassocfit[i]*suby[i]/suby_low;
        V3sube[i] = sqrt(V3e[i]*V3e[i] + V3e_low*Nassoc_lowfit/Nassocfit[i]*suby[i]/suby_low*V3e_low*Nassoc_lowfit/Nassocfit[i]*suby[i]/suby_low);
        
        v3sub[i] = sqrt(V3sub[i]);
        v3sube[i] = sqrt(V3sub[i])*(V3sube[i]/V3sub[i])/2;
        
	if(i==0){
	v3sub[i] = -0.001;
	v3sube[i] = 0;
	}

	if(V3[i]<0){
                v3[i]=-0.001;
                v3e[i]=0;
		continue;
        }
   
        v3[i] = sqrt(V3[i]);
        v3e[i] = sqrt(V3[i])*(V3e[i]/V3[i])/2;
    }

    //return;
    
    TFile ofile("v2sub1020_CMS_pt033_addLinearV2_slope3e05.root","RECREATE");
    
    /*signal_low->Divide(background_low);
    signal_low->Scale(Bz_low/nEvent_low/BW2D);
    
    for(unsigned int i=0;i<12;i++)
    {
        signal[i]->Divide(background[i]);
        signal[i]->Scale(Bz[i]/nEvent[i]/BW2D);
        
        signal[i]->Add(signal_low,-(sry[i]/sry_low));
        
        signal[i]->SetName(Form("corr_sub_%d",i));
        
        signal[i]->Write();
    }*/
    
double N2[11];
for(int j=0;j<9;j++)
{
	N2[j] = N[j+1];
}
double v2new[12];
double v2newe[12];
double V2new[12];
double V2newe[12];

double v3new[12];
double v3newe[12];
double V3new[12];
double V3newe[12];

v2new[0] = v2_low;
v2newe[0] = v2e_low;

V2new[0] = V2_low;
V2newe[0] = V2e_low;
v3new[0] = v3_low;
v3newe[0] = v3e_low;
V3new[0] = V3_low;
V3newe[0] = V3e_low;


for(int j=0;j<9;j++)
{
v2new[j+1] = v2[j];
v2newe[j+1] = v2e[j];
V2new[j+1] = V2[j];
V2newe[j+1] = V2e[j];
v3new[j+1] = v3[j];
v3newe[j+1] = v3e[j];
V3new[j+1] = V3[j];
V3newe[j+1] = V3e[j];

}
    
    TGraphErrors* pTplot = new TGraphErrors(8,N,pTmean,0,pTmeanerr);
    
    pTplot->SetName("meanPt");
    
    pTplot->Write();

    TGraphErrors* hadronv2 = new TGraphErrors(8,N,v2,0,v2e);
    TGraphErrors* hadronv2sub = new TGraphErrors(8,N,v2sub,0,v2sube);
    
    hadronv2->SetName("hadronv2_GplusPP");
    hadronv2sub->SetName("hadronv2sub_GplusPP");
    
    hadronv2->Write();
    hadronv2sub->Write();
    
    TGraphErrors* hadronv3 = new TGraphErrors(8,N,v3,0,v3e);
    TGraphErrors* hadronv3sub = new TGraphErrors(8,N,v3sub,0,v3sube);
    
    hadronv3->SetName("hadronv3_GplusPP");
    hadronv3sub->SetName("hadronv3sub_GplusPP");
    
    hadronv3->Write();
    hadronv3sub->Write();
    
    TGraphErrors* hadronV2 = new TGraphErrors(8,N,V2,0,V2e);
    TGraphErrors* hadronV2sub = new TGraphErrors(8,N,V2sub,0,V2sube);

    hadronV2->SetName("hadronV2_GplusPP");
    hadronV2sub->SetName("hadronV2sub_GplusPP");

    TGraphErrors* hadronV1 = new TGraphErrors(8,N,V1,0,V1e);
    TGraphErrors* hadronV1sub = new TGraphErrors(8,N,V1sub,0,V1sube);
    
    hadronV1->SetName("hadronV1_GplusPP");
    hadronV1sub->SetName("hadronV1sub_GplusPP");
    
    hadronV1->Write();
    hadronV1sub->Write();
    
    hadronV2->Write();
    hadronV2sub->Write();

    TGraphErrors* hadronV3 = new TGraphErrors(8,N,V3,0,V3e);
    TGraphErrors* hadronV3sub = new TGraphErrors(8,N,V3sub,0,V3sube);

    hadronV3->SetName("hadronV3_GplusPP");
    hadronV3sub->SetName("hadronV3sub_GplusPP");

    hadronV3->Write();
    hadronV3sub->Write();

TGraphErrors* lryplot = new TGraphErrors(8,N,lry,0,lrye);

lryplot->SetName("lry");
lryplot->Write();

 double NassocV2[11];
    double NassocV2e[11];
    double Nall[11];
    double jety[11];
    double jetye[11];

    NassocV2[0] = Nassoc_low*V2_low;
    NassocV2e[0] = Nassoc_low*V2e_low;
    Nall[0] = N_low;
    //Nall[0] = 10;
    jety[0] = suby_low;
    jetye[0] = subye_low;
    for(int j=0;j<8;j++)
        {
        Nall[j+1] = N[j];
        NassocV2[j+1] = Nassoc[j]*V2[j];
        NassocV2e[j+1] = Nassoc[j]*V2e[j];
        jety[j+1] = suby[j];
        jetye[j+1] = subye[j];
	//if(j==0) jety[j+1] = jety[j+1]+0.017;
	//if(j==1) jety[j+1] = jety[j+1]+0.003;
            //if(j==0) suby[j] = suby[j]+0.017;
            //if(j==1) suby[j] = suby[j]+0.003;
        }

        TGraphErrors* V2yield = new TGraphErrors(8,Nall,NassocV2,0,NassocV2e);
        TGraphErrors* jetyield = new TGraphErrors(8,N,suby,0,subye);

        V2yield->SetName("V2yield");
        jetyield->SetName("jetyield");

        V2yield->Write();
        jetyield->Write();

//    cout<<"v2low:"<<V2_low<<",error:"<<v2e_low<<endl;
//	cout<<Nassoc_low<<endl;
//cout<<suby_low<<endl;
    for(unsigned int i=0;i<8;i++)
    {
        cout<<"Nassoc:"<<Nassocfit[i]<<endl;
        cout<<"V2:"<<V2[i]<<",error:"<<V2e[i]<<endl;
	//cout<<Nassoc[i]<<endl;
	cout<<"suby: "<<suby[i]<<", lry: "<<lry[i]<<", sry: "<<sry[i]<<endl;
        cout<<"V1sub:"<<V1sub[i]<<",error:"<<V1sube[i]<<endl;
        cout<<"V2sub:"<<V2sub[i]<<",error:"<<V2sube[i]<<endl;
        cout<<"V3sub:"<<V3sub[i]<<",error:"<<V3sube[i]<<endl;
        cout<<"========="<<endl;
    }
    cout<<""<<endl;
    //cout<<"v3low:"<<v3_low<<",error:"<<v3e_low<<endl;
    for(unsigned int i=0;i<8;i++)
    {
        cout<<"N:"<<N[i]<<endl;
        cout<<"V3:"<<V3[i]<<",error:"<<V3e[i]<<endl;
        cout<<"V3sub:"<<V3sub[i]<<",error:"<<V3sube[i]<<endl;
    }

    for(unsigned int i=0;i<8;i++)
    {
	cout<<"suby:"<<suby[i]<<endl;
	cout<<"perisub factor:"<<NassocV2[i]/jety[i]<<endl;
    }
cout<<"subylow:"<<suby_low<<endl;
//TCanvas* tt = new TCanvas("tt","tt",900,600);
//tt->cd();
//alrs[1]->Draw();
}
