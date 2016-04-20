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

void perisubrefalltest()
{
    TH1::SetDefaultSumw2();

    TFile *_file0 = TFile::Open("../MB10_20/ppCorr_pt033.root");
    TFile *_file1 = TFile::Open("ppCorr_pt033.root");

    TFile *_file2 = TFile::Open("../MB10_20/ppCorr_all_pt033.root");
    TFile *_file3 = TFile::Open("ppCorr_all_pt033.root");
    
    TLatex* tex = new TLatex();
    tex->SetNDC();
    
    char Nplot[13][200] = {"0.2","0.4","0.6","0.8","1.0","1.4","1.8","2.2","2.8","3.6","4.6","6.0","9.0"};
 
    TCanvas* c1 = new TCanvas("lowNsr","lowNsr",1200,900);
    TCanvas* c2 = new TCanvas("lowNlr","lowNlr",1200,900);
    TCanvas* c3 = new TCanvas("highNsr","highNsr",1200,900);
    TCanvas* c4 = new TCanvas("highNlr","highNlr",1200,900);
    TCanvas* c5 = new TCanvas("lowFit","lowFit",1200,900);
    TCanvas* c6 = new TCanvas("highFit","highFit",1200,900);

    c1->Divide(4,3);
    c2->Divide(4,3);
    c3->Divide(4,3);
    c4->Divide(4,3);
    c5->Divide(4,3);
    c6->Divide(4,3);

    TCanvas* c = new TCanvas("c","c",400,400);
    
    c->cd();
   
    double BW2D = (9.9/33)*((2-1.0/16)*3.1416/31);
    double errfactor = sqrt(1.0);
    //V2ref high Y
    
    TH2D* signal_ref = _file1.Get("pp_MB10_GplusPP/signal");
    TH2D* background_ref = _file1.Get("pp_MB10_GplusPP/background");
    TH1D* mult_ref = _file1.Get("pp_MB10_GplusPP/mult");
    
    TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.6,2.2);
    fit->SetParameters(1,1,1);
    
    double nEvent_ref = mult_ref->Integral(3,10000);
    double Bz_ref = background_ref->GetBinContent(background_ref->FindBin(0,0));
    
    TH1D* srsks = signal_ref->ProjectionY("srsks",14,20);
    TH1D* srbks = background_ref->ProjectionY("srbks",14,20);
    srsks->Divide(srbks);
    srsks->Scale(Bz_ref/nEvent_ref/BW2D);
    
    srsks->Fit("fit","R");
    srsks->Fit("fit","R");
    srsks->Fit("fit","R");
    srsks->Fit("fit","R");
    srsks->Fit("fit","R");
    
    double srmks = fit->GetMinimum(0.6,2.2);
    double srmksx = fit->GetMinimumX(0.6,2.2);
    TF1* mfsrks = new TF1("mfsrks","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    mfsrks->SetParameter(0,-srmks);
    srsks->Add(mfsrks);
    double srye_ref;
    double sry_ref = srsks->IntegralAndError(srsks->FindBin(0.0),srsks->FindBin(srmksx),srye_ref,"width");
    double bin0yield = srsks->GetBinContent(srsks->FindBin(0.0))*0.19635;
    sry_ref = sry_ref*2 - bin0yield;
    
TF1* fit10 = new TF1("fit10","[0]*x^2+[1]*x+[2]",0,2.0);
    fit10->SetParameters(1,1,1);

    TH1D* alrs = signal_ref->ProjectionY("alrsks22",1,10);
    TH1D* alrs1 = signal_ref->ProjectionY("alrs1",24,33);
    TH1D* alrb = background_ref->ProjectionY("alrb",1,10);
    TH1D* alrb1 = background_ref->ProjectionY("alrb1",24,33);
    alrs->Add(alrs1);
    alrb->Add(alrb1);
    alrs->Divide(alrb);
    alrs->Scale(Bz_ref/nEvent_ref/BW2D);
    
    alrs->Fit("fit10","R");
    alrs->Fit("fit10","R");
    alrs->Fit("fit10","R");
    alrs->Fit("fit10","R");
    alrs->Fit("fit10","R");
    
    double lrm = fit10->GetMinimum(0,2.0);
    double lrmksx = fit10->GetMinimumX(0,2.0);
    TF1* mflr = new TF1("mflr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    mflr->SetParameter(0,-lrm);
    alrs->Add(mflr);
    double lrye_ref;
    double lry_ref = alrs->IntegralAndError(alrs->FindBin(0.0),alrs->FindBin(lrmksx),lrye_ref,"width");
    double bin0yield = alrs->GetBinContent(alrs->FindBin(0.0))*0.19635;
    lry_ref = lry_ref*2 - bin0yield;
    
    double suby_ref = sry_ref - lry_ref;
    double subye_ref = sqrt(srye_ref*srye_ref+lrye_ref*lrye_ref);
    
    TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetParNames("N","V1","V2","V3","V4");
    fit1->SetParameters(10,1,1,1,1);
    
    TH1D* lrs = signal_ref->ProjectionY("lrs",1,10);
    TH1D* lrs1 = signal_ref->ProjectionY("lrs1",24,33);
    TH1D* lrb = background_ref->ProjectionY("lrb",1,10);
    TH1D* lrb1 = background_ref->ProjectionY("lrb1",24,33);
    lrs->Add(lrs1);
    lrb->Add(lrb1);
    lrs->Divide(lrb);
lrs->Scale(Bz_ref/nEvent_ref/BW2D);
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    
    double V2_ref = fit1->GetParameter(2);
    double V2e_ref = fit1->GetParError(2)*errfactor;

    double V3_ref = fit1->GetParameter(3);
    double V3e_ref = fit1->GetParError(3)*errfactor;

    double V1_ref = fit1->GetParameter(1);
    double V1e_ref = fit1->GetParError(1)*errfactor;
    
    double v2_ref = sqrt(V2_ref);
    double v2e_ref = sqrt(V2_ref)*(V2e_ref/V2_ref)/2;
    double v3_ref = sqrt(V3_ref);
    double v3e_ref = sqrt(V3_ref)*(V3e_ref/V3_ref)/2;
    double v1_ref = sqrt(V1_ref);
    double v1e_ref = sqrt(V1_ref)*(V1e_ref/V1_ref)/2;
    
double Nassoc_ref_fit = fit1->GetParameter(0);   
    //V2ref low Y
    
    TH2D* signal_ref_low = _file0.Get("pp_MB10_GplusPP/signal");
    TH2D* background_ref_low = _file0.Get("pp_MB10_GplusPP/background");
    TH1D* mult_ref_low = _file0.Get("pp_MB10_GplusPP/mult");
    
    TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.6,2.2);
    fit->SetParameters(1,1,1);
    
    double nEvent_ref_low = mult_ref_low->Integral(3,10000);
    double Bz_ref_low = background_ref_low->GetBinContent(background_ref_low->FindBin(0,0));
    
    TH1D* srsks_low = signal_ref_low->ProjectionY("srsks_low",14,20);
    TH1D* srbks_low = background_ref_low->ProjectionY("srbks_low",14,20);
    srsks_low->Divide(srbks_low);
    srsks_low->Scale(Bz_ref_low/nEvent_ref_low/BW2D);
    
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    
    TCanvas* c22 = new TCanvas;
    c22->cd();
    srsks_low->Draw();
    c->cd();
    
    double srmks_low = fit->GetMinimum(0.6,2.2);
    double srmksx_low = fit->GetMinimumX(0.6,2.2);
    TF1* mfsrks_low = new TF1("mfsrks_low","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    mfsrks_low->SetParameter(0,-srmks_low);
    srsks_low->Add(mfsrks_low);
    double srye_ref_low;
    double sry_ref_low = srsks_low->IntegralAndError(srsks_low->FindBin(0.0),srsks_low->FindBin(srmksx_low),srye_ref_low,"width");
    double bin0yield = srsks_low->GetBinContent(srsks_low->FindBin(0.0))*0.19635;
    sry_ref_low = sry_ref_low*2 - bin0yield;
    
    TF1* fit10 = new TF1("fit10","[0]*x^2+[1]*x+[2]",0,2.0);
    fit10->SetParameters(1,1,1);
    
    TH1D* alrsks = signal_ref_low->ProjectionY("alrsks",1,10);
    TH1D* alrs1 = signal_ref_low->ProjectionY("alrs1",24,33);
    TH1D* alrb = background_ref_low->ProjectionY("alrb",1,10);
    TH1D* alrb1 = background_ref_low->ProjectionY("alrb1",24,33);
    alrsks->Add(alrs1);
    alrb->Add(alrb1);
    alrsks->Divide(alrb);
    alrsks->Scale(Bz_ref_low/nEvent_ref_low/BW2D);
    
    alrsks->Fit("fit10","R");
    alrsks->Fit("fit10","R");
    alrsks->Fit("fit10","R");
    alrsks->Fit("fit10","R");
    alrsks->Fit("fit10","R");
    
    double lrm = fit10->GetMinimum(0,2.0);
	double lrmksx_low = fit10->GetMinimumX(0,2.0);
    TF1* mflr = new TF1("mflr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    mflr->SetParameter(0,-lrm);
    alrsks->Add(mflr);
    double lrye_ref_low;
    double lry_ref_low = alrsks->IntegralAndError(alrsks->FindBin(0.0),alrsks->FindBin(lrmksx_low),lrye_ref_low,"width");
    double bin0yield = alrsks->GetBinContent(alrsks->FindBin(0.0))*0.19635;
    lry_ref_low = lry_ref_low*2 -bin0yield;
    
    double suby_ref_low = sry_ref_low - lry_ref_low;
    //double suby_ref_low = sry_ref_low;
    double subye_ref_low = sqrt(srye_ref_low*srye_ref_low+lrye_ref_low*lrye_ref_low);
    
    TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetParNames("N","V1","V2","V3","V4");
    fit1->SetParameters(10,1,1,1,1);
    
    TH1D* lrs = signal_ref_low->ProjectionY("lrs",1,10);
    TH1D* lrs1 = signal_ref_low->ProjectionY("lrs1",24,33);
    TH1D* lrb = background_ref_low->ProjectionY("lrb",1,10);
    TH1D* lrb1 = background_ref_low->ProjectionY("lrb1",24,33);
    lrs->Add(lrs1);
    lrb->Add(lrb1);
    lrs->Divide(lrb);
    lrs->Scale(Bz_ref_low/nEvent_ref_low/BW2D);
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    
    double V1_ref_low = fit1->GetParameter(1);
    double V1e_ref_low = fit1->GetParError(1)*errfactor;
    
    double v1_ref_low = sqrt(V1_ref_low);
    double v1e_ref_low = sqrt(V1_ref_low)*(V1e_ref_low/V1_ref_low)/2;
    
    double V2_ref_low = fit1->GetParameter(2);
    double V2e_ref_low = fit1->GetParError(2)*errfactor;
    
    double v2_ref_low = sqrt(V2_ref_low);
    double v2e_ref_low = sqrt(V2_ref_low)*(V2e_ref_low/V2_ref_low)/2;

    double V3_ref_low = fit1->GetParameter(3);
    double V3e_ref_low = fit1->GetParError(3)*errfactor;
    
    double v3_ref_low = sqrt(V3_ref_low);
    double v3e_ref_low = sqrt(V3_ref_low)*(V3e_ref_low/V3_ref_low)/2;

    double Nassoc_ref_fit_low = fit1->GetParameter(0);
    //LowMult Yield
    double Nassoc_ref_low;
    TH1D* mult_assoc_low_ref;
    _file0->GetObject("pp_MB10_GplusPP/mult_assoc",mult_assoc_low_ref);
    //mult_assoc_low_ref->GetXaxis()->SetRangeUser(2,300);
    Nassoc_ref_low = mult_assoc_low_ref->GetMean(1);
    //Nassoc_low = 11.68;
    
    int nEvent_low[12];
    double Bz_low[12];
    double Nassoc_low[12];
double Nassoc_fit_low[12];
    
    double srye_low[12];
    double sry_low[12];
    double lrye_low[12];
    double lry_low[12];
    double subye_low[12];
    double suby_low[12];
    
    double V2_low[12];
    double V2e_low[12];
    
    double V3_low[12];
    double V3e_low[12];

    double V1_low[12];
    double V1e_low[12];
    
    TH2D* signal_low[12];
    TH2D* background_low[12];
    
    TH1D* mult_low[12];
    TH1D* mult_assoc_low[12];
    
    TH1D* srs_low[12];
    TH1D* alrs_low[12];

    TH1D* lrb_low[12];
    TH1D* lrs_low[12];
    
    TH1D* proj_low;
    TH1D* proj_high;
    
    for(int i=0;i<12;i++){
        //sr Y
        _file2->GetObject(Form("pp_MB10_GplusPP/signal%d",i),signal_low[i]);
        _file2->GetObject(Form("pp_MB10_GplusPP/background%d",i),background_low[i]);
        _file2->GetObject(Form("pp_MB10_GplusPP/mult_good%d",i),mult_low[i]);
        _file2->GetObject(Form("pp_MB10_GplusPP/mult_assoc%d",i),mult_assoc_low[i]);
        Nassoc_low[i] = mult_assoc_low[i]->GetMean(1);

        TF1* fit2 = new TF1("fit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
        fit2->SetParameters(1,1,1);
        fit2->SetLineColor(2);
        TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.0,2.0);
        fit->SetParameters(1,1,1);
        fit->SetLineColor(2);

        nEvent_low[i] = mult_low[i]->Integral(2,10000);
        Bz_low[i] = background_low[i]->GetBinContent(background_low[i]->FindBin(0,0));
        c->cd();
        (TH1D*)srs_low[i] = signal_low[i]->ProjectionY(Form("srslow%d",i),14,20);
        TH1D* srb = background_low[i]->ProjectionY("srb",14,20);
        srs_low[i]->Divide(srb);
        srs_low[i]->Scale(Bz_low[i]/nEvent_low[i]/BW2D);
        
        srs_low[i]->Fit("fit2","R");
        srs_low[i]->Fit("fit2","R");
        srs_low[i]->Fit("fit2","R");
        srs_low[i]->Fit("fit2","R");
        srs_low[i]->Fit("fit2","R");
        
        double srm = fit2->GetMinimum(0.6,2.2);
	double srmx = fit2->GetMinimumX(0.6,2.2);
        TF1* mfsr = new TF1("mfsr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mfsr->SetParameter(0,-srm);
        TH1D* srslow = srs_low[i]->Clone();
        srslow->Add(mfsr);
        c1->cd(i+1);
        srs_low[i]->Draw();
        tex->DrawLatex(0.52,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.52,0.82,"10<N_{trk}^{offline}<20");
        tex->DrawLatex(0.52,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.52,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.52,0.59,"|#Delta#eta|<1");
        c->cd();
        sry_low[i] = srslow->IntegralAndError(srslow->FindBin(0.0),srslow->FindBin(srmx),srye_low[i],"width");
	double bin0yield = srslow->GetBinContent(srslow->FindBin(0.0))*0.19635;
	sry_low[i] = sry_low[i]*2 - bin0yield;
        //cout<<"111"<<endl;
        //lr Y
        (TH1D*)alrs_low[i] = signal_low[i]->ProjectionY(Form("alrslow%d",i),1,10);
        alrs1 = signal_low[i]->ProjectionY("alrs1",24,33);
        alrb = background_low[i]->ProjectionY("alrb",1,10);
        alrb1 = background_low[i]->ProjectionY("alrb1",24,33);
        alrs_low[i]->Add(alrs1);
        alrb->Add(alrb1);
        alrs_low[i]->Divide(alrb);
        alrs_low[i]->Scale(Bz_low[i]/nEvent_low[i]/BW2D);
        
        alrs_low[i]->Fit("fit","R");
        alrs_low[i]->Fit("fit","R");
        alrs_low[i]->Fit("fit","R");
        alrs_low[i]->Fit("fit","R");
        alrs_low[i]->Fit("fit","R");
        
        double lrm = fit->GetMinimum(0.0,2.0);
        double lrmx = fit->GetMinimumX(0.0,2.0);
        TF1* mflr = new TF1("mflr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mflr->SetParameter(0,-lrm);
        TH1D* alrslow = alrs_low[i]->Clone();
        alrslow->Add(mflr);
        if(i==11) proj_low = alrslow->Clone();
        c2->cd(i+1);
        alrs_low[i]->Draw();
        tex->DrawLatex(0.22,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.22,0.82,"10<N_{trk}^{offline}<20");
        tex->DrawLatex(0.22,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.22,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.22,0.59,"|#Delta#eta|>2");
        c->cd();
        lry_low[i] = alrslow->IntegralAndError(alrslow->FindBin(0.0),alrslow->FindBin(lrmx),lrye_low[i],"width");
	double bin0yield = alrslow->GetBinContent(alrslow->FindBin(0.0))*0.19635;
	lry_low[i] = lry_low[i]*2 - bin0yield;
        
        //sub Y
        suby_low[i] = sry_low[i] - lry_low[i];
        //suby_low[i] = sry_low[i];
        subye_low[i] = sqrt(srye_low[i]*srye_low[i]+lrye_low[i]*lrye_low[i]);
        c5->cd(i+1);

        TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3","V4");
        fit1->SetParameters(10,1,1,1,1);
        fit1->SetLineColor(2);
        
        //lr V2 obs
        lrs_low[i] = signal_low[i]->ProjectionY(Form("lrs%d",i),1,10);
        TH1D* lrsa = signal_low[i]->ProjectionY("lrsa",24,33);
        lrb_low[i] = background_low[i]->ProjectionY(Form("lrb%d",i),1,10);
        TH1D* lrba = background_low[i]->ProjectionY("lrba",24,33);
        lrs_low[i]->Add(lrsa);
        lrb_low[i]->Add(lrba);
        lrs_low[i]->Divide(lrb_low[i]);
        lrs_low[i]->Scale(Bz_low[i]/nEvent_low[i]/BW2D);
        //lrs_low[i]->Scale(1.0/4.2);
        lrs_low[i]->Fit("fit1","R");
        lrs_low[i]->Fit("fit1","R");
        lrs_low[i]->Fit("fit1","R");
        lrs_low[i]->Fit("fit1","R");
        
        V2_low[i] = fit1->GetParameter(2);
        V2e_low[i] = fit1->GetParError(2)*errfactor;
        V3_low[i] = fit1->GetParameter(3);
        V3e_low[i] = fit1->GetParError(3)*errfactor;
        V1_low[i] = fit1->GetParameter(1);
        V1e_low[i] = fit1->GetParError(1)*errfactor;

        Nassoc_fit_low[i] = fit1->GetParameter(0);
	c->cd();
    }
    
    //HighMult Yield
    double Nassoc_ref;
    TH1D* mult_assoc_ref;
    _file1->GetObject("pp_MB10_GplusPP/mult_assoc",mult_assoc_ref);
    //mult_assoc_ref->GetXaxis()->SetRangeUser(2,300);
    Nassoc_ref = mult_assoc_ref->GetMean(1);
    
    int nEvent[12];
    double Bz[12];
    double Nassoc[12];
double Nassoc_fit[12];
    
    double srye[12];
    double sry[12];
    double lrye[12];
    double lry[12];
    double subye[12];
    double suby[12];
    
    double V2[12];
    double V2e[12];

    double V3[12];
    double V3e[12];

    double V1[12];
    double V1e[12];

    TH2D* signal[12];
    TH2D* background[12];
    
    TH1D* mult[12];
    TH1D* mult_assoc[12];
    
    TH1D* srs_high[12];
    TH1D* alrs_high[12];
    
    TH1D* lrb_high[12];
    TH1D* lrs_high[12];

    for(int i=0;i<12;i++){
        //sr Y
        _file3->GetObject(Form("pp_MB10_GplusPP/signal%d",i),signal[i]);
        _file3->GetObject(Form("pp_MB10_GplusPP/background%d",i),background[i]);
        _file3->GetObject(Form("pp_MB10_GplusPP/mult_good%d",i),mult[i]);
        _file3->GetObject(Form("pp_MB10_GplusPP/mult_assoc%d",i),mult_assoc[i]);
	Nassoc[i] = mult_assoc[i]->GetMean(1);	

        TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0,2.0);
        fit->SetParameters(1,1,1);
        fit->SetLineColor(2);
        TF1* fit2 = new TF1("fit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
        fit2->SetParameters(1,1,1);
        fit2->SetLineColor(2);

        nEvent[i] = mult[i]->Integral(2,10000);
        Bz[i] = background[i]->GetBinContent(background[i]->FindBin(0,0));
        
        (TH1D*)srs_high[i] = signal[i]->ProjectionY(Form("srshigh%d",i),14,20);
        TH1D* srb = background[i]->ProjectionY("srb",14,20);
        srs_high[i]->Divide(srb);
        srs_high[i]->Scale(Bz[i]/nEvent[i]/BW2D);
        
        srs_high[i]->Fit("fit2","R");
        srs_high[i]->Fit("fit2","R");
        srs_high[i]->Fit("fit2","R");
        srs_high[i]->Fit("fit2","R");
        srs_high[i]->Fit("fit2","R");
        
        double srm = fit2->GetMinimum(0.6,2.2);
	double srmx = fit2->GetMinimumX(0.6,2.2);
        TF1* mfsr = new TF1("mfsr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mfsr->SetParameter(0,-srm);
        TH1D* srshigh = srs_high[i]->Clone();
        srshigh->Add(mfsr);
        c3->cd(i+1);
        srs_high[i]->Draw();
        tex->DrawLatex(0.52,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.52,0.82,"110<N_{trk}^{offline}<150");
        tex->DrawLatex(0.52,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.52,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.52,0.59,"|#Delta#eta|<1");
        c->cd();
        sry[i] = srshigh->IntegralAndError(srshigh->FindBin(0.0),srshigh->FindBin(srmx),srye[i],"width");
	double bin0yield = srshigh->GetBinContent(srshigh->FindBin(0.0))*0.19635;
	sry[i] = sry[i]*2 - bin0yield;
        
        //lr Y
        (TH1D*)alrs_high[i] = signal[i]->ProjectionY(Form("alrshigh%d",i),1,10);
        TH1D* alrs1 = signal[i]->ProjectionY("alrs1",24,33);
        TH1D* alrb = background[i]->ProjectionY("alrb",1,10);
        TH1D* alrb1 = background[i]->ProjectionY("alrb1",24,33);
        alrs_high[i]->Add(alrs1);
        alrb->Add(alrb1);
        alrs_high[i]->Divide(alrb);
        alrs_high[i]->Scale(Bz[i]/nEvent[i]/BW2D);
        
        alrs_high[i]->Fit("fit","R");
        alrs_high[i]->Fit("fit","R");
        alrs_high[i]->Fit("fit","R");
        alrs_high[i]->Fit("fit","R");
        alrs_high[i]->Fit("fit","R");
        
        double lrm = fit->GetMinimum(0,2.0);
	double lrmx = fit->GetMinimumX(0,2.0);
        TF1* mflr = new TF1("mflr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mflr->SetParameter(0,-lrm);
        TH1D* alrshigh = alrs_high[i]->Clone();
        alrshigh->Add(mflr);
        if(i==11) proj_high = alrshigh->Clone();
        c4->cd(i+1);
        alrs_high[i]->Draw();
        tex->DrawLatex(0.22,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.22,0.82,"110<N_{trk}^{offline}<150");
        tex->DrawLatex(0.22,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.22,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.22,0.59,"|#Delta#eta|>2");
        c->cd();
        lry[i] = alrshigh->IntegralAndError(alrshigh->FindBin(0.0),alrshigh->FindBin(lrmx),lrye[i],"width");
	double bin0yield = alrshigh->GetBinContent(alrshigh->FindBin(0.0))*0.19635;
	lry[i] = lry[i]*2 - bin0yield;
        
        //sub Y
        suby[i] = sry[i] - lry[i];
        subye[i] = sqrt(srye[i]*srye[i]+lrye[i]*lrye[i]);

        c6->cd(i+1);
        
        TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3","V4");
        fit1->SetParameters(10,1,1,1,1);
        fit1->SetLineColor(2);
        
        //lr V2 obs
        lrs_high[i] = signal[i]->ProjectionY(Form("lalrs%d",i),1,10);
        TH1D* lrsa = signal[i]->ProjectionY("lrsa",24,33);
        lrb_high[i] = background[i]->ProjectionY(Form("lalrb%d",i),1,10);
        TH1D* lrba = background[i]->ProjectionY("lrba",24,33);
	//lrs_high[i]->Scale(1.0/2.0);
	//lrb_high[i]->Scale(1.0/2.0);
	lrs_high[i]->Add(lrsa);
        lrb_high[i]->Add(lrba);
        lrs_high[i]->Divide(lrb_high[i]);
        lrs_high[i]->Scale(Bz[i]/nEvent[i]/BW2D);
        //lrs_high[i]->Scale(1.0/4.2);
        //lrs_high[i]->Scale(Bz[i]/nEvent[i]/BW2D);
        lrs_high[i]->Fit("fit1","R");
        lrs_high[i]->Fit("fit1","R");
        lrs_high[i]->Fit("fit1","R");
        lrs_high[i]->Fit("fit1","R");
        
        V2[i] = fit1->GetParameter(2);
        V2e[i] = fit1->GetParError(2)*errfactor;
        V3[i] = fit1->GetParameter(3);
        V3e[i] = fit1->GetParError(3)*errfactor;
        V1[i] = fit1->GetParameter(1);
        V1e[i] = fit1->GetParError(1)*errfactor;
        Nassoc_fit[i] = fit1->GetParameter(0);
	c->cd();
    }
    
    double V2sub_ref = V2_ref - V2_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low;
    double V2sube_ref = sqrt(V2e_ref*V2e_ref + V2e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low*V2e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low);
    
    double V3sub_ref = V3_ref - V3_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low;
    double V3sube_ref = sqrt(V3e_ref*V3e_ref + V3e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low*V3e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low);
    
    double V1sub_ref = V1_ref - V1_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low;
    double V1sube_ref = sqrt(V1e_ref*V1e_ref + V1e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low*V1e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low);
    
    double v2sub_ref = sqrt(V2sub_ref);
    double v2sube_ref = sqrt(V2sub_ref)*(V2sube_ref/V2sub_ref)/2;

    double v3sub_ref = sqrt(V3sub_ref);
    double v3sube_ref = sqrt(V3sub_ref)*(V3sube_ref/V3sub_ref)/2;
    
    double v1sub_ref = sqrt(V1sub_ref);
    double v1sube_ref = sqrt(V1sub_ref)*(V1sube_ref/V1sub_ref)/2;

    double V2sub[12];
    double V2sube[12];
    
    double V3sub[12];
    double V3sube[12];

    double V1sub[12];
    double V1sube[12];

    double jetYfactor[12];
    
    for(int i=0;i<12;i++)
    {
	Nassoc_low[i] = Nassoc_fit_low[i];
	Nassoc[i] = Nassoc_fit[i];
        V2sub[i] = V2[i] - V2_low[i]*Nassoc_low[i]/Nassoc[i]*suby[i]/suby_low[i];
        V2sube[i] = sqrt(V2e[i]*V2e[i] + sqrt(V2e_low[i]/V2_low[i]*V2e_low[i]/V2_low[i] + subye[i]/suby[i]*subye[i]/suby[i] + subye_low[i]/suby_low[i]*subye_low[i]/suby_low[i])*V2_low[i]*Nassoc_low[i]/Nassoc[i]*suby[i]/suby_low[i]*sqrt(V2e_low[i]/V2_low[i]*V2e_low[i]/V2_low[i] + subye[i]/suby[i]*subye[i]/suby[i] + subye_low[i]/suby_low[i]*subye_low[i]/suby_low[i])*V2_low[i]*Nassoc_low[i]/Nassoc[i]*suby[i]/suby_low[i]);
        
        V3sub[i] = V3[i] - V3_low[i]*Nassoc_low[i]/Nassoc[i]*suby[i]/suby_low[i];
        V3sube[i] = sqrt(V3e[i]*V3e[i] + sqrt(V3e_low[i]/V3_low[i]*V3e_low[i]/V3_low[i] + subye[i]/suby[i]*subye[i]/suby[i] + subye_low[i]/suby_low[i]*subye_low[i]/suby_low[i])*V3_low[i]*Nassoc_low[i]/Nassoc[i]*suby[i]/suby_low[i]*sqrt(V3e_low[i]/V3_low[i]*V3e_low[i]/V3_low[i] + subye[i]/suby[i]*subye[i]/suby[i] + subye_low[i]/suby_low[i]*subye_low[i]/suby_low[i])*V3_low[i]*Nassoc_low[i]/Nassoc[i]*suby[i]/suby_low[i]);

        V1sub[i] = V1[i] - V1_low[i]*Nassoc_low[i]/Nassoc[i]*suby[i]/suby_low[i];
        V1sube[i] = sqrt(V1e[i]*V1e[i] + sqrt(V1e_low[i]/V1_low[i]*V1e_low[i]/V1_low[i] + subye[i]/suby[i]*subye[i]/suby[i] + subye_low[i]/suby_low[i]*subye_low[i]/suby_low[i])*V1_low[i]*Nassoc_low[i]/Nassoc[i]*suby[i]/suby_low[i]*sqrt(V1e_low[i]/V1_low[i]*V1e_low[i]/V1_low[i] + subye[i]/suby[i]*subye[i]/suby[i] + subye_low[i]/suby_low[i]*subye_low[i]/suby_low[i])*V1_low[i]*Nassoc_low[i]/Nassoc[i]*suby[i]/suby_low[i]);
        
        jetYfactor[i] = suby[i]/suby_low[i];
        
    }
    
    
    
    double v2sub[12];
    double v2sube[12];
    double v2[12];
    double v2e[12];
    double v2_low[12];
    double v2e_low[12];
    
    double v3sub[12];
    double v3sube[12];
    double v3[12];
    double v3e[12];
    double v3_low[12];
    double v3e_low[12];
    
    double v1sub[12];
    double v1sube[12];
    double v1[12];
    double v1e[12];
    double v1_low[12];
    double v1e_low[12];
    
    for(int i=0;i<12;i++)
    {
        v2sub[i] = V2sub[i]/v2sub_ref;
        v2sube[i] = fabs(sqrt(V2sube[i]/V2sub[i]*V2sube[i]/V2sub[i] + v2sube_ref/v2sub_ref*v2sube_ref/v2sub_ref)*v2sub[i]);
        
        v2[i] = V2[i]/v2_ref;
        v2e[i] = sqrt(V2e[i]/V2[i]*V2e[i]/V2[i] + v2e_ref/v2_ref*v2e_ref/v2_ref)*v2[i];
        v2_low[i] = V2_low[i]/v2_ref_low;
        v2e_low[i] = sqrt(V2e_low[i]/V2_low[i]*V2e_low[i]/V2_low[i] + v2e_ref_low/v2_ref_low*v2e_ref_low/v2_ref_low)*v2_low[i];
        
        v3sub[i] = V3sub[i]/v3sub_ref;
        v3sube[i] = fabs(sqrt(V3sube[i]/V3sub[i]*V3sube[i]/V3sub[i] + v3sube_ref/v3sub_ref*v3sube_ref/v3sub_ref)*v3sub[i]);
        
        v3[i] = V3[i]/v3_ref;
        v3e[i] = sqrt(V3e[i]/V3[i]*V3e[i]/V3[i] + v3e_ref/v3_ref*v3e_ref/v3_ref)*v3[i];
        v3_low[i] = V3_low[i]/v3_ref_low;
        v3e_low[i] = sqrt(V3e_low[i]/V3_low[i]*V3e_low[i]/V3_low[i] + v3e_ref_low/v3_ref_low*v3e_ref_low/v3_ref_low)*v3_low[i];

        v1sub[i] = V1sub[i]/v1sub_ref;
        v1sube[i] = fabs(sqrt(V1sube[i]/V1sub[i]*V1sube[i]/V1sub[i] + v1sube_ref/v1sub_ref*v1sube_ref/v1sub_ref)*v1sub[i]);
        
        v1[i] = V1[i]/v1_ref;
        v1e[i] = sqrt(V1e[i]/V1[i]*V1e[i]/V1[i] + v1e_ref/v1_ref*v1e_ref/v1_ref)*v1[i];
        v1_low[i] = V1_low[i]/v1_ref_low;
        v1e_low[i] = sqrt(V1e_low[i]/V1_low[i]*V1e_low[i]/V1_low[i] + v1e_ref_low/v1_ref_low*v1e_ref_low/v1_ref_low)*v1_low[i];
    }
    
    TH1D* hPt[12];
    
    double pt[12];
    
    for(int i=0;i<12;i++)
    {
        _file3->GetObject(Form("pp_MB10_GplusPP/Pt%d",i),hPt[i]);

        pt[i] = hPt[i]->GetMean(1);
        
    }
    
	double perisubfactor[12];
	double perisubfactore[12];
	for(int i=0;i<12;i++)
	{
		perisubfactor[i] = V2_low[i]*Nassoc_low[i]/suby_low[i];
		perisubfactore[i] = sqrt((V2e_low[i]/V2_low[i])*(V2e_low[i]/V2_low[i]) + subye_low[i]/suby_low[i]*subye_low[i]/suby_low[i])*V2_low[i]*Nassoc_low[i]/suby_low[i];
	}


    TGraphErrors* V2plot = new TGraphErrors(12,pt,V2,0,V2e);
    
    TGraphErrors* V2_subplot = new TGraphErrors(12,pt,V2sub,0,V2sube);
    
    TGraphErrors* V3plot = new TGraphErrors(12,pt,V3,0,V3e);
    
    TGraphErrors* V3_subplot = new TGraphErrors(12,pt,V3sub,0,V3sube);

    TGraphErrors* V1plot = new TGraphErrors(12,pt,V1,0,V2e);
    
    TGraphErrors* V1_subplot = new TGraphErrors(12,pt,V1sub,0,V1sube);

    
    TFile ofile("v2_vspt_sub1020_NassFit_highpt_pt033.root","RECREATE");

    
    V2plot->SetName("hadronV2_GplusPP");
    V2_subplot->SetName("hadronV2sub_GplusPP");
    
    V2plot->Write();
    V2_subplot->Write();

    V3plot->SetName("hadronV3_GplusPP");
    V3_subplot->SetName("hadronV3sub_GplusPP");
    
    V3plot->Write();
    V3_subplot->Write();
    
    V1plot->SetName("hadronV1_GplusPP");
    V1_subplot->SetName("hadronV1sub_GplusPP");
    
    V1plot->Write();
    V1_subplot->Write();

    TGraphErrors* v2plot = new TGraphErrors(12,pt,v2,0,v2e);

    TGraphErrors* v2_subplot = new TGraphErrors(12,pt,v2sub,0,v2sube);
    
    TGraphErrors* v3plot = new TGraphErrors(12,pt,v3,0,v3e);
    
    TGraphErrors* v3_subplot = new TGraphErrors(12,pt,v3sub,0,v3sube);

    TGraphErrors* v1plot = new TGraphErrors(12,pt,v1,0,v1e);
    
    TGraphErrors* v1_subplot = new TGraphErrors(12,pt,v1sub,0,v1sube);

    v2plot->SetName("hadronv2_GplusPP");
    v2_subplot->SetName("hadronv2sub_GplusPP");

    v2plot->Write();
    v2_subplot->Write();
    
    v3plot->SetName("hadronv3_GplusPP");
    v3_subplot->SetName("hadronv3sub_GplusPP");
    
    v3plot->Write();
    v3_subplot->Write();

    v1plot->SetName("hadronv1_GplusPP");
    v1_subplot->SetName("hadronv1sub_GplusPP");
    
    v1plot->Write();
    v1_subplot->Write();

    TGraphErrors* sryplot = new TGraphErrors(12,pt,sry,0,srye);
    TGraphErrors* lryplot = new TGraphErrors(12,pt,lry,0,lrye);
    TGraphErrors* subyplot = new TGraphErrors(12,pt,suby,0,subye);
    TGraphErrors* sryplot_low = new TGraphErrors(12,pt,sry_low,0,srye_low);
    TGraphErrors* lryplot_low = new TGraphErrors(12,pt,lry_low,0,lrye_low);
    TGraphErrors* subyplot_low = new TGraphErrors(12,pt,suby_low,0,subye_low);

	TGraphErrors* perisubfactorplot = new TGraphErrors(11,pt,perisubfactor,0,perisubfactore);

	sryplot->SetName("sry_GplusPP");
	lryplot->SetName("lry_GplusPP");
	subyplot->SetName("suby_GplusPP");
	sryplot_low->SetName("srylow_GplusPP");
        lryplot_low->SetName("lrylow_GplusPP");
        subyplot_low->SetName("subylow_GplusPP");
	perisubfactorplot->SetName("perisubfactor_GplusPP");

	sryplot->Write();
	lryplot->Write();
	subyplot->Write();
	sryplot_low->Write();
        lryplot_low->Write();
        subyplot_low->Write();
	perisubfactorplot->Write();

    for(int i=0;i<12;i++)
    {
        cout<<"jet factor: "<<jetYfactor[i]<<endl;
    }
    
    cout<<"Nassoc: "<<Nassoc<<"  Nassoc_low: "<<Nassoc_low<<endl;
    
cout<<"jet yield: high:"<<suby_ref<<" low:"<<suby_ref_low<<endl;
    cout<<"sr yield: high:"<<sry_ref<<" low:"<<sry_ref_low<<endl;
    cout<<"lr yield: high:"<<lry_ref<<" low:"<<lry_ref_low<<endl;

    
    cout<<"Ref V2: before sub:"<<V2_ref<<" after sub:"<<V2sub_ref<<"  V2low:"<<V2_ref_low<<endl;
cout<<"Ref v2: before sub:"<<v2_ref<<" after sub:"<<v2sub_ref<<"  v2low:"<<v2_ref_low<<endl;
    cout<<"V2 vs pT before sub:"<<endl;
    for(int i=0;i<12;i++)
    {
        cout<<v2[i]<<"   "<<v2e[i]<<endl;
	//cout<<"sry:"<<sry[i]<<endl;
	//cout<<"sry_low:"<<sry_low[i]<<endl;
	//cout<<"lry:"<<lry[i]<<endl;
	//cout<<"lry_low:"<<lry_low[i]<<endl;
	cout<<"Nassfit:"<<Nassoc[i]<<endl;
	cout<<"Nassfit_low:"<<Nassoc_low[i]<<endl;
    }
    cout<<"V2 vs pT after sub:"<<endl;
    for(int i=0;i<12;i++)
    {
        cout<<v2sub[i]<<"   "<<v2sube[i]<<endl;
    }
    cout<<"V2_low:"<<endl;
    for(int i=0;i<12;i++)
    {
        cout<<v2_low[i]<<"   "<<v2e_low[i]<<endl;
    }

for(int i=0;i<12;i++)
{
	cout<<lry[i]<<endl;
}
    
    TCanvas* plot = new TCanvas("plot","plot",400,400);
    plot->cd();
    //proj_low->Scale(1.0/proj_low->Integral());
    //proj_high->Scale(1.0/proj_high->Integral());
    proj_low->Scale(suby[11]/suby_low[11]);
    proj_low->SetMarkerColor(2);
    proj_low->Draw();
    proj_high->Draw("SAME");
    
    TLegend* l = new TLegend(0.2,0.6,0.6,0.76);
    l->AddEntry(proj_low,"(10#leqN_{trk}^{offline}<20)#times(Y_{jet}^{high}/Y_{jet}^{low}) ","p");
    l->AddEntry(proj_high,"110 #leq N_{trk}^{offline} < 150 ","p");
    l->Draw();
    
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(tex->GetTextSize()*0.8);
    tex->DrawLatex(0.22,0.88,"CMS pp 13TeV");
    tex->DrawLatex(0.22,0.82,"0.3<p_{T}^{assoc}<3 GeV");
    tex->DrawLatex(0.22,0.76,"6<p_{T}^{trg}<9 GeV");
}
