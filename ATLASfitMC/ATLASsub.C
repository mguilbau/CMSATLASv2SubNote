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

void ATLASsub()
{
    TH1::SetDefaultSumw2();
    
    double BW2D = (9.9/33)*((2-1.0/16)*3.1416/31);
    double errfactor = sqrt(1.0);
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(tex->GetTextSize()*1.5);

    TFile* _file[20];
    TFile* _file_low;
    _file_low = TFile::Open("MB0_10/ppCorr_pt033.root");
    _file[0] = TFile::Open("MB10_20/ppCorr_pt033.root");
    _file[1] = TFile::Open("MB20_30/ppCorr_pt033.root");
    _file[2] = TFile::Open("MB30_40/ppCorr_pt033.root");
    _file[3] = TFile::Open("MB40_60/ppCorr_pt033.root");
    _file[4] = TFile::Open("MB60_85/ppCorr_pt033.root");
    _file[5] = TFile::Open("MB85_95/ppCorr_pt033.root");
    _file[6] = TFile::Open("MB95_110/ppCorr_pt033.root");
    _file[7] = TFile::Open("MB110_above/ppCorr_pt033.root");
    //_file[8] = TFile::Open("../HM115_125/ppCorr_pt033.root");
    //_file[9] = TFile::Open("../HM125_135/ppCorr_pt033.root");
    //_file[10] = TFile::Open("../HM135_150/ppCorr_pt033.root");
    //_file[11] = TFile::Open("../HM150_above/ppCorr_pt033.root");

    TH2D* signal_low;
    TH2D* background_low;
    TH1D* mult_low;
    double nEvent_low;
    double Bz_low;
    TH1D* alrs_low;
    
    TCanvas* c = new TCanvas("c","c",600,400);
    c->cd();
    
    _file_low->GetObject("pp_MB0_GplusPP/signal",signal_low);
    _file_low->GetObject("pp_MB0_GplusPP/background",background_low);
    _file_low->GetObject("pp_MB0_GplusPP/mult",mult_low);
    
    nEvent_low = mult_low->Integral(3,10000);
    Bz_low = background_low->GetBinContent(background_low->FindBin(0,0));
    double N_low = mult_low->GetMean(1);
    
    //lr Y
    (TH1D*)alrs_low = signal_low->ProjectionY("alrslow",1,10);
    alrs1 = signal_low->ProjectionY("alrs1",24,33);
    alrb = background_low->ProjectionY("alrb",1,10);
    alrb1 = background_low->ProjectionY("alrb1",24,33);
    alrs_low->Add(alrs1);
    alrb->Add(alrb1);
    alrs_low->Divide(alrb);
    alrs_low->Scale(Bz_low/nEvent_low/BW2D);
    
    TF1* fit2 = new TF1("fit2","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit2->SetLineColor(2);
    fit2->SetParNames("N","V1","V2","V3","V4");
    fit2->SetParameters(10,1,1,1,1);
    
    alrs_low->Fit("fit2","R");
    alrs_low->Fit("fit2","R");
    alrs_low->Fit("fit2","R");
    alrs_low->Fit("fit2","R");
    
    double Nlow = fit2->GetParameter(0);
    double V1low = fit2->GetParameter(1);
    double V2low = fit2->GetParameter(2);
    double V3low = fit2->GetParameter(3);
    double V4low = fit2->GetParameter(4);
    
    TF1* fit = new TF1("fit",Form("%f*(1.0+2.0*%f*cos(x)+2.0*%f*cos(2.0*x)+2.0*%f*cos(3.0*x)+2.0*%f*cos(4.0*x))*[0] + [1]*(1.0+2.0*[2]*cos(2.0*x))",Nlow,V1low,V2low,V3low,V4low),-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit->SetLineColor(2);
    fit->SetParNames("F","G","V2");
    fit->SetParameters(1,1,1);
    
    TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetLineColor(2);
    fit1->SetParNames("N","V1","V2","V3","V4");
    fit1->SetParameters(1,1,1,1,1);

    
    TCanvas* c1 = new TCanvas("StdFit","StdFit",1200,800);
    c1->Divide(4,3);
    
    TCanvas* c2 = new TCanvas("ATLASFit","ATLASFit",1200,800);
    c2->Divide(4,3);
    
    TH2D* signal[20];
    TH2D* background[20];
    TH1D* mult[20];
    TH1D* ntrk[20];
    int nEvent[20];
    double Bz[20];
    TH1D* alrs[20];
    TH1D* alrsATLAS[20];
    double V2[20];
    double V2sub[20];
    double V2e[20];
    double V2sube[20];
    double Ntrk[20];
    double Nassoc[20];
    double G[20];
    double Ge[20];
    double F[20];
    double Fe[20];
    
    char Nstring[20][200] = {"MB10","MB20","MB30","MB40","MB60","MB85","MB95","MB110","HM115","HM125","HM135","HM150"};
    
    for(int i=0;i<8;i++)
    {
        cout<<i<<endl;
        _file[i]->GetObject(Form("pp_%s_GplusPP/signal",Nstring[i]),signal[i]);
        _file[i]->GetObject(Form("pp_%s_GplusPP/background",Nstring[i]),background[i]);
        _file[i]->GetObject(Form("pp_%s_GplusPP/mult",Nstring[i]),mult[i]);
        _file[i]->GetObject(Form("pp_%s_GplusPP/mult_good",Nstring[i]),ntrk[i]);
    }
    
    for(int i=0;i<8;i++){
        //cout<<i<<endl;
        nEvent[i] = mult[i]->Integral(3,10000);
        Bz[i] = background[i]->GetBinContent(background[i]->FindBin(0,0));
        Ntrk[i] = ntrk[i]->GetMean(1);
        
        //lr Y
        (TH1D*)alrs[i] = signal[i]->ProjectionY(Form("alrs1%d",i),1,10);
        alrs1 = signal[i]->ProjectionY("alrs1",24,33);
        alrb = background[i]->ProjectionY("alrb",1,10);
        alrb1 = background[i]->ProjectionY("alrb1",24,33);
        alrs[i]->Add(alrs1);
        alrb->Add(alrb1);
        alrs[i]->Divide(alrb);
        alrs[i]->Scale(Bz[i]/nEvent[i]/BW2D);
        
        //Add constant V2 to PYTHIA
        alrs[i]->Fit("fit1","R");
        alrs[i]->Fit("fit1","R");
        alrs[i]->Fit("fit1","R");
        alrs[i]->Fit("fit1","R");
        Nassoc[i] = fit1->GetParameter(0);
        
        //TF1* f1 = new TF1("f1",Form("%f*2.0*%f*cos(2.0*x)",Nassoc[i],0.0),-0.5*TMath::Pi(),1.5*TMath::Pi());
        //alrs[i]->Add(f1);
        
        alrsATLAS[i] = (TH1D*)alrs[i]->Clone();
        
        c1->cd(i+1);
        
        alrs[i]->Fit("fit1","R");
        alrs[i]->Fit("fit1","R");
        alrs[i]->Fit("fit1","R");
        alrs[i]->Fit("fit1","R");
        V2[i] = fit1->GetParameter(2);
        V2e[i] = fit1->GetParError(2)*errfactor;
        
        c2->cd(i+1);
        
        alrsATLAS[i]->Fit("fit","R");
        alrsATLAS[i]->Fit("fit","R");
        alrsATLAS[i]->Fit("fit","R");
        alrsATLAS[i]->Fit("fit","R");
        
        tex->DrawLatex(0.2,0.8,Nstring[i]);

        V2sub[i] = fit->GetParameter(2);
        V2sube[i] = fit->GetParError(2)*errfactor;
        
        G[i] = fit->GetParameter(1);
        Ge[i] = fit->GetParError(1);
        F[i] = fit->GetParameter(0);
        Fe[i] = fit->GetParError(0);
        
        //if(i<=1) V2sub[i] = 99;
        //if(i<=1) V2sube[i] = 0;
    }
    
    double v2[20];
    double v2sub[20];
    double v2e[20];
    double v2sube[20];

    for(int i=0;i<8;i++)
    {
        v2[i] = sqrt(V2[i]);
        v2e[i] = sqrt(V2[i])*(V2e[i]/V2[i])/2;
        
        v2sub[i] = sqrt(V2sub[i]);
        v2sube[i] = sqrt(V2sub[i])*(V2sube[i]/V2sub[i])/2;
        cout<<"Ntrkoffline: "<<Ntrk[i]<<endl;
        cout<<"V2 from ATLAS: "<<V2sub[i]<<endl;
    }

    TGraphErrors* V2plot = new TGraphErrors(8,Ntrk,V2,0,V2e);
    TGraphErrors* V2subplot = new TGraphErrors(8,Ntrk,V2sub,0,V2sube);

    TGraphErrors* Gplot = new TGraphErrors(8,Ntrk,G,0,Ge);
    TGraphErrors* Fplot = new TGraphErrors(8,Ntrk,F,0,Fe);

    TGraphErrors* v2plot = new TGraphErrors(8,Ntrk,v2,0,v2e);
    TGraphErrors* v2subplot = new TGraphErrors(8,Ntrk,v2sub,0,v2sube);
    
    V2plot->SetName("hadronV2");
    V2subplot->SetName("hadronV2sub");

    Gplot->SetName("G");
    Fplot->SetName("F");
    
    v2plot->SetName("hadronv2");
    v2subplot->SetName("hadronv2sub");

    TFile ofile("v2sub010_ATLAS_pt033.root","RECREATE");
    V2plot->Write();
    V2subplot->Write();
    v2plot->Write();
    v2subplot->Write();
    Gplot->Write();
    Fplot->Write();
}
