void addN()
{
    TFile* file[12];
    file[0] = TFile::Open("../MB0_10/ppCorr_pt033.root");
    file[1] = TFile::Open("../MB10_20/ppCorr_pt033.root");

    TH2D* signal[11];
    TH2D* background[11];
    TH1D* mult[11];
    TH1D* ntrk[11];
    TH1D* hPt[11];
    
    signal[0] = (TH2D*)file[0]->Get("pp_MB0_GplusPP/signal");
    background[0] = (TH2D*)file[0]->Get("pp_MB0_GplusPP/background");
    mult[0] = (TH1D*)file[0]->Get("pp_MB0_GplusPP/mult");
    ntrk[0] = (TH1D*)file[0]->Get("pp_MB0_GplusPP/mult_good");
    hPt[0] = (TH1D*)file[0]->Get("pp_MB0_GplusPP/pT");
    signal[1] = (TH2D*)file[1]->Get("pp_MB10_GplusPP/signal");
    background[1] = (TH2D*)file[1]->Get("pp_MB10_GplusPP/background");
    mult[1] = (TH1D*)file[1]->Get("pp_MB10_GplusPP/mult");
    ntrk[1] = (TH1D*)file[1]->Get("pp_MB10_GplusPP/mult_good");
    hPt[1] = (TH1D*)file[1]->Get("pp_MB10_GplusPP/pT");
    
    signal[0]->Add(signal[1]);
    background[0]->Add(background[1]);
    mult[0]->Add(mult[1]);
    ntrk[0]->Add(ntrk[1]);
    hPt[0]->Add(hPt[1]);
    
    TFile ofile("ppCorr_pt033.root","RECREATE");

    signal[0]->Write();
    background[0]->Write();
    mult[0]->Write();
    ntrk[0]->Write();
    hPt[0]->Write();
}
