double mystepfunc(double* x, double* par);
void gendt(double jitter=0)//jitter in ns
{
  TFile fout(Form("fake_data_%ins.root",int(jitter)),"RECREATE");
  // Declaration of leaf types
  Int_t           run;
  Int_t           subrun;
  Int_t           event;
  Int_t           nhit;
  Double_t        pe_total;
  Double_t        dt;
  Double_t        t_width;
  Double_t        fy;
  Double_t        fy_err;
  Double_t        fz;
  Double_t        fz_err;
  Int_t           trig_word;

  TTree* fChain=new TTree("opflash_tree_filtered","filtered tree");
  fChain->Branch("run", &run,"run/I");
  fChain->Branch("subrun", &subrun,"subrun/I");
  fChain->Branch("event", &event,"event/I");
  fChain->Branch("nhit", &nhit,"nhit/I");
  fChain->Branch("pe_total", &pe_total,"pe_total/D");
  fChain->Branch("dt", &dt,"dt/D");
  fChain->Branch("t_width", &t_width,"t_width/D");
  fChain->Branch("fy", &fy,"fy/D");
  fChain->Branch("fy_err", &fy_err,"fy_err/D");
  fChain->Branch("fz", &fz,"fz/D");
  fChain->Branch("fz_err", &fz_err,"fz_err/D");
  fChain->Branch("trig_word", &trig_word,"trig_word/I");

  TF1* myf=new TF1("myf",mystepfunc,0,23,4);
  myf->SetParameter(0,3.2);
  myf->SetParameter(1,4.8);
  myf->SetParameter(2,1);
  myf->SetParameter(3,1.4);

  myf->Draw();
  
  TH1F* hsig=new TH1F("hsig","",4600,0,23);
  for (int i=0;i<1000000;i++) {
    dt=myf->GetRandom(0,23);
    dt+=gRandom->Uniform(-jitter/2000.,jitter/2000.);
    trig_word=2048;
    hsig->Fill(dt);
    fChain->Fill();
  }
  hsig->Scale(myf->Integral(0,23)/hsig->GetBinWidth(1)/hsig->GetEntries());
  hsig->Draw();

  TH1F* hbkg=new TH1F("hbkg","",4600,0,23);
  for (int i=0;i<1000000;i++) {
    dt=gRandom->Uniform(0,23);
    dt+=gRandom->Uniform(-jitter/2000.,jitter/2000.);
    trig_word=512;
    hbkg->Fill(dt);
    fChain->Fill();
  }
  hbkg->Scale(23./hbkg->GetBinWidth(1)/hbkg->GetEntries());
  hbkg->Draw();
  

  fChain->Write();
  fout.Close();
}


double mystepfunc(double* x, double* par)
{
  double val=0;
  if (x[0]<par[0]) val=par[2];
  else if (x[0]<par[1]) val=par[3];
  else val=par[2];

  return val;
}
