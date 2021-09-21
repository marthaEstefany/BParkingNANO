
void plot_macro() {

  // offline selections from BToKMuMu analysis
  const TCut offline_selections = "(BToKMuMu_fit_cos2D > 0.98 || BToKMuMu_fit_cos2D < -0.98) && ((BToKMuMu_mll_fullfit < 2.929) || (BToKMuMu_mll_fullfit > 3.263 && BToKMuMu_mll_fullfit < 3.574) || (BToKMuMu_mll_fullfit > 3.798)) && (BToKMuMu_fit_pt > 5.34) && (BToKMuMu_svprob > 0.083) && (BToKMuMu_l_xy/BToKMuMu_l_xy_unc > 4.5) && (BToKMuMu_fit_k_pt > 1.2)";

  // Displaced Lepton mumu selections
  const TCut mumu_selections = "nMuon >= 2 && abs(Muon_eta) < 2.4 && Muon_tightId  && Muon_pfRelIso04_custom < 0.10 && Muon_isGlobal && !(Muon_eta >= 0.3 && Muon_eta <= 1.2 && Muon_phi >= 0.4 && Muon_phi <= 0.8) && ((Muon_eta[0]-Muon_eta[1])*(Muon_eta[0]-Muon_eta[1]) + (Muon_phi[0]-Muon_phi[1])*(Muon_phi[0]-Muon_phi[1])) > 0.2 && ( Muon_pt[0]*cos(Muon_phi[0])*Muon_pt[1]*cos(Muon_phi[1]) +  Muon_pt[0]*sin(Muon_phi[0])*Muon_pt[1]*sin(Muon_phi[1]) + Muon_pt[0]/tan(Muon_phi[0])*Muon_pt[1]/tan(Muon_phi[1]) ) / ( sqrt(   Muon_pt[0]*cos(Muon_phi[0])*Muon_pt[0]*cos(Muon_phi[0]) + Muon_pt[1]*cos(Muon_phi[1])*Muon_pt[1]*cos(Muon_phi[1]) + Muon_pt[0]/tan(Muon_phi[0])*Muon_pt[0]/tan(Muon_phi[0]) ) * sqrt( Muon_pt[1]*cos(Muon_phi[1])*Muon_pt[1]*cos(Muon_phi[1]) + Muon_pt[1]*sin(Muon_phi[1])*Muon_pt[1]*sin(Muon_phi[1]) + Muon_pt[1]/tan(Muon_phi[1])*Muon_pt[1]/tan(Muon_phi[1])  )) > -0.99";

  const TCut onePrompt_oneDisplaced = "(10000*abs(Muon_dxy[0]) < 100 && 10000*abs(Muon_dxy[1]) > 100)  || (10000*abs(Muon_dxy[0]) > 100 &&  10000*abs(Muon_dxy[1]) < 100)";
  const TCut prompt_leading = "10000*abs(Muon_dxy[0]) < 100 && 10000*abs(Muon_dxy[1]) > 100";
  const TCut prompt_subleading = "10000*abs(Muon_dxy[1]) < 100)  && (10000*abs(Muon_dxy[0]) > 100";
  const TCut prompt_highPt = "10000*abs(Muon_dxy[0]) < 50)  && (10000*abs(Muon_dxy[1]) < 50";
  const TCut promp_lowPt = "10000*abs(Muon_dxy[0]) < 50)  && (10000*abs(Muon_dxy[1]) < 50 && Muon_pt[0] < 40 && Muon_pt[1] < 40";

  const TCut control_region1 = operator &&(mumu_selections, onePrompt_oneDisplaced);
  const TCut control_region2 = operator &&(mumu_selections, prompt_leading);
  const TCut control_region3 = operator &&(mumu_selections, prompt_subleading);
  const TCut control_prompt = operator &&(mumu_selections, prompt_highPt);

  TChain chain("Events");
  //TFileCollection fc("file", "", "data_2021July13/data_13July2021/data.txt");
  TFileCollection fc("file", "", "list.txt");
  chain.AddFileInfoList(fc.GetList());
  cout << "number of leaves: " << chain.GetNbranches() << endl;

  TCanvas *c30 = new TCanvas("c30","c30",500,500);
  TH1I* h30 = new TH1I("hist_30", "hist_30", 1000, 0, 2);
  chain.Project("hist_30", "Muon_pfRelIso04_all", control_prompt);
  h30->SetTitle("prompt control region");
  h30->GetXaxis()->SetTitle("default iso");
  h30->Draw("HIST EO");

  TCanvas *c31 = new TCanvas("c31","c31",500,500);
  TH1I* h31 = new TH1I("hist_31", "hist_31", 1000, 0, 2);
  chain.Project("hist_31", "Muon_pfRelIso04_custom", control_prompt);
  h31->SetTitle("prompt control region");
  h31->GetXaxis()->SetTitle("custom iso");
  h31->Draw("HIST EO");

  TCanvas *c22 = new TCanvas("c22","c22",500,500);
  TH1I* h22 = new TH1I("hist_22", "hist_22", 1000, 0, 1);
  chain.Project("hist_22", "((Muon_eta[0]-Muon_eta[1])*(Muon_eta[0]-Muon_eta[1]) + (Muon_phi[0]-Muon_phi[1])*(Muon_phi[0]-Muon_phi[1]))", control_prompt);
  h22->SetTitle("prompt control region");
  h22->GetXaxis()->SetTitle("deltaR");
  h22->Draw("HIST EO");

  TCanvas *c26 = new TCanvas("c26","c26",500,500);
  TH1I* h26 = new TH1I("hist_26", "hist_26", 400, -1, 1);
  chain.Project("hist_26", "( Muon_pt[0]*cos(Muon_phi[0])*Muon_pt[1]*cos(Muon_phi[1]) +  Muon_pt[0]*sin(Muon_phi[0])*Muon_pt[1]*sin(Muon_phi[1]) + Muon_pt[0]/tan(Muon_phi[0])*Muon_pt[1]/tan(Muon_phi[1]) ) / ( sqrt(   Muon_pt[0]*cos(Muon_phi[0])*Muon_pt[0]*cos(Muon_phi[0]) + Muon_pt[1]*cos(Muon_phi[1])*Muon_pt[1]*cos(Muon_phi[1]) + Muon_pt[0]/tan(Muon_phi[0])*Muon_pt[0]/tan(Muon_phi[0]) ) * sqrt( Muon_pt[1]*cos(Muon_phi[1])*Muon_pt[1]*cos(Muon_phi[1]) + Muon_pt[1]*sin(Muon_phi[1])*Muon_pt[1]*sin(Muon_phi[1]) + Muon_pt[1]/tan(Muon_phi[1])*Muon_pt[1]/tan(Muon_phi[1])  ))", control_prompt);
  h26->SetTitle("prompt control region");
  h26->GetXaxis()->SetTitle("cos(#alpha)");
  h26->Draw("HIST EO");

  TCanvas *c29 = new TCanvas("c29","c29",500,500);
  TH2I* h29 = new TH2I("hist_29", "hist_29", 50, 0, 50, 50, 0, 50);
  chain.Project("hist_29", "10000*abs(Muon_dxy[0]):10000*abs(Muon_dxy[1])", control_prompt);
  h29->SetTitle("prompt control region");
  h29->GetXaxis()->SetTitle("subleading #mu |d0|");
  h29->GetYaxis()->SetTitle("leading #mu |d0|");
  h29->Draw("COLZ");

}
