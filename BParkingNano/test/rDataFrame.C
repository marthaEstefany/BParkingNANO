#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"


using namespace ROOT::VecOps;


void rDataFrame() {

 // input files
 ROOT::EnableImplicitMT();
 ROOT::RDataFrame df("Events","/eos/uscms/store/user/manunezo/BParkingNANO_2021Aug23/ParkingBPH1/crab_data_Run2018A_part1/210823_144831/0000/BParkNANO_data_2021Jul13_*.root");
 
 // selections
 auto df_2mu = df.Filter("nMuon == 2", "Events with exactly two muons");
 auto df_os = df_2mu.Filter("Muon_charge[0] != Muon_charge[1]", "Muons with opposite charge");
 
 auto df_mass = df_os.Define("Dimuon_mass", InvariantMass<float>, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});
 auto h = df_mass.Histo1D({"Dimuon_mass", "Dimuon mass;m_{#mu#mu} (GeV);N_{Events}", 30000, 0.25, 300}, "Dimuon_mass");
 
 // Request cut-flow report
 auto report = df_mass.Report();
 
 // Produce plot
 gStyle->SetOptStat(0); gStyle->SetTextFont(42);
 auto c = new TCanvas("c", "", 800, 700);
 c->SetLogx(); c->SetLogy();
 
 h->GetXaxis()->SetTitleSize(0.04);
 h->GetYaxis()->SetTitleSize(0.04);
 h->DrawClone();

 TLatex label; label.SetNDC(true);
 label.DrawLatex(0.175, 0.740, "#eta");
 label.DrawLatex(0.205, 0.775, "#rho,#omega");
 label.DrawLatex(0.270, 0.740, "#phi");
 label.DrawLatex(0.400, 0.800, "J/#psi");
 label.DrawLatex(0.415, 0.670, "#psi'");
 label.DrawLatex(0.485, 0.700, "Y(1,2,3S)");
 label.DrawLatex(0.755, 0.680, "Z");
 label.SetTextSize(0.040); label.DrawLatex(0.100, 0.920, "#bf{CMS 2018 BParked Data}");
 label.SetTextSize(0.030); label.DrawLatex(0.630, 0.920, "#sqrt{s} = 13 TeV, L_{int} = __ fb^{-1}");

 c->SaveAs("dimuon_spectrum.pdf");

 // Print cut-flow report
 report->Print();

}
