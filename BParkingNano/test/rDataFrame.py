import ROOT 

ROOT.ROOT.EnableImplicitMT()

df = ROOT.RDataFrame("Events", "/eos/uscms/store/user/manunezo/BParkingNANO_2021Aug23/ParkingBPH1/crab_data_Run2018A_part1/210823_144831/0000/BParkNANO_data_2021Aug23_*.root")

df_2mu = df.Filter("nMuon >= 2", "Events with two or more muons")
df_eta = df_2mu.Filter("abs(Muon_eta[0]) < 2.4 && abs(Muon_eta[1]) < 2.4", "leading and subleading muon |eta| < 2.4")
df_tightId = df_eta.Filter("Muon_tightId[0] == 1 && Muon_tightId[1] == 1", "leading and subleading muons have tight id")
df_Iso = df_tightId.Filter("Muon_pfRelIso04_custom[0] < 0.10 && Muon_pfRelIso04_custom[1] < 0.10", "leading and subleading muon are isolated with custome isolation < 0.1")
df_global = df_Iso.Filter("Muon_isGlobal[0] == 1 && Muon_isGlobal[1] == 1", "leading and subleading muons are global")
df_deltaR = df_global.Filter("((Muon_eta[0]-Muon_eta[1])*(Muon_eta[0]-Muon_eta[1]) + (Muon_phi[0]-Muon_phi[1])*(Muon_phi[0]-Muon_phi[1])) > 0.2", "deltaR between leading and subleading muons > 0.2")
df_cosA = df_deltaR.Filter("( Muon_pt[0]*cos(Muon_phi[0])*Muon_pt[1]*cos(Muon_phi[1]) +  Muon_pt[0]*sin(Muon_phi[0])*Muon_pt[1]*sin(Muon_phi[1]) + Muon_pt[0]/tan(Muon_phi[0])*Muon_pt[1]/tan(Muon_phi[1]) ) / ( sqrt(Muon_pt[0]*cos(Muon_phi[0])*Muon_pt[0]*cos(Muon_phi[0]) + Muon_pt[1]*cos(Muon_phi[1])*Muon_pt[1]*cos(Muon_phi[1]) + Muon_pt[0]/tan(Muon_phi[0])*Muon_pt[0]/tan(Muon_phi[0]) ) * sqrt( Muon_pt[1]*cos(Muon_phi[1])*Muon_pt[1]*cos(Muon_phi[1]) + Muon_pt[1]*sin(Muon_phi[1])*Muon_pt[1]*sin(Muon_phi[1]) + Muon_pt[1]/tan(Muon_phi[1])*Muon_pt[1]/tan(Muon_phi[1])  )) > -0.99", "cos(alpha) btw leading and subleading muons > -0.99")
df_prompt = df_cosA.Filter("10000*abs(Muon_dxy[0]) < 50  && 10000*abs(Muon_dxy[1]) < 50", "leading and subleading muon |dxy| < 50 #mum")


# Perform the computation of the invariant mass in C++
ROOT.gInterpreter.Declare('''
using Vec_t = const ROOT::RVec<float>&;
float ComputeInvariantMass(Vec_t pt, Vec_t eta, Vec_t phi, Vec_t mass) {
    const ROOT::Math::PtEtaPhiMVector p1(pt[0], eta[0], phi[0], mass[0]);
    const ROOT::Math::PtEtaPhiMVector p2(pt[1], eta[1], phi[1], mass[1]);
    return (p1 + p2).M();
}
''')


# Add the result of the computation to the dataframe
df_mass = df_prompt.Define("Dimuon_mass", "ComputeInvariantMass(Muon_pt, Muon_eta, Muon_phi, Muon_mass)")
df_delR = df_mass.Define("DeltaR", "(Muon_eta[0]-Muon_eta[1])*(Muon_eta[0]-Muon_eta[1]) + (Muon_phi[0]-Muon_phi[1])*(Muon_phi[0]-Muon_phi[1])")
df_cosa = df_delR.Define("CosA", "( Muon_pt[0]*cos(Muon_phi[0])*Muon_pt[1]*cos(Muon_phi[1]) +  Muon_pt[0]*sin(Muon_phi[0])*Muon_pt[1]*sin(Muon_phi[1]) + Muon_pt[0]/tan(Muon_phi[0])*Muon_pt[1]/tan(Muon_phi[1]) ) / ( sqrt(Muon_pt[0]*cos(Muon_phi[0])*Muon_pt[0]*cos(Muon_phi[0]) + Muon_pt[1]*cos(Muon_phi[1])*Muon_pt[1]*cos(Muon_phi[1]) + Muon_pt[0]/tan(Muon_phi[0])*Muon_pt[0]/tan(Muon_phi[0]) ) * sqrt( Muon_pt[1]*cos(Muon_phi[1])*Muon_pt[1]*cos(Muon_phi[1]) + Muon_pt[1]*sin(Muon_phi[1])*Muon_pt[1]*sin(Muon_phi[1]) + Muon_pt[1]/tan(Muon_phi[1])*Muon_pt[1]/tan(Muon_phi[1])  ))")
df_dxySub = df_cosa.Define("Muon_dxy_um_sub", "10000*abs(Muon_dxy[1])")
df_dxyLead = df_dxySub.Define("Muon_dxy_um_lead", "10000*abs(Muon_dxy[0])")

# Book histogram of the dimuon mass spectrum (does not actually run the computation!)
h =  df_dxyLead.Histo1D(("Dimuon_mass", ";m_{#mu#mu} (GeV);N_{Events}", 30000, 0.25, 300), "Dimuon_mass")
h1 = df_dxyLead.Histo1D(("Muon_eta", ";#mu_{|eta|} ;N_{Events}", 30000, 0, 3), "Muon_eta")
h2 = df_dxyLead.Histo1D(("Muon_tightId", ";#mu tight Id ;N_{Events}", 2, 0, 2), "Muon_tightId")
h3 = df_dxyLead.Histo1D(("Muon_pfRelIso04_custom", ";#mu_{Iso} ;N_{Events}", 30000, 0, 3), "Muon_pfRelIso04_custom")
h4 = df_dxyLead.Histo1D(("DeltaR", ";#DeltaR ;N_{Events}", 30000, 0, 3), "DeltaR")
h5 = df_dxyLead.Histo1D(("CosA", ";Cos(#alpha) ;N_{Events}", 30000, 0, 3), "CosA")
h6 = df_dxyLead.Histo2D(("|d0|", "|d0|", 100, 0., 50., 100, 0., 50.), "Muon_dxy_um_sub", "Muon_dxy_um_lead")



# Request a cut-flow report (also does not run the computation yet!)
dfFinal = df_dxyLead
report = dfFinal.Report()

dfFinal.Snapshot("Events", "skim.root")

# Produce plot
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextFont(42)
c = ROOT.TCanvas("c", "", 800, 700)
# The contents of one of the dataframe results are accessed for the first time here:
# this is when all results will actually be produced!
#h.GetXaxis().SetTitleSize(0.04)
#h.GetYaxis().SetTitleSize(0.04)
#c.SetLogx(); c.SetLogy()
h.Draw()

c1 = ROOT.TCanvas("c1", "", 800, 700)
h1.Draw()

c2 = ROOT.TCanvas("c2", "", 800, 700)
h2.Draw()

c3 = ROOT.TCanvas("c3", "", 800, 700)
h3.Draw()

c4 = ROOT.TCanvas("c4", "", 800, 700)
h4.Draw()

c5 = ROOT.TCanvas("c5", "", 800, 700)
h5.Draw()

c6 = ROOT.TCanvas("c6", "", 800, 700)
h6.Draw("COLZ")

label = ROOT.TLatex()
label.SetNDC(True)
label.SetTextSize(0.040)
label.DrawLatex(0.100, 0.920, "#bf{CMS 2018 Parked Data}")
label.DrawLatex(0.550, 0.920, "#sqrt{s} = 13 TeV, L_{int} = __ fb^{-1}")

#Save as png file
c.SaveAs("dimuon_spectrum.png")
c1.SaveAs("eta.png")
c2.SaveAs("tightId.png")
c3.SaveAs("iso.png")
c4.SaveAs("deltaR.png")
c5.SaveAs("cosA.png")
c6.SaveAs("d0.png")


#Print cut-flow report
report.Print()

