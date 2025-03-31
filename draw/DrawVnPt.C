#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>
#include <iostream>
#include "../BootstrapProfile.h"

void DrawVnPt(TString deformation = "") {
  gSystem->Load("../libGGFW");

  TFile *f = TFile::Open(Form("output/AnalysisResult_xe_%s.root",deformation.Data()),"READ");

  OutputHandler* out = new OutputHandler();
  out->initialise(Form("output/AnalysisResult_xe_%s.root",deformation.Data()));
  //int samps = 10;

  //Flow-Pt Cumulants
  cout<<"--------------------------------C22pt--------------------------------"<<endl;
  TH1D* hC22pt = out->getHist(&Calculator::getCV22pt);
  // TH1D* hC22ptNorm = out->getNorm(2,2,1);


  cout<<"--------------------------------C22pt2--------------------------------"<<endl;
  TH1D* hC22pt2 = out->getHist(&Calculator::getCV22pt2);

  cout<<"--------------------------------C22pt3--------------------------------"<<endl;
  TH1D* hC22pt3 = out->getHist(&Calculator::getCV22pt3);

  cout<<"--------------------------------C22pt4--------------------------------"<<endl;
  TH1D* hC22pt4 = out->getHist(&Calculator::getCV22pt4);

  cout<<"--------------------------------C32pt--------------------------------"<<endl;
  TH1D* hC32pt = out->getHist(&Calculator::getCV32pt);

  cout<<"--------------------------------C42pt--------------------------------"<<endl;
  TH1D* hC42pt = out->getHist(&Calculator::getCV42pt);

  cout<<"--------------------------------C24pt--------------------------------"<<endl;
  TH1D* hC24pt = out->getHist(&Calculator::getCV24pt);

  cout<<"--------------------------------C24pt2--------------------------------"<<endl;
  TH1D* hC24pt2 = out->getHist(&Calculator::getCV24pt2);





  TFile *outf =  new TFile(Form("hists/HistsVnPt_%s_test.root",deformation.Data()), "RECREATE");
   // // cout<<"!!!!"<<endl;
  // ///////////////////////////////////////////////////


  TCanvas *C22pt = new TCanvas("C22pt","C22pt");
  hC22pt->SetStats(kFALSE);
  hC22pt->GetXaxis()->SetTitle("Centrality (%)");
  hC22pt->GetYaxis()->SetTitle("C(v_{2}^{2},#delta#it{p}_{T})");
  hC22pt->SetTitle("");
  hC22pt->SetAxisRange(0.,50,"X");
  //hC22pt->SetAxisRange(-2e-3,6e-3,"Y");
  hC22pt->SetLineColor(kBlack);
  hC22pt->SetMarkerStyle(20);
  hC22pt->Draw("SAME");

  TLine *line = new TLine(0,0,50,0);
  line->SetLineStyle(kDashed);
  line->Draw("SAME");

  auto label = new TLatex();
  label->SetTextAlign(23);
  label->SetTextSize(0.03);
  label->DrawLatex(15, 0.0047, "Xe-Xe, #sqrt{s_{NN}} = 5.44 TeV"); //#sqrt{s_{NN}} = 5.36 TeV
  label->DrawLatex(15, 0.00455, "0.2 < #it{p}_{T} < 3 GeV/#it{c}");
  label->DrawLatex(15, 0.0044, "|#eta| < 0.8");
  label->DrawLatex(40, 0.0029, "v_{n}: |#Delta#eta| > 0.8");
  label->DrawLatex(40, 0.00275, "[#it{p}^{(k)}]: |#eta| < 0.4");
  label->Draw("SAME"); 



  C22pt->Write();

  ///////////////////////////////////////////////////

  TCanvas *C22pt2 = new TCanvas("C22pt2","C22pt2");
  hC22pt2->SetStats(kFALSE);
  hC22pt2->GetXaxis()->SetTitle("Centrality (%)");
  hC22pt2->GetYaxis()->SetTitle("C(v_{2}^{2},#delta#it{p}_{T}^{(2)})");
  hC22pt2->SetTitle("");
  hC22pt2->SetAxisRange(0.,50,"X");
  hC22pt2->SetLineColor(kBlack);
  hC22pt2->SetMarkerStyle(20);
  hC22pt2->Draw();

  // auto legend4pc = new TLegend();
  // legend4pc->SetBorderSize(0);
  // legend4pc->AddEntry(hC22pt2,"C(v_{2}^{2},#delta#it{p}_{T}^{(2)})","ep");
  // legend4pc->Draw("SAME");

  auto label2 = new TLatex();
  label2->SetTextAlign(23);
  label2->SetTextSize(0.03);
  label2->DrawLatex(15, -1.16, "Xe-Xe, #sqrt{s_{NN}} = 5.44 TeV");
  label2->DrawLatex(15, -1.17, "0.2 < #it{p}_{T} < 3 GeV/#it{c}");
  label2->DrawLatex(15, -1.18, "|#eta| < 0.8");
  label2->DrawLatex(45, -1.06, "v_{n}: |#Delta#eta| > 0.8");
  label2->DrawLatex(45, -1.07, "[#it{p}^{(k)}]: |#eta| < 0.4");label2->Draw("SAME"); 

  C22pt2->Write();

  /////////////////////////////////////////////////

  TCanvas *C22pt3 = new TCanvas("C22pt3","C22pt3");
  hC22pt3->SetStats(kFALSE);
  hC22pt3->GetXaxis()->SetTitle("Centrality (%)");
  hC22pt3->GetYaxis()->SetTitle("C(v_{2}^{2},#delta#it{p}_{T}^{(3)})");
  hC22pt3->SetTitle("");
  hC22pt3->SetAxisRange(0.,50,"X");
  hC22pt3->SetLineColor(kBlack);
  hC22pt3->SetMarkerStyle(20);
  hC22pt3->Draw();
  // auto legend5pc = new TLegend();
  // legend5pc->SetBorderSize(0);
  // legend5pc->AddEntry(hC22pt3,"C(v_{2}^{2},#delta#it{p}_{T}^{(3)})","ep");
  // legend5pc->Draw("SAME");
  auto label3 = new TLatex();
  label3->SetTextAlign(23);
  label3->SetTextSize(0.03);
  label3->DrawLatex(15, 4.5, "Xe-Xe, #sqrt{s_{NN}} = 5.44 TeV");
  label3->DrawLatex(15, 4.45, "0.2 < #it{p}_{T} < 3 GeV/#it{c}");
  label3->DrawLatex(15, 4.4, "|#eta| < 0.8");
  label3->DrawLatex(45, 4.05, "v_{n}: |#Delta#eta| > 0.8");
  label3->DrawLatex(45, 4, "[#it{p}^{(k)}]: |#eta| < 0.4");label3->Draw("SAME"); 
  C22pt3->Write();

  // ///////////////////////////////////////////////////

  TCanvas *C22pt4 = new TCanvas("C22pt4","C22pt4");
  hC22pt4->SetStats(kFALSE);
  hC22pt4->GetXaxis()->SetTitle("Centrality (%)");
  hC22pt4->GetYaxis()->SetTitle("C(v_{2}^{2},#delta#it{p}_{T}^{4})");
  hC22pt4->SetTitle("");
  hC22pt4->SetAxisRange(0.,50,"X");
  hC22pt4->SetLineColor(kBlack);
  hC22pt4->SetMarkerStyle(20);
  hC22pt4->Draw();
  // auto legend6pc = new TLegend();
  // legend6pc->SetBorderSize(0);
  // legend6pc->AddEntry(hC22pt4,"C(v_{2}^{2},#delta#it{p}_{T}^{(4)})","ep");
  // legend6pc->Draw("SAME");
  auto label4 = new TLatex();
  label4->SetTextAlign(23);
  label4->SetTextSize(0.03);
  label4->DrawLatex(15, -4.4, "Xe-Xe, #sqrt{s_{NN}} = 5.44 TeV");
  label4->DrawLatex(15, -4.5, "0.2 < #it{p}_{T} < 3 GeV/#it{c}");
  label4->DrawLatex(15, -4.6, "|#eta| < 0.8");
  label4->DrawLatex(45, -3.8, "v_{n}: |#Delta#eta| > 0.8");
  label4->DrawLatex(45, -3.7, "[#it{p}^{(k)}]: |#eta| < 0.4");label4->Draw("SAME"); 
  C22pt4->Write();

  ///////////////////////////////////////////////////

  TCanvas *C24pt = new TCanvas("C24pt","C24pt");

  hC24pt->SetStats(kFALSE);
  hC24pt->GetXaxis()->SetTitle("Centrality (%)");
  hC24pt->GetYaxis()->SetTitle("C(v_{2}^{4},#delta#it{p}_{T})");
  hC24pt->SetTitle("");
  hC24pt->SetAxisRange(0.,50,"X");
  hC24pt->SetLineColor(kBlack);
  hC24pt->SetMarkerStyle(20);
  hC24pt->Draw();

  TLine *line7 = new TLine(0,0,50,0);
  line7->SetLineStyle(kDashed);
  line7->Draw("SAME");

  // auto legend5pc2 = new TLegend();
  // legend5pc2->SetBorderSize(0);
  // legend5pc2->AddEntry(hC24pt,"C(v_{2}^{4},#delta#it{p}_{T})","ep");
  // legend5pc2->Draw("SAME");

  auto label7 = new TLatex();
  label7->SetTextAlign(23);
  label7->SetTextSize(0.03);
  label7->DrawLatex(40, 0.004, "Xe-Xe, #sqrt{s_{NN}} = 5.44 TeV");
  label7->DrawLatex(40, 0.0033, "0.2 < #it{p}_{T} < 3 GeV/#it{c}");
  label7->DrawLatex(40, 0.0026, "|#eta| < 0.8");
  label7->DrawLatex(10, -0.005, "v_{n}: |#Delta#eta| > 0.8");
  label7->DrawLatex(10, -0.0057, "[#it{p}^{(k)}]: |#eta| < 0.4");
  label7->Draw("SAME"); 

  C24pt->Write();


  delete outf;

}