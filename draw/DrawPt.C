#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>
#include <iostream>
#include "../BootstrapProfile.h"

void DrawPt(TString deformation = "") {

    gSystem->Load("../libGGFW");
    TFile *f = TFile::Open(Form("output/AnalysisResult_xe_%s.root",deformation.Data()),"READ");

    //Pt
    OutputHandler* out = new OutputHandler();
    out->initialise(Form("output/AnalysisResult_xe_%s.root",deformation.Data()));


    TH1D* hMpt1 = out->getHist(&Calculator::getMpt);
    TH1D* hMpt2 = out->getHist(&Calculator::getMpt2);
    TH1D* hMpt3 = out->getHist(&Calculator::getMpt3);
    TH1D* hMpt4 = out->getHist(&Calculator::getMpt4);
    TH1D* hMpt5 = out->getHist(&Calculator::getMpt5);
    TH1D* hMpt6 = out->getHist(&Calculator::getMpt6);
    TH1D* hMpt7 = out->getHist(&Calculator::getMpt7);
    TH1D* hMpt8 = out->getHist(&Calculator::getMpt8);


    TFile *outf =  new TFile(Form("hists/HistsPt_%s.root", deformation.Data()),"RECREATE");

    ///////////////////////////////////////////////////

    TCanvas *Mpt1 = new TCanvas("Mpt1","Mpt1");
    hMpt1->SetStats(kFALSE);
    hMpt1->SetTitle("Mpt1 from calculator");
    //hC22pt->SetAxisRange(0.,50,"X");
    hMpt1->Draw();
    hMpt1->Write();

    ///////////////////////////////////////////////////

    TCanvas *Mpt2 = new TCanvas("Mpt2","Mpt2");
    hMpt2->SetStats(kFALSE);
    hMpt2->SetTitle("Mpt2 from calculator");
    //hC22pt->SetAxisRange(0.,50,"X");
    hMpt2->Draw();
    hMpt2->Write();

    ///////////////////////////////////////////////////

    TCanvas *Mpt3 = new TCanvas("Mpt3","Mpt3");
    hMpt3->SetStats(kFALSE);
    hMpt3->SetTitle("Mpt3 from calculator");
    //hC22pt->SetAxisRange(0.,50,"X");
    hMpt3->Draw();
    hMpt3->Write();

    ///////////////////////////////////////////////////

    TCanvas *Mpt4 = new TCanvas("Mpt4","Mpt4");
    hMpt4->SetStats(kFALSE);
    hMpt4->SetTitle("Mpt4 from calculator");
    //hC22pt->SetAxisRange(0.,50,"X");
    hMpt4->Draw();
    hMpt4->Write();

    ///////////////////////////////////////////////////

    TCanvas *Mpt5 = new TCanvas("Mpt5","Mpt5");
    hMpt5->SetStats(kFALSE);
    hMpt5->SetTitle("Mpt5 from calculator");
    //hC22pt->SetAxisRange(0.,50,"X");
    hMpt5->Draw();
    hMpt5->Write();

    ///////////////////////////////////////////////////

    TCanvas *Mpt6 = new TCanvas("Mpt6","Mpt6");
    hMpt6->SetStats(kFALSE);
    hMpt6->SetTitle("Mpt6 from calculator");
    //hC22pt->SetAxisRange(0.,50,"X");
    hMpt6->Draw();
    hMpt6->Write();

    ///////////////////////////////////////////////////

    TCanvas *Mpt7 = new TCanvas("Mpt7","Mpt7");
    hMpt7->SetStats(kFALSE);
    hMpt7->SetTitle("Mpt7 from calculator");
    //hC22pt->SetAxisRange(0.,50,"X");
    hMpt7->Draw();
    hMpt7->Write();

    ///////////////////////////////////////////////////

    TCanvas *Mpt8 = new TCanvas("Mpt8","Mpt8");
    hMpt8->SetStats(kFALSE);
    hMpt8->SetTitle("Mpt8 from calculator");
    //hC22pt->SetAxisRange(0.,50,"X");
    hMpt8->Draw();
    hMpt8->Write();

    ///////////////////////////////////////////////////

    TCanvas *MptAll = new TCanvas("MptAll","MptAll");

    hMpt1->GetYaxis()->SetTitle("[p_{T}^{(k)}]");
    hMpt1->GetXaxis()->SetTitle("Centrality (%)");
    hMpt1->SetStats(kFALSE);
    hMpt1->SetLineColor(kBlack);
    hMpt1->SetMarkerColor(kBlack);
    hMpt1->SetMarkerStyle(20);
    hMpt1->SetAxisRange(0.,1.1,"Y");
    hMpt1->SetAxisRange(0., 50.,"X");
    hMpt1->Draw();

    hMpt2->SetStats(kFALSE);
    hMpt2->SetLineColor(kRed);
    hMpt2->SetMarkerColor(kRed);
    hMpt2->SetMarkerStyle(21);
    hMpt2->Draw("SAME");

    hMpt3->SetStats(kFALSE);
    hMpt3->SetLineColor(kOrange);
    hMpt3->SetMarkerColor(kOrange);
    hMpt3->SetMarkerStyle(33);
    hMpt3->Draw("SAME");

    hMpt4->SetStats(kFALSE);
    hMpt4->SetLineColor(kGreen);
    hMpt4->SetMarkerColor(kGreen);
    hMpt4->SetMarkerStyle(34);
    hMpt4->Draw("SAME");

    hMpt5->SetStats(kFALSE);
    hMpt5->SetLineColor(kAzure+8);
    hMpt5->SetMarkerColor(kAzure+8);
    hMpt5->SetMarkerStyle(43);
    hMpt5->Draw("SAME");

    hMpt6->SetStats(kFALSE);
    hMpt6->SetLineColor(kBlue);
    hMpt6->SetMarkerColor(kBlue);
    hMpt6->SetMarkerStyle(47);
    hMpt6->Draw("SAME");

    hMpt7->SetStats(kFALSE);
    hMpt7->SetLineColor(kViolet+2);
    hMpt7->SetMarkerColor(kViolet+2);
    hMpt7->SetMarkerStyle(29);
    hMpt7->Draw("SAME");

    hMpt8->SetLineColor(kViolet);
    hMpt8->SetMarkerColor(kViolet);
    hMpt8->SetStats(kFALSE);
    hMpt8->SetMarkerStyle(22);
    hMpt8->Draw("SAME");

    auto legend = new TLegend();
    legend->SetNColumns(4);
    legend->AddEntry(hMpt1,"[p_{T}^{(1)}]","ep");
    legend->AddEntry(hMpt2,"[p_{T}^{(2)}]","ep");
    legend->AddEntry(hMpt3,"[p_{T}^{(3)}]","ep");
    legend->AddEntry(hMpt4,"[p_{T}^{(4)}]","ep");
    legend->AddEntry(hMpt5,"[p_{T}^{(5)}]","ep");
    legend->AddEntry(hMpt6,"[p_{T}^{(6)}]","ep");
    legend->AddEntry(hMpt7,"[p_{T}^{(7)}]","ep");
    legend->AddEntry(hMpt8,"[p_{T}^{(8)}]","ep");
    legend->Draw("SAME");

    MptAll->Write();


    delete outf;

}