#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include "TString.h"

#include <cmath>
#include <iostream>
#include "../BootstrapProfile.h"
#include "../Calculator.h"

void DrawVn(TString deformation = "") {
    gSystem->Load("../libGGFW");
    TFile *f = TFile::Open(Form("output/AnalysisResult_xe_%s.root",deformation.Data()),"READ");

    //Flow 
    FlowContainer* fc = (FlowContainer*)f->Get("FlowContainer");
    fc->SetIDName("ChGap");
    fc->SetPropagateErrors(kTRUE);

    TH1D* hV22 = fc->GetVN2VsMulti();
    TH1D* hV32 = fc->GetVN2VsMulti(3);
    TH1D* hV42 = fc->GetVN2VsMulti(4);
    TH1D* hV24 = fc->GetVN4VsMulti();



    TFile *outf =  new TFile(Form("hists/HistsVn_%s.root", deformation.Data()),"RECREATE");
    // // cout<<"!!!!"<<endl;
    // // ///////////////////////////////////////////////////

    TCanvas *V22 = new TCanvas("V22","V22");
    hV22->SetStats(kFALSE);
    hV22->SetTitle("V22 from FlowContainer");
    //hC22pt->SetAxisRange(0.,50,"X");
    hV22->Draw();
    hV22->Write();

    // // ///////////////////////////////////////////////////

    TCanvas *Vall = new TCanvas("Vall","Vall");


    hV22->SetAxisRange(0.,60,"X");
    hV22->SetAxisRange(0.,0.14,"Y");
    hV22->GetYaxis()->SetTitle("v_{n} {m}");
    hV22->GetXaxis()->SetTitle("Centrality (%)");

    hV22->SetStats(kFALSE);
    hV22->SetLineColor(kBlack);
    hV22->SetMarkerStyle(20);
    hV22->Draw();

    hV32->SetStats(kFALSE);
    hV32->SetLineColor(kGreen+1);
    hV32->SetMarkerColor(kGreen+1);
    hV32->SetMarkerStyle(21);
    hV32->Draw("SAME");

    hV42->SetStats(kFALSE);
    hV42->SetLineColor(kBlue-4);
    hV42->SetMarkerColor(kBlue-4);
    hV42->SetMarkerStyle(34);
    hV42->Draw("SAME");

    hV24->SetStats(kFALSE);
    hV24->SetLineColor(kRed);
    hV24->SetMarkerColor(kRed);
    hV24->SetMarkerStyle(24);
    hV24->Draw("SAME");


    auto legend = new TLegend();
    legend->AddEntry(hV22,"v_{2} {2,|#Delta#eta|>0.8}","ep");
    legend->AddEntry(hV32,"v_{3} {2,|#Delta#eta|>0.8}","ep");
    legend->AddEntry(hV42,"v_{4} {2,|#Delta#eta|>0.8}","ep");
    legend->AddEntry(hV24,"v_{2} {4,|#Delta#eta|>0.8}","ep");
    legend->Draw("SAME");

    // auto label = new TLatex();
    // label->DrawLatex(0.0, 0.0, "Xe-Xe");
    // label->DrawLatex(0.0, 0.0, "");
    // label->Draw("SAME"); 

    Vall->Write();

    delete outf;

}