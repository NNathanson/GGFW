#include "Calculator.h"

namespace CalculatorFunctions {
  Calculator::ObsFunc k22 = &Calculator::get22;
  Calculator::ObsFunc k32 = &Calculator::get32;
  Calculator::ObsFunc k42 = &Calculator::get42;
  Calculator::ObsFunc k24 = &Calculator::get24;

  Calculator::ObsFunc kV22 = &Calculator::getV22;
  Calculator::ObsFunc kV32 = &Calculator::getV32;
  Calculator::ObsFunc kV42 = &Calculator::getV42;
  Calculator::ObsFunc kV24 = &Calculator::getV24;
  Calculator::ObsFunc kV26 = &Calculator::getV26;
  Calculator::ObsFunc kV28 = &Calculator::getV28;


  Calculator::PtDiffFunc kV22pt = &Calculator::getV22pt;
  Calculator::PtDiffFunc kV32pt = &Calculator::getV32pt;
  Calculator::PtDiffFunc kV42pt = &Calculator::getV42pt;
  Calculator::PtDiffFunc kV24pt = &Calculator::getV24pt;

  Calculator::ObsFunc kC22 = &Calculator::getC22;
  Calculator::ObsFunc kC32 = &Calculator::getC32;
  Calculator::ObsFunc kC42 = &Calculator::getC42;
  Calculator::ObsFunc kC24 = &Calculator::getC24;
  Calculator::ObsFunc kC26 = &Calculator::getC26;
  Calculator::ObsFunc kC28 = &Calculator::getC28;

  Calculator::ObsFunc kCV22pt = &Calculator::getCV22pt;
  Calculator::ObsFunc kCV32pt = &Calculator::getCV32pt;
  Calculator::ObsFunc kCV42pt = &Calculator::getCV42pt;
  Calculator::ObsFunc kCV22pt2 = &Calculator::getCV22pt2;
  Calculator::ObsFunc kCV22pt3 = &Calculator::getCV22pt3;
  Calculator::ObsFunc kCV22pt4 = &Calculator::getCV22pt4;
  Calculator::ObsFunc kCV24pt = &Calculator::getCV24pt;
  Calculator::ObsFunc kCV24pt2 = &Calculator::getCV24pt2;

  Calculator::ObsFunc kV22pt_Mpt0 = &Calculator::getV22pt_Mpt0;
  Calculator::ObsFunc kV22pt_Mpt1 = &Calculator::getV22pt_Mpt1;
  Calculator::ObsFunc kV24pt2_Mpt0 = &Calculator::getV24pt2_Mpt0;
  Calculator::ObsFunc kV24pt2_Mpt1 = &Calculator::getV24pt2_Mpt1;
  Calculator::ObsFunc kV24pt2_Mpt2 = &Calculator::getV24pt2_Mpt2;

  Calculator::ObsFunc kMpt = &Calculator::getMpt;
  Calculator::ObsFunc kMpt2 = &Calculator::getMpt2;
  Calculator::ObsFunc kMpt3 = &Calculator::getMpt3;
  Calculator::ObsFunc kMpt4 = &Calculator::getMpt4;
  Calculator::ObsFunc kK2 = &Calculator::getK2;
  Calculator::ObsFunc kK3 = &Calculator::getK3;
  Calculator::ObsFunc kK4 = &Calculator::getK4;
  Calculator::ObsFunc kMu2 = &Calculator::getMu2;
  Calculator::ObsFunc kMu3 = &Calculator::getMu3;
  Calculator::ObsFunc kMu4 = &Calculator::getMu4;
}
Calculator::Calculator() : flowcont{0}, flowptcont{0}, idname(""), current_suffix("")
{
}
Calculator::Calculator(TFile* infile) : flowcont{0}, flowptcont{0}, idname(""), current_suffix("")
{
  readData(infile);
}
Calculator::~Calculator()
{
  delete flowcont;
  delete flowptcont;
};

void Calculator::readData(TFile* infile)
{
  // Get flow container
  flowcont = dynamic_cast<FlowContainer*>(infile->Get(Form("FlowContainer")));
  if (!flowcont) {
    std::cout << "Warning: could not find FlowContainer in " << "generic-framework" << current_suffix.c_str() << std::endl;
  }
  // Get pt container
  flowptcont = dynamic_cast<FlowPtContainer*>(infile->Get(Form("FlowPtContainer")));
  if (!flowptcont) {
    std::cout << "Warning: could not find FlowPtContainer in " << "generic-framework" << current_suffix.c_str() <<  std::endl;
  }
  else {
    // Get the correlator profiles
    corr_list = flowptcont->getCovList();
    std::cout<<"Uses gap? "<<flowptcont->usesGap()<<std::endl;
    idname = (flowptcont->usesGap())?"Gap":"Full";
    std::cout<<"idname? "<<idname<<std::endl;
    idname = "Gap";

    std::cout<<"test: mpar? "<<flowptcont->getMpar()<<std::endl;
   // flowptcont->setMpar(8);
    //flowptcont->setUseCentralMoments(kFALSE);
    std::cout<<"test: cm? "<<flowptcont->usesCentralMoments()<<std::endl;
    
  }

  return;
}
TH1D* Calculator::calculate(ObsFunc func, int i)
{
  return (this->*func)(i);
}
TH1D* Calculator::calculate(PtDiffFunc func, int i, double min, double max)
{
  return (this->*func)(i, min, max);
}


TH1D* Calculator::getNormNC(int n, int m, int k, int i){
    auto denom = getNN(n, m, Form("Ch%s",idname.c_str()), i);
    denom->Multiply(flowptcont->getCentralMomentHist(i, k));
    return denom;
}

TH1D* Calculator::getV22pt_Mpt0(int i){
    auto vn2pt_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt1_Mpt0", idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(vn2pt_Mpt0_profile->getHist(i));
    reth->SetNameTitle("V22pt_mpt0","#LT v_{2}^{2}[#it{p}_{T}]#GT");
    return reth;
}
TH1D* Calculator::getV22pt_Mpt1(int i){
    auto vn2pt_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt1_Mpt1", idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(vn2pt_Mpt1_profile->getHist(i));
    reth->SetNameTitle("V22pt_mpt1","#LT v_{2}^{2}#GT");
    return reth;
}

TH1D* Calculator::getV32pt_Mpt0(int i){
    auto vn2pt_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s32pt1_Mpt0", idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(vn2pt_Mpt0_profile->getHist(i));
    reth->SetNameTitle("V32pt_mpt0","#LT v_{3}^{2}[#it{p}_{T}]#GT");
    return reth;
}

TH1D* Calculator::getV32pt_Mpt1(int i){
    auto vn2pt_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s32pt1_Mpt1", idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(vn2pt_Mpt1_profile->getHist(i));
    reth->SetNameTitle("V32pt_mpt1","#LT v_{3}^{2}#GT");
    return reth;
}

TH1D* Calculator::getV42pt_Mpt0(int i){
    auto vn2pt_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s42pt1_Mpt0", idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(vn2pt_Mpt0_profile->getHist(i));
    reth->SetNameTitle("V42pt_mpt0","#LT v_{4}^{2}[#it{p}_{T}]#GT");
    return reth;
}

TH1D* Calculator::getV42pt_Mpt1(int i){
    auto vn2pt_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s42pt1_Mpt1", idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(vn2pt_Mpt1_profile->getHist(i));
    reth->SetNameTitle("V42pt_mpt1","#LT v_{4}^{2}#GT");
    return reth;
}

TH1D* Calculator::getV24pt2_Mpt0(int i){
    auto v24pt2_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt2_Mpt0", idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(v24pt2_Mpt0_profile->getHist(i));
    reth->SetNameTitle("V24pt2_mpt0","#LT v_{2}^{4}[#it{p}_{T}^{(2)}]#GT");
    return reth;
}
TH1D* Calculator::getV24pt2_Mpt1(int i){
    auto v24pt2_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt2_Mpt1", idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(v24pt2_Mpt1_profile->getHist(i));
    reth->SetNameTitle("V24pt2_mpt1","#LT v_{2}^{4}[#it{p}_{T}]#GT");
    return reth;
}
TH1D* Calculator::getV24pt2_Mpt2(int i){
    auto v24pt2_Mpt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt2_Mpt2", idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(v24pt2_Mpt2_profile->getHist(i));
    reth->SetNameTitle("V24pt2_mpt2","#LT v_{2}^{4}#GT");
    return reth;
}
TH1D* Calculator::getCVN2pt(int nhar, int i)
{
  if (flowptcont->usesCentralMoments()) {
    auto vn2pt_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s%i2pt1_Mpt0", idname.c_str(), nhar)));
    if(!vn2pt_Mpt0_profile) std::cout << Form("Ch%s%i2pt1_Mpt0", idname.c_str(), nhar) << " not found in list" << std::endl;
    auto vn2pt_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s%i2pt1_Mpt1", idname.c_str(), nhar)));
    if(!vn2pt_Mpt1_profile) std::cout << Form("Ch%s%i2pt1_Mpt1", idname.c_str(), nhar) << " not found in list" << std::endl;
    TH1D* reth = dynamic_cast<TH1D*>(vn2pt_Mpt0_profile->getHist(i));
    TH1D* v22_pt = dynamic_cast<TH1D*>(vn2pt_Mpt1_profile->getHist(i));
    auto mpt_hist(getMpt(i));
    v22_pt->Multiply(mpt_hist);
    reth->Add(v22_pt, -1);
    reth->SetNameTitle(Form("CV%i2dpt",nhar),Form("C(v_{%i}^{2},#delta#it{p}_{T})",nhar));
    return reth;
  }
  auto vn2pt_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s%i2pt1", idname.c_str(), nhar)));
  if(!vn2pt_profile) std::cout << Form("Ch%s%i2pt1", idname.c_str(), nhar) << std::endl;
  TH1D* reth = dynamic_cast<TH1D*>(vn2pt_profile->getHist(i));
  auto mpt = getMpt(i);
  auto vn2_mpt = getNN(nhar, 2, Form("Ch%s",idname.c_str()), i);
  vn2_mpt->Multiply(mpt);
  reth->Add(vn2_mpt,-1);
  reth->SetNameTitle(Form("CV%i2pt",nhar),Form("C(v_{%i}^{2},[#it{p}_{T}])",nhar));
  return reth;
}
TH1D* Calculator::getCV22pt2(int i)
{
  if (flowptcont->usesCentralMoments()) {
    auto v22pt2_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt2_Mpt0",idname.c_str())));
    // First get the <v2^2*deltapt^2> term
    TH1D* reth = dynamic_cast<TH1D*>(v22pt2_Mpt0_profile->getHist(i));
    auto mpt_hist = getMpt(i);
    auto v22pt2_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt2_Mpt1",idname.c_str())));
    auto v22pt_pt = dynamic_cast<TH1D*>(v22pt2_Mpt1_profile->getHist(i));
    v22pt_pt->Multiply(mpt_hist);
    reth->Add(v22pt_pt);
    auto v22pt2_Mpt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt2_Mpt2",idname.c_str())));
    auto v22_pt2 = dynamic_cast<TH1D*>(v22pt2_Mpt2_profile->getHist(i));
    auto mpt_sq_hist = dynamic_cast<TH1D*>(mpt_hist->Clone("mpt_sq_hist"));
    mpt_sq_hist->Multiply(mpt_hist);
    v22_pt2->Multiply(mpt_sq_hist);
    reth->Add(v22_pt2);
    // Then subtract <v2^2><deltapt^2>
    auto v22_cm2 = getNN(2, 2, Form("Ch%s",idname.c_str()), i);
    v22_cm2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22_cm2, -1);

    reth->SetNameTitle("CV22dpt2","C(v_{2}^{2},#delta#it{p}_{T}^{2})");
    return reth;
  }
  auto v22pt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt2",idname.c_str())));
  TH1D* reth = dynamic_cast<TH1D*>(v22pt2_profile->getHist(i));

  //subtract 2 x <v2^2*[pt]>*<[pt]>
  auto v22pt_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt1",idname.c_str())));
  auto v22pt_mpt = dynamic_cast<TH1D*>(v22pt_profile->getHist(i));
  auto mpt = getMpt(i);
  v22pt_mpt->Multiply(mpt);
  reth->Add(v22pt_mpt,-2);

  //Add 2 x <v2^2>*<[pt]>^2
  auto v22_mpt_mpt = getNN(2,2,Form("Ch%s",idname.c_str()),i);
  v22_mpt_mpt->Multiply(mpt);
  v22_mpt_mpt->Multiply(mpt);
  reth->Add(v22_mpt_mpt,2);

  //Subtract <v2^2>*<[pt^2]>
  auto v22_mpt2 = getNN(2,2,Form("Ch%s",idname.c_str()),i);
  v22_mpt2->Multiply(flowptcont->getCorrHist(i,2));
  reth->Add(v22_mpt2,-1);

  reth->SetNameTitle("CV22dpt2","C(v_{2}^{2},#delta#it{p}_{T}^{(2)})");
  return reth;
}
TH1D* Calculator::getCV22pt3(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^2*deltapt^3> term
    auto v22pt3_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt3_Mpt0",idname.c_str())));
    auto v22pt3_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt3_Mpt1",idname.c_str())));
    auto v22pt3_Mpt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt3_Mpt2",idname.c_str())));
    auto v22pt3_Mpt3_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt3_Mpt3",idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(v22pt3_Mpt0_profile->getHist(i));
    auto mpt_hist = getMpt(i);
    auto v22pt2_pt = dynamic_cast<TH1D*>(v22pt3_Mpt1_profile->getHist(i));
    v22pt2_pt->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22pt_pt2 = dynamic_cast<TH1D*>(v22pt3_Mpt2_profile->getHist(i));
    v22pt_pt2->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22_pt3 = dynamic_cast<TH1D*>(v22pt3_Mpt3_profile->getHist(i));
    v22_pt3->Multiply(mpt_hist);
    reth->Add(v22pt2_pt);
    reth->Add(v22pt_pt2);
    reth->Add(v22_pt3);
    // Then subtract the 3 x <v2^2*deltapt>*<deltapt^2> term
    auto v22deltapt_deltapt2 = getCV22pt(i);
    v22deltapt_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22deltapt_deltapt2, -3);
    // And then subtract the <v2^2>*<deltapt^3> term
    auto v22_deltapt3 = getNN(2, 2, Form("Ch%s",idname.c_str()), i);
    v22_deltapt3->Multiply(flowptcont->getCentralMomentHist(i, 3));
    reth->Add(v22_deltapt3, -1);
    reth->SetNameTitle("CV22dpt3","C(v_{2}^{2},#delta#it{p}_{T}^{(3)})");
    return reth;
  }
  return nullptr;
}
TH1D* Calculator::getCV22pt4(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^2*deltapt^4> term
    auto v22pt4_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt4_Mpt0",idname.c_str())));
    auto v22pt4_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt4_Mpt1",idname.c_str())));
    auto v22pt4_Mpt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt4_Mpt2",idname.c_str())));
    auto v22pt4_Mpt3_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt4_Mpt3",idname.c_str())));
    auto v22pt4_Mpt4_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt4_Mpt4",idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(v22pt4_Mpt0_profile->getHist(i));
    auto mpt_hist = getMpt(i);
    auto v22pt3_pt = dynamic_cast<TH1D*>(v22pt4_Mpt1_profile->getHist(i));
    v22pt3_pt->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22pt2_pt2 = dynamic_cast<TH1D*>(v22pt4_Mpt2_profile->getHist(i));
    v22pt2_pt2->Multiply(mpt_hist);
    auto v22pt_pt3 = dynamic_cast<TH1D*>(v22pt4_Mpt3_profile->getHist(i));
    v22pt_pt3->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22_pt4 = dynamic_cast<TH1D*>(v22pt4_Mpt4_profile->getHist(i));
    v22_pt4->Multiply(mpt_hist);
    reth->Add(v22pt3_pt);
    reth->Add(v22pt2_pt2);
    reth->Add(v22pt_pt3);
    reth->Add(v22_pt4);
    // Then subtract the 6 x <v2^2*deltapt2>*<deltapt^2> term
    auto v22deltapt2_deltapt2 = getCV22pt2(i);
    v22deltapt2_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22deltapt2_deltapt2, -6);
    // Then subtract the 4 x <v2^2*deltapt>*<deltapt^3> term
    auto v22deltapt_deltapt3 = getCV22pt(i);
    v22deltapt_deltapt3->Multiply(flowptcont->getCentralMomentHist(i, 3));
    reth->Add(v22deltapt_deltapt3, -4);
    // Then add the 6 x <v2^2>*<deltapt^2>*<deltapt^2> term
    auto v22_deltapt2_deltapt2 = getNN(2, 2, Form("Ch%s",idname.c_str()), i);
    v22_deltapt2_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    v22_deltapt2_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22_deltapt2_deltapt2, -4);
    // Then subtract the <v2^2>*<deltapt^4> term
    auto v22_deltapt4 = getNN(2, 2, Form("Ch%s",idname.c_str()), i);
    v22_deltapt4->Multiply(flowptcont->getCentralMomentHist(i, 4));
    reth->Add(v22_deltapt4, 1);
    reth->SetNameTitle("CV22dpt4","C(v_{2}^{2},#delta#it{p}_{T}^{(4)})");
    return reth;
  }
  return nullptr;
}
TH1D* Calculator::getCV24pt(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^4*deltapt> term
    BootstrapProfile* v24pt_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt1_Mpt0",idname.c_str())));
    BootstrapProfile* v24pt_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt1_Mpt1",idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(v24pt_Mpt0_profile->getHist(i));
    auto v24_pt = dynamic_cast<TH1D*>(v24pt_Mpt1_profile->getHist(i));
    auto mpt_hist = getMpt(i);
    v24_pt->Multiply(mpt_hist);
    reth->Add(v24_pt, -1);
    // Then subtract the 4 x <v2^2*deltapt>*<v2^2> term
    auto v22pt_v22 = getCV22pt(i);
    v22pt_v22->Multiply(getNN(2, 2, Form("Ch%s",idname.c_str()), i));
    reth->Add(v22pt_v22, -4);
    reth->SetNameTitle("CV24dpt","C(v_{2}^{4},#delta#it{p}_{T})");
    return reth;
  }
  //Get the <v2^4*[pt]> term
  BootstrapProfile* v24pt_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt1",idname.c_str())));
  TH1D* reth = dynamic_cast<TH1D*>(v24pt_profile->getHist(i));
  //Subtract the <v24><[pt]> term
  auto v24_mpt = getNN(2, 4, Form("Ch%s",idname.c_str()), i);
  auto mpt = getMpt(i);
  v24_mpt->Multiply(mpt);
  reth->Add(v24_mpt, -1);
  //Subtract the 4 x <v2^2><v2^2*[pt]> term
  BootstrapProfile* v22pt_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt1",idname.c_str())));
  TH1D* v22pt_v22 = dynamic_cast<TH1D*>(v22pt_profile->getHist(i));
  v22pt_v22->Multiply(getNN(2, 2, Form("Ch%s",idname.c_str()), i));
  reth->Add(v22pt_v22,-4);
  //Add the 4 x <v2^2>*<v2^2>*<[pt]> term
  auto v22_v22_mpt = getNN(2, 2, Form("Ch%s",idname.c_str()), i);
  v22_v22_mpt->Multiply(v22_v22_mpt);
  v22_v22_mpt->Multiply(mpt);
  reth->Add(v22_v22_mpt,4);

  reth->SetNameTitle("CV24pt","C(v_{2}^{4},[#it{p}_{T}])");
  return reth;
}
TH1D* Calculator::getCV24pt2(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^4*deltapt^2> term
    BootstrapProfile* v24pt2_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt2_Mpt0",idname.c_str())));
    BootstrapProfile* v24pt2_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt2_Mpt1",idname.c_str())));
    BootstrapProfile* v24pt2_Mpt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt2_Mpt2",idname.c_str())));
    TH1D* reth = dynamic_cast<TH1D*>(v24pt2_Mpt0_profile->getHist(i));
    auto v24pt_pt = dynamic_cast<TH1D*>(v24pt2_Mpt1_profile->getHist(i));
    auto mpt = getMpt(i);
    v24pt_pt->Multiply(mpt);
    auto v24_pt2 = dynamic_cast<TH1D*>(v24pt2_Mpt2_profile->getHist(i));
    mpt->Multiply(mpt);
    v24_pt2->Multiply(mpt);
    reth->Add(v24pt_pt);
    reth->Add(v24_pt2);

    // Then subtract the 4 x <v2^2>*<vn^2*deltapt^2> term
    BootstrapProfile* v22pt2_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt2_Mpt0",idname.c_str())));
    BootstrapProfile* v22pt2_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt2_Mpt1",idname.c_str())));
    BootstrapProfile* v22pt2_Mpt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt2_Mpt2",idname.c_str())));
    auto v22deltapt2_v22 = dynamic_cast<TH1D*>(v22pt2_Mpt0_profile->getHist(i));
    auto v22pt_mpt = dynamic_cast<TH1D*>(v22pt2_Mpt1_profile->getHist(i));
    v22pt_mpt->Multiply(mpt);
    auto v22_mpt_mpt = dynamic_cast<TH1D*>(v22pt2_Mpt2_profile->getHist(i));
    v22_mpt_mpt->Multiply(mpt);
    v22_mpt_mpt->Multiply(mpt);
    v22deltapt2_v22->Add(v22pt_mpt,-2);
    v22deltapt2_v22->Add(v22_mpt_mpt);
    v22deltapt2_v22->Multiply(getNN(2, 2, Form("Ch%s",idname.c_str()), i));
    reth->Add(v22deltapt2_v22,-4);

    // Then subtract the 4 x <v2^2*deltapt>*<v2^2*deltapt> term
    BootstrapProfile* v22pt_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt1_Mpt0",idname.c_str())));
    BootstrapProfile* v22pt_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt1_Mpt1",idname.c_str())));
    auto v22deltapt_v22deltapt = dynamic_cast<TH1D*>(v22pt_Mpt0_profile->getHist(i));
    auto v22_mpt = dynamic_cast<TH1D*>(v22pt_Mpt1_profile->getHist(i));
    v22_mpt->Multiply(mpt);
    v22deltapt_v22deltapt->Add(v22_mpt,-1);
    v22deltapt_v22deltapt->Multiply(v22deltapt_v22deltapt);
    reth->Add(v22deltapt_v22deltapt,-4);

    // Then add the 4 x <v2^2>*<v2^2>*<deltapt^2> term
    auto v22_v22_deltapt2 = getNN(2, 2, Form("Ch%s",idname.c_str()), i);
    v22_v22_deltapt2->Multiply(v22_v22_deltapt2);
    v22_v22_deltapt2->Multiply(flowptcont->getCentralMomentHist(i,2));
    reth->Add(v22_v22_deltapt2,4);

    // Then subtract the <v2^4>*<deltapt^2> term
    auto v24_deltapt2 = getNN(2,2,Form("Ch%s",idname.c_str()),i);
    v24_deltapt2->Multiply(flowptcont->getCentralMomentHist(i,2));
    reth->Add(v24_deltapt2,-1);

    reth->SetNameTitle("CV24dpt2","C(v_{2}^{4},#delta#it{p}_{T}^{2})");
    return reth;
  }
  // Get the <v2^4*[pt]> term
  BootstrapProfile* v24pt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt2",idname.c_str())));
  auto reth = dynamic_cast<TH1D*>(v24pt2_profile->getHist(i));

  // Subtract -2 x <v2^4*[pt]>*<[pt]>
  BootstrapProfile* v24pt_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s24pt1",idname.c_str())));
  auto v24pt_mpt = dynamic_cast<TH1D*>(v24pt_profile->getHist(i));
  auto mpt = getMpt(i);
  v24pt_mpt->Multiply(mpt);
  reth->Add(v24pt_mpt,-2);

  // Add 2 x <v2^4><[pt]>^2
  auto v24_mpt_mpt = getNN(2,4,Form("Ch%s",idname.c_str()),i);
  v24_mpt_mpt->Multiply(mpt);
  v24_mpt_mpt->Multiply(mpt);
  reth->Add(v24_mpt_mpt,2);

  // Subtract <v2^4><[pt^2]>
  auto v24_mpt2 = getNN(2,4,Form("Ch%s",idname.c_str()),i);
  v24_mpt2->Multiply(flowptcont->getCorrHist(i,2));
  reth->Add(v24_mpt2,-1);

  // Subtract 4 x <v2^2*[pt^2]>*<v2^2>
  BootstrapProfile* v22pt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt2",idname.c_str())));
  auto v22pt2_v22 = dynamic_cast<TH1D*>(v22pt2_profile->getHist(i));
  v22pt2_v22->Multiply(getNN(2,2,Form("Ch%s",idname.c_str()),i));
  reth->Add(v22pt2_v22,-4);

  // Add 16 x <v2^2*[pt]>*<v2^2>*<[pt]>
  BootstrapProfile* v22pt_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("Ch%s22pt1",idname.c_str())));
  auto v22pt_v22_mpt = dynamic_cast<TH1D*>(v22pt_profile->getHist(i));
  v22pt_v22_mpt->Multiply(getNN(2,2,Form("Ch%s",idname.c_str()),i));
  v22pt_v22_mpt->Multiply(mpt);
  reth->Add(v22pt_v22_mpt,16);

  // Add 4 x <v2^2>^2*<[pt^2]>
  auto v22_v22_mpt2 = getNN(2,2,Form("Ch%s",idname.c_str()),i);
  v22_v22_mpt2->Multiply(v22_v22_mpt2);
  v22_v22_mpt2->Multiply(flowptcont->getCorrHist(i,2));
  reth->Add(v22_v22_mpt2,4);

  // Subtract 12 x <v2^2>^2*<[pt]>^2
  auto v22_v22_mpt_mpt = getNN(2,2,Form("Ch%s",idname.c_str()),i);
  v22_v22_mpt_mpt->Multiply(v22_v22_mpt_mpt);
  v22_v22_mpt_mpt->Multiply(mpt);
  v22_v22_mpt_mpt->Multiply(mpt);
  reth->Add(v22_v22_mpt_mpt,-12);

  reth->SetNameTitle("CV24pt2","C(v_{2}^{4},[#it{p}_{T}]^{(2)})");
  return reth;
}
TH1D* Calculator::getVNN(int nhar, int npar, std::string id, int i, bool onPt, double min, double max)
{
  TH1D* reth = nullptr;
  //auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  auto clone = dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
  };
  if(id_changed) id = idname;
  clone->SetIDName(id.c_str());
  switch (npar) {
    case 2:
      reth = dynamic_cast<TH1D*>(clone->GetVN2VsX(nhar, onPt, min, max));
      break;
    case 4:
      reth = dynamic_cast<TH1D*>(clone->GetVN4VsX(nhar, onPt, min, max));
      break;
    case 6:
      reth = dynamic_cast<TH1D*>(clone->GetVN6VsX(nhar, onPt, min, max));
      break;
    case 8:
      reth = dynamic_cast<TH1D*>(clone->GetVN8VsX(nhar, onPt, min, max));
      break;
    default:
      printf("Number of particles in correlation should be multiple of 2 and only up to 8 (for now)\n");
      reth = nullptr;
      break;
  }
  reth->SetNameTitle(Form("v%i%i",nhar,npar),Form("v_{%i}{%i}",nhar,npar));
  return reth;
}
TH1D* Calculator::getCNN(int nhar, int npar, std::string id, int i)
{
  TH1D* reth = nullptr;
  auto clone = dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  //auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
  };
  if(id_changed) id = idname;
  clone->SetIDName(id.c_str());
  switch (npar) {
    case 2:
      reth = dynamic_cast<TH1D*>(clone->GetCN2VsX(nhar, kFALSE, 0));
      break;
    case 4:
      reth = dynamic_cast<TH1D*>(clone->GetCN4VsX(nhar, kFALSE, 0));
      break;
    case 6:
      reth = dynamic_cast<TH1D*>(clone->GetCN6VsX(nhar, kFALSE, 0));
      break;
    case 8:
      reth = dynamic_cast<TH1D*>(clone->GetCN8VsX(nhar, kFALSE, 0));
      break;
    default:
      printf("Number of particles in correlation should be multiple of 2 and only up to 8 (for now)\n");
      reth = nullptr;
      break;
  }
  reth->SetNameTitle(Form("c%i%i",nhar,npar),Form("c_{%i}{%i}",nhar,npar));
  return reth;
}
TH1D* Calculator::getNN(int nhar, int npar, std::string id, int i)
{
  auto clone = dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  //auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
  };
  if(id_changed) id = idname;
  clone->SetIDName(id.c_str());
  TH1D* reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti(Form("%i%u", nhar, npar), 0));
  reth->SetNameTitle(Form("%i%i",nhar,npar),Form("%i%i",nhar,npar));
  return reth;
}
TH1D* Calculator::getPtCorr(int i, int n){
  TH1D* reth = dynamic_cast<TH1D*>(flowptcont->getCorrHist(i, n));
  reth->SetNameTitle(Form("mpt%i",n),Form("mpt%i",n));
  return reth;
}


TH1D* Calculator::getKN(int i, int n){
  TH1D* reth = dynamic_cast<TH1D*>(flowptcont->getCumulantHist(i, n));
  reth->SetNameTitle(Form("k%i",n),Form("k%i",n));
  return reth;
}
TH1D* Calculator::getMuN(int i, int n){
  TH1D* reth = dynamic_cast<TH1D*>(flowptcont->getCentralMomentHist(i, n));
  reth->SetNameTitle(Form("mu%i",n),Form("mu%i",n));
  return reth;
}