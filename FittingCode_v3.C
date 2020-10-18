//N. Nepal (nepal1n@cmich.edu) for combined fit

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <string>


Double_t funcA(double *x, double *p){
  double a = 0;
  a = p[0]*exp(-p[1]*x[0]);
  return a*p[1];
}

Double_t funcB(double *x, double *p){
  double a = 0;
  double a1 = exp(-p[1]*x[0]);
  double a4 = exp(-p[4]*x[0]);
  a = p[0]*(1-p[2]-p[3])*p[1]*(a1/(p[4]-p[1])+a4/(p[1]-p[4]));
  return a*p[4];
}

Double_t funcC(double *x, double *p){
  double a = 0;
  double a1 = exp(-p[1]*x[0]);
  double a4 = exp(-p[4]*x[0]);
  double a7 = exp(-p[7]*x[0]);
  a = p[0]*(1-p[2]-p[3])*(1-p[5]-p[6])*p[1]*p[4]*(a1/(p[7]-p[1])/(p[4]-p[1])+a4/(p[1]-p[4])/(p[7]-p[4])+a7/(p[1]-p[7])/(p[4]-p[7]));
  return a*p[7];
}

Double_t funcE(double *x, double *p){
  double a = 0;
  double a1 = exp(-p[1]*x[0]);
  double a10 = exp(-p[10]*x[0]);
  a = p[0]*p[2]*p[1]*(a1/(p[10]-p[1])+a10/(p[1]-p[10]));
  return a*p[10];
}

Double_t funcF(double *x, double *p){
  double a = 0;
  double a1 = exp(-p[1]*x[0]);
  double a4 = exp(-p[4]*x[0]);
  double a10 = exp(-p[10]*x[0]);
  double a13 = exp(-p[13]*x[0]);
  double A_B_F = p[0]*(1-p[2]-p[3])*p[1]*p[4]*p[5]*(a1/(p[4]-p[1])/(p[13]-p[1])+a4/(p[1]-p[4])/(p[13]-p[4])+a13/(p[1]-p[13])/(p[4]-p[13]));
  double A_E_F = p[0]*p[2]*(1-p[11])*p[1]*p[10]*(a1/(p[10]-p[1])/(p[13]-p[1])+a10/(p[1]-p[10])/(p[13]-p[10])+a13/(p[1]-p[13])/(p[10]-p[13]));
  a = A_B_F + A_E_F;
  return a*p[13];
}

Double_t funcH(double *x, double *p){
  double a = 0;
  double a1 = exp(-p[1]*x[0]);
  double a16 = exp(-p[16]*x[0]);
  a = p[0]*p[3]*p[1]*(a1/(p[16]-p[1])+a16/(p[1]-p[16]));
  return a*p[16];
}

Double_t funcI(double *x, double *p){
  double a = 0;
  double a1 = exp(-p[1]*x[0]);
  double a4 = exp(-p[4]*x[0]);
  double a10 = exp(-p[10]*x[0]);
  double a16 = exp(-p[16]*x[0]);
  double a19 = exp(-p[19]*x[0]);  
  double A_B_I = p[0]*(1-p[2]-p[3])*p[6]*p[1]*p[4]*(a1/(p[4]-p[1])/(p[19]-p[1])+a4/(p[1]-p[4])/(p[19]-p[4])+a19/(p[1]-p[19])/(p[4]-p[19]));
  double A_E_I = p[0]*p[2]*p[11]*p[1]*p[10]*(a1/(p[10]-p[1])/(p[19]-p[1])+a10/(p[1]-p[10])/(p[19]-p[10])+a19/(p[1]-p[19])/(p[10]-p[19]));
  double A_H_I = p[0]*p[3]*p[1]*p[16]*(a1/(p[16]-p[1])/(p[19]-p[1])+a16/(p[1]-p[16])/(p[19]-p[16])+a19/(p[1]-p[19])/(p[16]-p[19]));
  a = A_B_I + A_E_I + A_H_I;
  return a*p[19];
}

//defining the global function
Double_t FT(double *x, double *p){
  Double_t returnValue = p[20];
  if(x[0]>-5 && x[0]<5){
    TF1::RejectPoint();
    return returnValue;
  }
  Double_t  func = funcA(x,p)+funcB(x,p)+funcC(x,p)+funcE(x,p)+funcF(x,p)+funcH(x,p)+funcI(x,p);
  if(x[0]>0.0) returnValue += func;
  return returnValue;
}

Double_t F1NT(double *x, double *p){
  Double_t returnValue = p[21];
  if(x[0]>-5 && x[0]<5){
    TF1::RejectPoint();
    return returnValue;
  }
  Double_t  func = (funcA(x,p)*p[2]+funcB(x,p)*p[5]+funcC(x,p)*p[8]+funcE(x,p)*p[11]+funcF(x,p)*p[14])*p[23]+2*p[23]*(1-p[23])*(funcA(x,p)*p[3]+funcB(x,p)*p[6]+funcC(x,p)*p[9]);
  if(x[0]>0.0) returnValue += func;
  return returnValue;
}

Double_t F2NT(double *x, double *p){
  Double_t returnValue = p[22];
  if(x[0]>-5 && x[0]<5){
    TF1::RejectPoint();
    return returnValue;
  }
  Double_t  func = p[23]*p[23]*(funcA(x,p)*p[3]+funcB(x,p)*p[6]+funcC(x,p)*p[9]);
  if(x[0]>0.0) returnValue += func;
  return returnValue;
}

//defination of shared parameter
//total decay function
int iparB[24] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
//gate on p1n function
int iparSB[24] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
//gate on p2n function
int iparSB1[24] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};

struct GlobalChi2 {
  GlobalChi2( ROOT::Math::IMultiGenFunction &f1,
	      ROOT::Math::IMultiGenFunction &f2,
	      ROOT::Math::IMultiGenFunction &f3) :
    fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}
  double operator() (const double *par) const{
    double p1[24];
    for(int i = 0; i<24; ++i) p1[i] = par[iparB[i]];
    double p2[24];
    for(int i = 0; i<24; ++i) p2[i] = par[iparSB[i]];
    double p3[24];
    for(int i = 0; i<24; ++i) p3[i] = par[iparSB1[i]];    
    return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
  }
  const ROOT::Math::IMultiGenFunction *fChi2_1;
  const ROOT::Math::IMultiGenFunction *fChi2_2;
  const ROOT::Math::IMultiGenFunction *fChi2_3;
};

void FittingIsotope(TString isoName, Int_t AtNo, Int_t NeuNo, char *decayFile, char *bkgFile, double st, double et, Int_t rebi=1) {
  
#if defined(__CINT__) && !defined(__MAKECINT__)
  cout << "ERROR: This tutorial can run only using ACliC, you must run it by doing: " << endl;
  cout << "\t .x $ROOTSYS/tutorials/fit/combinedFit.C+" << endl;
  return;
#endif
  
  cout<<"I am going to fit the isotope which has"<<endl;
  cout<<"z = "<<AtNo<<endl;
  cout<<"n = "<<NeuNo<<endl;
  Int_t k1=9999, k2=9999, k3=9999, k4=9999, k5=9999, k6=9999, k7=9999; //the crossponding isotope line number

  Double_t a, b, c, d, e, f, g, h, l, n, o;
  Int_t m = 157; 
  Double_t Z[m], N[m], T[m], TL[m], TU[m], Pn[m], PnL[m], PnU[m], Pnn[m], PnnL[m], PnnU[m];

  ifstream infile;
  infile.open("parameters.txt"); //open nuclear data file or parameters
  //to check wheather the file is properly open or not
  if( !infile ){
    cout<<"Unable to open the file!!! check the file name and try again"<<endl;
    exit(1); //terminate with error
  }

  Int_t i = 0;
  while (!infile.eof()) {
    infile>>a>>b>>c>>d>>e>>f>>g>>h>>l>>n>>o;
    Z[i] = a;
    N[i] = b;
    T[i] = c;
    TL[i] = d;
    TU[i] = e;
    Pn[i] = f;
    PnL[i] = g;
    PnU[i] = h;
    Pnn[i] = l;
    PnnL[i] = n;
    PnnU[i] = o;
    i++;
  }

  //Let's print one random row
  // cout<<Z[4]<<"\t"<<N[4]<<"\t"<<T[4]<<"\t"<<TL[4]<<"\t"<<TU[4]<<"\t"<<Pn[4]<<"\t"<<PnL[4]<<"\t"<<PnU[4]<<"\t"<<Pnn[4]<<"\t"<<PnnL[4]<<"\t"<<PnnU[4]<<endl;


  //parent nuclei I
  for(Int_t j=0; j<m; j++){
    if( Z[j] == AtNo && N[j] == NeuNo ) k1 = j;
  }
  Double_t act1 = log(2)/T[k1]; //activity of parent nuclei
  Double_t act1L = (log(2)/T[k1])+(log(2)/TL[k1]);
  Double_t act1U = (log(2)/T[k1])+(log(2)/TU[k1]);
  Double_t Pn1 = Pn[k1]/100;
  Double_t Pnn1 = Pnn[k1]/100;

  //direct beta-decay daughter
  for(Int_t j=0; j<m; j++){
    if( Z[j] == Int_t(AtNo+1) && N[j] == Int_t(NeuNo-1)) k2 = j;
  }
  Double_t act2 = log(2)/T[k2]; //activity of daughter nuclei
  Double_t act2L = (log(2)/T[k2])+(log(2)/TL[k2]);
  Double_t act2U = (log(2)/T[k2])+(log(2)/TU[k2]);
  Double_t Pn2 = Pn[k2]/100;
  Double_t Pn2L = Pn2+(PnL[k2]/100);
  Double_t Pn2U = Pn2+(PnU[k2]/100);
  Double_t Pnn2 = Pnn[k2]/100;
  Double_t Pnn2L = Pnn2+(PnnL[k2]/100);
  Double_t Pnn2U = Pnn2+(PnnU[k2]/100);

  //direct beta-decay grand daughter
  for(Int_t j=0; j<m; j++){
    if( Z[j] == Int_t(AtNo+2) && N[j] == Int_t(NeuNo-2)) k3 = j;
  } 
  Double_t act3 = log(2)/T[k3];
  Double_t act3L = (log(2)/T[k3])+(log(2)/TL[k3]);
  Double_t act3U = (log(2)/T[k3])+(log(2)/TU[k3]);
  Double_t Pn3 = Pn[k3]/100;
  Double_t Pn3L = Pn3+(PnL[k3]/100);
  Double_t Pn3U = Pn3+(PnU[k3]/100);
  Double_t Pnn3 = Pnn[k3]/100;
  Double_t Pnn3L = Pnn3+(PnnL[k3]/100);
  Double_t Pnn3U = Pnn3+(PnnU[k3]/100);

  //daughter by beta-decay 1n emission
  for(Int_t j=0; j<m; j++){
    if( Z[j] == Int_t(AtNo+1) && N[j] == Int_t(NeuNo-2)) k4 = j;
  }
  Double_t act4 = log(2)/T[k4]; 
  Double_t act4L = (log(2)/T[k4])+(log(2)/TL[k4]);
  Double_t act4U = (log(2)/T[k4])+(log(2)/TU[k4]);
  Double_t Pn4 = Pn[k4]/100;
  Double_t Pn4L = Pn4+(PnL[k4]/100);
  Double_t Pn4U = Pn4+(PnU[k4]/100);
  Double_t Pnn4 = Pnn[k4]/100;
  Double_t Pnn4L = Pnn4+(PnnL[k4]/100);
  Double_t Pnn4U = Pnn4+(PnnU[k4]/100);

  //daughter by beta-decay 1n emission and then beta-decay
  for(Int_t j=0; j<m; j++){
    if( Z[j] == Int_t(AtNo+2) && N[j] == Int_t(NeuNo-3)) k5 = j;
  } 
  Double_t act5 = log(2)/T[k5]; 
  Double_t act5L = (log(2)/T[k5])+(log(2)/TL[k5]);
  Double_t act5U = (log(2)/T[k5])+(log(2)/TU[k5]);
  Double_t Pn5 = Pn[k5]/100;
  Double_t Pn5L = Pn5+(PnL[k5]/100);
  Double_t Pn5U = Pn5+(PnU[k5]/100);
  Double_t Pnn5 = Pnn[k5]/100;
  Double_t Pnn5L = Pnn5+(PnnL[k5]/100);
  Double_t Pnn5U = Pnn5+(PnnU[k5]/100);

  //daughter by beta-decay 2n emission
  for(Int_t j=0; j<m; j++){
    if( Z[j] == Int_t(AtNo+1) && N[j] == Int_t(NeuNo-3)) k6 = j;
  }
  Double_t act6 = log(2)/T[k6]; 
  Double_t act6L = (log(2)/T[k6])+(log(2)/TL[k6]);
  Double_t act6U = (log(2)/T[k6])+(log(2)/TU[k6]);
  Double_t Pn6 = Pn[k6]/100;
  Double_t Pn6L = Pn6+(PnL[k6]/100);
  Double_t Pn6U = Pn6+(PnU[k6]/100);
  Double_t Pnn6 = Pnn[k6]/100;
  Double_t Pnn6L = Pnn6+(PnnL[k6]/100);
  Double_t Pnn6U = Pnn6+(PnnU[k6]/100);

  //daughter by beta-decay 2n emission and then beta decay
  for(Int_t j=0; j<m; j++){
    if( Z[j] == Int_t(AtNo+2) && N[j] == Int_t(NeuNo-4)) k7 = j;
  }
  Double_t act7 = log(2)/T[k7];
  Double_t act7L = (log(2)/T[k7])+(log(2)/TL[k7]);
  Double_t act7U = (log(2)/T[k7])+(log(2)/TU[k7]);

  //Initializing the backgrounds, neutron efficiency etc ..
  Double_t rN0, rN0L, rN0H; 
  Double_t rbkg1, rbkg1L, rbkg1H;
  Double_t rbkg2, rbkg2L, rbkg2H;
  Double_t rbkg3, rbkg3L, rbkg3H; 
  Double_t neuEff;

  //open text file which has information of histogram backgrounds and neutron efficeincy
  cout<<"Hey I am going to opent background txt file"<<endl;

  std::ifstream ifs(bkgFile);
  if(ifs.is_open()){
    cout<<"I am opening txt file, which has information of histogram backgrounds"<<endl;
    ifs>>rN0>>rN0L>>rN0H;
    ifs>>rbkg1>>rbkg1L>>rbkg1H;
    ifs>>rbkg2>>rbkg2L>>rbkg2H;
    ifs>>rbkg3>>rbkg3L>>rbkg3H;
    ifs>>neuEff;
  }
  else cout<<"Can't open background file"<<endl; 

  //for the rebining
  Double_t N0 = rN0*rebi;
  Double_t N0L = rN0L*rebi;
  Double_t N0H = rN0H*rebi;
  Double_t bkg1 = rbkg1*rebi;
  Double_t bkg1L = rbkg1L*rebi;
  Double_t bkg1H = rbkg1H*rebi;
  Double_t bkg2 = rbkg2*rebi;
  Double_t bkg2L = rbkg2L*rebi;
  Double_t bkg2H = rbkg2H*rebi;
  Double_t bkg3 = rbkg3*rebi;
  Double_t bkg3L = rbkg3L*rebi;
  Double_t bkg3H = rbkg3H*rebi;

  /*
  cout<<"Total number of activity = "<<N0<<endl;
  cout<<"background of ion-beta histogram = "<<bkg1<<endl;
  cout<<"background of ion-beta histogram when one neutron is detected = "<<bkg2<<endl;
  cout<<"background of ion-beta histogram when two neutron is detected = "<<bkg3<<endl;
  cout<<"neutron efficiency = "<<neuEff<<endl;
  */

  TFile *in = new TFile(decayFile);
  TH1F *decaycurve1 = (TH1F*)in->Get(isoName+"Decay");
  TH1F *n1decaycurve1 = (TH1F*)in->Get(isoName+"Decay_1n_corrtd");
  TH1F *n2decaycurve1 = (TH1F*)in->Get(isoName+"Decay_2n_corrtd");

  //avoiding bins which has negative counts (replaced by ZERO)
  Double_t n_bins = n1decaycurve1->GetXaxis()->GetNbins();
  cout<<"#####TOTAL NUMBER OF BINS IS = "<<n_bins<<endl;
  //create new histogram for 1n and 2n decay curve
  TH1F *n1decaycurve = new TH1F("n1decaycurve","decay curve with 1n emission", n1decaycurve1->GetNbinsX(), n1decaycurve1->GetXaxis()->GetXmin(), n1decaycurve1->GetXaxis()->GetXmax());
  for(Int_t i=1; i<n_bins+1; i++){
    Double_t xy = 0;
    xy = n1decaycurve1->GetBinContent(i);
    if(xy<0) n1decaycurve->SetBinContent(i,0);
    else if(xy==0 || xy>0) n1decaycurve->SetBinContent(i,xy);
  }
  TH1F *n2decaycurve = new TH1F("n2decaycurve","decay curve with 2n emission", n2decaycurve1->GetNbinsX(), n2decaycurve1->GetXaxis()->GetXmin(), n2decaycurve1->GetXaxis()->GetXmax());
  //  for(Int_t i=1; i<n2decaycurve1->GetNbinsX(); i++){
  for(Int_t i=1; i<n_bins+1; i++){
    Double_t xy = 0;
    xy = n2decaycurve1->GetBinContent(i);
    if(xy<0) n2decaycurve->SetBinContent(i,0);
    else if(xy==0 || xy>0) n2decaycurve->SetBinContent(i,xy);
  }
  TH1F *decaycurve = new TH1F("decaycurve","decay curve", decaycurve1->GetNbinsX(), decaycurve1->GetXaxis()->GetXmin(), decaycurve1->GetXaxis()->GetXmax());
  //  for(Int_t i=1; i<decaycurve1->GetNbinsX(); i++){
  for(Int_t i=1; i<n_bins+1; i++){
    Double_t xy = 0;
    xy = decaycurve1->GetBinContent(i);
    if(xy<0) decaycurve->SetBinContent(i,0);
    else if(xy==0 || xy>0) decaycurve->SetBinContent(i,xy);
  }

  decaycurve->Rebin(rebi);
  n1decaycurve->Rebin(rebi);
  n2decaycurve->Rebin(rebi);

  TF1 *ftotal = new TF1("ftotal",FT,st,et,24);
  TF1 *f1n = new TF1("f1n",F1NT,st,et,24);
  TF1 *f2n = new TF1("f2n",F2NT,st,et,24);

  //perform now global fit
  ROOT::Math::WrappedMultiTF1 wfB(*ftotal,1);
  ROOT::Math::WrappedMultiTF1 wfSB(*f1n, 1);
  ROOT::Math::WrappedMultiTF1 wfSB1(*f2n,1);
  
  ROOT::Fit::DataOptions opt;
  //set the data range
  //total decay
  ROOT::Fit::DataRange rangeB;
  rangeB.SetRange(st,et);
  ROOT::Fit::BinData dataB(opt, rangeB);
  ROOT::Fit::FillData(dataB, decaycurve);
  //1n decay
  ROOT::Fit::DataRange rangeSB;
  rangeSB.SetRange(st,et);
  ROOT::Fit::BinData dataSB(opt, rangeSB);
  ROOT::Fit::FillData(dataSB, n1decaycurve);
  //2n decay
  ROOT::Fit::DataRange rangeSB1;
  rangeSB1.SetRange(st,et);
  ROOT::Fit::BinData dataSB1(opt, rangeSB1);
  ROOT::Fit::FillData(dataSB1, n2decaycurve);
  
  ROOT::Fit::Chi2Function chi2_B(dataB, wfB);
  ROOT::Fit::Chi2Function chi2_SB(dataSB, wfSB);
  ROOT::Fit::Chi2Function chi2_SB1(dataSB1, wfSB1);

  GlobalChi2 GlobalChi2(chi2_B, chi2_SB, chi2_SB1);

  ROOT::Fit::Fitter fitter;

  const int Npar = 24;

   //initializing the fit parameters!!!  
  double par0[Npar] = {N0, act1, Pn1, Pnn1, act2, Pn2, Pnn2, act3, Pn3, Pnn3, act4, Pn4, Pnn4, act5, Pn5, Pnn5, act6, Pn6, Pnn6, act7, bkg1, bkg2, bkg3, neuEff};
  
  fitter.Config().SetParamsSettings(24, par0);
  
  //FIX or SETLIMITS CONDITIONS
  ifstream parSet;
  parSet.open("parSetting.txt"); //open parameter setting text file
  if( !parSet ){
    cout<<"Unable to open parSetting.txt, where is this file?"<<endl;
    exit(1); //terminate with error
  }
  Double_t parNum[25], parS[25]; 
  Double_t pn, pv;
  Int_t ij = 0;
  while(!parSet.eof()){
    parSet>>pn>>pv;
    parNum[ij] = pn;
    parS[ij] = pv;
    ij++;
   }
  //set limits or fixing  on parameters 
  if(parS[0]==1)   fitter.Config().ParSettings(0).SetLimits(N0L,N0H);
  if(parS[0]==0) fitter.Config().ParSettings(0).Fix();
  if(parS[1]==1) fitter.Config().ParSettings(1).SetLimits(0,360);
  if(parS[1]==0) fitter.Config().ParSettings(1).Fix();
  if(parS[2]==1) fitter.Config().ParSettings(2).SetLimits(0,100);
  if(parS[2]==0) fitter.Config().ParSettings(2).Fix();
  if(parS[3]==1) fitter.Config().ParSettings(3).SetLimits(0,100);
  if(parS[3]==0) fitter.Config().ParSettings(3).Fix();
  
  if(parS[4]==1) fitter.Config().ParSettings(4).SetLimits(act2L,act2U);
  if(parS[5]==1) fitter.Config().ParSettings(5).SetLimits(Pn2L,Pn2U);
  if(parS[6]==1) fitter.Config().ParSettings(6).SetLimits(Pnn2L,Pnn2U);
  if(parS[7]==1) fitter.Config().ParSettings(7).SetLimits(act3L,act3U);
  if(parS[8]==1) fitter.Config().ParSettings(8).SetLimits(Pn3L,Pn3U);
  if(parS[9]==1) fitter.Config().ParSettings(9).SetLimits(Pnn3L,Pnn3U);
  if(parS[10]==1) fitter.Config().ParSettings(10).SetLimits(act4L,act4U);
  if(parS[11]==1) fitter.Config().ParSettings(11).SetLimits(Pn4L,Pn4U);
  if(parS[12]==1) fitter.Config().ParSettings(12).SetLimits(Pnn4L,Pnn4U);
  if(parS[13]==1) fitter.Config().ParSettings(13).SetLimits(act5L,act5U);
  if(parS[14]==1) fitter.Config().ParSettings(14).SetLimits(Pn5L,Pn5U);
  if(parS[15]==1) fitter.Config().ParSettings(15).SetLimits(Pnn5L,Pnn5U);
  if(parS[16]==1) fitter.Config().ParSettings(16).SetLimits(act6L,act6U);
  if(parS[17]==1) fitter.Config().ParSettings(17).SetLimits(Pn6L,Pn6U);
  if(parS[18]==1) fitter.Config().ParSettings(18).SetLimits(Pnn6L,Pnn6U);
  if(parS[19]==1) fitter.Config().ParSettings(19).SetLimits(act7L,act7U);
  
  if(parS[4]==0) fitter.Config().ParSettings(4).Fix();
  if(parS[5]==0) fitter.Config().ParSettings(5).Fix();
  if(parS[6]==0) fitter.Config().ParSettings(6).Fix();
  if(parS[7]==0) fitter.Config().ParSettings(7).Fix();
  if(parS[8]==0) fitter.Config().ParSettings(8).Fix();
  if(parS[9]==0) fitter.Config().ParSettings(9).Fix();
  if(parS[10]==0) fitter.Config().ParSettings(10).Fix();
  if(parS[11]==0) fitter.Config().ParSettings(11).Fix();
  if(parS[12]==0) fitter.Config().ParSettings(12).Fix();
  if(parS[13]==0) fitter.Config().ParSettings(13).Fix();
  if(parS[14]==0) fitter.Config().ParSettings(14).Fix();
  if(parS[15]==0) fitter.Config().ParSettings(15).Fix();
  if(parS[16]==0) fitter.Config().ParSettings(16).Fix();
  if(parS[17]==0) fitter.Config().ParSettings(17).Fix();
  if(parS[18]==0) fitter.Config().ParSettings(18).Fix();
  if(parS[19]==0) fitter.Config().ParSettings(19).Fix();

  if(parS[20]==1) fitter.Config().ParSettings(20).SetLimits(bkg1L,bkg1H);
  if(parS[21]==1) fitter.Config().ParSettings(21).SetLimits(bkg2L,bkg2H);
  if(parS[22]==1) fitter.Config().ParSettings(22).SetLimits(bkg3L,bkg3H);

  if(parS[20]==0) fitter.Config().ParSettings(20).Fix();
  if(parS[21]==0) fitter.Config().ParSettings(21).Fix();
  if(parS[22]==0) fitter.Config().ParSettings(22).Fix();
  //neutron efficiency always fixed!!!
  if(parS[23]==0) fitter.Config().ParSettings(23).Fix();

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");

  //fit FCN function directly
  //specify optionally data size and flag to indicate that is a chi2 fit
  fitter.FitFCN(24,GlobalChi2,0,dataB.Size()+dataSB.Size()+dataSB1.Size(),true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  //newly added to store result as a root files
  TFile *fOut;
  TString foutname ="fittingResults/" + isoName+"_fitResults.root";
  fOut = new TFile(foutname,"recreate");
  
  TCanvas *c1 = new TCanvas("c1","c1");
  gStyle->SetOptFit(0);
  ftotal->SetFitResult(result, iparB);
  ftotal->SetRange(rangeB().first, rangeB().second);
  ftotal->SetLineColor(kRed);
  decaycurve->GetListOfFunctions()->Add(ftotal);
  decaycurve->Draw();
  c1->Write();

  TCanvas *c2 = new TCanvas("c2","c2");
  f1n->SetFitResult(result, iparSB);
  f1n->SetRange(rangeSB().first,rangeSB().second);
  f1n->SetLineColor(kGreen);
  n1decaycurve->GetListOfFunctions()->Add(f1n);
  n1decaycurve->Draw();
  c2->Write();

  TCanvas *c3 = new TCanvas("c3","c3");
  f2n->SetFitResult(result, iparSB1);
  f1n->SetRange(rangeSB1().first, rangeSB1().second);
  f2n->SetLineColor(kYellow);
  n2decaycurve->GetListOfFunctions()->Add(f2n);
  n2decaycurve->Draw();
  c3->Write();

  double halflife = log(2)/ftotal->GetParameter(1);
  double halflife_error = log(2)/ftotal->GetParameter(1)/ftotal->GetParameter(1)*ftotal->GetParError(1);
  double P1n_value = ftotal->GetParameter(2)*100;
  double P1n_value_error = ftotal->GetParError(2)*100;
  double P2n_value = ftotal->GetParameter(3)*100;
  double P2n_value_error = ftotal->GetParError(3)*100;

  cout<<"Half-life "<<halflife<<" ms +-"<< halflife_error <<endl;
  cout<<"P1n value "<<P1n_value<<"% +-"<<P1n_value_error <<endl;
  cout<<"P2n value "<<P2n_value<<"% +_"<<P2n_value_error <<endl;
   
  TH1F *histo1 = new TH1F("histo1","Error from Decay Hist",40,-50,50);
  for(Int_t i=1; i<n_bins+1;i++){
    double residual = (decaycurve->GetBinContent(i+1) - ftotal->Eval(decaycurve->GetBinCenter(i+1)))/decaycurve->GetBinError(i+1);
    histo1->Fill(residual);
  }
  histo1->Write();
  
  TH1F *histo2 = new TH1F("histo2","Error from Decay Hist 1n emi",40,-50,50);
  for(Int_t i=1; i<n_bins+1;i++){
    double residual = (n1decaycurve->GetBinContent(i+1) - f1n->Eval(n1decaycurve->GetBinCenter(i+1)));
    histo2->Fill(residual);
  }
  histo2->Write();
  
  TH1F *histo3 = new TH1F("histo3","Error from Decay Hist 2n emi",40,-50,50);
  for(Int_t i=1; i<n_bins+1;i++){
    double residual = (n2decaycurve->GetBinContent(i+1) - f2n->Eval(n2decaycurve->GetBinCenter(i+1)));
    histo3->Fill(residual);
  }
  histo3->Write();

   //THE MAIN PROBLEM HERE IS HAVING GETBINERROR VALUE EQUALS TO ZERO!!!
    
   TH1F *hist1 = new TH1F("hist1","residual decay",decaycurve->GetNbinsX(),decaycurve->GetXaxis()->GetXmin(),decaycurve->GetXaxis()->GetXmax());
   for(Int_t i=1;i<n_bins+1;i++){
     //  double residual = (decaycurve->GetBinContent(i+1) - ftotal->Eval(decaycurve->GetBinCenter(i+1) )) / decaycurve->GetBinError(i+1);
     double residual = (decaycurve->GetBinContent(i+1) - ftotal->Eval(decaycurve->GetBinCenter(i+1) ));
    hist1->SetBinContent(i,residual);
    hist1->SetBinError(i,1);
   }

   TH1F *hist2 = new TH1F("hist2","residual decay P1n",n1decaycurve->GetNbinsX(),n1decaycurve->GetXaxis()->GetXmin(),n2decaycurve->GetXaxis()->GetXmax());
   for(Int_t i=1;i<n_bins+1;i++){
     //  double residual = (n1decaycurve->GetBinContent(i+1) - f1n->Eval(n1decaycurve->GetBinCenter(i+1) )) / n1decaycurve->GetBinError(i+1);
     double residual = (n1decaycurve->GetBinContent(i+1) - f1n->Eval(n1decaycurve->GetBinCenter(i+1) ));
    hist2->SetBinContent(i,residual);
    hist2->SetBinError(i,1);
   }

   TH1F *hist3 = new TH1F("hist3","residual decay P2n",n2decaycurve->GetNbinsX(),n2decaycurve->GetXaxis()->GetXmin(),n2decaycurve->GetXaxis()->GetXmax());
   for(Int_t i=1;i<n_bins+1;i++){
     //     double residual = (n2decaycurve->GetBinContent(i+1) - f2n->Eval(n2decaycurve->GetBinCenter(i+1) )) / n2decaycurve->GetBinError(i+1);
     double residual = (n2decaycurve->GetBinContent(i+1) - f2n->Eval(n2decaycurve->GetBinCenter(i+1) ));
    hist3->SetBinContent(i,residual);
     hist3->SetBinError(i,1);
     //    cout<<"Bin Error at "<<i<<" is = "<<n2decaycurve->GetBinError(i+1)<<endl;
   }
   hist1->Write("e1c");
   hist2->Write("e1c");
   hist3->Write("e1c");
   //these histograms are not produced because the value of bin error is 0 when there is zero count
   /*/CHECK THIS AGAINNNN    
   TCanvas *cc2 = new TCanvas("cc2","cc2");
   cc2->Divide(1,3);
   cc2->cd(1);
   hist1->Draw("e1c");
   cc2->cd(2);
   hist2->Draw("e1c");
   cc2->cd(3);
   hist3->Draw("e1c");
   */


  fOut->Close();

  //now to save the result and fitting parameters in separate file
  ofstream myfile;
  myfile.open("fittingResults/"+isoName+"_result.txt");
  myfile<<"Half Life of "<<isoName<<" is = "<<halflife<<" ms +- "<<halflife_error<<endl;
  myfile<<"P1n value of "<<isoName<<" is = "<<P1n_value<<"% +-"<<P1n_value_error <<endl;
  myfile<<"P2n value of "<<isoName<<" is = "<<P2n_value<<"% +-"<<P2n_value_error <<endl;
  myfile<<"decay time range of fitting is "<<st<<" (ms) to "<<et<<" (ms) "<<endl;
  myfile<<"\n\n ---- PARAMETERS USED IN THIS FIT ARE ---- "<<endl;
  myfile<<"Note: parameters as input (in terms of half life and neu emission percentage\n"<<endl;
  myfile<<"N0 = "<<N0<<endl;
  myfile<<"act1 = "<<log(2)/act1<<endl;
  myfile<<"Pn1 = "<<Pn1*100<<endl;
  myfile<<"Pnn1 = "<<Pnn1*100<<endl;
  myfile<<"act2 = "<<log(2)/act2<<endl;
  myfile<<"Pn2 = "<<Pn2*100<<endl;
  myfile<<"Pnn2 = "<<Pnn2*100<<endl;
  myfile<<"act3 = "<<log(2)/act3<<endl;
  myfile<<"Pn3 = "<<Pn3*100<<endl;
  myfile<<"Pnn3 = "<<Pnn3*100<<endl;
  myfile<<"act4 = "<<log(2)/act4<<endl;
  myfile<<"Pn4 = "<<Pn4*100<<endl;
  myfile<<"Pnn4 = "<<Pnn4*100<<endl;
  myfile<<"act5 = "<<log(2)/act5<<endl;
  myfile<<"Pn5 = "<<Pn5*100<<endl;
  myfile<<"Pnn5 = "<<Pnn5*100<<endl;
  myfile<<"act6 = "<<log(2)/act6<<endl;
  myfile<<"Pn6 = "<<Pn6*100<<endl;
  myfile<<"Pnn6 = "<<Pnn6*100<<endl;
  myfile<<"act7 = "<<log(2)/act7<<endl;
  myfile<<"bkg1 = "<<bkg1<<endl;
  myfile<<"bkg2 = "<<bkg2<<endl;
  myfile<<"bkg3 = "<<bkg3<<endl;
  myfile<<"neuEff = "<<neuEff<<endl;
  myfile<<"rebin = "<<rebi<<endl;
  myfile.close();
  
}

