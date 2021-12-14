#include <iostream>
#include <vector>
using namespace std;

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"

#include "SetStyle.h"
#include "PlotHelper4.h"
#include "CommandLine.h"
#include "ProgressBar.h"
#include "TauHelperFunctions2.h"

int main(int argc, char *argv[]);
double GetCorrection(int hiBin, double eta);
void HistogramXYDivide(TH2D &H);
void HistogramShiftCopy(TH2D &H, TH2D &H2);
void PlotDifference(PdfFileHelper &PdfFile, TH1D &H1, TH1D &H2, string Name);

int main(int argc, char *argv[])
{
   SetThesisStyle();

   CommandLine CL(argc, argv);

   string InputFileName = CL.Get("Input", "vJetTrkSkim_pbpb_2018_data_zmm.root");
   string OutputFileName = CL.Get("Output", "Plots");
   double ZPtCut = CL.GetDouble("ZPtCut",20);
   double TrkEtaCut = CL.GetDouble("TrkEtaCut",-10);
   cout <<ZPtCut<<" "<<TrkEtaCut<<endl;
   bool DoEE = CL.GetBool("EE", false);


   TFile File(InputFileName.c_str());

   TTree *HiTree = (TTree *)File.Get("HiTree");
//   TTree *MixEventTree = (TTree *)File.Get("mixEventSkim");
   TTree *EventTree = (TTree *)File.Get("EventTree");
   TTree *TrackTree = (TTree *)File.Get("trackSkim");

   int hiBin;
   HiTree->SetBranchAddress("hiBin", &hiBin);

   int nEle;
   vector<int> *eleCharge = nullptr;
   vector<float> *elePt = nullptr, *eleEta = nullptr, *elePhi = nullptr;
   vector<float> *eleSCEn = nullptr, *eleSCEta = nullptr, *eleSCPhi = nullptr;
   int nMu;
   vector<int> *muCharge = nullptr;
   vector<int> *muIDTight = nullptr;
   vector<float> *muPt = nullptr, *muEta = nullptr, *muPhi = nullptr;
   EventTree->SetBranchAddress("nEle", &nEle);
   EventTree->SetBranchAddress("eleCharge", &eleCharge);
   EventTree->SetBranchAddress("elePt", &elePt);
   EventTree->SetBranchAddress("eleEta", &eleEta);
   EventTree->SetBranchAddress("elePhi", &elePhi);
   EventTree->SetBranchAddress("eleSCEn", &eleSCEn);
   EventTree->SetBranchAddress("eleSCEta", &eleSCEta);
   EventTree->SetBranchAddress("eleSCPhi", &eleSCPhi);
   EventTree->SetBranchAddress("nMu", &nMu);
   EventTree->SetBranchAddress("muCharge", &muCharge);
   EventTree->SetBranchAddress("muPt", &muPt);
   EventTree->SetBranchAddress("muEta", &muEta);
   EventTree->SetBranchAddress("muPhi", &muPhi);
   EventTree->SetBranchAddress("muIDTight", &muIDTight);

   TFile OutputFile((OutputFileName + ".root").c_str(), "RECREATE");

   TNtuple *ntZ = new TNtuple ("ntZ","","mass:pt:ch1:ch2:hiBin");


   double ZCount = 0;
   double ZCount0030 = 0, ZCount3090 = 0;

   int EntryCount = EventTree->GetEntries();
   ProgressBar Bar(cout, EntryCount);

   for(int iE = 0; iE < EntryCount; iE++)
   {
      Bar.Update(iE);
      Bar.Print();

      HiTree->GetEntry(iE);
      EventTree->GetEntry(iE);

      int nLepton = (DoEE ? nEle : nMu);


      FourVector P1, P2;
      
      map<double, FourVector> Leptons;
      for(int i = 0; i < nLepton; i++)
      {
         double Correction = (DoEE ? GetCorrection(hiBin, (*eleSCEta)[i]) : 1);

         if(DoEE)
            if((*eleSCEta)[i] < -1.392 && (*eleSCPhi)[i] < -0.872665 && (*eleSCPhi)[i] > -1.570796)
               continue;
         if (!DoEE) {
	    if ((*muIDTight)[i]==0) continue;
	 }
         FourVector P;
         if(DoEE)
            P.SetPtEtaPhiMass((*elePt)[i] * Correction, (*eleEta)[i], (*elePhi)[i], 0.000511);
         else
            P.SetPtEtaPhiMass((*muPt)[i] * Correction, (*muEta)[i], (*muPhi)[i], 0.105658);

         if(DoEE == true && (*eleCharge)[i] < 0)
            P[0] = -P[0];
         if(DoEE == false && (*muCharge)[i] < 0)
            P[0] = -P[0];

         Leptons.insert(pair<double, FourVector>(-P.GetPT(), P));
      }

      if(Leptons.size() < 2)
         continue;

      map<double, FourVector>::iterator iter = Leptons.begin();
      P1 = iter->second;
      iter++;
      P2 = iter->second;
      if(GetDR(P1, P2) < 0.2 && iter != Leptons.end())
      {
         iter++;
         if(iter != Leptons.end())
            P2 = iter->second;
         else
            P2.SetPtEtaPhiMass(0, 0, 0, 0);
      }

      int Charge1 = (P1[0] > 0) ? 1 : -1;
      int Charge2 = (P2[0] > 0) ? 1 : -1;

      P1[0] = fabs(P1[0]);
      P2[0] = fabs(P2[0]);

      if(P1.GetPT() < 20)
         continue;
      if(P2.GetPT() < 20)
         continue;

      if(P1.GetAbsEta() > 2.1)
         continue;
      if(P2.GetAbsEta() > 2.1)
         continue;

      FourVector LL = P1 + P2;

      ntZ->Fill(LL.GetMass(),LL.GetPT(),Charge1,Charge2,hiBin);
      
   }

   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();


   ntZ->Write();

   OutputFile.Close();

   return 0;
}

double GetCorrection(int hiBin, double eta)
{
   eta = fabs(eta);

   if(eta < 1.442)
   {
      if(hiBin < 20)   return 0.990;
      if(hiBin < 60)   return 1.006;
      return 1.016;
   }
   if(eta > 1.566 && eta < 2.5)
   {
      if(hiBin < 20)   return 0.976;
      if(hiBin < 60)   return 1.015;
      return 1.052;
   }

   return 0;
}

void HistogramXYDivide(TH2D &H)
{
   TH1D *H1 = (TH1D *)H.ProjectionX();
   TH1D *H2 = (TH1D *)H.ProjectionY();

   for(int iX = 1; iX <= H1->GetNbinsX(); iX++)
      for(int iY = 1; iY <= H2->GetNbinsX(); iY++)
         H.SetBinContent(iX, iY, H.GetBinContent(iX, iY) / H1->GetBinContent(iX) / H2->GetBinContent(iY));
}

void HistogramShiftCopy(TH2D &H1, TH2D &H2)
{
   int NX = H1.GetNbinsX();
   int NY = H1.GetNbinsY();

   for(int iX = 1; iX <= NX; iX++)
      for(int iY = 1; iY <= NY; iY++)
         H2.SetBinContent(iX, iY, H1.GetBinContent((iX + NX / 2 + iY - 2) % NX + 1, iY));
}

void PlotDifference(PdfFileHelper &PdfFile, TH1D &H1, TH1D &H2, string Name)
{
   TH1D *H = (TH1D *)H1.Clone(Name.c_str());
   H->Add(&H2, -1);
   H->SetMinimum(0);
   H->SetMaximum(25);

   H->Write();

   PdfFile.AddPlot(H);
}




