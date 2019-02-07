#include "Pulse.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "Riostream.h"
#include <iostream>     
#include <algorithm>

TChain *tr;
unsigned short nPU;
unsigned short nhits;
unsigned short hashedIndex[61200];
float  eXtal[61200];
float  tXtal[61200];
float  eAPD1[61200];
float  tAPD1[61200];
float  eAPD2[61200];
float  tAPD2[61200];
TBranch *b_nPU;
TBranch *b_nhits;
TBranch *b_hashedIndex;
TBranch *b_eXtal;
TBranch *b_tXtal;
TBranch *b_eAPD1;
TBranch *b_tAPD1;
TBranch *b_eAPD2;
TBranch *b_tAPD2;

TChain *trSig;
unsigned short nhitsSig;
unsigned short hashedIndexSig[61200];
float  eSig[61200];
TBranch *b_nhitsSig;
TBranch *b_hashedIndexSig;
TBranch *b_eSig;

int h2eta[61200];
int h2phi[61200];
double sinTheta[85];
double tZero[85];

const int NCHMAX = 61200;
const int NBXMAX = 3564;
float ene[NBXMAX][NCHMAX];

const int NBXFILLED = 2835;
int indexBXFilled[NBXFILLED];

// readout phase
double rophase = 1.0;

TRandom rnd;



void MakeLHCFillingScheme()
{
  // 3564 = 12 x 297 = ((81b + 8e) x 3 + 30e) x 11 + ((81b + 8e) x 2 + 119e)

  int idFilled = 0;
  int idAll    = 0;
  for(int itrain=0; itrain<11; itrain++){
    for(int igroup=0; igroup<3; igroup++){
      for(int ib=0; ib<81; ib++){
	indexBXFilled[idFilled] = idAll;
	idFilled++;
	idAll++;
      }
      for(int ib=0; ib<8; ib++){
	idAll++;
      }
    }
    for(int ib=0; ib<30; ib++){
      idAll++;
    }
  }
  for(int igroup=0; igroup<2; igroup++){
    for(int ib=0; ib<81; ib++){
      indexBXFilled[idFilled] = idAll;
      idFilled++;
      idAll++;
    }
    for(int ib=0; ib<8; ib++){
      idAll++;
    }
  }
  for(int ib=0; ib<119; ib++){
    idAll++;
  }
}



void MapChannels()
{
  for(int ich=0; ich<61200; ich++){
    int ieta = (ich / 360) - 85;
    if( ieta >=0 ) ieta++;
    int iphi = ich % 360 + 1;
    h2eta[ich] = ieta;
    h2phi[ich] = iphi;
  }
  for(int i=0; i<85; i++){
    double eta = 1.5 / 85. * (i + 0.5);
    sinTheta[i] = sin( atan( exp( -1. * eta ) ) * 2.0 );
    tZero[i] = 4.22839e+00/sinTheta[i] + 3.48572e-01;
  }  
}



void Load()
{
  tr = new TChain("T");
  tr->Add("pu200_minbias_xtals_01.root");
  tr->Add("pu200_minbias_xtals_02.root");
  tr->Add("pu200_minbias_xtals_03.root");
  tr->Add("pu200_minbias_xtals_04.root");

  tr->SetBranchAddress("nPU",&nPU,&b_nPU);
  tr->SetBranchAddress("nhits",&nhits,&b_nhits);
  tr->SetBranchAddress("hashedIndex",hashedIndex,&b_hashedIndex);
  tr->SetBranchAddress("eXtal",eXtal,&b_eXtal);
  tr->SetBranchAddress("tXtal",tXtal,&b_tXtal);
  tr->SetBranchAddress("eAPD1",eAPD1,&b_eAPD1);
  tr->SetBranchAddress("tAPD1",tAPD1,&b_tAPD1);
  tr->SetBranchAddress("eAPD2",eAPD2,&b_eAPD2);
  tr->SetBranchAddress("tAPD2",tAPD2,&b_tAPD2);
}



void LoadDecay(TString fname)
{
  trSig = new TChain("T");
  trSig->Add(fname.Data());

  trSig->SetBranchAddress("nhits",&nhitsSig,&b_nhitsSig);
  trSig->SetBranchAddress("hashedIndex",hashedIndexSig,&b_hashedIndexSig);
  trSig->SetBranchAddress("eXtal",eSig,&b_eSig);
}



void MakeTreeWaveforms(double tstep=6.25,
		       double timeMax =  1000.0,
		       int ietaMin = -85,
		       int ietaMax = +85,
		       int iphiMin =  1,
		       int iphiMax =  360,
		       int nOrbs   =  1,
		       TString fname = "output.root")
{
  MapChannels();
  MakeLHCFillingScheme();

  // Mask channels to process
  bool isGood[61200];
  for(int i=0; i<61200; i++){
    if(h2eta[i]>=ietaMin && h2eta[i]<=ietaMax && h2phi[i]>=iphiMin && h2phi[i]<=iphiMax){
      isGood[i] = true;
    }else{
      isGood[i] = false;
    }
  }

  Load();
  int nentries = (int)tr->GetEntries();
  
  Pulse *pulse = new Pulse(2);
  
  double timeMaxOrbit =  25.0 * NBXMAX;
  
  TFile *fout = new TFile(fname.Data(), "recreate");
  TTree *tr_out = new TTree("WF", "");
  double amp[61200];
  double timeNow;
  tr_out->Branch("timeNow", &timeNow, "timeNow/D");
  tr_out->Branch("amp",     amp,      "amp[61200]/D");
  
  int ientry = 0;
  for(int iorb=0; iorb<nOrbs; iorb++){

    for(int ibx=0; ibx<NBXMAX; ibx++){
      for(int ich=0; ich<NCHMAX; ich++){
	ene[ibx][ich] = 0;
      }
    }
    
    // Randomly populate BXs in this orbit with MinBias events from
    // the file
    
    bool isEmpty[NBXFILLED];
    for(int i=0; i<NBXFILLED; i++){
      isEmpty[i] = true;
    }
    int nBXEmpty = NBXFILLED;
    int debug[NBXMAX];
    for(int i=0; i<NBXMAX; i++){
      debug[i] = -1;
    }
    
    while(nBXEmpty>0){

      if(ientry>=nentries)
	ientry = 0;
      tr->GetEntry(ientry);
      ientry++;
      
      int id = (int) nBXEmpty * rnd.Rndm();
      int counter = 0;
      for(int i=0; i<NBXFILLED; i++){
	if(isEmpty[i]){
	  if(counter==id){
	    for(int ih=0; ih<nhits; ih++){
	      ene[indexBXFilled[i]][hashedIndex[ih]] = eXtal[ih];
	    }
	    debug[indexBXFilled[i]] = ientry;
	    nBXEmpty--;
	    isEmpty[i] = false;
	    break;
	    if(nBXEmpty<=0) break;
	  }
	  counter++;
	}
      }
    }
    cout << " Finished orbit " << iorb << endl;
    
    // Start making waveforms

    timeNow = -10.0;
    while(timeNow< min(timeMax, timeMaxOrbit) ){
    
      for(int ich=0; ich<NCHMAX; ich++){
	amp[ich] = 0;
	if(isGood[ich]){
	  for(int ibx=0; ibx<NBXMAX; ibx++){
	    double timeRelative = timeNow - 25.0 * ibx;
	    if(timeRelative > timeMax) continue;
	    if(timeRelative < 0 ) break;
	    amp[ich] += ene[ibx][ich] * pulse->Value(timeRelative);
	  }
	}
      }
      cout << " Time within orbit " << iorb+1  << " of " << nOrbs << " is " << timeNow << " ns" << endl;
      tr_out->Fill();      
      timeNow += tstep;
    }
  }
  fout->cd();
  tr_out->Write();
  fout->Close();
}



void MakeTreeWaveformsWithSignal(double tstep=6.25,
				 double timeMax =  1000.0,
				 int ietaMin = -85,
				 int ietaMax = +85,
				 int iphiMin =  1,
				 int iphiMax =  360,
				 int nOrbs   =  1,
				 double eneSignal = 0.,
				 int firstBXsignal = 10,
				 int stepBXsignal  = 20,
				 TString fname = "output.root")
{
  MapChannels();
  MakeLHCFillingScheme();

  // Mask channels to process
  bool isGood[61200];
  for(int i=0; i<61200; i++){
    if(h2eta[i]>=ietaMin && h2eta[i]<=ietaMax && h2phi[i]>=iphiMin && h2phi[i]<=iphiMax){
      isGood[i] = true;
    }else{
      isGood[i] = false;
    }
  }

  Load();
  int nentries = (int)tr->GetEntries();
  cout << " nentries " << nentries << endl;
  
  Pulse *pulse = new Pulse(2);
  
  double timeMaxOrbit =  25.0 * NBXMAX;
  
  TFile *fout = new TFile(fname.Data(), "recreate");
  TTree *tr_out = new TTree("WF", "");
  double amp[61200];
  double timeNow = -10.0;
  tr_out->Branch("timeNow", &timeNow, "timeNow/D");
  tr_out->Branch("amp",     amp,      "amp[61200]/D");

  
  int ientry = 0;
  for(int iorb=0; iorb<nOrbs; iorb++){

    for(int ibx=0; ibx<NBXMAX; ibx++){
      for(int ich=0; ich<NCHMAX; ich++){
	ene[ibx][ich] = 0;
      }
    }
    
    // Fill with signal energies
    
    for(int i=0; i<NBXFILLED; i++){
      if( i>=firstBXsignal && (i-firstBXsignal)%stepBXsignal==0 ){
	for(int ich=0; ich<NCHMAX; ich++){
	  if(isGood[ich]){
	    ene[indexBXFilled[i]][ich] = eneSignal;
	  }
	}
      }
    }
    
    // Randomly populate BXs in this orbit with MinBias events from
    // the file
    
    bool isEmpty[NBXFILLED];
    for(int i=0; i<NBXFILLED; i++){
      isEmpty[i] = true;
    }
    int nBXEmpty = NBXFILLED;
    int debug[NBXMAX];
    for(int i=0; i<NBXMAX; i++){
      debug[i] = -1;
    }
    
    while(nBXEmpty>0){

      if(ientry>=nentries)
	ientry = 0;
      tr->GetEntry(ientry);
      ientry++;

      int id = (int) nBXEmpty * rnd.Rndm();
      int counter = 0;
      for(int i=0; i<NBXFILLED; i++){
	if(isEmpty[i]){
	  if(counter==id){
	    for(int ih=0; ih<nhits; ih++){
	      ene[indexBXFilled[i]][hashedIndex[ih]] += eXtal[ih];
	    }
	    debug[indexBXFilled[i]] = ientry;
	    nBXEmpty--;
	    isEmpty[i] = false;
	    break;
	    if(nBXEmpty<=0) break;
	  }
	  counter++;
	}
      }
    }
    cout << " Finished orbit " << iorb << endl;

    // Start making waveforms

    while(timeNow< min(timeMax, (iorb+1) * timeMaxOrbit) ){
    
      for(int ich=0; ich<NCHMAX; ich++){
	amp[ich] = 0;
	if(isGood[ich]){
	  for(int ibx=0; ibx<NBXMAX; ibx++){
	    double timeRelative = timeNow - timeMaxOrbit * iorb - 25.0 * ibx;
	    if(timeRelative > timeMax) continue;
	    if(timeRelative < 0 ) break;
	    if(timeRelative > 3000.) continue;  // hard cut on 120 bx, ignore long tales
	    amp[ich] += ene[ibx][ich] * pulse->Value(timeRelative + rophase);
	  }
	}
      }
      cout << " Time within orbit " << iorb+1  << " of " << nOrbs << " is " << timeNow << " ns" << endl;
      tr_out->Fill();      
      timeNow += tstep;
    }
  }
  fout->cd();
  tr_out->Write();
  fout->Close();
}



void PlotWaveform(TString fname="output.root", int ieta=85, int iphi=1)
{
  TChain *trwf = new TChain("WF");
  trwf->Add(fname.Data());
  double amp[61200];
  double timeNow;
  TBranch *b_amp = 0;
  TBranch *b_timeNow = 0;
  trwf->SetBranchAddress("timeNow",&timeNow,&b_timeNow);
  trwf->SetBranchAddress("amp",amp,&b_amp);
  int nentries = (int)trwf->GetEntries();
  cout << " N entries " << nentries << endl;
  
  if(ieta<-85 || ieta==0 || ieta>85 || iphi<=0 || iphi>360){
    cout << " Wrong range of ieta or iphi " << ieta << " " << iphi << endl;
    return;
  }
  
  int ietaTmp = ieta;
  if(ietaTmp>0) ietaTmp--;
  ietaTmp += 85;
  int hashedIndex = 360 * ietaTmp + iphi - 1;
  cout << " plotting channel " << hashedIndex << endl;
  
  TGraph *gr = new TGraph();
  for(int ievt=0; ievt<nentries; ievt++){
    trwf->GetEntry(ievt);
    gr->SetPoint(ievt, timeNow, amp[hashedIndex]);
  }

  gr->Draw("APL");
    
}



void PlotPulse(int iopt=2)
{
  Pulse *pulse = new Pulse(iopt);
  TGraph *gr = new TGraph();
  for(int i=0; i<1000; i++){
    double t = 0.1 * i;
    double y = pulse->Value(t);
    gr->SetPoint( i, t, y);
  }

  gr->Draw("AL");
}



void PlotCellEnergy()
{
  MapChannels();
  Load();
  int nentries = (int)tr->GetEntries();

  TH1D *h = new TH1D("h", "", 1000, -2, 3);

  for(int ievt=0; ievt<nentries; ievt++){
    tr->GetEntry(ievt);
    for(int ih=0; ih<nhits; ih++){
      if(eXtal[ih]>0){
	h->Fill( log10(eXtal[ih]) );
      }
    }
  }

  h->Draw();  
}


void PlotNPU()
{
  MapChannels();
  Load();
  int nentries = (int)tr->GetEntries();

  TH1D *h = new TH1D("h", "", 500, -0.5, 499.7);

  for(int ievt=0; ievt<nentries; ievt++){
    tr->GetEntry(ievt);
    h->Fill(nPU);
  }

  h->Draw();  
}



void PlotFiveWaveform(TString fname="output.root", int ieta=85, int iphiMin=1, int iphiMax=5)
{
  TChain *trwf = new TChain("WF");
  trwf->Add(fname.Data());
  double amp[61200];
  double timeNow;
  TBranch *b_amp = 0;
  TBranch *b_timeNow = 0;
  trwf->SetBranchAddress("timeNow",&timeNow,&b_timeNow);
  trwf->SetBranchAddress("amp",amp,&b_amp);
  int nentries = (int)trwf->GetEntries();
  cout << " N entries " << nentries << endl;
  
  if(ieta<-85 || ieta==0 || ieta>85 || iphiMin<=0 || iphiMin>360 || iphiMax<=0 || iphiMax>360){
    cout << " Wrong range of ieta or iphi " << endl;
    return;
  }

  const int ngr = iphiMax - iphiMin + 1;
  TGraph *gr[ngr];
  for(int igr=0; igr<ngr; igr++){
    gr[igr] = new TGraph();
    gr[igr]->SetLineColor(igr+1);
  }
  
  int ietaTmp = ieta;
  if(ietaTmp>0) ietaTmp--;
  ietaTmp += 85;
  
  for(int ievt=0; ievt<nentries; ievt++){
    trwf->GetEntry(ievt);
    for(int iphi = iphiMin; iphi<=iphiMax; iphi++){
      int hashedIndex = 360 * ietaTmp + iphi - 1;	
      gr[iphi-iphiMin]->SetPoint(ievt, timeNow, amp[hashedIndex]);
    }
  }

  gr[0]->SetMinimum(-1.0);
  gr[0]->SetMaximum(10.0);
  gr[0]->Draw("AL");
  for(int igr=1; igr<ngr; igr++){
    gr[igr]->Draw("L");
  }
  TLegend *leg = new TLegend( 0.1, 0.7, 0.48, 0.9 );
  leg->SetFillStyle(1001);
  leg->SetTextFont(42);
  leg->AddEntry(gr[0],"ieta=85, iphi=1","l");
  leg->AddEntry(gr[1],"ieta=85, iphi=2","l");
  leg->AddEntry(gr[2],"ieta=85, iphi=3","l");
  leg->AddEntry(gr[3],"ieta=85, iphi=4","l");
  leg->AddEntry(gr[4],"ieta=85, iphi=5","l");
  leg->Draw(); 
   
    
}



void MakeTreeWaveformsDecayMode(double tstep=6.25,
				double timeMax =  1000.0,
				int nOrbs   =  1,
				int firstBXsignal = 10,
				int stepBXsignal  = 20,
				TString fnameDecay = "hgg_1000evt_pu200.root",
				TString fnameOut = "output.root")
{
  
  MapChannels();
  MakeLHCFillingScheme();

  // Mask channels to process: ALL channels should be processed for Signal Decay mode
  bool isGood[61200];
  for(int i=0; i<61200; i++){
      isGood[i] = true;
  }
  /*
  int ietaMin = -85;
  int ietaMax = +85;
  int iphiMin =  1;
  int iphiMax =  360;
  for(int i=0; i<61200; i++){
    if(h2eta[i]>=ietaMin && h2eta[i]<=ietaMax && h2phi[i]>=iphiMin && h2phi[i]<=iphiMax){
      isGood[i] = true;
    }else{
      isGood[i] = false;
    }
  }
  */
  
  Load();
  int nentries = (int)tr->GetEntries();
  cout << " PU nentries " << nentries << endl;

  LoadDecay(fnameDecay);
  int nentriesSig = (int)trSig->GetEntries();
  cout << " Signal nentries " << nentriesSig << endl;
  int counterSig = 0;
  
  Pulse *pulse = new Pulse(2);
  
  double timeMaxOrbit =  25.0 * NBXMAX;
  
  TFile *fout = new TFile(fnameOut.Data(), "recreate");
  TTree *tr_out = new TTree("WF", "");
  double amp[61200];
  double timeNow = -10.0;
  tr_out->Branch("timeNow", &timeNow, "timeNow/D");
  tr_out->Branch("amp",     amp,      "amp[61200]/D");

  
  int ientry = 0;
  for(int iorb=0; iorb<nOrbs; iorb++){

    cout << " Initializing energies for orbit " << iorb+1 << " of " << nOrbs << endl;
    for(int ibx=0; ibx<NBXMAX; ibx++){
      for(int ich=0; ich<NCHMAX; ich++){
	ene[ibx][ich] = 0;
      }
    }
    
    
    // Randomly populate BXs in this orbit with MinBias events from
    // the file
    
    bool isEmpty[NBXFILLED];
    for(int i=0; i<NBXFILLED; i++){
      isEmpty[i] = true;
    }
    int nBXEmpty = NBXFILLED;
    int debug[NBXMAX];
    for(int i=0; i<NBXMAX; i++){
      debug[i] = -1;
    }
    
    while(nBXEmpty>0){

      if(ientry>=nentries)
	ientry = 0;
      tr->GetEntry(ientry);
      ientry++;

      int id = (int) nBXEmpty * rnd.Rndm();
      int counter = 0;
      for(int i=0; i<NBXFILLED; i++){
	if(isEmpty[i]){
	  if(counter==id){

	    cout << " Filling PU energies in BX " << indexBXFilled[i]
		 << " Remaining empty BXs: " << nBXEmpty << endl;
	    
	    for(int ih=0; ih<nhits; ih++){
	      ene[indexBXFilled[i]][hashedIndex[ih]] += eXtal[ih];
	    }
	    debug[indexBXFilled[i]] = ientry;
	    nBXEmpty--;
	    isEmpty[i] = false;
	    break;
	    if(nBXEmpty<=0) break;
	  }
	  counter++;
	}
      }
    }
    
    // Fill with signal energies
    
    for(int i=0; i<NBXFILLED; i++){
      if( i>=firstBXsignal && (i-firstBXsignal)%stepBXsignal==0 ){
	if(counterSig >= nentriesSig){
	  counterSig = 0;
	}
	trSig->GetEntry(counterSig);
	counterSig++;
	   
	cout << " Filling Signal energies in BX " << indexBXFilled[i] << " ";

	double emax = 0;
	int    idmax = 0;
	for(int ih=0; ih<nhitsSig; ih++){
	  if(isGood[hashedIndexSig[ih]]){	    
	    // Overwrite (not add) if signal MC includes PU
	    ene[indexBXFilled[i]][hashedIndexSig[ih]] = eSig[ih];
	    if(eSig[ih] >= emax){
	      emax = eSig[ih];
	      idmax = hashedIndexSig[ih];
	    }
	  }
	}
	cout << "  Emax hit: hashedIndex=" << idmax << " energy=" << emax << endl;
      }
    }

    cout << " Finished orbit " << iorb << endl;

    // Start making waveforms

    while(timeNow< min(timeMax, (iorb+1) * timeMaxOrbit) ){
    
      for(int ich=0; ich<NCHMAX; ich++){
	amp[ich] = 0;
	if(isGood[ich]){
	  for(int ibx=0; ibx<NBXMAX; ibx++){
	    double timeRelative = timeNow - timeMaxOrbit * iorb - 25.0 * ibx;
	    if(timeRelative > timeMax) continue;
	    if(timeRelative < 0 ) break;
	    if(timeRelative > 3000.) continue;  // hard cut on 120 bx, ignore long tales
	    amp[ich] += ene[ibx][ich] * pulse->Value(timeRelative + rophase);
	  }
	}
      }
      cout << " Time within orbit " << iorb+1  << " of " << nOrbs << " is " << timeNow << " ns" << endl;
      tr_out->Fill();      
      timeNow += tstep;
    }
  }
  fout->cd();
  tr_out->Write();
  fout->Close();
}
