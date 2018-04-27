# PileUpWaveForms
Generate Waveforms for PU events

## Installation
To install ROOT macros
````
cd NewDirectory
git init
git clone https://github.com/cms-eb-upgrade/PileUpWaveForms.git
cd PileUpWaveForms
root -l
.L WaveformsPU200.C+
.q
````
There should be no errors during compile.
Also, you need four files with MinBias info from CMSSW MC
````
eos ls /eos/cms/store/user/ledovsk/minbias200
pu200_minbias_xtals_01.root
pu200_minbias_xtals_02.root
pu200_minbias_xtals_03.root
pu200_minbias_xtals_04.root
````
Each file is about 160 MB

## Signal Pulse shape

The default pulse shape is from Marc, obtained from run 9862 at TB-H4 in October 2017. 
The hit crystal is C3@150 GeV

To plot it
````
root -l
.L WaveformsPU200.C+
PlotPulse()
````

## LHC filling scheme

The default LHC filling scheme is taken from CERN-AC-95-05-LHC
````
((81b + 8e) x 3 + 30e) x 11 + ((81b + 8e) x 2 + 119e)
````
It can be modified in ````MakeLHCFillingScheme()````

## Generate PileUP waveform in one ECAL channel

This is a waveform from pileup only with ````npu=200````, no signal injected.
````
root -l
.L WaveformPU200.C+
MakeTreeWaveforms(0.2, 1000.0, 85, 85, 1, 1, 1, "180410_wf_01.root")
````
Will generate a sequence of amplitude samples with
````0.2```` ns steps, during first 
````1000.0```` ns of LHC orbit, for channels with ieta from
 ````85```` to ````85````, and iphi from 
 ````1```` to ````1````, in number of LHC orbits ````1````, 
 and save TTree in output file ````180410_wf_01.root````

One can plot this waveform with
````
PlotWaveform("180410_wf_01.root",85, 1)
````

## Generate Signal on top of Pileup

One can inject a signal with fixed amplituded periodically with chosen frequency. For example,
````
MakeTreeWaveformsWithSignal( 0.2, 1000.0, 85, 85, 1, 1, 1, 100, 10, 20, "180410_wf_01.root")
````
will add ````100```` GeV pulses in every chosen channel, starting with BX ````10````, 
every ````20```` filled BX 

## Complete orbits

One generage stream of samples during several LHC orbits.
For example, to generage 160 MHz samples for 500 GeV signal pulses,
on top of PU=200,
in one channel ieta=85, iphi=1, 
during 2 orbits,
starting with BX=10, every 20 BXs
````
MakeTreeWaveformsWithSignal( 6.25, 1e+9, 85, 85, 1, 1, 2, 500, 10, 20, "180410_wf_01.root")
````

