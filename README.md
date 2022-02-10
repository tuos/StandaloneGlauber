# Standalone Glauber Model from arXiv:0805.4411

## Setting up ROOT at ACCRE

```
We will use ROOT in CMSSW at ACCRE.
Login to ACCRE and go to a directory that you want to work on the Glauber model:
[tuos@gw344 class2022]$ mkdir glauber
[tuos@gw344 class2022]$ cd glauber/
[tuos@gw344 glauber]$ source /cvmfs/cms.cern.ch/cmsset_default.sh
[tuos@gw344 glauber]$ cmsrel CMSSW_10_3_0
[tuos@gw345 glauber]$ cd CMSSW_10_3_0/src/
[tuos@gw345 src]$ cmsenv
[tuos@gw345 src]$ root -l
root [0] .q
```
## Getting the code

```
git clone https://github.com/tuos/StandaloneGlauber.git 
```

## Running the code

```
[tuos@gw344 src]$ cd StandaloneGlauber/
[tuos@gw344 StandaloneGlauber]$ root -l
root [0] .L runglauber_v1.5.C+
Info in <TUnixSystem::ACLiC>: creating shared library /gpfs23/scratch/tuos/PbPb2015/glauber/glauber_v1p5/./runglauber_v1.5_C.so
root [1] runAndSaveNtuple(10000,"Pb","Pb",67.6,0.4,"glauber_PbPb_default_v1p5_10k.root")
Setting up nucleus Pb
Setting up nucleus Pb
Event # 9950 x-sect = 7.69905 +- 0.0771798 b        
Done!
root [2] .q
The output file "glauber_PbPb_default_v1p5_10k.root" is created.
```

### Analyzing the output
```
[tuos@gw344 StandaloneGlauber]$ cd macros/
Read the output ROOT file and draw figures:
[tuos@gw345 macros]$ root -l readFileAndDraw.C 
root [0] 
Processing readFileAndDraw.C...
<Npart> = 113.804
# of entries = 10000
Info in <TCanvas::Print>: file fig_npart_b.png has been created
root [1]
```
![Npart vs. B](https://github.com/tuos/StandaloneGlauber/blob/master/macros/fig_npart_b.png)
```
Read the tree inside the file event by event and fill histograms:
[tuos@gw345 macros]$ root -l analyzeEventTree.C 
root [0] 
Processing analyzeEventTree.C...
Have run 0 out of 10000 events; 
Have run 1000 out of 10000 events; 
Have run 2000 out of 10000 events; 
Have run 3000 out of 10000 events; 
Have run 4000 out of 10000 events; 
Have run 5000 out of 10000 events; 
Have run 6000 out of 10000 events; 
Have run 7000 out of 10000 events; 
Have run 8000 out of 10000 events; 
Have run 9000 out of 10000 events; 
<Npart> = 113.804
# of entries = 10000
root [1]

We can produce the 2D correlation figure for Npart vs. b in PbPb at 5.02 TeV.
b is the impact parameter for each collision.
```

### Changing collision energy

```
Find the input value of nucleon-nucleon inelastic cross section (67.6 mb in the above example) in table 2 of the paper:
 
https://arxiv.org/abs/1710.07098
```

### Homework

```
1. Produce the Npart vs. b 2D figure for Oxygen-Oxygen collisions at 7 TeV
2. Produce the b distribution for proton-Oxygen at 8.8 TeV
```



