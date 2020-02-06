# Standalone Glauber Model from arXiv:0805.4411

## Setting up ROOT at ACCRE

```
We will use ROOT in CMSSW at ACCRE.
In any directory at ACCRE, do the following (find available release by doing 'scram list'):
[tuos@gw345 glauber]$ cmsrel CMSSW_10_3_0
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
[tuos@gw343 glauber_v1p5]$ root -l
root [0] .L runglauber_v1.5.C+
Info in <TUnixSystem::ACLiC>: creating shared library /gpfs23/scratch/tuos/PbPb2015/glauber/glauber_v1p5/./runglauber_v1.5_C.so
root [1] runAndSaveNtuple(10000,"Pb","Pb",67.6,0.4,"glauber_PbPb_default_v1p5_10k.root")
Setting up nucleus Pb
Setting up nucleus Pb
Event # 9950 x-sect = 7.69905 +- 0.0771798 b        
Done!
root [2] .q
```

### Analyzing the output
```
cd macros/
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
```

### Changing collision energy

```
Find the input value of nucleon-nucleon inelastic cross section (67.6 mb in the above example) in table 2 of the paper:
 
https://arxiv.org/abs/1710.07098
```


