$rungauss
%chk=checkpoint.chk
%Mem=4Gb
%Nproc=1 
#P BHandHLYP/6-31++g** SCF=XQC guess=(Read,TCheck)

comment

1 2
 O -5.2187669866104937e-02 0.0000000000000000e+00 -4.0773472597237397e-02
 H -5.2753809927809404e-02 0.0000000000000000e+00 9.2125442922104372e-01
 H 8.8112705653754120e-01 0.0000000000000000e+00 -2.7405837778623499e-01
 
--Link1--
%chk=checkpoint.chk
%Mem=4Gb
%Nproc=1 
#P BHandHLYP/6-31++g** SCF=XQC FORCE TD=(nstates=2,root=2,conver=5) Geom=Check Guess=Read

comment

1 2

