$rungauss
%chk=checkpoint.chk
%Mem=4Gb
%Nproc=1 
#P BHandHLYP/6-31g** SCF=XQC guess=(Read,TCheck)

comment

0 1
 O -5.2727618512262896e-02 0.0000000000000000e+00 -4.1195326882842283e-02
 H -5.3049954482923364e-02 0.0000000000000000e+00 9.3046644065869999e-01
 H 8.8999378901443715e-01 0.0000000000000000e+00 -2.7657430938243982e-01
 
--Link1--
%chk=checkpoint.chk
%Mem=4Gb
%Nproc=1 
#P BHandHLYP/6-31g** SCF=XQC FORCE TD=(nstates=2,root=2,conver=5) Geom=Check Guess=Read

comment

0 1

