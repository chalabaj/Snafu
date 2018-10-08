$rungauss
%chk=checkpoint.chk
%Mem=Mb
%Nproc=1 
#P BLYP/6-31g** SCF=XQC guess=(Read,TCheck)

comment

0 1
 
--Link1--
%chk=checkpoint.chk
%Mem=Mb
%Nproc=1 
#P BLYP/6-31g** SCF=XQC FORCE TD=(nstates=2,root=2) Geom=Check Guess=Read

comment

0 1

