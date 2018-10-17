$rungauss
%chk=checkpoint.chk
%Mem=4Gb
%Nproc=1 
#P BHandHLYP/6-31++g** SCF=XQC guess=(Read,TCheck)

comment

1 2
 O -5.2664343331159727e-02 0.0000000000000000e+00 -4.1145890875282726e-02
 H -5.0446454659512967e-02 0.0000000000000000e+00 9.2609900208748297e-01
 H 8.8638592422637674e-01 0.0000000000000000e+00 -2.7299156701873128e-01
 
--Link1--
%chk=checkpoint.chk
%Mem=4Gb
%Nproc=1 
#P BHandHLYP/6-31++g** SCF=XQC FORCE TD=(nstates=2,root=2,conver=5) Geom=Check Guess=Read

comment

1 2

