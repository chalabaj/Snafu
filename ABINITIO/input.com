$rungauss
%chk=/home/chalabaj/snafu/ABINITIO/checkpoint.chk
%Mem=500Mb
%Nproc=2
#BLYP/6-31g* Force  guess=(Read,TCheck)

comment

0 1
 O -5.219047819e-02 -5.206365822e-02 -5.216273858e-02
 H -5.297947680e-02 -5.206365822e-02 9.102021387e-01
 H 8.813972998e-01 -5.206365822e-02 -2.857454041e-01
 
