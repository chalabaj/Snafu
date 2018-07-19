$rungauss
%chk=/home/chalabaj/snafu/ABINITIO/checkpoint.chk
%Mem=500Mb
%Nproc=2
#BLYP/6-31g* Force  guess=(Read,TCheck)

comment

0 1
 O -5.242174926e-02 -5.206365822e-02 -5.234342741e-02
 H -5.272937415e-02 -5.206365822e-02 9.136653615e-01
 H 8.848181556e-01 -5.206365822e-02 -2.863405586e-01
 
