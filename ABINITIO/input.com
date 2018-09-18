%maxcore 8000          # memory per core
%pal 
   nprocs            # number of cores
end

%basis
Basis "LANL2DZ"        # The orbital expansion basis set
end

! ! ENERGY B3LYP RIJCOSX           # for small systems increase accuracy by: Grid5 FinalGrid6 
! AUTOSTART         # try to read from previous step
* xyz -3 6 
    3
FINAL HEAT OF FORMATION =   -76.328992
 O    -0.000064     0.000000    -0.000050
 H    -0.000979     0.000000     0.961319
 H     0.932527     0.000000    -0.233514
*

! ! ENGRAD B3LYP RIJCOSX
! moread
%basis
 Basis "LANL2DZ"        # The orbital expansion basis set
end
%tddft
  maxdim 5
  nroots 4
  maxcore 8000
  iroot 1  
end
* xyz -3 6 
    3
FINAL HEAT OF FORMATION =   -76.328992
 O    -0.000064     0.000000    -0.000050
 H    -0.000979     0.000000     0.961319
 H     0.932527     0.000000    -0.233514
