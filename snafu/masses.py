# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 11:29:42 2018

@author: chalabaj
"""

    
def assign_masses(at_names):
    """
    Atomic weights of the elements 2013 (IUPAC Technical Report) Meija, Juris; et al.
    http://doi.org/10.1515/pac-2015-0305, Pure and Applied Chemistry, 88, 3, pp. 265-291, 2016-02-24
    """
    
    atomic_masses = {    'H'    :  1.00794,	
                         'D'		:	2.01410177785,                   	
                         'He'		:	4.002602	,
                         'Li'		:	6.93800,	
                         'Be'		:	9.0121831	,
                         'B'		:	10.811	,
                         'C'		:	12.011	,
                         'N'		:	14.00674	,
                         'O'		:	15.999	,
                        'F'		:	18.9984032	,
                        'Ne'		:	20.1797	,
                        'Na'		:	22.989768	,
                        'Mg'		:	24.305	,
                        'Al'		:	26.981539	,
                        'Si'		:	28.0855	,
                        'P'		:	30.973762	,
                        'S'		:	32.066	,
                        'Cl'		:	35.4527	,
                        'Ar'		:	39.948	,
                        'K'		:	39.0983	,
                        'Ca'		:	40.078	,
                        'Sc'		:	44.95591	,
                        'Ti'		:	47.88	,
                        'V'		:	50.9415	,
                        'Cr'		:	51.9961	,
                        'Mn'		:	54.93805	,
                        'Fe'		:	55.847	,
                        'Co'		:	58.9332	,
                        'Ni'		:	58.6934	,
                        'Cu'		:	63.546	,
                        'Zn'		:	65.39	,
                        'Ga'		:	69.723	,
                        'Ge'		:	72.61	,
                        'As'		:	74.92159	,
                        'Se'		:	78.96	,
                        'Br'		:	79.904	,
                        'Kr'		:	83.8	,
                        'Rb'		:	85.4678	,
                        'Sr'		:	87.62	,
                        'Y'		  :	88.90585	,
                        'Zr'		:	91.224	,
                        'Nb'		:	92.90638	,
                        'Mo'		:	95.94	,
                        'Tc'		:	-97.9072	,
                        'Ru'		:	101.07	,
                        'Rh'		:	102.9055	,
                        'Pd'	   :	106.42	,
                        'Ag'	   :	107.8682	,
                        'Cd'	   :	112.411	,
                        'In'   	:	114.818	,
                        'Sn'   	:	118.71	,
                        'Sb'	   :	121.757	,
                        'Te'		:	127.6	,
                        'I'      :	126.90447	,
                        'Xe'		:	131.29	,
                        'Cs'		:	132.90543	,
                        'Ba'		:	137.327	}

    mass = [ 0.0 for x in range(0,len(at_names))]  
    
    for iat in range(0,len(at_names)):
       mass[iat] = (atomic_masses[at_names[iat]]) 
    
    return(mass)