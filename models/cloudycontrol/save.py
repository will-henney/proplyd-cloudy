"""
Save commands for cloudy model
"""

def default():
    return """\
* Output options
set save hash "return"
set save flush
save last overview ".ovr"
save last physical conditions ".phy"
save last grain charge ".grc"
save last grain drift velocity ".grv"
save last grain temperature ".grt"
save last PDR ".pdr"
save last element oxygen ".ion_o"
save last element nitrogen ".ion_n"
save last element sulphur ".ion_s"
save last element carbon ".ion_c"
save last element neon ".ion_ne"
save last element argon ".ion_ar"
save last element chlorine ".ion_cl"
save last element silicon ".ion_si"
save last dr ".dr"
save last continuum ".cont"
save last pressure ".pre"
save last heat ".heat" 
save last cool ".cool" 
save last lines, emissivity ".em"
* Hydrogen Balmer lines
H  1  4861A  // H beta
H  1  6563A  // H alpha
Ca B  3704A  // H16 3703.86
Ca B  3712A  // H15 3711.97
Ca B  3734A  // H13 3734.37
Ca B  3750A  // H12 3750.15
Ca B  3771A  // H11 3770.63
H  1  3835A  // H9  3835.39
H  1  4102A  // H6  4101.74
H  1  4340A  // H gama 4340.47
Ca B  8467A  // P17 8467.25
Ca B  8502A  // P16 8502.48
Ca B  8665A  // P13 8665.02
Ca B  8750A  // P12 8750.47
Ca B  8863A  // P11 8862.79
H  1  9015A  // P10 9014.91
H  1  9229A  // P9  9229.01
H  1  9546A  // P8  9545.97
H  1 1.005m  // P7  10049.4
* Helium recombination lines (WARNING: some may be optically thick)
He 1  7281A  // 7281.35
He 1  7065A  // 7065.28
He 1  6678A  // 6678.15
He 1  5876A  // 5875.64
He 1  5016A  // 5015.68
He 1  4922A  // 4958.91
He 1  4471A
He 1  3889A
He 1  3188A
He 1  2945A  // 2946 in STIS identification
He 1  2829A  // 2830 in STIS identification
* Collisionally excited optical diagnostic lines
N  2  6548A  // [N II] 6548.03
N  2  6584A  // [N II] nebular line 6583.41
N  2  5755A  // [N II] auroral line 5754.64
TOTL  4363A  // [O III] auroral line
O  3  4959A  // [O III] nebular line
O  3  5007A  // ditto
S II  6731A  // [S II] nebular line 
S II  6716A  // ditto
S II  4070A  // [S II] auroral line (4068A in Mesa-Delgado)
S II  4078A  // [S II] auroral line (4076A in Mesa-Delgado)
S  3  3722A  // [S III] 3721.83 same upper level as 6312
S  3  6312A  // [S III] 6312.10
S  3  9069A  // [SIII] 9068.90
S  3  9532A  // [SIII] 9530.60
O  1  6300A  // [O I] 6300.3
O  1  5577A  // [O I]
O II  3729A  // [O II] 
O II  3726A  // [O II]
O II  7323A  // [O II] (7320 in Tsamis)
O II  7332A  // [O II] (7330 in Tsamis)
Ne 3  3869A  // [Ne III] 
Ne 3  3968A  // [Ne III]
Ar 3  5192A  // [Ar III]
Ar 3  7135A  // [Ar III] 7135.78
Ar 4  4740A  // [Ar IV]
Ar 4  4711A  // [Ar IV]
Cl 3  5518A  // [Cl III]
Cl 3  5538A  // [Cl III] given as 5539A in Tsamis et al 2011
* Metal optical recombination lines
6lev  8446A  // O I 8446 from six level atom
O 2r  4651A   // O II 4651 total recombination, 4638.86-4696.35 (8 lines) 
C  2  4267A   // C II recombination line
TOTL  6580A   // C II recombination 6578 + 6580 A
C 2r  6580A   // C II recombination part 6580A
C 2p  6580A   // C II pump part 6580A
C  2  6580A   // C II we do not exactly who Cloudy calculate this
* Fluorescent lines
TOTL  5199A  // [N I] line
* Iron lines for comparison with Mesa-Delgado
Fe 3  5271A
Fe 3  4988A
Fe 3  4881A  // 4881.00
Fe 3  4755A  // 4754.83
Fe 3  4734A  // 4733.93
Fe 3  4702A  // 4701.62
Fe 3  4659A  // 4658.10 Strongest line
Fe 3  4608A  // 4607.13 
* NUV lines
TOTL  1750A  // N III] 
N  2  2141A  // N 2 intercombination line (as Cloudy said). Not sure if it is the same as N III 2143 + 2144
TOTL  2326A  // C II] Semi-forbidden 2324 - 2329 multiplet (5 lines)
C  3 1910A   // C III] Semi-forbidden
C  3 1907A   // [C III] forbidden
O II  2471A  // [O II] forbidden (2470 in Tsamis)
Mg 2  2796A  //
Mg 2  2803A  // 
* FIR lines, just in case
S II 1.029m
S II 1.032m
S II 1.034m
Ne 2 12.81m
Ne 3 15.55m
end of lines
"""
