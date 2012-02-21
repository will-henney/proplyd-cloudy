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
save last dr ".dr"
save last continuum ".cont"
save last pressure ".pre"
save last heat ".heat" 
save last cool ".cool" 
save last lines, emissivity ".em"
* Hydrogen Balmer lines
H  1  4861A  // H beta
H  1  6563A  // H alpha
* Helium recombination lines (WARNING: some may be optically thick)
He 1  6678A  // 6678.15
He 1  5876A  // 5875.64
He 1  4471A
* Collisionally excited optical diagnostic lines
N  2  6584A  // [N II] nebular line 6583.41
N  2  5755A  // [N II] auroral line 5754.64
TOTL  4363A  // [O III] auroral line
O  3  4959A  // [O III] nebular line
O  3  5007A  // ditto
S II  6731A  // [S II] nebular line 
S II  6716A  // ditto
S II  4070A  // [S II] auroral line (4068A in Mesa-Delgado)
S II  4078A  // [S II] auroral line (4076A in Mesa-Delgado)
S  3  6312A  // [S III] 6312.10
O  1  6300A  // [O I] 6300.3
O  1  5577A  // [O I]
O II  3729A  // [O II] 
O II  3726A  // [O II]
O II  7323A  // [O II] (7320 in Tsamis)
O II  7332A  // [O II] (7330 in Tsamis)
Ne 3  3869A  // [Ne III] 
Ar 3  5192A  // [Ar III]
Ar 3  7135A  // [Ar III] 7135.78
Ar 4  4740A  // [Ar IV]
Ar 4  4740A  // [Ar IV]
Ar 4  4711A  // [Ar IV]
Cl 3  5518A  // [Cl III]
Cl 3  5538A  // [Cl III] given as 5539A in Tsamis et al 2011
* Metal optical recombination lines
O 2r  4651A   // O II 4651 total recombination, 4638.86-4696.35 (8 lines) 
C  2  4267A   // C II recombination line
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
TOTL  2326A  // C II] Semi-forbidden 2324 - 2329 multiplet (5 lines)
C  3 1910A   // C III] Semi-forbidden
C  3 1907A   // [C III] forbidden
O II  2471A  // [O II] forbidden (2470 in Tsamis)
* FIR lines, just in case
Ne 2 12.81m
Ne 3 15.55m
end of lines
"""
