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
save last element oxygen ".ion_O"
save last element nitrogen ".ion_N"
save last element sulphur ".ion_S"
save last element carbon ".ion_C"
save last dr ".dr"
save last continuum ".cont"
save last pressure ".pre"
save last heat ".heat" 
save last cool ".cool" 
save last lines, emissivity ".em"
H  1  4861A 
H  1  6563A 
TOTr  5199A 
N  2  6584A 
O  1  6300A 
O II  3726A 
TOTL  4363A 
O  3  4959A 
O  3  5007A 
S II  6731A 
S II  6716A 
S  3  6312A
He 1  5876A
Ne 2 12.81m
Ne 3 15.55m
O  2 4651A
C  3 1910A
C  3 1907A
end of lines
"""
