import numpy as np

K_list=[5,6,7,8,9,10,11,12,13,14]

NK = len(K_list)

for k in range(NK):
    kneg = max(0, k - 1)
    kpos = min(NK-1, k + 1)
    print "K_list[kpos], K_list[kneg]", K_list[kpos], K_list[kneg]
    DPhi = (1./2.)*(K_list[kpos]-K_list[kneg])
    print "Dphi = ", DPhi
