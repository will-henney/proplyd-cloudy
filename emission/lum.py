import numpy as np

K_list=[5,6,7,8,9,10,11,12,13,14]

NK = len(K_list)

for k in range(len(K_list)):
    print k, K_list[k], k+1, K_list[k+1];
    print K_list[k+1], K_list[k-1];
    DPhi = (1./2.)*(K_list[k+1]-K_list[k-1]);
    print DPhi
