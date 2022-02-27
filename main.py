import numpy as np
import p_minhash_0
from collections import Counter


def prob_minhash_0(W, upper,z=[100000000]*64, m=int(64)):
    W = np.array(W, dtype=np.int32, order='F')
    z = np.array(z, dtype=np.int32, order='F')
    n = int(W.shape[1])
    p_minhash_0.p_minhash_0(W, n=n, m=m, upper=upper, z=z)
    sig_counter = Counter()
    sig_counter.update(z)
    return z, sig_counter.most_common(2)[-1][0]

if __name__ == "__main__":
    W1 = np.array([[123,321,888,9,100],[1,3,10,2,1]])
    z, sig = prob_minhash_0(W=W1, upper=3000, z=[3000]*64)
    print(z)
    print(sig)
