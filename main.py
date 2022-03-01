import numpy as np
import p_minhash
#import p_minhash_0
#import p_minhash_1
#import p_minhash_2
#import p_minhash_3
#import p_minhash_4
from collections import Counter
from functools import reduce


def prob_minhash_0(W, upper,z=[100000000]*64, m=int(64)):
    W = np.array(W, dtype=np.int32, order='F')
    z = np.array(z, dtype=np.int32, order='F')
    n = int(W.shape[1])
    p_minhash.p_minhash_0(W, n=n, m=m, upper=upper, z=z)
    #sig_counter = Counter()
    #sig_counter.update(z)
    return z, reduce(lambda x,y: x^y, z)


def prob_minhash_1(W, upper,z=[100000000]*64, m=int(64)):
    W = np.array(W, dtype=np.int32, order='F')
    z = np.array(z, dtype=np.int32, order='F')
    n = int(W.shape[1])
    p_minhash.p_minhash_1(W, n=n, m=m, upper=upper, z=z)
    #sig_counter = Counter()
    #sig_counter.update(z)
    return z, reduce(lambda x,y: x^y, z)


def prob_minhash_2(W, upper,z=[100000000]*64, m=int(64)):
    W = np.array(W, dtype=np.int32, order='F')
    z = np.array(z, dtype=np.int32, order='F')
    n = int(W.shape[1])
    p_minhash.p_minhash_2(W, n=n, m=m, upper=upper, z=z)
    #sig_counter = Counter()
    #sig_counter.update(z)
    return z, reduce(lambda x,y: x^y, z)


def prob_minhash_3(W, upper,z=[100000000]*64, m=int(64)):
    W = np.array(W, dtype=np.int32, order='F')
    z = np.array(z, dtype=np.int32, order='F')
    n = int(W.shape[1])
    p_minhash.p_minhash_3(W, n=n, m=m, upper=upper, z=z)
    #sig_counter = Counter()
    #sig_counter.update(z)
    return z, reduce(lambda x,y: x^y, z)


def prob_minhash_4(W, upper,z=[100000000]*64, m=int(64)):
    W = np.array(W, dtype=np.int32, order='F')
    z = np.array(z, dtype=np.int32, order='F')
    n = int(W.shape[1])
    p_minhash.p_minhash_4(W, n=n, m=m, upper=upper, z=z)
    #sig_counter = Counter()
    #sig_counter.update(z)
    return z, reduce(lambda x,y: x^y, z)



if __name__ == "__main__":
    for i in range(3):
      W1 = np.array([[123,321,888,9,100],[1,3,10,2,1]])
      z, sig = prob_minhash_0(W=W1, upper=3000, z=[3000]*400, m=int(400))
      print(sig)
      W1 = np.array([[123,321,888,9,100],[1,3,10,2,1]])
      z, sig = prob_minhash_1(W=W1, upper=3000, z=[3000]*400, m=int(400))
      print(sig)
      W1 = np.array([[123,321,888,9,100],[1,3,10,2,1]])
      z, sig = prob_minhash_2(W=W1, upper=3000, z=[3000]*400, m=int(400))
      print(sig)
      W1 = np.array([[123,321,888,9,100],[1,3,10,2,1]])
      z, sig = prob_minhash_3(W=W1, upper=3000, z=[3000]*400, m=int(400))
      print(sig)
      W1 = np.array([[123,321,888,9,100],[1,3,10,2,1]])
      z, sig = prob_minhash_3(W=W1, upper=3000, z=[3000]*400, m=int(400))
      print(sig)
      W1 = np.array([[123,321,888,9,100],[1,3,10,2,1]])
      z, sig = prob_minhash_4(W=W1, upper=3000, z=[3000]*400, m=int(400))
      print(sig)


