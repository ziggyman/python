import arrayfire as af
import numpy as np
n = 2000
A = np.random.rand(n,n.astype(np.float32)
A_af = af.np_to_af_array(A)