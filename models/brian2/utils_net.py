import numpy as np
# try to import numba
# or define dummy decorator
try:
    from numba import autojit
except:
    def autojit(func):
        return func

# util functions for network simulation
def smooth_trace(trace, scale):
    scale = int(scale)
    if scale == 1 or scale == 0:
        return trace
    slen = int(len(trace) / scale)
    if slen == 0:
        return trace
    return np.array([np.mean(trace[i*scale:(i+1)*scale]) for i in xrange(slen)])

@autojit
def choose_k_from_n(n, k):
    # use vaguely estimated metric of when sorting random numbers is better
    if float(k) / float(n) > 0.125:
        ans = np.argsort(np.random.rand(n))[:k]
        return ans
    nums = range(n)
    swaps = (np.random.rand(k) * xrange(n, n - k, -1)).astype('int') + xrange(k)
    for i in xrange(k):
        # swap with some random element from here to end - these swap positions precalculated
        nums[i], nums[swaps[i]] = nums[swaps[i]], nums[i]
    ans = nums[:k]
    return ans

def fixed_connectivity(n, k):
    prelist = np.zeros(k * n, dtype = int)
    postlist = np.zeros_like(prelist)
    for j in xrange(n):
        presynapses = choose_k_from_n(n, k)
        prelist[j * k:(j + 1) * k] = presynapses
        postlist[j * k:(j + 1) * k] = j * np.ones(k, dtype = int)
    return prelist, postlist