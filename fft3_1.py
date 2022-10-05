import numpy as np

def FFT(data):
    kpower = np.fft.fftn(data)

    kpower_2 = np.sum(kpower, axis=2)
    kpower_21 = np.sum(kpower_2, axis=1)
    kpower0 = np.abs(kpower_21)

    kpower_0 = np.sum(kpower, axis=0)
    kpower_01 = np.sum(kpower_0, axis=1)
    kpower1 = np.abs(kpower_01)

    kpower_00 = np.sum(kpower_0, axis=0)
    kpower2 = np.abs(kpower_00)

    kpower1dim = (kpower0 + kpower1 + kpower2) / 3
    return kpower1dim
#%%
