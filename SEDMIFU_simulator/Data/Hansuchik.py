import numpy as np

hc = 1.98644521e-8 # erg angstrom

# From Hanuschik
uves_sky = np.array([
    (5730 , -0.2857),
    (6546 , -0.2857),
    (6815 , -0.9870),
    (6860 , -1.064),
    (6893 , -1.084),
    (6960 , -1.084),
    (7005 , -1.084),
    (7061 , -0.9090),
    (7095 , -0.9090),
    (7151 , -0.8311),
    (7195 , -0.9870),
    (7262 , -1.006),
    (7330 , -1.064),
    (7386 , -1.064),
    (7497 , -0.7532),
    (7598 , -0.6363),
    (7632 , -0.7532),
    (7699 , -0.8311),
    (7755 , -0.8896),
    (7811 , -0.8506),
    (7855 , -0.8701),
    (7911 , -0.9870),
    (7967 , -0.9675),
    (8001 , -0.9090),
    (8012 , -0.6363),
    (8079 , -0.4999),
    (8169 , -0.4805),
    (8225 , -0.5389),
    (8292 , -0.5389),
    (8325 , -0.7337),
    (8393 , -0.8506),
    (8437 , -1.006),
    (8482 , -1.084),
    (8549 , -1.123),
    (8583 , -0.5194),
    (8650 , -0.4999),
    (8695 , -0.3051),
    (8751 , -0.3831),
    (8806 , -0.5779),
    (8851 , -0.7532),
    (8918 , -0.9675),
    (8952 , -0.4999),
    (9019 , -0.4025),
    (9075 , -1.006),
    (9097 , -0.4025),
    (9109 , -0.2272),
    (9165 , -0.2662),
    (9232 , -0.4025),
    (9276 , -0.3246),
    (9332 , -0.4805),
    (9400 , -0.5974),
    (9444 , -1.045),
    (9489 , -1.103),
    (9545 , -1.084),
    (9601 , -1.084),
    (9657 , -1.103),
    (9702 , -0.07142),
    (9746 , -0.1883),
    (9813 , -0.2662),
    (9869 , -0.2857),
    (9925 , -0.3441),
    (9970 , -0.4999),
    (10030 , -0.6948),
    (10070 , -0.4025),
    (10120 , -0.4025),
    (10180 , -0.1883),
    (10230 , -0.5779),
    (10290 , -0.2662),
    (10320 , 0.02597),
    (10390 , -0.1493),
    (10440 , 0.006493),
    (10490 , -0.1688),
    (10550 , -0.1688),
    (10600 , 0.1038),
    (10660 , -0.01298),
    (10700 , 0.2402),
    (10750 , -0.2272)
])

# Fix teh wavelength value of the UVES sky spectrum.
uves_sky[:,1] = 10**(uves_sky[:,1])*1e-17 # erg/sec/cm^2/Ang/Arc^2
uves_sky[:,1] /= (hc/uves_sky[:,0]) # #/sec/cm^2/Ang/Arc^2