import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_med_fn(L, T, repeat, shi=1.0, nm='sqT'):
    fname = "MED"
    fname += "_L" + str(L) 
    fname += "_T%.1f" % T
    fname += "_shi%.1f" % shi
    fname += "_sHi%.2f" % np.sqrt(T)
    fname += "_r%d" % repeat
    fname += "_" + nm + ".csv"
    return(fname)
    


