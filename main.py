# MD code in microcanonical ensemble
# force field = LJ with eps = alpha * kb * T and r0 = 15 angstrom
# the treatment of the discontinuity at the cutoff is not implemented, i.e. the potential does not go smoothly to 0 at r = cutoff
#import numpy as np
import settings
import initialize
import force
import update
import debug
import sys
import time
import misc
from numba import njit, prange
from Simulation import SimulationSheet7a
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit 


def main(subtask):

    
    settings.init()

    if subtask == "a":
        print('a')

        SimulationSheet7a(True, '7a_', thermostat_type='berendsen')

if __name__ == "__main__":    
    main(subtask="a")
        

