import os
import numpy as np
import ptarcade.models_utils as aux

name = 'stable_n'

smbhb = False

parameters = {
    'log10_Gmu' : aux.prior("Uniform", -14, -6),
    }

group = ['log10_Gmu']   

cwd = os.path.dirname(os.path.abspath(__file__))
log_spectrum = aux.spec_importer(cwd +'/models_data/stable_n.h5')

def spectrum(f, log10_Gmu):
    return 10**log_spectrum(np.log10(f), log10_Gmu=log10_Gmu)
