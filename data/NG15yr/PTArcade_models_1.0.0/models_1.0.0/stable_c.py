import os
import numpy as np
import ptarcade.models_utils as aux

name = 'stable_c'

smbhb = False

parameters = {
    'log10_Gmu' : aux.prior("Uniform",-14, -6),
    }

cwd = os.path.dirname(os.path.abspath(__file__))
log_spectrum = aux.spec_importer(cwd +'/models_data/stable_c.h5')

def spectrum(f, log10_Gmu):
    return 10**log_spectrum(np.log10(f), log10_Gmu=log10_Gmu)
