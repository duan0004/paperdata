import os
import numpy as np
import ptarcade.models_utils as aux

name = 'meta_l' 

smbhb = False

parameters = {
    'log10_Gmu' : aux.prior("Uniform", -14, -1.5),
    'sqrt_kappa' : aux.prior("Uniform", 7.0, 9.5)
    }

cwd = os.path.dirname(os.path.abspath(__file__))
log_spectrum = aux.spec_importer(cwd +'/models_data/meta_l.h5')

def spectrum(f, log10_Gmu, sqrt_kappa):
    return 10**log_spectrum(np.log10(f), log10_Gmu=log10_Gmu, sqrt_kappa=sqrt_kappa)
