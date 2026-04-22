import os
import numpy as np
import ptarcade.models_utils as aux

name = 'super'

smbhb = False

parameters = {
    'log10_Gmu' : aux.prior("Uniform", -14, -6),
    'log10_P': aux.prior("Uniform", -3, 0)
    }

group = ['log10_Gmu', 'log10_P']   

cwd = os.path.dirname(os.path.abspath(__file__))
log_spectrum = aux.spec_importer(cwd +'/models_data/super.h5')

def spectrum(f, log10_Gmu, log10_P):
    return 10**log_spectrum(np.log10(f), log10_Gmu=log10_Gmu, log10_P=log10_P)
