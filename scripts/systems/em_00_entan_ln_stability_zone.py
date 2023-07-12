# dependencies
import numpy as np
import os
import sys

# qom modules
from qom.utils.loopers import run_loopers_in_parallel, wrap_looper
from qom.utils.solvers import get_func_quantum_correlation_measures, get_func_stability_zone

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('.')))
# import system
from systems.EndMirror import EM_00

# all parameters
params = {
    'looper': {
        'show_progress': True,
        'file_path_prefix': 'data/systems/em_00_en_sz',
        'X': {
            'var': 'Delta_0_norm',
            'min': -1.5,
            'max': -0.5,
            'dim': 1001
        },
        'Y': {
            'var': 'g_0_norm',
            'min': 0.002,
            'max': 0.005,
            'dim': 901
        }
    },
    'solver': {
        'show_progress': False,
        'method': 'eig',
        'measure_codes': ['entan_ln'],
        'indices': (0, 1),
        'system_measure_name': 'A'
    },
    'system': {
        'A_l_norm': 10.0,
        'Delta_0_norm': -1.0, 
        'g_0_norm': 0.005,
        'gamma_norm': 0.005,
        'kappa_norm': 0.15,
        'T_norm': 0.0
    },
    'plotter': {
        'type': 'contourf',
        'palette': 'Blues',
        'x_label': '$\\Delta_{0} / \\omega_{m}$',
        'x_ticks': [i * 0.25 - 1.5 for i in range(5)],
        'x_ticks_minor': [i * 0.125 - 1.5 for i in range(9)],
        'y_label': '$g_{0} / \\omega_{m} \\times 10^{3}$',
        'y_tick_labels': list(range(2, 6)),
        'y_ticks': [i * 0.001 + 0.002 for i in range(4)],
        'y_ticks_minor': [i * 0.0005 + 0.002 for i in range(7)],
        'show_cbar': True,
        'cbar_title': '$10^{2} \\times E_{N}$',
        'cbar_tick_labels': [0, 3, 6, 9],
        'cbar_ticks': [0.0, 0.03, 0.06, 0.09],
        'annotations': [{
            'text': 'multi\n-stable',
            'xy': (0.67, 0.73),
            'orientation': 'vertical'
        }]
    }
}

# function to obtain the steady state entanglement in stable zones
def func(system_params):
    # get stability zone
    sz = get_func_stability_zone(
        SystemClass=EM_00,
        params=params['solver'],
        steady_state=True,
        use_rhc=False
    )(system_params)[0]
    # get entanglement if stable
    if sz == 1:
        return get_func_quantum_correlation_measures(
            SystemClass=EM_00,
            params=params['solver'],
            steady_state=True
        )(system_params)[0, 0]
    
    return np.NaN

if __name__ == '__main__':
    # wrap loopers and plot
    looper = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func,
        params=params['looper'],
        params_system=params['system'],
        plot=True,
        params_plotter=params['plotter']
    )