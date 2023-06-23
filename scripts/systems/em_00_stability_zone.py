# dependencies
import os
import sys

# qom modules
from qom.utils.looper import run_loopers_in_parallel, wrap_looper
from qom.utils.solver import get_func_stability_zone

# add path to local libraries
sys.path.append(os.path.abspath(os.path.join('..', 'qom-examples')))
# import system
from systems.EndMirror import EM_00

# frequently used parameters
use_rhc=False

# all parameters
params = {
    'looper'    : {
        'show_progress'     : True,
        'file_path_prefix'  : 'data/em_00/sz_' + ('rhc' if use_rhc else 'eig'),
        'X' : {
            'var'   : 'Delta_0_norm',
            'min'   : -0.5,
            'max'   : 0.5,
            'dim'   : 1001,
        },
        'Y' : {
            'var'   : 'A_l_norm',
            'min'   : 0.0,
            'max'   : 20.0,
            'dim'   : 1001
        }
    },
    'solver'    : {
        'show_progress'         : False,
        'method'                : 'eig',
        'system_measure_code'   : 'A'
    },
    'system'    : {
        'A_l_norm'      : 5.0,
        'Delta_0_norm'  : -1.0, 
        'g_0_norm'      : 0.005,
        'gamma_norm'    : 0.005,
        'kappa_norm'    : 0.15,
        'T_norm'        : 0.0
    },
    'plotter'   : {
        'type'              : 'pcolormesh',
        'palette'           : 'Reds',
        'x_label'           : '$\\Delta_{0} / \\omega_{m}$',
        'x_ticks'           : [i * 0.25 - 0.5 for i in range(5)],
        'x_ticks_minor'     : [i * 0.125 - 0.5 for i in range(9)],
        'y_label'           : '$A_{l} / \\omega_{m}$',
        'y_ticks'           : list(range(0, 21, 5)),
        'y_ticks_minor'     : [i * 2.5 for i in range(9)],
        'show_cbar'         : True,
        'cbar_position'     : 'top',
        'cbar_tick_labels'  : ['0S1U', '1S0U', '0S3U', '1S2U', '2S1U', '3S0U'],
        'cbar_ticks'        : [0, 1, 2, 3, 4, 5],
        'width'             : 5.0,
        'height'            : 5.0
    }
}

# wrapper function for parallel processes
def func(system_params):
    return get_func_stability_zone(
        SystemClass=EM_00,
        params=params['solver'],
        steady_state=True,
        use_rhc=use_rhc
    )(system_params)

if __name__ == '__main__':
    # wrap looper and plot
    looper = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=func,
        params=params['looper'],
        params_system=params['system'],
        plot=True,
        params_plotter=params['plotter']
    )