from qom.utils.looper import wrap_looper
import numpy as np

# function to obtain the square of the values
def func_square(system_params):
    return system_params['x']**2
    
# loop and plot
looper = wrap_looper(
    looper_name='XLooper',
    func=func_square,
    params={
        'X' : {
            'var'   : 'x',
            'val'   : np.linspace(0, 10, 101),
        }
    },
    params_system={},
    plot=True,
    params_plotter={
        'type'      : 'lines',
        'x_label'   : '$x$',
        'x_ticks'   : list(range(0, 11, 2)),
        'v_label'   : '$x^{2}$',
        'v_ticks'   : list(range(0, 101, 20))
    }
)