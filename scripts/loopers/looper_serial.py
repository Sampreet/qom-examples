from qom.loopers import XLooper
import numpy as np

# function to obtain the square of the values
def func_square(system_params, val, logger, res):
    res.append((val, val**2))
    
# initialize looper
looper = XLooper(func=func_square, params={
    'looper': {
        'X': np.linspace(0, 10, 101),
    },
    'plotter': {
        'x_label': '$x$',
        'x_ticks': list(range(0, 11, 2)),
        'v_label': '$x^{2}$',
        'v_ticks': list(range(0, 101, 20))
    }
})
# loop and plot
looper.wrap(file_path_prefix=None, plot=True)