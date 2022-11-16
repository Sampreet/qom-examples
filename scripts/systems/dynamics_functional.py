####################################################################################                              DEPENDENCIES                              ####
################################################################################

# dependencies
import logging
import numpy as np  

# qom modules
from qom.systems import SOSMSystem
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

####################################################################################                               FUNCTIONS                                ####
################################################################################

# function to obtain the normalized drift matrix
def func_A(modes, params, t=None):
    # extract frequently used variables
    A_l_norm, Delta_0_norm, g_0_norm, gamma_norm, kappa_norm = params
    alpha, beta = modes

    # effective values
    Delta_norm = Delta_0_norm + 2 * g_0_norm * np.real(beta)
    g_norm = g_0_norm * alpha

    # drift matrix
    A_norm = np.zeros((4, 4), dtype=np.float_)
    # X quadratures
    A_norm[0][0] = - kappa_norm / 2
    A_norm[0][1] = - Delta_norm
    A_norm[0][2] = - 2 * np.imag(g_norm)
    # Y quadratures
    A_norm[1][0] = Delta_norm
    A_norm[1][1] = - kappa_norm / 2
    A_norm[1][2] = 2 * np.real(g_norm)
    # Q quadratures
    A_norm[2][2] = - gamma_norm / 2
    A_norm[2][3] = 1.0
    # P quadratures
    A_norm[3][0] = 2 * np.real(g_norm)
    A_norm[3][1] = 2 * np.imag(g_norm)
    A_norm[3][2] = - 1.0
    A_norm[3][3] = - gamma_norm / 2

    return A_norm

# function to obtain the initial values and constants required for the IVP
def func_ivc():
    # extract frequently used variables
    A_l_norm       = 25.0
    Delta_0_norm   = -1.0
    g_0_norm       = 0.005
    gamma_norm     = 0.005
    kappa_norm     = 0.15
    n_th            = 0.0

    # initial mode values as 1D list
    modes_0 = np.zeros(2, dtype=np.complex_).tolist()
    
    # normalized parameters
    c_params = [A_l_norm, Delta_0_norm, g_0_norm, gamma_norm, kappa_norm]

    # # uncomment this block to simulate only classical mode amplitudes
    # return modes_0, c_params

    # initial quadrature correlations
    corrs_0 = np.zeros((4, 4), dtype=np.float_)
    corrs_0[0][0] = 0.5
    corrs_0[1][1] = 0.5
    corrs_0[2][2] = n_th + 0.5
    corrs_0[3][3] = n_th + 0.5

    # convert to 1D list and concatenate all variables
    iv = modes_0 + [np.complex_(element) for element in corrs_0.flatten()]

    # normalized noise correlation matrix
    D_norm = np.zeros((4, 4), dtype=np.float_)
    D_norm[0][0] = kappa_norm / 2
    D_norm[1][1] = kappa_norm / 2
    D_norm[2][2] = gamma_norm * (2 * n_th + 1) / 2
    D_norm[3][3] = gamma_norm * (2 * n_th + 1) / 2

    # all constants
    c = D_norm.flatten().tolist() + c_params

    return iv, c

# function to obtain the rates of the optical and mechanical modes.
def func_mode_rates(modes, params, t=None):
    # extract frequently used variables
    A_l_norm, Delta_0_norm, g_0_norm, gamma_norm, kappa_norm = params
    alpha, beta = modes

    # effective detuning
    Delta_norm = Delta_0_norm + 2 * g_0_norm * np.real(beta)

    # calculate mode rates
    dalpha_dt = - kappa_norm / 2 * alpha + 1j * Delta_norm * alpha + A_l_norm
    dbeta_dt = 1j * g_0_norm * np.conjugate(alpha) * alpha - gamma_norm / 2 * beta - 1j * beta
    # normalize
    mode_rates = [dalpha_dt, dbeta_dt]

    return mode_rates

# function to obtain the required parameters to calculate the optical steady state
def func_oss_args(params):
    # extract frequently used variables
    A_l_norm, Delta_0_norm, g_0_norm, gamma_norm, kappa_norm = params

    # coefficient of occupancy
    C_norm = 2 * g_0_norm**2 / (gamma_norm**2 / 4 + 1)
    
    return A_l_norm, Delta_0_norm, kappa_norm, C_norm

####################################################################################                                SYSTEM                                  ####
################################################################################

# initialize logger
init_log()
logger = logging.getLogger('qom.examples')

# initialize system
system = SOSMSystem()
# set the functions to appropriate methods
system.get_A = func_A
system.get_ivc = func_ivc
system.get_mode_rates = func_mode_rates

####################################################################################                       MEAN OPTICAL OCCUPANCIES                         ####
################################################################################

# add function to system method
system.get_oss_args = func_oss_args

# obtain and display mean optical occupancies
N_os, _ = system.get_mean_optical_occupancies()
logger.info('------------------------------------------Number of real roots: {}\n'.format(len(N_os)))
logger.info('N_os: {}\n'.format(N_os))

####################################################################################                       DYNAMICAL VALUES (MODES)                         ####
################################################################################

# parameters for the solver
solver_params = {
    'show_progress': True,
    'method': 'zvode',
    'cache': False,
    'measure_type': 'mode_amp',
    'idx_e': [0, 1],
    't_min': 0.0,
    't_max': 200.0,
    't_dim': 2001
}
# calculate dynamical values
Modes, T = system.get_measure_dynamics(solver_params=solver_params)
# # uncomment for alternative expression
# Modes, _, T = system.get_modes_corrs_dynamics(solver_params=solver_params)
# extract values
x_d = (np.sqrt(2) * np.real(np.transpose(Modes)[0])).tolist()
y_d = (np.sqrt(2) * np.imag(np.transpose(Modes)[0])).tolist()
q_d = (np.sqrt(2) * np.real(np.transpose(Modes)[1])).tolist()
p_d = (np.sqrt(2) * np.imag(np.transpose(Modes)[1])).tolist()
# initialize plotter and display dynamical values
plotter = MPLPlotter(axes={
    'X': T,
    'Y': list(range(4))
}, params={
    'type': 'lines',
    'x_label': '$\\omega_{m} t$',
    'x_ticks': list(range(0, 201, 50)),
    'y_colors': ['r', 'r', 'b', 'b'],
    'y_styles': ['-', '--', '-', '--'],
    'y_legend': ['$x$', '$y$', '$q$', '$p$'],
    'v_ticks': list(range(-80, 41, 20)),
    'show_legend': True,
    'width': 5.0
})
plotter.update(xs=T, vs=[x_d, y_d, q_d, p_d])
plotter.show(hold=True)

####################################################################################                      STATIONARY VALUES (MODES)                         ####
################################################################################

# parameters for the solver
solver_params = {
    'measure_type': 'mode_amp',
    'idx_e': [0, 1]
}
# calculate stationary values
modes = system.get_measure_stationary(solver_params=solver_params)
# # uncomment for alternative method
# modes, _ = system.get_modes_corrs_stationary(solver_params=solver_params)
# extract values
x_s = np.sqrt(2) * np.real(modes[0])
y_s = np.sqrt(2) * np.imag(modes[0])
q_s = np.sqrt(2) * np.real(modes[1])
p_s = np.sqrt(2) * np.imag(modes[1])
# display stationary values
logger.info('------------------------------------Stationary value of x: {:6.2f}\n'.format(x_s))
logger.info('------------------------------------Stationary value of y: {:6.2f}\n'.format(y_s))
logger.info('------------------------------------Stationary value of q: {:6.2f}\n'.format(q_s))
logger.info('------------------------------------Stationary value of p: {:6.2f}\n'.format(p_s))

####################################################################################                    DYNAMICAL VALUES (CORRELATIONS)                     ####
################################################################################

# parameters for the solver
solver_params = {
    'show_progress': True,
    'method': 'zvode',
    'cache': False,
    'measure_type': 'corr_ele',
    'idx_e': [(2, 2), (3, 3)],
    't_min': 0.0,
    't_max': 200.0,
    't_dim': 2001
}
# calculate dynamic values
Corrs, T = system.get_measure_dynamics(solver_params=solver_params)
# extract values
Q_2_expect_d = np.transpose(Corrs)[0].tolist()
P_2_expect_d = np.transpose(Corrs)[1].tolist()
# # uncomment for alternative expression
# _, Corrs, T = system.get_modes_corrs_dynamics(solver_params=solver_params)
# # extract values
# Q_2_expect_d = [corrs[2][2] for corrs in Corrs]
# P_2_expect_d = [corrs[3][3] for corrs in Corrs]
# initialize plotter and display dynamical values
plotter = MPLPlotter(axes={
    'X': T,
    'Y': list(range(4))
}, params={
    'type': 'lines',
    'x_label': '$\\omega_{m} t$',
    'x_ticks': list(range(0, 201, 50)),
    'y_colors': ['b', 'r'],
    'y_styles': ['-', '--'],
    'y_legend': ['$\\langle Q^{2} \\rangle$', '$\\langle P^{2} \\rangle$'],
    'v_ticks': [0.4, 0.5, 0.6, 0.7],
    'show_legend': True,
    'width': 5.0
})
plotter.update(xs=T, vs=[Q_2_expect_d, P_2_expect_d])
plotter.show(hold=True)

####################################################################################                    STATIONARY VALUES (CORRELATIONS)                    ####
################################################################################

# parameters for the solver
solver_params = {
    'measure_type': 'corr_ele',
    'idx_e': [(2, 2), (3, 3)]
}
# calculate stationary values
corrs = system.get_measure_stationary(solver_params=solver_params)
# extract values
Q_2_expect_s = corrs[0]
P_2_expect_s = corrs[1]
# # uncomment for alternative method
# _, corrs = system.get_modes_corrs_stationary(solver_params=solver_params)
# # extract values
# Q_2_expect_s = corrs[2][2]
# P_2_expect_s = corrs[3][3]
# initialize logger and display stationary values
logger = logging.getLogger('qom.examples.systems')
logger.info('-------------------Stationary value of Q_2_expect: {:6.2f}\n'.format(Q_2_expect_s))
logger.info('-------------------Stationary value of P_2_expect: {:6.2f}\n'.format(P_2_expect_s))

####################################################################################                    DYNAMICAL VALUES (ENTANGLEMENT)                     ####
################################################################################

# parameters for the solver
solver_params = {
    'measure_type': 'entan_ln',
    'idx_e': (0, 1)
}
# calculate and display dynamical value
Entan_ln, T = system.get_measure_dynamics(solver_params=solver_params)
plotter = MPLPlotter(axes={}, params={
    'type': 'lines',
    'x_label': '$\\omega_{m} t$',
    'x_ticks': list(range(0, 201, 50)),
    'y_colors': ['b'],
    'v_label': '$E_{N}$',
    'v_ticks': [0.0, 0.1, 0.2, 0.3, 0.4],
    'width': 5.0
})
plotter.update(xs=T, vs=np.transpose(Entan_ln)[0])
plotter.show(hold=True)

####################################################################################                    DYNAMICAL VALUES (ENTANGLEMENT)                     ####
################################################################################

# calculate and display stationary value
entan_ln = system.get_measure_stationary(solver_params=solver_params)
logger.info('---------------------Stationary value of entan_ln: {:6.2f}\n'.format(entan_ln[0]))