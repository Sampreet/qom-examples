#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
"""Class to simulate an optomechanical system with a moveable end-mirror demonstrated in arXiv:2211.02596."""

__authors__ = ['Sampreet Kalita']
__version__ = '1.0.0'
__created__ = '2021-01-01'
__updated__ = '2023-06-23'

# dependencies
import numpy as np
import scipy.constants as sc

# qom modules
from qom.systems import BaseSystem

class EM_00(BaseSystem):
    r"""Class to simulate an optomechanical system with a moveable end-mirror demonstrated in arXiv:2211.02596.
    
    Parameters
    ----------
    params : dict
        Parameters for the system. The system parameters are:
        ============    ====================================================
        key             meaning
        ============    ====================================================
        A_l_norm        (float) normalized amplitude of the laser :math:`A_{l} / \omega_{m}`. Default is :math:`25.0`.
        Delta_0_norm    (float) normalized detuning of the laser from the cavity :math:`( \omega_{l} - \omega_{o} ) / \omega_{m}`. Default is :math:`-1.0`.
        g_0_norm        (float) normalized optomechanical coupling strength :math:`g_{0} / \omega_{m}`. Default is :math:`5 \times 10^{-3}`.
        gamma_norm      (float) normalized mechanical damping rate :math:`\gamma / \omega_{m}`. Default is :math:`5 \times 10^{-3}`.
        kappa_norm      (float) normalized optical decay rate :math:`\kappa / \omega_{m}`. Default is :math:`0.15`.
        T_norm          (float) normalized bath temperature :math:`T / \omega_{m}`. Default is :math:`0.0` K/Hz.
        ============    ====================================================
    cb_update : callable, optional
        Callback function to update status and progress, formatted as ``cb_update(status, progress, reset)``, where ``status`` is a string, ``progress`` is an integer and ``reset`` is a boolean.

    References
    ----------

    .. [1] A. K. Sarma and S. Kalita, *Tutorial: Cavity Quantum Optomechanics*, arXiv:2211.02596 (2022).
    """

    # default looping parameters
    looper_defaults = {
        'A_l_norm'  : {
            'var'   : 'A_l_norm',
            'min'   : 0.0,
            'max'   : 100.0,
            'dim'   : 1001,
            'scale' : 'linear'
        },
        'Delta_0_norm'  : {
            'var'   : 'Delta_0_norm',
            'min'   : -2.0,
            'max'   : 2.0,
            'dim'   : 1001,
            'scale' : 'linear'
        },
        'g_0_norm'  : {
            'var'   : 'g_0_norm',
            'min'   : 1e-6,
            'max'   : 1e-2,
            'dim'   : 1001,
            'scale' : 'log'
        },
        'gamma_norm': {
            'var'   : 'gamma_norm',
            'min'   : 1e-8,
            'max'   : 1e-4,
            'dim'   : 1001,
            'scale' : 'log'
        },
        'kappa_norm': {
            'var'   : 'kappa_norm',
            'min'   : 1e-2,
            'max'   : 1e2,
            'dim'   : 1001,
            'scale' : 'log'
        },
        'T_norm'    : {
            'var'   : 'T_norm',
            'min'   : 1e-8,
            'max'   : 1e-4,
            'dim'   : 1001,
            'scale' : 'log'
        }
    }

    # default system parameters
    system_defaults = {
        'A_l_norm'      : 25.0,
        'Delta_0_norm'  : -1.0,
        'g_0_norm'      : 0.005,
        'gamma_norm'    : 0.005,
        'kappa_norm'    : 0.15,
        'T_norm'        : 0.0
    }

    def __init__(self, params={}, cb_update=None):
        """Class constructor for EM_00."""
        
        # initialize super class
        super().__init__(
            params=params,
            name='EM_00',
            desc='Moveable End-mirror System',
            num_modes=2,
            cb_update=cb_update
        )

    def get_A(self, modes, c, t):
        """Method to obtain the drift matrix.

        Parameters
        ----------
        modes : numpy.ndarray
            Classical modes.
        c : numpy.ndarray
            Derived constants and controls.
        t : float
            Time at which the drift matrix is calculated.
        
        Returns
        -------
        A : numpy.ndarray
            Drift matrix.
        """

        # extract frequently used variables
        alpha, beta = modes

        # effective values
        Delta_norm  = self.params['Delta_0_norm'] + 2 * self.params['g_0_norm'] * np.real(beta)
        g_norm      = self.params['g_0_norm'] * alpha

        # X quadratures
        self.A[0][0]    = - self.params['kappa_norm'] / 2.0
        self.A[0][1]    = - Delta_norm
        self.A[0][2]    = - 2.0 * np.imag(g_norm)
        # Y quadratures
        self.A[1][0]    = Delta_norm
        self.A[1][1]    = - self.params['kappa_norm'] / 2.0
        self.A[1][2]    = 2.0 * np.real(g_norm)
        # Q quadratures
        self.A[2][2]    = - self.params['gamma_norm'] / 2.0
        self.A[2][3]    = 1.0
        # P quadratures
        self.A[3][0]    = 2.0 * np.real(g_norm)
        self.A[3][1]    = 2.0 * np.imag(g_norm)
        self.A[3][2]    = - 1.0
        self.A[3][3]    = - self.params['gamma_norm'] / 2.0

        return self.A
    
    def get_D(self, modes, corrs, c, t):
        """Method to obtain the noise matrix.
        
        Parameters
        ----------
        modes : numpy.ndarray
            Classical modes.
        corrs : numpy.ndarray
            Quantum correlations.
        params : numpy.ndarray
            Derived constants and controls.
        t : float
            Time at which the drift matrix is calculated.
        
        Returns
        -------
        D : numpy.ndarray
            Noise matrix.
        """

        # thermal phonon number
        n_th = 0.0 if self.params['T_norm'] == 0.0 else 1.0 / (np.exp(sc.hbar / sc.k / self.params['T_norm']) - 1.0)
        
        # update drift matrix
        self.D[0][0]    = self.params['kappa_norm'] / 2.0
        self.D[1][1]    = self.params['kappa_norm'] / 2.0
        self.D[2][2]    = self.params['gamma_norm'] * (n_th + 0.5)
        self.D[3][3]    = self.params['gamma_norm'] * (n_th + 0.5)

        return self.D

    def get_ivc(self):
        """Method to obtain the initial values of the modes, correlations and derived constants and controls.
        
        Returns
        -------
        iv_modes : numpy.ndarray
            Initial values of the classical modes.
        iv_corrs : numpy.ndarray
            Initial values of the quantum correlations.
        c : numpy.ndarray
            Derived constants and controls.
        """

        # extract frequently used variables
        n_th = 0.0 if self.params['T_norm'] == 0.0 else 1.0 / (np.exp(sc.hbar / sc.k / self.params['T_norm']) - 1.0)
 
        # initial mode values
        iv_modes = np.zeros(self.num_modes, dtype=np.complex_)

        # initial quadrature correlations
        iv_corrs        = np.zeros(self.dim_corrs, dtype=np.float_)
        iv_corrs[0][0]  = 0.5
        iv_corrs[1][1]  = 0.5
        iv_corrs[2][2]  = n_th + 0.5
        iv_corrs[3][3]  = n_th + 0.5

        return iv_modes, iv_corrs, np.empty(0)

    def get_mode_rates(self, modes, c, t):
        """Method to obtain the rates of the classical modes.

        Parameters
        ----------
        modes : numpy.ndarray
            Classical modes.
        c : numpy.ndarray
            Derived constants and controls.
        t : float
            Time at which the drift matrix is calculated.
        
        Returns
        -------
        mode_rates : numpy.ndarray
            Normalized rates for each mode.
        """

        # extract frequently used variables
        alpha, beta = modes

        # effective values
        Delta_norm  = self.params['Delta_0_norm'] + 2 * self.params['g_0_norm'] * np.real(beta)
        g_norm      = self.params['g_0_norm'] * alpha

        # calculate mode rates
        dalpha_dt   = - self.params['kappa_norm'] / 2 * alpha + 1.0j * Delta_norm * alpha + self.params['A_l_norm']
        dbeta_dt    = 1.0j * g_norm * np.conjugate(alpha) - self.params['gamma_norm'] / 2.0 * beta - 1.0j * beta

        return np.array([dalpha_dt, dbeta_dt], dtype=np.complex_)

    def get_optomechanical_damping_rate_norm(self, c):
        r"""Method to obtain the normalized optomechanical damping rate :math:`\gamma_{om} \omega_{m} / g_{s}^{2}`.

        Refer [1]_ for the simplified expression.

        Parameters
        ----------
        c : numpy.ndarray
            Derived constants and controls.
        
        Returns
        -------
        gamma_om_norm : numpy.float_
            Normalized optomechanical damping rate.
        """

        # calculate using optomechanical self energy
        gamma_om_norm = np.imag(self.get_optomechanical_self_energy_norm(c))
        # # calculate using simplified experession
        # gamma_om_norm = self.params['kappa_norm'] * (1.0 / (self.params['kappa_norm']**2 / 4.0 + (1.0 + self.params['Delta_0_norm'])**2) - 1 / (self.params['kappa_norm']**2 / 4.0 + (1.0 - self.params['Delta_0_norm'])**2))

        return gamma_om_norm

    def get_optomechanical_frequency_shift_norm(self, c):
        r"""Method to obtain the normalized optomechanical frequency shift :math:`d \omega_{m} \omega_{m} / g_{s}^{2}`.

        Refer [1]_ for the simplified expression.

        Parameters
        ----------
        c : numpy.ndarray
            Derived constants and controls.
        
        Returns
        -------
        d_omega_m_norm : numpy.float_
            Normalized optomechanical frequency shift.
        """

        # calculate using optomechanical self energy
        d_omega_m_norm = - np.real(self.get_optomechanical_self_energy_norm(c)) / 2.0
        # # calculate using simplified experession
        # d_omega_m_norm = (1.0 + self.params['Delta_0_norm']) / (self.params['kappa_norm']**2 / 4.0 + (1.0 + self.params['Delta_0_norm'])**2) - (1.0 - self.params['Delta_0_norm']) / (self.params['kappa_norm']**2 / 4.0 + (1.0 - self.params['Delta_0_norm'])**2)

        return d_omega_m_norm

    def get_optomechanical_self_energy_norm(self, c):
        r"""Method to obtain the normalized optomechanical self energy :math:`\Sigma (\omega_{m}) / (m g_{s}^{2})`.

        Refer [1]_ for the compact expression.

        Parameters
        ----------
        c : numpy.ndarray
            Derived constants and controls.
        
        Returns
        -------
        Sigma_norm : numpy.complex_
            Normalized optomechanical self energy.
        """

        # extract frequently used variables
        chi_o_omega_norm = lambda omega_norm: 1.0 / (self.params['kappa_norm'] / 2.0 - 1.0j * (self.params['Delta_0_norm'] + omega_norm))

        # calculate using compact expression
        Sigma_norm = 2.0j * 1.0 * (chi_o_omega_norm(1.0) - np.conjugate(chi_o_omega_norm(- 1.0)))

        return Sigma_norm

    def get_params_steady_state(self, c):
        r"""Method to obtain the parameters required to calculate the optical steady states.
        
        Parameters
        ----------
        c : numpy.ndarray
            Derived constants and controls.
        
        Returns
        -------
        A_l_norm : float
            Normalized amplitude of the laser.
        Delta_0_norm : float
            Normalized detuning of the laser.
        kappa_norm : float
            Normalized optical decay rate.
        C : float
            Coefficient of :math:`|\alpha_{s}|^{2}`.
        """

        # coefficient of the mean optical occupancies
        C = 2.0 * self.params['g_0_norm']**2 * 1.0 / (self.params['gamma_norm']**2 / 4.0 + 1.0**2)
        
        return self.params['A_l_norm'], self.params['Delta_0_norm'], self.params['kappa_norm'], C

    def get_modes_steady_state(self, c):
        """Method to obtain the steady state mode apmlitudes from the mean optical occupancies.
        
        Parameters
        ----------
        c : numpy.ndarray
            Derived constants and controls.
        
        Returns 
        -------
        Modes : numpy.ndarray
            Optical and mechanical mode amplitudes for each mean optical occupancy.
        """

        # frequently used variables
        A_l_norm, Delta_0_norm, kappa_norm, C = self.get_params_steady_state(c)

        # initialize lists
        Modes = list()

        # get mean optical occupancies
        N_os, _ = self.get_mean_optical_occupancies(
            method='cubic'
        )
        # for each mean optical occupancy
        for N_o in N_os:
            # calculate mode amplitudes
            alpha   = A_l_norm / (kappa_norm / 2.0 - 1.0j * (Delta_0_norm + C * N_o))
            beta    = 1.0j * self.params['g_0_norm'] * N_o * (self.params['gamma_norm'] / 2.0 + 1.0j * 1.0) / (self.params['gamma_norm']**2 / 4.0 + 1.0**2)

            # append to list
            Modes.append([alpha, beta])

        return np.array(Modes, dtype=np.complex_)