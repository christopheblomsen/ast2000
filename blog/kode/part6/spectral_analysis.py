"""Egen kode."""

import numpy as np
import matplotlib.pyplot as plt
import ast2000tools.constants as c


class SpectralAnalysis:
    """Class for analysing spectral data from the atmosphere."""
    def __init__(self, loaddata=True):
        """Initialize and load data."""
        if(loaddata is True):
            print('Loading flux observations')
            self.spectrum = np.loadtxt('spectrum_seed82_600nm_3000nm.txt')
            print('Loading lambda deviations')
            self.std_dev = np.loadtxt('sigma_noise.txt')
        else:
            print('Data loading skipped')
        # atomic mass unit in kg
        self.amu = 1.66054e-27
        # Calculate gass particle masses
        self.m_O2 = self.amu*16*2
        print(f'O2 weight is {self.m_O2}')
        self.m_H2O = self.amu*(2+16)
        self.m_CO2 = self.amu*(12+(16*2))
        self.m_CH4 = self.amu*(12+4)
        self.m_CO = self.amu*(12+16)
        self.m_N20 = self.amu*(14*2+16)

    def get_array_range(self, arr, lam, width):
        """Gets the range of the array centered around lam, with width samples."""
        idx = (np.abs(arr[:, 0] - lam)).argmin()
        if(idx < width//2):
            start = 0
        else:
            start = idx - width//2

        if(idx > arr.shape[0] - width//2):
            end = arr.shape[0]
        else:
            end = idx + width//2
        return arr[start:end, :]

    def plot_spectrum(self, lam):
        """Plot the spectrum measurements."""
        spectrum = self.get_array_range(self.spectrum, lam, 500)

        plt.xlabel('Bølgelengde i nm')
        plt.ylabel('Fluks')
        plt.plot(spectrum[:, 0], spectrum[:, 1], label='Observert')

    def plot_deviation(self, lam):
        """Plots the standard deviation."""
        deviation = self.get_array_range(self.std_dev, lam, 500)
        plt.plot(deviation[:, 0], deviation[:, 1])

    def plot_bestfit(self, lam0, lam, T):
        """Plot the best fit curve."""
        # Get the wavelengths to plot
        lam_range = self.get_array_range(self.spectrum, lam0, 500)[:, 0]
        fmod = self.Fmod(lam_range, lam, T, self.m_O2)
        plt.plot(lam_range, 1-fmod, label='Best fit')

    def Fmod(self, lam, lam0, T, m):
        """Model for the flux data.

        lam is a numpy array of wave lengths
        lam0 is the median wave length we are looking at
        T is the temperature of the gas in Kelvin
        m is the molecular weight of the gas
        """
        stdd = np.sqrt((2*c.k_B*T)/m)
        # print(f'Standard deviation is {c.k_B}*{T}/{m}={stdd}')
        return 1/(np.sqrt(2*np.pi)*stdd)*np.exp(-0.5*(lam-lam0)**2/stdd**2)

    def chi_square(self, Fobs, Fmod, deviation):
        """Calculate the Chi square.

        Fobs is the observed flux
        Fmod is the modelled flux
        deviation is the standard deviation of the observed flux
        """
        return ((Fobs - Fmod)/deviation)**2

    def doppler_shift(self, lam, v):
        """Calculates the range of wavelengths to check with max speed v.

        Returns min and max values
        """
        dlam = (v*lam)/c.c
        return lam - dlam, lam + dlam

    def best_fit(self, lam, m):
        """Find best fit for the wavelength lambda."""
        # This is the range we need to check for lam0
        lam_min, lam_max = self.doppler_shift(lam, 10000)  # Get doppler shift wavelength range with max speed 10 km/s

        # The temperature range we are looking at
        T_min = 150
        T_max = 450

        N = 100  # Number of steps to take for the lambda 0
        NT = 300  # Number of steps to take for the different temperatures
        chi2min = np.inf
        lam0min = 0
        Tmin = 0

        for lam0 in np.linspace(lam_min, lam_max, N):
            obs = self.get_array_range(self.spectrum, lam0, 500)
            dev = self.get_array_range(self.std_dev, lam0, 500)
            lam_range = obs[:, 0]  # get the wavelengths for the model
            # print(f'Working on lambda 0 = {lam0}')
            for T in np.linspace(T_min, T_max, NT):
                fmod = self.Fmod(lam_range, lam0, T, m)
                chi2 = self.chi_square(obs[:, 1], fmod, dev[:, 1])
                chi2sum = chi2.sum()
                # print(f'Chi^2 sum is {chi2sum} for T = {T}')
                if(chi2sum < chi2min):
                    chi2min = chi2sum
                    lam0min = lam0
                    Tmin = T
                    # print(f'Current chi2min is {chi2min}')
        return lam0min, Tmin


if __name__ == '__main__':
    analysis = SpectralAnalysis(True)
    #lamb = np.linspace(631, 633, 1000)
    #mod = analysis.Fmod(lamb, 632, 300, analysis.m_O2)
    #plt.plot(lamb, 1-mod)

    gasses = {
        "02": [analysis.m_O2, 632, 690, 760],
        "H2O": [analysis.m_H2O, 720, 820, 940],
        "CO2": [analysis.m_CO2, 1400, 1600],
        "CH4": [analysis.m_CH4, 1660, 2200],
        "CO": [analysis.m_CO, 2340],
        "N2O": [analysis.m_N20, 2870]
    }

    for gas, data in gasses.items():
        for wl in range(1, len(data)):
            lam, temp = analysis.best_fit(data[wl], data[0])
            print(f'Found best fit for {gas} lambda 0 = {lam} and T = {temp}')
            analysis.plot_spectrum(data[wl])
            analysis.plot_bestfit(data[wl], lam, temp)

            # analysis.plot_spectrum(690)
            # analysis.plot_spectrum(720)
            # analysis.plot_deviation()
            plt.legend()
            plt.show()
