"""sncosmo Source implementing the "dm15" model from Sako et al (2008).

Example::

    >>> import sncosmo
    >>> import dm15source
    >>> model = sncosmo.Model(source='dm15')
    >>> model.set(z=0.5, t0=55000., amplitude=1e-10, dm15=0.9)
    >>> model.bandmag('sdssi', 'ab', 55011.5)
    24.661108429127694
    >>> model.flux(55011.5, [4000., 5000.])
    array([  2.62194736e-20,   1.40864925e-19])
"""

import numpy as np
from scipy.interpolate import RectBivariateSpline as Spline2d
import sncosmo


class DM15Source(sncosmo.Source):
    """dm15 source from Sako et al (2008).

    Rest frame flux is given by:

        F(t, \lambda) = A * f_0(t * (\tau / 15), \lambda) * f_1(\lambda, dm15)

    where:

    f_0(t, lambda) = Hsiao template
    f_1(lambda, dm15) = 10^(-0.4*(a*x+b*x*x))

    and:

        x = dm15 - 1.1

    if lambda < 12000:
        a = 1.248 - 1.045e-4 * lambda
        b = 0.633
    else:
        a = 0
        b = 0

    tau is from Equation C4:

        tau = 3.455 + 13.719*dm15 - 3.601*dm15^2 + 0.946*dm15^3
    """

    _param_names = ['amplitude', 'dm15']
    param_names_latex = ['A', r'\Delta m_{15}']
    
    def __init__(self, phase, wave, flux, name=None, version=None):

        self.name = name
        self.version = version
        self._phase = phase
        self._wave = wave
        self._parameters = np.array([1.0, 1.1])
        self._model_flux = Spline2d(phase, wave, flux, kx=3, ky=3)

    def _flux(self, phase, wave):
        dm15 = self._parameters[1]

        # f0 is base template with time axis stretched
        tau = 3.455 + 13.719*dm15 - 3.601*dm15**2 + 0.946*dm15**3
        basephase = phase * (tau / 15.)
        f0 = self._model_flux(basephase, wave)

        # ft depends on dm15 and wave (coefficients both zero above 12000 AA)
        a = 1.248 - 1.045e-4 * wave
        a[wave >= 12000.] = 0.0
        b = 0.633 * (wave < 12000.)
        x = dm15 - 1.1
        f1 = 10**(-0.4 * (a + b * x) * x)

        return self._parameters[0] * (f0 * f1)


# register a loader for the model with Hsiao base template
def load_dm15(name=None, version=None):
    # hack to get hsiao template arrays
    hsiao = sncosmo.get_source('hsiao')
    phase = hsiao._phase
    wave = hsiao._wave
    flux = hsiao.flux(phase, wave)
    return DM15Source(phase, wave, flux, name=name, version=version)


sncosmo.register_loader(sncosmo.Source, 'dm15', load_dm15, args=())
