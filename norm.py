# Normalizations
# Normalizing constants to convert between
# code units and cgs units. Conversion
# routines that convert in either direction
# given a variable are provided.
# Requires python >= 2.7 because of OrderedDict

import numpy as np
import numpy.linalg as la
from collections import OrderedDict


class PhysNorm:
    """
    Class that stores units (dimensions) and scaling factors for a selection of physical
    quantities which have dimensions are derivable from the set of five SI base quantities,
    length, mass, time, current, and temperature.

    Note, this class does not make any assumption of an equation of state, from
    which one could obtain the scaling factor of one unknown (e.g. temperature).
    The class also doesn't assume between which two systems the scaling obtains.

    """

    def __init__(self, **kwargs):
        """
        kwargs          Varname - scaling value pairs. Currently supported
                        varnames are those in the keys of self.defs, namely:
                        x, m, t, curr, temp, v, dens, n, pres, pmom, pdot, force,
                        acc, pot, pflx, ener, epwr, eflx, eint, edot, cool,
                        cooln, coolm, mdot, area, volume, newton
        """

        # Independent SI base dimensions
        self.dimdefs = ['length', 'mass', 'time', 'current', 'temperature']
        self.ndims = len(self.dimdefs)

        # The tuple of dimensions (fundamental, or base quantities) is
        # length, mass, time, current, and temperature
        # (L, M, T, A, K)
        # The number in the tuple determines the power of the dimension.
        #
        #           Dimensions [0]        Description [2]
        self.defs = OrderedDict([
            ('x'      , (( 1,  0,  0,  0,  0), 'position or displacement')),
            ('m'      , (( 0,  1,  0,  0,  0), 'mass')),
            ('t'      , (( 0,  0,  1,  0,  0), 'time')),
            ('curr'   , (( 0,  0,  0,  1,  0), 'electric current')),
            ('temp'   , (( 0,  0,  0,  0,  1), 'temperature')),
            ('v'      , (( 1,  0, -1,  0,  0), 'speed')),
            ('dens'   , ((-3,  1,  0,  0,  0), 'mass density')),
            ('n'      , ((-3,  0,  0,  0,  0), 'number density')),
            ('pres'   , ((-1,  1, -2,  0,  0), 'pressure, or energy density')),
            ('pmom'   , (( 1,  1, -1,  0,  0), 'linear momentum')),
            ('pdot'   , (( 1,  1, -2,  0,  0), 'force, rate of change of linear momentum')),
            ('force'  , (( 1,  1, -2,  0,  0), 'force, rate of change of linear momentum')),
            ('acc'    , (( 1,  0, -2,  0,  0), 'grav field (or any force field)')),
            ('pot'    , (( 2,  0, -2,  0,  0), 'grav potential (or any force potential)')),
            ('pflx'   , ((-1,  1, -2,  0,  0), 'linear momentum flux (pressure)')),
            ('ener'   , (( 2,  1, -2,  0,  0), 'energy')),
            ('epwr'   , (( 2,  1, -3,  0,  0), 'power or luminosity (energy per unit time)')),
            ('eflx'   , (( 0,  1, -3,  0,  0), 'energy flux')),
            ('eint'   , (( 2,  0, -2,  0,  0), 'specific (internal) energy, energy per unit mass')),
            ('edot'   , (( 2,  0, -3,  0,  0), 'rate of change of specific internal energy density')),
            ('cool'   , ((-1,  1, -3,  0,  0), 'rate of change of internal energy density')),
            ('coolm'  , (( 5, -1, -3,  0,  0), 'rate of change of internal energy density per unit mass density^2 ')),
            ('cooln'  , (( 5,  1, -3,  0,  0), 'rate of change of internal energy density per unit number density^2')),
            ('mdot'   , (( 0,  1, -1,  0,  0), 'mass outflow/accretion/loading/etc rate')),
            ('area'   , (( 2,  0,  0,  0,  0), 'area')),
            ('volume' , (( 3,  0,  0,  0,  0), 'volume')),
            ('newton' , (( 3, -1, -2,  0,  0), 'Newtons gravitational constant')),
            ('none'   , (( 0,  0,  0,  0,  0), 'dimensionless quantity'))
        ])

        # Test if all keys are known.
        err = 0
        for k in kwargs:
            if k not in self.defs:
                print('Error, Unknown key '+k+'.')
                exit()
        if err == 1:
            exit()

        # Create ordered dictionary of kwargs
        kwargs_od = OrderedDict.fromkeys(self.defs)
        for k, v in list(kwargs_od.items()):
            if k in kwargs:
                kwargs_od[k] = kwargs[k]
            else:
                kwargs_od.pop(k)
        self.kwargs = kwargs_od
        dims = [self.defs[k][0] for k in kwargs_od]

        # Test if one of the dimension powers is zero.
        # If so, issue error message and exit
        cs = map(lambda l: np.sum(map(abs, l)), zip(*dims))
        err = 0
        for i, s in enumerate(cs):
            if s == 0:
                print('This set of scalings is incomplete. No finite dimension\
                      for '+self.dimdefs[i]+'.')
                err = 1
        if err == 1:
            exit()

        # Coefficient matrix
        cm = np.array(dims)

        # Calculate powers to construct all scalings
        scalings = OrderedDict()
        powers = OrderedDict.fromkeys(self.defs)
        for ip in powers:
            powers[ip] = la.solve(cm.T, np.array(self.defs[ip][0]))
            scalings[ip] = np.prod(np.array(list(kwargs_od.values())) ** powers[ip])

        self.powers = powers
        self.scalings = scalings

        # Create attribute for this class of each variable
        for k, v in list(scalings.items()):
            setattr(self, k, v)

        return

    def print_scalings(self):
        """
        Output a two column table of dvar name and scaling factor for all
        variables in this class.
        """
        for k, v in list(self.scalings.items()):
            print(format(k, '16s') + format(v, '>16.8e'))
