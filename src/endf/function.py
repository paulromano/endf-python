# SPDX-FileCopyrightText: 2023-2025 OpenMC contributors and Paul Romano
# SPDX-License-Identifier: MIT

from collections.abc import Iterable

import numpy as np

from .data import EV_PER_MEV


class Tabulated1D:
    """A one-dimensional tabulated function.

    This class mirrors the TAB1 type from the ENDF-6 format. A tabulated
    function is specified by tabulated (x,y) pairs along with interpolation
    rules that determine the values between tabulated pairs.

    Once an object has been created, it can be used as though it were an actual
    function, e.g.:

    >>> f = Tabulated1D([0, 10], [4, 5])
    >>> [f(xi) for xi in numpy.linspace(0, 10, 5)]
    [4.0, 4.25, 4.5, 4.75, 5.0]

    Parameters
    ----------
    x : Iterable of float
        Independent variable
    y : Iterable of float
        Dependent variable
    breakpoints : Iterable of int
        Breakpoints for interpolation regions
    interpolation : Iterable of int
        Interpolation scheme identification number, e.g., 3 means y is linear in
        ln(x).

    Attributes
    ----------
    x : Iterable of float
        Independent variable
    y : Iterable of float
        Dependent variable
    breakpoints : Iterable of int
        Breakpoints for interpolation regions
    interpolation : Iterable of int
        Interpolation scheme identification number, e.g., 3 means y is linear in
        ln(x).
    n_regions : int
        Number of interpolation regions
    n_pairs : int
        Number of tabulated (x,y) pairs

    """

    def __init__(self, x, y, breakpoints=None, interpolation=None):
        if breakpoints is None or interpolation is None:
            # Single linear-linear interpolation region by default
            self.breakpoints = np.array([len(x)])
            self.interpolation = np.array([2])
        else:
            self.breakpoints = np.asarray(breakpoints, dtype=int)
            self.interpolation = np.asarray(interpolation, dtype=int)

        self.x = np.asarray(x)
        self.y = np.asarray(y)

    def __repr__(self):
        return f"<Tabulated1D: {self.x.size} points, {self.breakpoints.size} regions>"

    def __call__(self, x):
        # Check if input is scalar
        if not isinstance(x, Iterable):
            return self._interpolate_scalar(x)

        x = np.array(x)

        # Create output array
        y = np.zeros_like(x)

        # Get indices for interpolation
        idx = np.searchsorted(self.x, x, side='right') - 1

        # Loop over interpolation regions
        for k in range(len(self.breakpoints)):
            # Get indices for the begining and ending of this region
            i_begin = self.breakpoints[k-1] - 1 if k > 0 else 0
            i_end = self.breakpoints[k] - 1

            # Figure out which idx values lie within this region
            contained = (idx >= i_begin) & (idx < i_end)

            xk = x[contained]                 # x values in this region
            xi = self.x[idx[contained]]       # low edge of corresponding bins
            xi1 = self.x[idx[contained] + 1]  # high edge of corresponding bins
            yi = self.y[idx[contained]]
            yi1 = self.y[idx[contained] + 1]
            
            p = self.interpolation[k]
            y[contained] = self._interpolate(p, xk, xi, yi, xi1, yi1)

        # In some cases, x values might be outside the tabulated region due only
        # to precision, so we check if they're close and set them equal if so.
        y[np.isclose(x, self.x[0], atol=1e-14)] = self.y[0]
        y[np.isclose(x, self.x[-1], atol=1e-14)] = self.y[-1]

        return y

    def _interpolate_scalar(self, x):
        if x <= self._x[0]:
            return self._y[0]
        elif x >= self._x[-1]:
            return self._y[-1]

        # Get the index for interpolation
        idx = np.searchsorted(self._x, x, side='right') - 1

        # Loop over interpolation regions
        for b, p in zip(self.breakpoints, self.interpolation):
            if idx < b - 1:
                break

        xi = self._x[idx]       # low edge of the corresponding bin
        xi1 = self._x[idx + 1]  # high edge of the corresponding bin
        yi = self._y[idx]
        yi1 = self._y[idx + 1]

        return self._interpolate(p, x, xi, yi, xi1, yi1)

    @staticmethod
    def _interpolate(p, x, xi, yi, xi1, yi1):
        if p == 1:
            # Histogram
            return yi

        elif p == 2:
            # Linear-linear
            return yi + (x - xi)/(xi1 - xi)*(yi1 - yi)

        elif p == 3:
            # Linear-log
            return yi + np.log(x/xi)/np.log(xi1/xi)*(yi1 - yi)

        elif p == 4:
            # Log-linear
            return yi*np.exp((x - xi)/(xi1 - xi)*np.log(yi1/yi))

        elif p == 5:
            # Log-log
            return yi*np.exp(np.log(x/xi)/np.log(xi1/xi)*np.log(yi1/yi))

        else:
            raise ValueError(f"Unknown interpolation rule {p}")

    def __len__(self):
        return len(self.x)

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def breakpoints(self):
        return self._breakpoints

    @property
    def interpolation(self):
        return self._interpolation

    @property
    def n_pairs(self):
        return len(self.x)

    @property
    def n_regions(self):
        return len(self.breakpoints)

    @x.setter
    def x(self, x):
        self._x = x

    @y.setter
    def y(self, y):
        self._y = y

    @breakpoints.setter
    def breakpoints(self, breakpoints):
        self._breakpoints = breakpoints

    @interpolation.setter
    def interpolation(self, interpolation):
        self._interpolation = interpolation

    def integral(self):
        """Integral of the tabulated function over its tabulated range.

        Returns
        -------
        numpy.ndarray
            Array of same length as the tabulated data that represents partial
            integrals from the bottom of the range to each tabulated point.

        """

        # Create output array
        partial_sum = np.zeros(len(self.x) - 1)

        i_low = 0
        for k in range(len(self.breakpoints)):
            # Determine which x values are within this interpolation range
            i_high = self.breakpoints[k] - 1

            # Get x values and bounding (x,y) pairs
            x0 = self.x[i_low:i_high]
            x1 = self.x[i_low + 1:i_high + 1]
            y0 = self.y[i_low:i_high]
            y1 = self.y[i_low + 1:i_high + 1]
            
            p = self.interpolation[k]
            partial_sum[i_low:i_high] = self._integrate(p, x0, y0, x1, y1)

            i_low = i_high

        return np.concatenate(([0.], np.cumsum(partial_sum)))

    def integrate(self, a: float, b: float, clip_bounds: bool = True) -> float:
        """Performs a definite integral over a portion of the function.

        Performs the definite integral over the range [a,b] for the tabulated
        function. If clip_bounds is True (default), then the integration bounds
        will be truncated to be within the tabulated domain. This is equivalent
        to the function being zero outside the tabulation. If clip_bounds is
        False, then a ValueError is raised if a or b are outside the domain.

        Parameters
        ----------
        a : float
            Lower bound of integration.
        b : float
            Upper bound of integration.
        clip_bounds : bool, default True
            If True, the integration bounds are clipped to cover the tabulated
            domain of the function. If False, a ValueError will be raised if
            either integration bound is outside the tabulated domain.

        Returns
        -------
        float
            Value of the definite integral.

        Raises
        ------
        ValueError
            If clip_bounds is False and either of the integration bounds a or b
            are outside the domain of the tabulated function.
        """
        # If the integration range doesn't go from low to high, flip the order
        flipped = False
        if a > b:
            flipped = True
            tmp = a
            a = b
            b = tmp
        
        # Next we check that the integration range is valid. If the bounds go
        # off the grid, we clip them if clib_bounds is True. Otherwise, we
        # raise an exception.
        if not clip_bounds:
            if a < self._x[0] or self._x[-1] < b:
                raise ValueError("Integration bounds are outside of the "
                                 "tabulated function domain.")
        else:
            # Clip the bounds if necessary
            if a < self._x[0]:
                a = self._x[0]
            if self._x[-1] < b:
                b = self._x[-1]
        
        # Check for this special case
        if a == b:
            return 0.

        # This function finds the interpolation rule for a given index
        def get_interpolation(indx: int) -> int:
            # Loop over interpolation regions
            for b, p in zip(self.breakpoints, self.interpolation):
                if indx < b - 1:
                    return p
            # Should never get here
            return self.interpolation[-1]

        # Now we can start to perform the real integral
        integral = 0.
        x_lower_bound = a
        x_upper_bound = b

        # Get the first index for interpolation
        idx = np.searchsorted(self._x, x_lower_bound, side='right') - 1

        while idx < self._x.size-1:
            # Get the interpolation rule
            interp = get_interpolation(idx)

            # Get tabulated values
            xi = self._x[idx]       # low edge of the corresponding bin
            xi1 = self._x[idx + 1]  # high edge of the corresponding bin
            yi = self._y[idx]
            yi1 = self._y[idx + 1]
            
            # If we are at one of the end points, perform the necessary interpolation
            if xi < x_lower_bound:
                yi = self._interpolate(interp, x_lower_bound, xi, yi, xi1, yi1)
                xi = x_lower_bound

            if x_upper_bound < xi1:
                yi1 = self._interpolate(interp, x_upper_bound, xi, yi, xi1, yi1)
                xi1 = x_upper_bound
                idx = self._x.size # This makes sure we exit the loop

            # Contribute to the integral
            integral += self._integrate(interp, xi, yi, xi1, yi1)
            
            # Prepare for next iteration
            idx += 1
            x_lower_bound = xi1
        
        # If we had to flip integration bounds, multiply by -1
        if flipped:
            integral = -integral

        return integral

    @staticmethod
    def _integrate(p, x0, y0, x1, y1):
        if p == 1:
            # Histogram
            return y0*(x1 - x0)

        elif p == 2:
            # Linear-linear
            m = (y1 - y0)/(x1 - x0)
            return (y0 - m*x0)*(x1 - x0) + 0.5*m*(x1**2 - x0**2)

        elif p == 3:
            # Linear-log
            logx = np.log(x1/x0)
            m = (y1 - y0)/logx
            return y0 + m*(x1*(logx - 1) + x0)

        elif p == 4:
            # Log-linear
            m = np.log(y1/y0)/(x1 - x0)
            return y0/m*(np.exp(m*(x1 - x0)) - 1)

        elif p == 5:
            # Log-log
            m = np.log(y1/y0)/np.log(x1/x0)
            return y0/((m + 1)*x0**m)*(x1**(m + 1) - x0**(m + 1))
        
        else:
            raise ValueError(f"Unknown interpolation rule {p}")

    @classmethod
    def from_ace(cls, ace, idx=0, convert_units=True):
        """Create a Tabulated1D object from an ACE table.

        Parameters
        ----------
        ace : openmc.data.ace.Table
            An ACE table
        idx : int
            Offset to read from in XSS array (default of zero)
        convert_units : bool
            If the abscissa represents energy, indicate whether to convert MeV
            to eV.

        Returns
        -------
        openmc.data.Tabulated1D
            Tabulated data object

        """

        # Get number of regions and pairs
        n_regions = int(ace.xss[idx])
        n_pairs = int(ace.xss[idx + 1 + 2*n_regions])

        # Get interpolation information
        idx += 1
        if n_regions > 0:
            breakpoints = ace.xss[idx:idx + n_regions].astype(int)
            interpolation = ace.xss[idx + n_regions:idx + 2*n_regions].astype(int)
        else:
            # 0 regions implies linear-linear interpolation by default
            breakpoints = np.array([n_pairs])
            interpolation = np.array([2])

        # Get (x,y) pairs
        idx += 2*n_regions + 1
        x = ace.xss[idx:idx + n_pairs].copy()
        y = ace.xss[idx + n_pairs:idx + 2*n_pairs].copy()

        if convert_units:
            x *= EV_PER_MEV

        return Tabulated1D(x, y, breakpoints, interpolation)


class Tabulated2D:
    """Metadata for a two-dimensional function.

    This is a dummy class that is not really used other than to store the
    interpolation information for a two-dimensional function. Once we refactor
    to adopt GNDS-like data containers, this will probably be removed or
    extended.

    Parameters
    ----------
    breakpoints : Iterable of int
        Breakpoints for interpolation regions
    interpolation : Iterable of int
        Interpolation scheme identification number, e.g., 3 means y is linear in
        ln(x).

    """
    def __init__(self, breakpoints, interpolation):
        self.breakpoints = breakpoints
        self.interpolation = interpolation
