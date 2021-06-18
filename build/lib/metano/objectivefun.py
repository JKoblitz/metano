"""
This module defines a number of objective functions for nonlinear optimization.


This file is part of metano.
Copyright (C) 2010-2019 Alexander Riemer, Julia Helmecke
Braunschweig University of Technology,
Dept. of Bioinformatics and Biochemistry

metano is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

metano is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with metano.  If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import division

from builtins import range
from past.utils import old_div
from numpy import dot

CUTOFF = 1E-20
LAX_CUTOFF = 1E-8


def biomass_per_flux(flux_vec, biomass_index):
    """ nonlinear objective function biomass per flux unit
    """
    denominator = (flux_vec**2).sum()
    if denominator < CUTOFF:
        return 0.   # Punish all-zero flux (in this case, the biomass flux is
                    # very low as well)
    return old_div(flux_vec[biomass_index],denominator)


def neg_grad_biomass_per_flux(flux_vec, biomass_index):
    """ compute gradient of biomass_per_flux function
    """
    denominator = (flux_vec**2).sum()
    sqdenom = denominator*denominator
    if sqdenom < CUTOFF:
        return 0.   # Punish all-zero flux (in this case, the numerator is very
                    # low as well)
    factor = 2.*flux_vec[biomass_index]/sqdenom
    solution = flux_vec*factor
    solution[biomass_index] -= 1./denominator
    return solution


def linear_to_lambda_func(vec):
    """ transform the linear function given by the coefficient vector into a
        function object
    """
    nonzero_indices = []
    for i in range(len(vec)):
        if abs(vec[i]) >= LAX_CUTOFF:
            nonzero_indices.append(i)

    # Special case: only one non-zero index
    if len(nonzero_indices) == 1:
        i = nonzero_indices[0]
        return lambda x : x[i]*vec[i]

    return lambda x : dot(x, vec)
