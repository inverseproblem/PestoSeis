
#------------------------------------------------------------------------
#
#    PestoSeis, a numerical laboratory to learn about seismology, written
#    in the Python language.
#    Copyright (C) 2021  Andrea Zunino 
#
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------

"""
 Various functions to calculate travel times and ray in a rectilinear grid.

"""

######################################################################

from .traveltime_rays import rollmod,unrollmod,setupgrid,lininv,traveltime,traceallrays,traceall_straight_rays,buildtomomat,plotgrid,plotrays,plotvelmod,plotttimemod


from .rayhorlayers import tracerayhorlay


# __all__ = ["__version__","ttimerays.traveltime_rays"]
# __all__.extend(traveltime_rays.__all__)

