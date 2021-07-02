""" Functions to calculate acoustic and elastic wave propagation in 2D

"""

from .acousticwaveprop2D import solveacoustic2D

from .elasticwaveprop2D import solveelastic2D

from .sourcetimefuncs import gaussource, rickersource

from .animatewaves import animateacousticwaves,animateelasticwaves
