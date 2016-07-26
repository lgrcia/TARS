import math
import random
import bisect

import numpy as np
from scipy import interpolate

import pythagor.bench.lgarcia.TARS.lib.util as util


class particle:
    def __init__(self, particle_specs_p, events_catalog, CCD_part):
        self.particle_specs = particle_specs_p

        self.CCD = CCD_part

        self.part_energy, self.angle = self.random_generation()

        self.coordinates = 0

        self.init_coordinates = 0

        self.energy = 0

        self.index = np.size(events_catalog) - 1

        self.dd = math.sin(self.angle[0])
        self.dx = math.cos(self.angle[0]) * math.cos(self.angle[1])
        self.dy = math.cos(self.angle[0]) * math.sin(self.angle[1])

        self.depth = 0

        if self.dd < 0:
            self.depth = self.CCD.dep + self.CCD.ff + self.CCD.sub

        self.electrons = 0

        self.minmax = np.array([100000, -1, 100000, -1])

    def random_generation(self):

        part_energy = 2500
        angles = np.array([2. * math.pi * 0.25, 2. * math.pi * 0.25])

        return part_energy, angles


class particle_specs:
    def __init__(self, particule_type='none'):
        self.particule_type = particule_type

        self.spectre = 0

        self.spectre_function = 0

        self.CDF = 0

    def add_spectre_from_file(self, n_value, n_data_type, file_name):
        self.spectre, spectre_function = util.read_doc(n_value, n_data_type, file_name)

        rand = np.arange(np.min(self.spectre[:, 0]), np.max(self.spectre[:, 0]) - 0.1, 0.1)

        cdf_function = interpolate.interp1d(self.spectre[:, 0],
                                            np.cumsum(self.spectre[:, 1]) / np.max(np.cumsum(self.spectre[:, 1])),
                                            kind='cubic')

        self.CDF = cdf_function(rand + 0.1)

        return 0
