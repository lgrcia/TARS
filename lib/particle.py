#   --------------------------------------------------------------------------
#   Copyright 2015 SRE-F, ESA (European Space Agency)
#       Lionel Garcia <lionel_garcia@live.fr>
#
#   This is restricted software and is only to be used with permission
#   from the author, or from ESA.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#   --------------------------------------------------------------------------
#
# Fully documented
# Not fully commented

import math
import random
import bisect

import numpy as np

import pythagor.bench.lgarcia.TARS.lib.util as util


class particle:
    """
     particle class define a particle together with its characteristics
    """

    def __init__(self, particle_specs_p,
                 events_catalog,
                 CCD_part,
                 in_energy='random',
                 in_alpha='random',
                 in_beta='random'):
        """
        Creation of a particle according to some parameters
        :param particle_specs particle_specs_p: specs of the particle generated
        :param events_catalog events_catalog: events catalog in which the particle is going to be referenced
        :param CCD_part: CCD in which the particle is going to interact
        :param  float or 'random' in_energy: initial energy of the incident particle
        :param  float or 'random' in_alpha: alpha incident angle of the particle in the ccd
        :param float or 'random' in_beta: beta incident angle of the particle in the ccd
        """

        self.particle_specs = particle_specs_p

        self.CCD = CCD_part

        if in_energy is 'random':
            u = np.random.random()
            self.part_energy = self.particle_specs.CDF[0][bisect.bisect(self.particle_specs.CDF[1], u) - 1]
        else:
            self.part_energy = in_energy

        if in_alpha is 'random':
            alpha = 2. * math.pi * random.random()
        else:
            alpha = in_alpha

        if in_beta is 'random':
            beta = 2. * math.pi * random.random()
        else:
            beta = in_beta

        self.angle = np.array([alpha, beta])

        self.actual_position = 0

        self.initial_position = 0

        self.energy = 0

        self.index = np.size(events_catalog) - 1

        self.dd = math.sin(self.angle[0])
        self.dx = math.cos(self.angle[0]) * math.cos(self.angle[1])
        self.dy = math.cos(self.angle[0]) * math.sin(self.angle[1])

        self.depth = 0

        if self.dd < 0:
            self.depth = self.CCD.ccd_cfg.total_thickness

        self.electrons = 0

        self.minmax = np.array([100000, -1, 100000, -1])


class particle_specs:
    """
    specs for particle generation
    """
    def __init__(self, particule_type='none'):
        """
        Initialisation of particle specs
        :param string particule_type: type of the particle
        """
        self.particule_type = particule_type

        self.spectre = 0

        self.spectre_function = 0

        self.CDF = 0

    def add_spectre_from_file(self, file_name):
        """
        Setting up the particle specs according to a spectrum
        :param string file_name: path of the file containing the spectrum
        """
        self.spectre, spectre_function = util.read_doc(file_name)

        self.spectre[:, 1] *= 4 * math.pi * 1e-4

        rand = np.arange(np.min(self.spectre[:, 0]), np.max(self.spectre[:, 0]), 0.001)

        spectre = spectre_function(rand)

        Csum = np.cumsum(spectre)
        Csum /= np.max(Csum)

        self.CDF = (rand, Csum)

        return 0
