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
# Fully Commented

import math
import random
import time
import scipy
import configparser


import pythagor.bench.lgarcia.TARS.lib.ccd_events_catalog as CCDeCat
import pythagor.bench.lgarcia.TARS.lib.nice_figure as nice_figure
import numpy as np

from pythagor.bench.lgarcia.TARS.lib.particle import particle_specs
from pythagor.bench.lgarcia.TARS.lib.particle import particle
from pythagor.bench.lgarcia.TARS.lib.ccd import CCD
from tqdm import tqdm


class simulation_config:
    """
    simulation config is a class making the interface between a configuration file containing all the simulator parameters
    and the simulation program
    """

    def __init__(self, config_file):
        """
        :param string config_file: path of the configuration file
        """

        random = "random"

        # reading of the simulation specs
        config = configparser.ConfigParser()
        config.read_file(open(config_file))
        self.cfg_file = config
        self.real_cfg_file = config_file
        self.number_of_particles = eval(config.get('SIMULATION CONFIG', 'number_of_particles'))
        self.energy = eval(config.get('SIMULATION CONFIG', 'energy'))
        self.position_x = eval(config.get('SIMULATION CONFIG', 'position_x'))
        self.position_y = eval(config.get('SIMULATION CONFIG', 'position_y'))
        self.alpha_angle = eval(config.get('SIMULATION CONFIG', 'alpha_angle'))
        self.beta_angle = eval(config.get('SIMULATION CONFIG', 'beta_angle'))
        self.spreading_step = eval(config.get('SIMULATION CONFIG', 'spreading_step'))

        # reading of the simulation inputs specs
        if 'INPUTS FILE' in config.sections():
            input_positions_file = config.get('INPUTS FILE', 'positions')
            input_energies_file = config.get('INPUTS FILE', 'energies')
            # reading of the simulation inputs choice : 'yes' for outside inputs
            self.input_file_choice = config.get('INPUTS FILE', 'input_file')

            # if outside inputs are chosen then the simulator take them as inputs
            if self.input_file_choice is 'yes':
                self.input_positions = np.loadtxt("../Data/Inputs/Inputs from external source/" + input_positions_file)
                self.input_energies = np.loadtxt("../Data/Inputs/Inputs from external source/" + input_energies_file)

        # set date and processing time. This are empty (or 0) in the input_simulation_config but nit in
        # the simulation_report
        if 'SIMULATION REPORT' in config.sections():
            self.date = (config.get('SIMULATION REPORT', 'date'))
            self.processing_time = eval(config.get('SIMULATION REPORT', 'processing_time'))

        # set date and processing time. This are empty (or 0) in the input_simulation_config
        # but nit in the simulation_report
        if self.energy is not 'random':
            self.actualise()

    def actualise(self):
        """
        This function actualise the config file options (used after a set for exemple)
        :return:
        """
        config_file = open(self.real_cfg_file, 'w')
        self.cfg_file.write(config_file)

    def set_value(self, cfg_option, value):
        """
        This function set the option to a choosen value
        :param string cfg_option: option to be set
        :param value: new value taken by the option
        :return:
        """
        self.cfg_file.set('SIMULATION CONFIG', cfg_option, str(value))
        x = 'self.' + cfg_option
        # option value as to be a string that will be re-interpreted in the simulation_config class
        if not isinstance(value, str):
            exec("%s = %d" % (x, value))
        else:
            exec("%s = %s" % (x, value))
        self.actualise()


class SimResults:
    """
    SimResults is a class containing all the useful results produces by the simulation
    """

    def __init__(self, events_catalog, n_CCD):
        """
        Initialisation of SimResults
        :param events_catalog events_catalog: Events catalog (events_catalog object from CCD_events_catalog library) wich will contain
        all the CR events
        :param CCD n_CCD: CCD object(from CCD library) containing all the simulated CCD specs
        """

        #   Here is a catalog of all the CR events - see CCD_events_catalog
        self.events_catalog = events_catalog

        #   Here is an image of all the simulated CRs events on the CCD
        self.pcmap = np.zeros((n_CCD.ccd_cfg.number_of_pixels_ac, n_CCD.ccd_cfg.number_of_pixels_al))

        #   Here is an image of all the last simulated CRs events on the CCD
        self.pcmap_clean = np.zeros((n_CCD.ccd_cfg.number_of_pixels_ac, n_CCD.ccd_cfg.number_of_pixels_al))

        #   Here is all the energy of the simulated particles arriving on the CCD
        self.en = []

        #   Here is all the alpha incidence angles of the simulated particles arriving on the CCD
        self.alpha = []

        #   Here is all the beta incidence angles of the simulated particles arriving on the CCD
        self.beta = []


class Simulation:
    """
    Main class of the program, Simulation contain all the methods to set and run a simulation
    """

    def __init__(self, CCD_sim, simulation_name, sim_cfg_file, result_folder):
        """
        Initialisation of the simulation

        :param CCD CCD_sim: CCD object(from CCD library) containing all the simulated CCD specs
        :param string simulation_name: Name of the simulation
        :param string sim_cfg_file: Path of the simulation configuration file
        :param string result_folder: Name of the folder containing all the simulation results
        """

        self.simulation_name = simulation_name

        self.CCD = CCD_sim

        self.particle_spec = particle_specs(simulation_name + "_particle_spec")

        self.Results = SimResults(CCDeCat.events_catalog(simulation_name + "_events_catalog"), self.CCD)

        self.processing_time = 0

        self.sim_cfg = simulation_config(sim_cfg_file)

        self.result_folder = result_folder

    def event_generation(self, particle_energy, particle_angle_alpha, particle_angle_beta):
        """
        Generation of an event on the CCD due to an incident particle taken according to the simulation configuration
        file
        :param float or 'random' particle_energy: incident particle energy
        ('random' for a random energy taken in the spectrum)
        :param float or 'random' particle_angle_alpha: alpha incidence angle of the particle
        :param float or 'random' particle_angle_beta: betaincidence angle of the particle
        :return:
        """

        self.Results.pcmap_clean[:, :] = 0

        part_energy_test = 100000.0

        #   part_energy_test, energy of the particle, has to be < 10000 to be take in account in our simulation
        # (exclusion of high energy events)
        while part_energy_test > 10000.0:

            p = particle(self.particle_spec, self.Results.events_catalog, self.CCD, particle_energy,
                         particle_angle_alpha, particle_angle_beta)

            part_energy_test = p.part_energy

            p.actual_position = np.array([self.CCD.ccd_cfg.ccd_dimension_ac * random.random(), self.CCD.ccd_cfg.ccd_dimension_al * random.random()])
            p.initial_position = np.copy(p.actual_position)

        #   starting to fill the simulation results
        self.Results.en.append(p.part_energy)
        self.Results.alpha.append(p.angle[0])
        self.Results.beta.append(p.angle[1])

        #   single particle visualisation _ only for unique n_event
        #   if number_of_particles is 1 we keep all the geometric parameters of the particle to make a nice plot
        if self.sim_cfg.number_of_particles is 1:

            #   particle position
            x = []
            y = []
            z = []

            #   particle energy
            en = []
            #   number of electrons generated by the particle at each step
            e = []
            #   dimension of the e- cloud generated at each step
            p_sig = []

            #   start taking particle information
            z.append(p.depth)
            x.append(p.actual_position[0])
            y.append(p.actual_position[1])
            en.append(p.part_energy)
            p_sig.append(0)
            e.append(0)

        # main loop : electrons generation and collection at each step while the particle is in the CCD and
        # have enough energy to spread
        if 0 <= p.depth <= (self.CCD.ccd_cfg.depletion_zone_thickness + self.CCD.ccd_cfg.field_free_zone_thickness + self.CCD.ccd_cfg.substrate_thickness) \
                and 0.0 < p.actual_position[0] < self.CCD.ccd_cfg.ccd_dimension_ac \
                and 0.0 < p.actual_position[1] < self.CCD.ccd_cfg.ccd_dimension_al \
                and 0.001 < p.part_energy < 10000.0:

            while 0 <= p.depth <= (self.CCD.ccd_cfg.depletion_zone_thickness + self.CCD.ccd_cfg.field_free_zone_thickness + self.CCD.ccd_cfg.substrate_thickness) \
                 and 0.0 < p.actual_position[0] < self.CCD.ccd_cfg.ccd_dimension_ac \
                 and 0.0 < p.actual_position[1] < self.CCD.ccd_cfg.ccd_dimension_al \
                 and 0.001 < p.part_energy < 10000.0:

                p.actual_position[0] += p.dx
                p.actual_position[1] += p.dy
                p.depth += p.dd * self.sim_cfg.spreading_step

                s = self.CCD.stopping_power_function(p.part_energy)

                #   p.energy in MeV
                p.energy = 100 * s * 2.32 * 1e-6 * self.sim_cfg.spreading_step

                if p.energy >= p.part_energy:
                    p.energy = p.part_energy
                    p.part_energy = 0
                else:
                    p.part_energy -= p.energy

                p.electrons = p.energy * 1e6 / 3.65

                if p.energy < 0:
                    print("energy", p.energy, "part_energy", p.part_energy, "s", s)

                sig = self.spread(p)
                elec = self.mat2pix(p, sig, sig)

                # single particle visualisation _ only for unique n_event
                if self.sim_cfg.number_of_particles is 1:
                    z.append(p.depth)
                    x.append(p.actual_position[0])
                    y.append(p.actual_position[1])
                    en.append(p.part_energy)
                    p_sig.append(sig)
                    e.append(elec)

            # single particle visualisation _ only for unique n_event
            if self.sim_cfg.number_of_particles is 1:
                f_e = np.sum(e)
                nice_figure.make_trace(x, y, z, en, e, f_e, p_sig, self, True, True)

        new_map = np.copy(self.Results.pcmap_clean[p.minmax[0]:(p.minmax[1] + 1), p.minmax[2]:(p.minmax[3] + 1)])

        self.Results.events_catalog.add_event(new_map, p.initial_position)

    def event_from_input(self, particle_index):
        """
        Generation of an event on the CCD due to an incident particle taken from a particles input file
        :param int particle_index: index of the particle to be generated in the particles input file
        """

        self.Results.pcmap_clean[:, :] = 0

        p = particle(self.particle_spec, self.Results.events_catalog, self.CCD, self.sim_cfg.input_energies[0][particle_index], 0, 0)

        p.actual_position = self.sim_cfg.input_positions[0][:]
        p.initial_position = np.copy(p.actual_position)

        #   single particle visualisation _ only for unique n_event
        #   if number_of_particles is 1 we keep all the geometric parameters of the particle to make a nice plot
        if self.sim_cfg.number_of_particles is 1:

            #   particle energy
            en = []
            #   number of electrons generated by the particle at each step
            e = []
            #   dimension of the e- cloud generated at each step
            p_sig = []

            en.append(p.part_energy)
            p_sig.append(0)
            e.append(0)

        for i in range(1, np.size(self.sim_cfg.input_positions, axis=0)):

            p.actual_position = self.sim_cfg.input_positions[:][i]

            s = self.CCD.stopping_power_function(self.sim_cfg.input_energies[i][particle_index])

            displacement_vector = self.sim_cfg.input_positions[i][particle_index] - self.sim_cfg.input_positions[i-1][particle_index]

            #   p.energy in MeV
            p.energy = 100 * s * 2.32 * 1e-6 * np.linalg.norm(displacement_vector)

            if p.energy >= p.part_energy:
                p.energy = p.part_energy
                p.part_energy = 0
            else:
                p.part_energy -= p.energy

            p.electrons = p.energy * 1e6 / 3.65

            if p.energy < 0:
                print("energy", p.energy, "part_energy", p.part_energy, "s", s)

            sig = self.spread(p)
            elec = self.mat2pix(p, sig, sig)

            # single particle visualisation _ only for unique n_event
            if self.sim_cfg.number_of_particles is 1:
                en.append(p.part_energy)
                p_sig.append(sig)
                e.append(elec)

        # single particle visualisation _ only for unique n_event
        if self.sim_cfg.number_of_particles is 1:
            f_e = np.sum(e)
            nice_figure.make_trace(self.sim_cfg.input_positions[:, 0],
                                   self.sim_cfg.input_positions[:, 1],
                                   self.sim_cfg.input_positions[:, 2], en, e, f_e, p_sig, self, True, True)

        new_map = np.copy(self.Results.pcmap_clean[p.minmax[0]:(p.minmax[1] + 1), p.minmax[2]:(p.minmax[3] + 1)])

        self.Results.events_catalog.add_event(new_map, p.initial_position)

    def spread(self, particle_s):
        """
        spread the particle into the material and compute the density and size of the electronic cloud generated
        at each step

        :param particle particle_s: particle to spread
        :return: float sigma : diameter of the electronic cloud at the generation point (um)

        """

        #   specify temperature in Kelvin
        temperature = self.CCD.ccd_cfg.temperature
        #     specify na in /m3 for evaluation of con in SI units
        na = 1e19
        #     specify diffusion length in um (field free region)
        l1 = 1000.
        #     depletion/field free boundary parameter
        bound = 2.
        #     constant includes factor of 1.d6 for conversion of m to um
        con = 1e6 * math.sqrt((2. * 1.38e-23 * temperature * 11.8 * 8.85e-12) / (na * 1.6e-19 ** 2))

        #     electron velocity saturation parameter
        sat = 1.6e-19 * na * self.CCD.ccd_cfg.depletion_zone_thickness / 11.8 / 8.85e-12 / temperature ** 1.55 / 1.01e8

        #     calculate initial 1 sigma cloud size in um (many refs)
        ci = 0.0044 * ((particle_s.electrons * 3.65 / 1000.) ** 1.75)

        #     spreading across entire depletion region
        cfr = con * math.sqrt(sat + bound)

        sig = 0

        if 0 <= particle_s.depth < self.CCD.ccd_cfg.depletion_zone_thickness:
            cf = con * math.sqrt(
                sat * particle_s.depth / self.CCD.ccd_cfg.depletion_zone_thickness + math.log(self.CCD.ccd_cfg.depletion_zone_thickness / (self.CCD.ccd_cfg.depletion_zone_thickness - particle_s.depth)))
            if cf > cfr:
                cf = cfr
            sig = math.sqrt(ci ** 2 + cf ** 2)
            hh = 1.0

        elif self.CCD.ccd_cfg.depletion_zone_thickness <= particle_s.depth < self.CCD.ccd_cfg.depletion_zone_thickness + self.CCD.ccd_cfg.field_free_zone_thickness:
            d = particle_s.depth - self.CCD.ccd_cfg.depletion_zone_thickness
            hh = (math.exp(self.CCD.ccd_cfg.field_free_zone_thickness / l1 - d / l1) + math.exp(d / l1 - self.CCD.ccd_cfg.field_free_zone_thickness / l1)) / (
                math.exp(self.CCD.ccd_cfg.field_free_zone_thickness / l1) + math.exp(-self.CCD.ccd_cfg.field_free_zone_thickness / l1))
            cff = self.CCD.ccd_cfg.field_free_zone_thickness / 1.0 * math.sqrt(1 - ((self.CCD.ccd_cfg.field_free_zone_thickness - d) / self.CCD.ccd_cfg.field_free_zone_thickness) ** 2)
            sig = math.sqrt(ci ** 2 + cfr ** 2 + cff ** 2)

        elif self.CCD.ccd_cfg.depletion_zone_thickness + self.CCD.ccd_cfg.field_free_zone_thickness <= particle_s.depth <= self.CCD.ccd_cfg.depletion_zone_thickness + self.CCD.ccd_cfg.field_free_zone_thickness + self.CCD.ccd_cfg.substrate_thickness:
            d = particle_s.depth - self.CCD.ccd_cfg.field_free_zone_thickness - self.CCD.ccd_cfg.depletion_zone_thickness
            cff = self.CCD.ccd_cfg.field_free_zone_thickness / 1.0
            hhsub = math.sinh((self.CCD.ccd_cfg.substrate_thickness - d) / 10.) / math.sinh(self.CCD.ccd_cfg.substrate_thickness / 10.)
            hhff = 2. / (math.exp(self.CCD.ccd_cfg.field_free_zone_thickness / l1) + math.exp(-self.CCD.ccd_cfg.field_free_zone_thickness / l1))
            hh = hhsub * hhff
            cfsub = 0.5 * self.CCD.ccd_cfg.substrate_thickness * math.sqrt(1 - ((self.CCD.ccd_cfg.substrate_thickness - d) / self.CCD.ccd_cfg.substrate_thickness) ** 2)
            sig = math.sqrt(ci ** 2 + cfr ** 2 + cfsub ** 2 + cff ** 2)

        else:
            hh = 0

        particle_s.electrons *= hh

        return sig

    def mat2pix(self, particle_m2p, sigac, sigal):
        """
        Compute the charge collection function to determine the number of electron collected by each pixel based on the
        generated electronic cloud shape
        :param particle particle_m2p: particle responsible of the electronic cloud
        :param float sigac: diameter of the resulting electronic cloud in the AC (across scan, vertical) dimension
        :param float sigal: diameter of the resulting electronic cloud in the AL (along scan, horizontal) dimension
        """

        px = []
        py = []

        dx = particle_m2p.actual_position[0] - self.CCD.ccd_cfg.pixel_ac_size * int(particle_m2p.actual_position[0] / self.CCD.ccd_cfg.pixel_ac_size)
        dy = particle_m2p.actual_position[1] - self.CCD.ccd_cfg.pixel_al_size * int(particle_m2p.actual_position[1] / self.CCD.ccd_cfg.pixel_al_size)

        try:
            int(4 * sigac / self.CCD.ccd_cfg.pixel_ac_size)

        except ValueError:
            print(sigac, particle_m2p.electrons)

        xsteps = int(4 * sigac / self.CCD.ccd_cfg.pixel_ac_size)

        if xsteps > 49:
            xsteps = 49
        if xsteps < 1:
            xsteps = 1

        ysteps = int(4 * sigal / self.CCD.ccd_cfg.pixel_al_size)
        if ysteps > 49:
            ysteps = 49
        if ysteps < 1:
            ysteps = 1

        for xi in np.arange(-(xsteps * self.CCD.ccd_cfg.pixel_ac_size + dx), ((xsteps + 1) * self.CCD.ccd_cfg.pixel_ac_size - dx),
                            self.CCD.ccd_cfg.pixel_ac_size):

            if sigac is not 0:
                case1 = (xi + self.CCD.ccd_cfg.pixel_ac_size) / 1.41 / sigac
                case2 = xi / 1.41 / sigac
            else:
                case1 = 0
                case2 = 0

            px.append((scipy.special.erf(case1) - scipy.special.erf(case2)) / 2)

        for yi in np.arange(-(ysteps * self.CCD.ccd_cfg.pixel_al_size + dy), ((ysteps + 1) * self.CCD.ccd_cfg.pixel_al_size - dy),
                            self.CCD.ccd_cfg.pixel_al_size):

            if sigac is not 0:
                case1 = (yi + self.CCD.ccd_cfg.pixel_al_size) / 1.41 / sigal
                case2 = yi / 1.41 / sigal

            else:
                case1 = 0
                case2 = 0

            py.append((scipy.special.erf(case1) - scipy.special.erf(case2)) / 2)

        cx = 0

        for ix in range(int(particle_m2p.actual_position[0] / self.CCD.ccd_cfg.pixel_ac_size) - xsteps,
                        int(particle_m2p.actual_position[0] / self.CCD.ccd_cfg.pixel_ac_size) + xsteps + 1, 1):

            if ix < particle_m2p.minmax[0]:
                particle_m2p.minmax[0] = ix

            if ix > particle_m2p.minmax[1]:
                particle_m2p.minmax[1] = ix

            cy = 0

            for iy in range(int(particle_m2p.actual_position[1] / self.CCD.ccd_cfg.pixel_al_size) - ysteps,
                            int(particle_m2p.actual_position[1] / self.CCD.ccd_cfg.pixel_al_size) + ysteps + 1, 1):

                if iy < particle_m2p.minmax[2]:
                    particle_m2p.minmax[2] = iy

                if iy > particle_m2p.minmax[3]:
                    particle_m2p.minmax[3] = iy

                if 0 <= ix < self.CCD.ccd_cfg.number_of_pixels_ac and 0 <= iy < self.CCD.ccd_cfg.number_of_pixels_al:
                    self.Results.pcmap[ix, iy] += px[cx] * py[cy] * particle_m2p.electrons
                    self.Results.pcmap_clean[ix, iy] += px[cx] * py[cy] * particle_m2p.electrons
                cy += 1

            cx += 1

            if particle_m2p.minmax[0] < 0:
                particle_m2p.minmax[0] = 0

            if particle_m2p.minmax[2] < 0:
                particle_m2p.minmax[2] = 0

        return particle_m2p.electrons

    def run(self):
        """
        Main running of the simulation based on the number of iteration specified in the main input file
        """
        #   computation of the processing time
        start_time = time.time()
        print("Simulation processing...")

        # simulation processing depending on many case:
        #   - if the particle energy is taken randomly from an input spectrum
        #   - if the particle energy is fixed and unique or following a range (linear or logarithmic one)
        #   - if the particle energy  and position are taken from a particles input file

        if self.sim_cfg.input_file_choice is 'yes':

            for i in range(0, np.size(self.sim_cfg.input_energies)):
                self.event_from_input(i)

        else:

            if self.sim_cfg.energy is 'random':

                for i in tqdm(range(0, self.sim_cfg.number_of_particles)):
                    self.event_generation(self.sim_cfg.energy, self.sim_cfg.alpha_angle, self.sim_cfg.beta_angle)

            elif isinstance(self.sim_cfg.energy, int) or isinstance(self.sim_cfg.energy, float):

                for i in tqdm(range(0, self.sim_cfg.number_of_particles)):
                    self.event_generation(self.sim_cfg.energy, self.sim_cfg.alpha_angle, self.sim_cfg.beta_angle)

            else:

                if self.sim_cfg.number_of_particles > np.size(self.sim_cfg.energy):
                    iteration_by_p = int(self.sim_cfg.number_of_particles / np.size(self.sim_cfg.energy))

                else:
                    iteration_by_p = 1

                print(self.sim_cfg.energy)

                for u in tqdm(self.sim_cfg.energy):
                    for j in range(0, iteration_by_p):
                        self.event_generation(u, self.sim_cfg.alpha_angle, self.sim_cfg.beta_angle)

        proc_t = (time.time() - start_time)
        print("--- %s seconds ---" % proc_t)

        #   writing of all the results and data into a folder for further consultation
        nice_figure.create_result_folder(self.result_folder, self)

        self.processing_time = proc_t
