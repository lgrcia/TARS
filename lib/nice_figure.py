import math
import os
import time
import datetime
import configparser

import matplotlib.patches as patches
import numpy as np

from scipy.stats import gaussian_kde
from matplotlib import pyplot as plt
from shutil import copyfile
from matplotlib import cm

#   bokeh importing

from bokeh.models.annotations import Label, Arrow, BoxAnnotation
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.models.arrow_heads import NormalHead
from bokeh.plotting import figure, show
from bokeh.io import gridplot


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
                self.input_positions = np.loadtxt("../Data/Inputs/" + input_positions_file)
                self.input_energies = np.loadtxt("../Data/Inputs/" + input_energies_file)

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


def plot_electrons(folder_name, bokeh=False):

    os.chdir("../Data/Results/")

    energy = np.load(folder_name + "/Particles_energies.npy")
    electrons = np.load(folder_name + "/Electrons_generated.npy")

    os.chdir("../../Script")

    if bokeh is True:

        el_distr_fig = figure(width=400,
                              plot_height=300,
                              title='Histogram of simulated CRs',
                              tools="pan,wheel_zoom,save,reset",
                              x_axis_label='Electrons',
                              y_axis_label="Occurrence", )

        hist, edges = np.histogram(electrons, bins=1000, range=(0, 1e5))
        hist_style(el_distr_fig)
        el_distr_fig.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color="#578BB5", line_color="#578BB5")

        density_el_scatter = figure(width=1200,
                                    plot_height=500,
                                    title="Density scatter map of Particle_energy vs. Number_of_electrons_generated",
                                    x_axis_type="log",
                                    tools="pan,wheel_zoom,save,reset,hover",
                                    x_axis_label='Incident particle energy (MeV)',
                                    y_axis_label="Number of electrons generated", )

        hist_style(density_el_scatter)
        xy = np.vstack([energy, electrons])
        z = gaussian_kde(xy)(xy)
        idx = z.argsort()
        energy, electrons, z = energy[idx], electrons[idx], z[idx]

        z /= np.max(z)

        source = ColumnDataSource(data=dict(energy=energy, electrons=electrons))

        colors = ["#%02x%02x%02x" % (int(200 * cm.jet(i)[0]),
                                     int(200 * cm.jet(i)[1]),
                                     int(200 * cm.jet(i)[2])) for i in z]

        density_el_scatter.circle('energy', 'electrons', energy, electrons, source=source, color=colors, size=5)

        density_el_scatter.select_one(HoverTool).tooltips = [('incoming energy', '@energy'),
                                                             ('electrons generated', '@electrons'),]

        return el_distr_fig,  density_el_scatter

    else:

        fig = plt.figure(facecolor = 'white')
        plt.style.use('ggplot')

        ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
        xy = np.vstack([energy, electrons])
        z = gaussian_kde(xy)(xy)
        idx = z.argsort()
        energy, electrons, z = energy[idx], electrons[idx], z[idx]
        ax3.scatter(energy, electrons, c=z, s=50, edgecolor='')
        ax3.set_title('Density scatter map of Particle_energy vs. Number_of_electrons_generated', fontsize=10)
        ax3.set_xlabel('Incident particle energy (MeV)', fontsize=9)
        ax3.set_ylabel('Number of electrons generated', fontsize=9)
        ax3.set_xscale('log')

        ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
        plt.hist(electrons.flatten(), bins=1000, histtype='step', range=(0, 1e5))
        ax1.set_title('Histogram of simulated CRs', fontsize=10)
        ax1.set_xlabel('Electrons', fontsize=9)
        ax1.set_ylabel("Occurrence", fontsize=9)

        fig.suptitle(folder_name + " Electrons", fontsize=12)

        plt.show()

        return 0, 0


def plot_energy(folder_name, bokeh=False, size=(400, 300)):

    os.chdir("../Data/Results/")

    energy = np.load(folder_name + "/Particles_energies.npy")
    spectre_x = np.load(folder_name + "/spectre_x.npy")
    spectre_y = np.load(folder_name + "/spectre_y.npy")
    sim_configuration = simulation_config(folder_name + "/simulator_config_report.cfg")

    os.chdir("../../Script")

    if bokeh is True:

        en_distr_fig = figure(width=size[0],
                              plot_height=size[1],
                              title='Histogram of particles energies',
                              tools="pan,wheel_zoom, save,reset",
                              x_axis_label='Energy (MeV)',
                              y_axis_label="Occurrence", )

        hist, edges = np.histogram(energy, density=True, bins=50, range=(0, 10000))
        hist_style(en_distr_fig)
        en_distr_fig.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color="#578BB5", line_color="#578BB5")

        if sim_configuration.energy is "random":

            spectre_fig = figure(width=size[0],
                                 plot_height=size[1],
                                 title='Incident particles energy spectrum',
                                 tools="pan,wheel_zoom,save, reset",
                                 x_axis_type="log",
                                 x_axis_label='Energy (MeV)',
                                 y_axis_label='Differential flux', )

            hist_style(spectre_fig)
            spectre_fig.line(spectre_x, spectre_y)

        else:

            spectre_fig = figure(width=400,
                                 plot_height=300,
                                 title='Incident particles energy spectrum',
                                 tools="pan,wheel_zoom,save, reset",
                                 x_axis_type="log",
                                 x_axis_label='Energy (MeV)',
                                 y_axis_label='Differential flux', )

        hist_style(spectre_fig)
        spectre_fig.line(edges[0:np.size(edges)], hist)

        # spectre_fig.logo = None
        # spectre_fig.toolbar_location = None
        # en_distr_fig.logo = None
        # en_distr_fig.toolbar_location = None

        return spectre_fig, en_distr_fig

    else:

        fig = plt.figure(facecolor = 'white')
        plt.style.use('ggplot')

        ax = plt.subplot2grid((1, 2), (0, 0))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title('Incident particles energy spectrum', fontsize=10)
        ax.set_xlabel('Energy (' + r'$MeV.nuc^{-1}$' + ')', fontsize=9)
        ax.set_ylabel('Differential flux (' + r'$m^{-2}.s^{-1}.sr^{-1}.MeV^{-1}.nuc$' + ')', fontsize=9)
        plt.plot(spectre_x, spectre_y)

        ax2 = plt.subplot2grid((1, 2), (0, 1))
        ax2.set_title('Histogram of particles energies', fontsize=10)
        ax2.set_xlabel('Energy(MeV)', fontsize=9)
        ax2.set_ylabel("Occurrence", fontsize=9)
        plt.hist(energy, bins=300, histtype='step', range=(0, 1500))

        fig.suptitle(folder_name + " Energy", fontsize=12)

        plt.show()

        return 0, 0


def plot_angles(folder_name, bokeh=False):

    os.chdir("../Data/Results/")

    alpha = np.load(folder_name + "/alpha.npy")
    beta = np.load(folder_name + "/beta.npy")

    os.chdir("../../Script")

    fig = plt.figure(facecolor = 'white')
    plt.style.use('ggplot')

    ax2 = plt.subplot2grid((2, 3), (0, 0))
    ax2.set_title('Histogram of particles incident angle alpha', fontsize=10)
    ax2.set_xlabel('Incident angle(rad)', fontsize=9)
    ax2.set_ylabel("Occurrence", fontsize=9)
    plt.hist(alpha, bins=100, histtype='step', range=(0, 2*math.pi))

    ax1 = plt.subplot2grid((2, 3), (0, 2))
    ax1.set_title('Histogram of particles incident angle beta', fontsize=10)
    ax1.set_xlabel('Incident angle(rad)', fontsize=9)
    ax1.set_ylabel("Occurrence", fontsize=9)
    plt.hist(beta, bins=100, histtype='step', range=(0, 2*math.pi))

    fig.suptitle(folder_name + " Angles", fontsize=12)

    dd = []
    dx = []
    dy = []

    for i in range(np.size(alpha, axis=0)):
        dd.append(math.sin(alpha[i]))
        dx.append(math.cos(alpha[i]) * math.cos(beta[i]))
        dy.append(math.cos(alpha[i]) * math.sin(beta[i]))

    ax3 = plt.subplot2grid((2, 3), (1, 0))
    ax3.set_title('Histogram of particles dz', fontsize=10)
    ax3.set_xlabel('dz(um)', fontsize=9)
    ax3.set_ylabel("Occurrence", fontsize=9)
    plt.hist(dd, bins=100, histtype='step', range=(-1, 1))

    ax4 = plt.subplot2grid((2, 3), (1, 1))
    ax4.set_title('Histogram of particles dx', fontsize=10)
    ax4.set_xlabel('dx(um)', fontsize=9)
    ax4.set_ylabel("Occurrence", fontsize=9)
    plt.hist(dx, bins=100, histtype='step', range=(-1, 1))

    ax5 = plt.subplot2grid((2, 3), (1, 2))
    ax5.set_title('Histogram of particles dy', fontsize=10)
    ax5.set_xlabel('dy(um)', fontsize=9)
    ax5.set_ylabel("Occurrence", fontsize=9)
    plt.hist(dy, bins=100, histtype='step', range=(-1, 1))

    plt.show()


def display_simulation_report(simulation_folder_name, bokeh=False):

    os.chdir("../Data/Results/")
    CFG = simulation_config(simulation_folder_name + "/simulator_config_report.cfg")
    os.chdir("../../Script")

    spectre_fig, en_distr_fig = plot_energy(simulation_folder_name, bokeh)
    el_distr_fig, density_el_scatter = plot_electrons(simulation_folder_name, bokeh)

    # n_part = create_info_cfg_box("NUMBER OF PARTICLES", CFG.number_of_particles)
    # spread_step = create_info_cfg_box("SPREADING STEPS", CFG.spreading_step)
    # alpha = create_info_cfg_box("ALPHA", CFG.alpha_angle)

    if bokeh is True:
        p = gridplot([[spectre_fig, en_distr_fig, el_distr_fig], [density_el_scatter]])
        show(p)


def create_info_cfg_box(title, value, bokeh=False):

    if bokeh is True:

        cfg_box = figure(width=150, plot_height=150, title=title, toolbar_location=None, tools="")

        if isinstance(value, int or float):
            cfg_box.text(500, 500, text=[str('{:,}'.format(value).replace(',', ' '))], text_align="center",
                         text_font_size="15pt", text_font_style='bold', text_color='white', text_baseline='middle')

        elif isinstance(value, float):
            cfg_box.text(500, 500, text=[str(value) + ' um'], text_align="center", text_font_size="15pt",
                         text_font_style='bold', text_color='white', text_baseline='middle')

        elif isinstance(value, str):
            cfg_box.text(500, 500, text=[value], text_align="center", text_font_size="15pt", text_font_style='bold',
                         text_color='white', text_baseline='middle')

        cfg_box.xgrid.grid_line_color = None
        cfg_box.ygrid.grid_line_color = None
        cfg_box.title.text_font_size = '7pt'
        cfg_box.background_fill_color = "#7F7F82"
        # cfg_box.background_fill_color = random.choice(("#7F8A8A", "#607D7F", "#87A8B5", "#567173", "#303F40"))
        cfg_box.axis.visible = None

        return cfg_box


def hist_style(bokeh_plot):

    bokeh_plot.xaxis.minor_tick_line_color = None
    bokeh_plot.yaxis.minor_tick_line_color = None
    bokeh_plot.outline_line_color = None
    bokeh_plot.xaxis.axis_line_color = None
    bokeh_plot.yaxis.axis_line_color = None
    bokeh_plot.xaxis.major_tick_line_color = None
    bokeh_plot.yaxis.major_tick_line_color = None
    bokeh_plot.title.text_font_size = '8pt'
    bokeh_plot.xaxis.axis_label_text_font = '8pt'
    bokeh_plot.yaxis.axis_label_text_font = '8pt'
    bokeh_plot.background_fill_color = "#EBEBF0"
    bokeh_plot.xgrid.grid_line_color = "white"
    bokeh_plot.ygrid.grid_line_color = "white"


def create_result_folder(new_folder_name, Simulation):
    os.chdir("../Data/Results/")

    now_time = datetime.datetime.now()

    new_folder_name += "(" + "%s-%s-%s" % (str(now_time.day).zfill(2), str(now_time.month).zfill(2),
                                           str(now_time.year).zfill(2)) + "_" + "%sh%s" \
                                                                                % (str(now_time.hour).zfill(2),
                                                                                   str(now_time.minute).zfill(
                                                                                       2)) + ")"

    if not os.path.exists(new_folder_name):
        os.makedirs(new_folder_name)

    os.chdir(new_folder_name)

    copyfile("../../Inputs/Simulation config/input_simulation_config.cfg", "simulator_config_report.cfg")
    copyfile("../../Inputs/CCD specs/ccd_configuration_file.cfg", "ccd_configuration_report.cfg")

    np.save("Particles_energies.npy", Simulation.Results.en)
    np.save("Electrons_generated.npy", Simulation.Results.events_catalog.electrons_generated)
    np.save("spectre_x.npy", Simulation.particle_spec.spectre[:, 0])
    np.save("spectre_y.npy", Simulation.particle_spec.spectre[:, 1])
    np.save("alpha.npy", Simulation.Results.alpha)
    np.save("beta.npy", Simulation.Results.beta)

    CFG = Simulation.sim_cfg
    CFG.date = time.strftime("%x") + " at " + time.strftime("%X")
    CFG.processing_time = str(Simulation.processing_time)

    os.chdir("../../../Script")

    return new_folder_name


def make_trace(x_p, y_p, z_p, en_p, e_p, f_e_p, sig, Sim, values=True, bokeh=False):

    if bokeh is False:

        plt.style.use('ggplot')
        fig1 = plt.figure(figsize=(20, 10))
        ax1 = fig1.add_subplot(111, aspect='equal')
        # ax2 = fig1.add_subplot(212, aspect='equal')
        ax1.add_patch(patches.Rectangle((0, 0),  Sim.CCD.ccd_cfg.ccd_dimension_ac,  Sim.CCD.ccd_cfg.depletion_zone_thickness, alpha=0.3, edgecolor='none'))
        ax1.add_patch(patches.Rectangle((0, Sim.CCD.ccd_cfg.depletion_zone_thickness), Sim.CCD.ccd_cfg.ccd_dimension_ac, Sim.CCD.ccd_cfg.field_free_zone_thickness, alpha=0.45, edgecolor='none'))
        ax1.add_patch(patches.Rectangle((0, Sim.CCD.ccd_cfg.depletion_zone_thickness + Sim.CCD.ccd_cfg.field_free_zone_thickness), Sim.CCD.ccd_cfg.ccd_dimension_ac, Sim.CCD.ccd_cfg.substrate_thickness, alpha=0.7, edgecolor='none'))
        # ax2.add_patch(patches.Rectangle((0, 0), Sim.CCD.ccd_cfg.ccd_dimension_ac , Sim.CCD.ccd_cfg.ccd_dimension_al, alpha=0.3, edgecolor='none'))

        for i in range(0, np.size(x_p)):
            ax1.add_patch(patches.Circle((x_p[i], z_p[i]), sig[i], alpha=0.1, color='#2EB4BF'))
            ax1.add_patch(patches.Circle((x_p[i], z_p[i]), 0.2, color='k', alpha=(en_p[i]/np.max(en_p))))
            if values is True:
                ax1.text(x_p[i], z_p[i] + 0.3, str((en_p[i])) + ' MeV', color='k', fontsize=10, alpha=0.2)
                ax1.text(x_p[i], z_p[i] + 0.5, str((e_p[i])) + ' e-', color='#209FBF', fontsize=8, alpha=0.7)
            ax1.set_xlim([np.min(x_p) - 5, np.max(x_p) + 5])
            ax1.set_ylim([np.max(z_p) + 5, np.min(z_p) - 5])
            ax1.set_xlabel("Along-scan dimension")
            ax1.set_ylabel("Depth")
            ax1.set_title("Particle spreading within CCD")

        if values is True:
            ax1.text(x_p[i] - 3, z_p[i] - sig[i], 'Total : ' + str(int(f_e_p)) + ' e-')


            # ax2.plot(x_p[i], y_p[i], 'k', color='b')
            # ax2.add_patch(patches.Circle((x_p[i], y_p[i]), 0.3, alpha=(en_p[i] / np.max(en_p))))
            # ax2.text(x_p[i], y_p[i] + 0.3, str(int(en_p[i])) + ' MeV')

        for i in range(1, np.size(x_p)):
            ax1.arrow(x_p[i], z_p[i], x_p[i-1] - x_p[i], z_p[i-1] - z_p[i], head_width=0.05, head_length=0.1, fc='w', ec='w')


        plt.show()

    else:

        unique_trace= figure(width=1500,
                             plot_height=800,
                             title='Particle spreading within CCD',
                             tools="pan,wheel_zoom,save,reset",
                             x_axis_label='Along-scan dimension',
                             y_axis_label="Depth")

        b1 = BoxAnnotation(top=Sim.CCD.ccd_cfg.depletion_zone_thickness+Sim.CCD.ccd_cfg.field_free_zone_thickness+Sim.CCD.ccd_cfg.substrate_thickness,
                           bottom=Sim.CCD.ccd_cfg.field_free_zone_thickness + Sim.CCD.ccd_cfg.substrate_thickness,
                           left=0,
                           right=Sim.CCD.ccd_cfg.ccd_dimension_ac ,
                           fill_alpha=0.3,
                           fill_color='#578BB5')

        b2 = BoxAnnotation(top=Sim.CCD.ccd_cfg.field_free_zone_thickness + Sim.CCD.ccd_cfg.substrate_thickness,
                           bottom=Sim.CCD.ccd_cfg.substrate_thickness,
                           left=0,
                           right=Sim.CCD.ccd_cfg.ccd_dimension_ac ,
                           fill_alpha=0.50,
                           fill_color='#578BB5')

        b3 = BoxAnnotation(top=Sim.CCD.ccd_cfg.substrate_thickness,
                           bottom=0,
                           left=0,
                           right=Sim.CCD.ccd_cfg.ccd_dimension_ac ,
                           fill_alpha=0.60,
                           fill_color='#578BB5')

        unique_trace.add_layout(b1)
        unique_trace.add_layout(b2)
        unique_trace.add_layout(b3)

        for i in range(0, np.size(x_p)):
            unique_trace.circle(x_p[i], z_p[i], size=sig[i], color='#2EB4BF', alpha=0.2)
            unique_trace.circle(x_p[i], z_p[i], size=12, color="black", alpha=(en_p[i] / np.max(en_p)))
            label = Label(x=x_p[i], y=z_p[i], x_offset=12, text=(str(en_p[i]) + ' MeV'), text_baseline="middle", text_font_size='10px')
            unique_trace.add_layout(label)
            label = Label(x=x_p[i], y=z_p[i], x_offset=12, y_offset=12, text=(str(e_p[i]) + ' e-'),text_baseline="middle", text_font_size='10px', text_color="#2EB4BF")
            unique_trace.add_layout(label)

        for i in range(1, np.size(x_p)):
            unique_trace.add_layout(Arrow(end=NormalHead(fill_color="white", size=5, line_color="white"),
                                          x_start=x_p[i-1], y_start=z_p[i-1], x_end=x_p[i], y_end=z_p[i], line_width=1,
                                          line_color="white"))

        show(unique_trace)
