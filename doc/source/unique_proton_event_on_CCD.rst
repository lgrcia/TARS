.. _unique_proton_event_on_CCD.rst:

Unique random proton event on CCD
==========================

.. raw:: html

    <h2>Initialisation</h2>
    <br>

First we have to import all the useful tools and class provided by TARS


.. code-block:: python

   import TARS.lib.util

   from TARS.lib.simulation import Simulation
   from TARS.lib.ccd import CCD
   from TARS.lib import nice_figure

.. raw:: html

    <h2>Inputs</h2>
    <br>

Then we can create all the important objects which will constitute our simulation. For each of this object we specify the proper configuration file path.

.. code-block:: python

   #  Specify the input files paths
   new_sim_folder_name = "test_unique_proton_on_CCD"

   simulation_config = '../Data/Inputs/Simulation config/input_simulation_config.cfg'

   ccd_config = '../Data/Inputs/CCD specs/ccd_configuration_file.cfg'

   stopping_power_file = "../Data/Inputs/CCD specs/stopping_in.dat"

   spectre_file = "../Data/Inputs/Proton spectres/proton_L2_solarMax_NoShielding.txt"

   #  CCD creation and setup
   Gaia_CCD = CCD(ccd_config)
   Gaia_CCD.stopping_power_tab, Gaia_CCD.stopping_power_function = util.read_doc(stopping_power_file)

   #  Simulation creation and setup
   CR_Sim = Simulation(Gaia_CCD, "cr_simulation_test", simulation_config, new_sim_folder_name)
   CR_Sim.particle_spec.add_spectre_from_file(spectre_file)

   # Simulation run
   CR_Sim.run()



Indeed, all the simulation parameters and inputs are situated in configuration files (.cfg) or text files (.txt) that you can directly edit to configure your simulation.

Here are all the kind of inputs we need for this example :

.. toctree::
   :maxdepth: 1

   ccd_config
   simulation_config
   stopping_power
   input_spectrum
   api
   tutorial


Here is the simulation and ccd configuration we are going to use for our single particle simulation (CCD specs taken from Gaia BAM *CCD91-72* ) :


.. raw:: html

    <br>
    <div class="container-fluid">
        <div class="col-md-6">
            <h5><b>ccd_configuration_file.cfg</b><h6><b>examples</b> &#9658; <b>Data</b> &#9658; <b>Inputs</b> &#9658; </h6></h5>
            <br>
             <table class="table table-striped table-hover ">
               <tr>
                  <th>CCD CONFIG</th>
                   <th> </th>
               </tr>
               <tr>
                  <td>number_of_pixels_ac</td>
                   <td>4200</td>
               </tr>
               <tr>
                  <td>number_of_pixels_al</td>
                   <td>1956</td>
               </tr>
               <tr>
                  <td>temperature</td>
                   <td>163.</td>
               </tr>
               <tr>
                  <td>depletion_zone_thickness</td>
                   <td>38.00</td>
               </tr>
               <tr>
                  <td>field-free_zone_thickness</td>
                   <td>2.00</td>
               </tr>
               <tr>
                  <td>substrate_thickness</td>
                   <td>0.0000001</td>
               </tr>
               <tr>
                  <td>pixel_ac_size</td>
                   <td>30.0</td>
               </tr>
               <tr>
                  <td>pixel_al_size</td>
                   <td>10.0</td>
               </tr>
               <tr>
                  <td>electrons_saturation</td>
                   <td>350000.</td>
               </tr>
            </table>
        </div>

         <div class="col-md-6">
            <h5><b>input_simulation_config.cfg</b><h6><b>examples</b> &#9658; <b>Data</b> &#9658; <b>Inputs</b> &#9658; </h6></h5>
            <br>
         <table class="table table-striped table-hover ">
            <tr>
               <th>SIMULATION CONFIG</th>
                <th> </th>
            </tr>
            <tr>
               <th class="info">
               number_of_particles</th>
                <th class="info">1</th>
            </tr>
            <tr>
               <td>energy</td>
                <td>random</td>
            </tr>
            <tr>
               <td>position_x</td>
                <td>random</td>
            </tr>
            <tr>
               <td>position_y</td>
                <td>random</td>
            </tr>
            <tr>
               <td>alpha_angle</td>
                <td>random</td>
            </tr>
            <tr>
               <td>beta_angle</td>
                <td>random</td>
            </tr>
            <tr>
               <td>spreading_step</td>
                <td>0.5</td>
            </tr>
            <tr>
               <th>INPUTS FILE</th>
                <th> </th>
            </tr>
            <tr>
               <td>input_file</td>
                <td>no</td>
            </tr>
            <tr>
               <td>positions</td>
                <td>test_inputs_positions.txt</td>
            </tr>
            <tr>
               <td>energies</td>
                <td>test_inputs_energy.txt</td>
            </tr>
            <tr>
               <th>SIMULATION REPORT</th>
                <th> </th>
            </tr>
            <tr>
               <td>date</td>
                <td>0</td>
            </tr>
            <tr>
               <td>processing_time</td>
                <td>0</td>
            </tr>
         </table>
    </div>


.. raw:: html

    <h2>Running and results</h2>
    <br>

We can then run the simulation

.. code-block:: python

   # Simulation run
   CR_Sim.run()

   #  Creation of a result folder
   nice_figure.create_result_folder(new_sim_folder_name, CR_Sim)

A plot of the trace of your unique particle in the CCD will then automatically be plotted (on your web browser via bokeh or on a new window via matplotlib)

The results will be then saved in a folder. Here is a way to see the results of your simulation at any time without running again the simulation:

.. code-block:: python

   #  Complete simulation results and plots
   nice_figure.display_simulation_report(complete_path_of_you_folder, bokeh=True)

For more details about the simulation report folder see


.. raw:: html

    <h2>Complete code</h2>
    <br>

.. code-block:: python

   import TARS.lib.util as util
   from TARS.lib.simulation import Simulation
   from TARS.lib.ccd import CCD
   from TARS.lib import nice_figure


   new_sim_folder_name = "test_before_committing"

   simulation_config = '../Data/Inputs/Simulation config/input_simulation_config.cfg'

   ccd_config = '../Data/Inputs/CCD specs/ccd_configuration_file.cfg'

   stopping_power_file = "../Data/Inputs/CCD specs/stopping_in.dat"

   spectre_file = "../Data/Inputs/Proton spectres/proton_L2_solarMax_NoShielding.txt"

   Gaia_CCD = CCD(ccd_config)
   Gaia_CCD.stopping_power_tab, Gaia_CCD.stopping_power_function = util.read_doc(stopping_power_file)

   CR_Sim = Simulation(Gaia_CCD, "cr_simulation_test", simulation_config, new_sim_folder_name)
   CR_Sim.particle_spec.add_spectre_from_file(spectre_file)

   CR_Sim.run()

   nice_figure.create_result_folder(new_sim_folder_name, CR_Sim)
