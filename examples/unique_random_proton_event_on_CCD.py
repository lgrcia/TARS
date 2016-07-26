import pythagor.bench.lgarcia.TARS.lib.util as util
from pythagor.bench.lgarcia.TARS.lib.simulation import Simulation
from pythagor.bench.lgarcia.TARS.lib.ccd import CCD
from pythagor.bench.lgarcia.TARS.lib import nice_figure


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

end = True
