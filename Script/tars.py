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
# Not fully commented

import pythagor.bench.lgarcia.TARS.lib.util as util
from pythagor.bench.lgarcia.TARS.lib.simulation import Simulation
from pythagor.bench.lgarcia.TARS.lib.ccd import CCD
from pythagor.bench.lgarcia.TARS.lib import nice_figure


new_sim_folder_name = "simulator_profile"

simulation_config = '../Data/Inputs/Simulation config/input_simulation_config.cfg'

ccd_config = '../Data/Inputs/CCD specs/ccd_configuration_file.cfg'

stopping_power_file = "../Data/Inputs/CCD specs/stopping_in.dat"

spectre_file = "../Data/Inputs/Proton spectres/_2ideal_spectre.txt"

Gaia_CCD = CCD(ccd_config)
Gaia_CCD.stopping_power_tab, Gaia_CCD.stopping_power_function = util.read_doc(stopping_power_file)

CR_Sim = Simulation(Gaia_CCD, "cr_simulation_test", simulation_config, new_sim_folder_name)
CR_Sim.particle_spec.add_spectre_from_file(spectre_file)

CR_Sim.run()

nice_figure.create_result_folder(new_sim_folder_name, CR_Sim)

end = True
