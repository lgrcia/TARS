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

import configparser


class CCD_config:
    """
    CCD config is a class making the interface between a configuration file containing all the CCD parameters
    and the simulation program
    """

    def __init__(self, config_file):
        """
        :param string config_file: path of the configuration file
        """

        # reading of the simulation specs
        config = configparser.ConfigParser()
        config.read_file(open(config_file))
        self.cfg_file = config
        self.real_cfg_file = config_file
        self.number_of_pixels_ac = eval(config.get('CCD CONFIG', 'number_of_pixels_ac'))
        self.number_of_pixels_al = eval(config.get('CCD CONFIG', 'number_of_pixels_al'))
        self.temperature = eval(config.get('CCD CONFIG', 'temperature'))
        self.depletion_zone_thickness = eval(config.get('CCD CONFIG', 'depletion_zone_thickness'))
        self.field_free_zone_thickness = eval(config.get('CCD CONFIG', 'field-free_zone_thickness'))
        self.substrate_thickness = eval(config.get('CCD CONFIG', 'substrate_thickness'))
        self.pixel_ac_size = eval(config.get('CCD CONFIG', 'pixel_ac_size'))
        self.pixel_al_size = eval(config.get('CCD CONFIG', 'pixel_al_size'))
        self.electrons_saturation = eval(config.get('CCD CONFIG', 'electrons_saturation'))

        self.total_thickness = self.depletion_zone_thickness + self.field_free_zone_thickness + self.substrate_thickness

        self.ccd_dimension_ac = self.number_of_pixels_ac * self.pixel_ac_size
        self.ccd_dimension_al = self.number_of_pixels_al * self.pixel_al_size

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


class CCD:
    """
    Creation of a simulated CCD following its configuration file
    """

    def __init__(self, CCD_cfg_file):
        """
        Creation of a simulated CCD following a configuration file containing all the specs
        :param string CCD_cfg_file: path of the specifications file of the CCD (configuration file)
        """

        self.ccd_cfg = CCD_config(CCD_cfg_file)

        self.stopping_power_tab = 0

        self.stopping_power_function = 0

        self.binning = 0