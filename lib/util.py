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
# Not fully documented
# not fully Commented

import os
import numpy as np

from matplotlib import pyplot as plt
from scipy import interpolate


def save(path, ext='png', close=True, verbose=True):
    """Save a figure from pyplot.
    Parameters
    ----------
    path : string
        The path (and filename, without the extension) to save the
        figure to.
    ext : string (default='png')
        The file extension. This must be supported by the active
        matplotlib backend (see matplotlib.backends module).  Most
        backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.
    close : boolean (default=True)
        Whether to close the figure after saving.  If you want to save
        the figure multiple times (e.g., to multiple formats), you
        should NOT close it in between saves or you will have to
        re-plot it.
    verbose : boolean (default=True)
        Whether to print information about when and where the image
        has been saved.
    """

    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
        directory = '.'

    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)

    # The final path to save to
    savepath = os.path.join(directory, filename)

    if verbose:
        print("Saving figure to '%s'..." % savepath),

    # Actually save the figure
    plt.savefig(savepath)

    # Close it
    if close:
        plt.close()

    if verbose:
        print("Done")


def get_float(text_file):
    """Get float separated by spaces and \n from a file.
        Parameter
        ----------
        text_file : file
            file opened with file=open(...)

    """

    l = ''
    s = ''

    test = False

    while not test:
        s = text_file.read(1)
        if str.isdigit(s):
            test = True

    while test:
        l = l + s
        s = text_file.read(1)

        if s is ' ':
            test = False

        elif s is 'e':
            l = l + s
            s = text_file.read(1)
            if s is '-':
                l = l + s
                s = text_file.read(1)

        if not str.isdigit(s) and s is not '.':
            test = False

    return float(l)


def binning(image, dim, bin_fact):

    if dim is 1:
        binned_image = np.zeros((np.size(image, 0), int(np.size(image, 1) / bin_fact)))

        for j in range(0, int(np.size(image, 1) / bin_fact) - bin_fact, 1):
            binned_image[:, j] = np.sum(image[:, j * bin_fact:(j + 1) * bin_fact], 1)

    if dim is 0:
        binned_image = np.zeros((np.size(image, 0) / bin_fact, int(np.size(image, 1))))

        for j in range(0, int(np.size(image, 0) / bin_fact) - bin_fact, 1):
            binned_image[j, :] = np.sum(image[j * bin_fact:(j + 1) * bin_fact, :], 0)

    return binned_image


def read_doc(file_name):

    data = np.loadtxt(file_name, 'float', '#')

    data_function = interpolate.interp1d(data[:, 0], data[:, 1], kind='linear')

    return data, data_function
