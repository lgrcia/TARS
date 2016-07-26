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

import numpy as np
import random
import glob

from matplotlib import pylab
from matplotlib import pyplot as plt
from scipy import ndimage

import pythagor.bench.lgarcia.TARS.lib.util as util
import pythagor.bench.lgarcia.TARS.lib.cosmics as cosmics
import matplotlib.patches as patches


class events_catalog:
    """
    events_catalog is a catalog containing all the events information

    Attributes
    ----------
    name : string
        name of the catalog

    events_images : list of np.ndarray
        list of events images

    number_of_events : integer
        number of images (called shapes) in the catalog

    events_dimensions : list of np.array([*, *])
        x,y dimensions of each shape
        example: shape_dim[[ dim x of shape 1 , dim y of shape 1 ],
                           [ dim x of shape 2 , dim y of shape 2 ],
                            ...]]

    events_positions_in_ccd : list of np.array([*, *])
        i,j position of each shapes into their original image
        example: shapes_position[[ position i of shape 1 , position j of shape 1 ],
                                 [ position i  of shape 2 , position j shape 2 ],
                                 ...]]

    electrons_generated : list of integer
        list of electrons generated

    CCD_image : np.ndarray
        image of the CCD containing all the events

    """

    def __init__(self, name):
        """
        Initialisation of the catalog
        :param string name: name of the catalog
        """

        self.name = name

        self.events_images = []

        self.number_of_events = 0

        self.events_dimensions = []

        self.events_positions_in_ccd = []

        self.electrons_generated = []

        self.CCD_image = 0

    def __getattr__(self, item):
        """
        get function of the catalog
        """
        return self.number_of_events

    def save_cat_as_image(self, files_name, how_much, outline=True):
        """
        saving all the events images in png files
        :param string files_name: prefix of the files containing the images
        :param int or 'all' how_much: number of events images to save
        :param bool outline: outline around the images
        """

        if how_much is 'all':
            how_much == self.number_of_events

        if outline is True:
            new_cat = [None] * how_much

            for i in range(0, how_much):

                if new_cat[i] is None:
                    new_cat[i] = np.zeros(((self.events_dimensions[i][0] + 8), (self.events_dimensions[i][1] + 8)))

                _i = int(self.events_dimensions[i][0] + 4)
                _j = int(self.events_dimensions[i][1] + 4)
                new_cat[i][4:_i, 4:_j] = self.events_images[i][:]

        if how_much != 'all':
            for i in range(0, how_much):
                name = files_name + str(i).zfill(4)
                pylab.imshow(np.rot90(new_cat[i]), cmap='Greys_r', interpolation='None')
                plt.colorbar()
                util.save(name, ext="png", close=True, verbose=False)

    def add_event(self, added_shape, position=None):
        """
        Add an event to the catalog
        :param np.ndarray(2*2 image matrix) added_shape: event_image to add
        :param np.array([*, *]) position: position of the event in the ccd
        """

        self.number_of_events += 1
        self.events_images.append(added_shape)
        self.events_dimensions.append([np.size(added_shape, 0), np.size(added_shape, 1)])
        if position is not None:
            self.events_positions_in_ccd.append(position)
        _k = added_shape.reshape((1, np.size(added_shape, 0) * np.size(added_shape, 1)))
        self.electrons_generated.append(np.sum(_k[np.where(_k > 0)]))

    def cat_histogram(self, how_much='all'):
        """
        Make an histogram of number of electrons repartition from catalog
        :param int or 'all' how_much: number of events to take in account in the histogram
        """

        if how_much is 'all':
            how_much = self.number_of_events

        fig = plt.figure()
        plt.style.use('ggplot')
        plt.hist(self.electrons_generated, bins=1000, histtype='step', label="cr_catS_sim")
        fig.suptitle('Particle Catalog Histogram (Simulated)', fontsize=14, fontweight='bold')
        plt.xlabel('number of e-')
        plt.ylabel('number of cr_cat')
        plt.show()

    def shape_minus_local_mean(self, array):
        """
        shape_minus_local_mean is a local background extraction method to clean all the extracted CCDs events from
        their background in the original image (usefull when CCDs events are extracted from an original image)

        :param np.ndarray (2*2 image matrix) array: image to clean
        :return:
        """
        """
            Inputs
            --------
            array:
                original image
            mask : np.ndarray (binary 2*2 image matrix containing 0 or 1)
            mask of shape detection

        """
        _i = int(np.size(array, 0) + 4)
        _j = int(np.size(array, 1) + 4)

        new_array = np.random.randint(1490, high=1510, size=(np.size(array, 0) + 8, np.size(array, 1) + 8))

        new_array[4:_i, 4:_j] = np.copy(array)

        for i in range(0, self.number_of_events):
            _i = int(self.events_dimensions[i][0] + 2)
            _j = int(self.events_dimensions[i][1] + 2)

            shape_extended_local = np.copy(new_array[(self.events_positions_in_ccd[i][0] + 2):(self.events_positions_in_ccd[i][0] + _i + 4),
                                           (self.events_positions_in_ccd[i][1] + 2):(self.events_positions_in_ccd[i][1] + _j + 4)])

            shape_extended_local[2:_i, 2:_j] -= self.events_images[i].astype('uint16')

            new_shape = self.events_images[i] - np.mean(shape_extended_local[shape_extended_local > 0])

            new_shape[np.where(new_shape < 0)] = 0

            self.events_images[i] = new_shape


def from_image_to_cat(image, mask, name, local_mean=True):
    """
    fromImageToCat extract a catalog from an image and an associated mask
    :param np.ndarray image: original image containing the shapes to extract
    :param np.ndarray mask: binarize version of image showing the most clearly the shape delimitation
    :param string name: Name of the created catalog
    :param bool local_mean: subtraction of the local mean of the event
    :return events_catalog : catalog from image
    """

    s_cat = events_catalog(name)

    labeled_image, n_shape = ndimage.label(mask)

    shapes = [None] * n_shape

    minmax = np.zeros((n_shape, 4))
    minmax[:, :] = np.array([100000, -1, 100000, -1])

    for i in range(0, np.size(mask, 0)):

        for j in range(0, np.size(mask, 1)):

            if mask[i, j] != 0:

                idx = labeled_image[i, j] - 1

                if shapes[idx] is None:
                    shapes[idx] = []

                shapes[idx].append([i, j])

                if i < minmax[idx][0]:
                    minmax[idx][0] = i
                if i > minmax[idx][1]:
                    minmax[idx][1] = i
                if j < minmax[idx][2]:
                    minmax[idx, 2] = j
                if j > minmax[idx][3]:
                    minmax[idx][3] = j

    for k in range(0, n_shape):

        shape_dim = np.array([(minmax[k][1] - minmax[k][0] + 1), (minmax[k][3] - minmax[k][2] + 1)])

        positions = np.array([minmax[k][0], minmax[k][2]])

        out = np.zeros((int(shape_dim[0]), int(shape_dim[1])))

        for i in range(0, int(shape_dim[0])):
            for j in range(0, int(shape_dim[1])):

                if [positions[0] + i, positions[1] + j] in shapes[k]:
                    out[i, j] = 1

        _s = np.multiply(out, image[positions[0]:(positions[0] + shape_dim[0]),
                                    positions[1]:(positions[1] + shape_dim[1])])

        s_cat.add_event(_s, positions)

    if local_mean is True:
        s_cat.shape_minus_local_mean(image)

    return s_cat


def make_outline(image):

    new_image = np.zeros(((np.size(image, axis=0) + 8), (np.size(image, axis=1)+ 8)))

    _i = int(np.size(image, axis=0) + 4)
    _j = int(np.size(image, axis=1) + 4)
    new_image[4:_i, 4:_j] = image[:]

    return new_image


def from_file_to_cat(fileName, how_much):
    """
    fromFileToCat extract a catalog from a file
    :param string fileName: path of the Folder/file containing all the data to extract
    :param int or string how_much: how_much = n : Extraction of n shapes from the file
            how_much = 'all' : = Extraction of all the shapes from the file
    :return: events_catalog : catalog from file
    """

    s_cat = events_catalog(fileName)

    cat_file = open(fileName, 'r')

    dim = np.zeros(2, dtype=int)

    s_dim = np.zeros(2, dtype=int)

    if how_much != 'all':

        for k in range(0, how_much):

            util.get_float(cat_file)

            dim[0] = int(util.get_float(cat_file))

            dim[1] = int(util.get_float(cat_file))

            s_dim[0] = dim[0]
            s_dim[1] = dim[1]

            cat = np.zeros((dim[0], dim[1]))

            for i in range(0, dim[0]):
                for j in range(0, dim[1]):
                    cat[i, j] = util.get_float(cat_file)

            s_cat.add_event(cat)

    cat_file.close()

    return s_cat


def from_folder_to_cat(Folder):
    """
    extraction of a catalog from fits in a folder
    :param string Folder: path of the folder containing the fits
    :return: events_catalog : catalog
    """

    crs = events_catalog("CRS")

    for file in glob.glob(Folder + "/*.fits"):

        print(glob.glob(Folder + "/*.fits"))

        array, header = cosmics.fromfits(file)

        c = cosmics.cosmicsimage(array, gain=3.853, readnoise=10.0, sigclip=5.0, sigfrac=0.3, objlim=1.0, verbose=False)

        c.run(maxiter=4)

        new_catalog = from_image_to_cat(array, c.mask, "CR_extraction", local_mean=True)

        for j in range(0, new_catalog.n_shapes):
            crs.add_event(new_catalog.events_images[j], new_catalog.events_positions_in_ccd[j])

    return crs


def make_nice_cat_plot(array, cat):

    fig = plt.figure(facecolor='white', figsize=(30, 10))
    plt.style.use('ggplot')

    ax = plt.subplot(211)
    ax.set_title('Incident particles energy spectrum', fontsize=10)
    ax.set_xlabel('AC dimension (' + r'$\mu$' + ')', fontsize=10)
    ax.set_ylabel('AL dimension (' + r'$\mu$' + ')', fontsize=10)
    plt.imshow(array, cmap='Greys_r', interpolation='None', aspect='auto', vmin=np.mean(array), vmax=np.mean(array)+200)

    for i in range(0, cat.number_of_events):
        ax.add_patch(patches.Rectangle((cat.events_positions_in_ccd[i][1]-0.5, cat.events_positions_in_ccd[i][0]-0.5),
                                        cat.events_dimensions[i][1], cat.events_dimensions[i][0],
                                       fill=None, edgecolor='r', linewidth=1))

    rand_idx = int(random.random()*(cat.number_of_events-1))

    ax.add_patch(patches.Rectangle((cat.events_positions_in_ccd[rand_idx][1] - 0.5, cat.events_positions_in_ccd[rand_idx][0] - 0.5),
                                   cat.events_dimensions[rand_idx][1], cat.events_dimensions[rand_idx][0],
                                   fill=None, edgecolor='g', linewidth=1))
    ax1 = plt.subplot(212)
    ax1.set_title("Example of extracted CR", fontsize=10)
    plt.imshow(make_outline(cat.events_images[rand_idx]), cmap='Greys_r',
               interpolation='None', aspect='auto', vmin=0, vmax=np.mean(cat.events_images[5]))
    plt.show()


