from pythagor.bench.lgarcia.TARS.lib import ccd_events_catalog
from pythagor.bench.lgarcia.TARS.lib import cosmics

array, h = cosmics.fromfits("Data/Gaia images from bam/CR_000017380011161851100_2014-05-18.fits")

c = cosmics.cosmicsimage(array, gain=3.853, readnoise=10.0, sigclip=5.0, sigfrac=0.3, objlim=1.0, verbose=False)

c.run(maxiter=4)

new_catalog = ccd_events_catalog.from_image_to_cat(array, c.mask, "CR_extraction", local_mean=True)

ccd_events_catalog.make_nice_cat_plot(array, new_catalog)

end = True