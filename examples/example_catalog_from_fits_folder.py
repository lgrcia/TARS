from pythagor.bench.lgarcia.TARS.lib import ccd_events_catalog as ccd_events_catalog

#   Creation of a catalog from folder
CR_cat_Gaia = ccd_events_catalog.from_folder_to_cat("../../../../../Gaia_bam_data")

#   plot of electrons repartition histogram of the events catalog
CR_cat_Gaia.cat_histogram('all')

end = True