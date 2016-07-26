.. _ccd_config.rst:

CCD configuration file
==========================

.. raw:: html

    <h5><b>ccd_configuration_file.cfg</b><h6><b>TARS</b> &#9658; <b>Data</b> &#9658; <b>Inputs</b> &#9658; <b>CCD specs</b></h6></h5><br>

ccd_configuration_file.cfg is a configuration files (.cfg) that you can directly change to configure your ccd specs for your simulation

Here is the content of this configuration file :

.. raw:: html

    <script type="text/x-mathjax-config">
    MathJax.Hub.Config({
    tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
    });
    </script>
    <script type="text/javascript" async src="path-to-mathjax/MathJax.js?config=TeX-AMS_CHTML"></script>

    <br>
    <div class="row">
        <div class="col-md-6">
            <table class="table table-striped table-hover ">
                <tr>
                    <th>CCD CONFIG</th>
                    <th> </th>
                </tr>
                <tr>
                    <td>number_of_pixels_ac</td>
                    <td>Number of pixels in the across-scan direction </td>
               </tr>
               <tr>
                  <td>number_of_pixels_al</td>
                   <td>Number of pixels in the along-scan direction  </td>
               </tr>
               <tr>
                  <td>temperature</td>
                   <td> CCD temperature ($K$) </td>
               </tr>
               <tr>
                  <td>depletion_zone_thickness</td>
                   <td> Depletion zone thickness ($\mu m $)</td>
               </tr>
               <tr>
                  <td>field-free_zone_thickness</td>
                   <td> Field free zone thickness ($\mu m $) </td>
               </tr>
               <tr>
                  <td>substrate_thickness</td>
                   <td> Substrate thickness ($\mu m $) </td>
               </tr>
               <tr>
                  <td>pixel_ac_size</td>
                   <td> Pixel across-scan size ($\mu m $) </td>
               </tr>
               <tr>
                  <td>pixel_al_size</td>
                   <td> Pixel along-scan size ($\mu m $) </td>
               </tr>
               <tr>
                  <td>electrons_saturation</td>
                   <td> Electrons saturation (number of $e-$)  </td>
               </tr>
            </table>
        </div>
    </div>
    <br>

.. raw:: html

    <div class="row">
        <div class="col-md-6">

Here is an example of some simulation and ccd configuration (CCD specs taken from Gaia BAM *CCD91-72* ) :

.. raw:: html

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
    </div>


