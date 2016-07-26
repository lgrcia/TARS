.. _simulation_config.rst:

Simulation configuration file
==========================

.. raw:: html

   <h5><b>input_simulation_config.cfg</b><h6><b>TARS</b> &#9658; <b>Data</b> &#9658; <b>Inputs</b> &#9658; <b>Simulation config</b></h6></h5><br>

input_simulation_config.cfg is a configuration files (.cfg) that you can directly change to configure your simulation

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
        <div class="col-md-8">
            <table class="table table-striped table-hover ">
               <tr>
                  <th>SIMULATION CONFIG</th>
                   <th> </th>
                   <th> </th>
               </tr>
               <tr>
                  <td>number_of_particles</td>
                  <td> <i>integer </i> </td>
                   <td> Total number of incident particle on CCD </td>
               </tr>
               <tr>
                    <td>energy</td>
                    <td> <i>"random"</i> <br> <i>float</i> <br> <i>range()</i> <br> <i>np.logspace()</i> </td>
                    <td> Energies of the particles </td>
               </tr>
               <tr>
                  <td>position_x</td>
                  <td> <i> float </i> </td>
                   <td> across_scan positions of the particles ($\mu m $) </td>
               </tr>
               <tr>
                  <td>position_y</td>
                  <td><i> float </i></td>
                   <td> along_scan positions of the particles ($\mu m $) </td>
               </tr>
               <tr>
                  <td>alpha_angle</td>
                  <td><i> float </i></td>
                   <td> alpha incident angle of the particles ($ rad $) </td>
               </tr>
               <tr>
                  <td>beta_angle</td>
                  <td><i> float </i></td>
                   <td> beta incident angle of the particles ($ rad $) </td>
               </tr>
               <tr>
                  <td>spreading_step</td>
                  <td><i> float </i></td>
                   <td> spreading step in CCD </td>
               </tr>
               <tr>
                  <th>INPUTS FILE</th>
                   <th> </th>
                   <th> </th>
               </tr>
               <tr>
                  <td>input_file</td>
                   <td><i>['yes'/'no']</i></td>
                   <td> input file containing a pre-computed spreading of the particule <td>
               </tr>
               <tr>
                  <td>positions</td>
                  <td><i> file_path (string) </i></td>
                   <td> file of positions pre-computed</td>
               </tr>
               <tr>
                  <td>energies</td>
                  <td><i> file_path (string) </i></td>
                   <td> file of energies pre-computed</td>
               </tr>
               <tr>
                  <th>SIMULATION REPORT</th>
                   <th> </th>
                   <th> </th>
               </tr>
               <tr>
                  <td>date</td>
                   <td><i> year-month-day </i></td>
                   <td> date of the simulation - written after simulation </td>
               </tr>
               <tr>
                  <td>processing_time</td>
                  <td><i> hours-minutes-seconds </i></td>
                   <td> processing time of the simulation - written after simulation </td>
               </tr>
            </table>
        </div>

        <div class="col-md-4">

If *energy* is:

``float``   : the simulator will generate *number_of_particles* particles with energy *energy*

``range`` or ``logspace``   : the simulator will generate particles with energy *energy* following the **range** or **logspace**

``'random'``  : the simulator will generate *number_of_particles* particles with random *energy* following the :ref:`input spectrum <input_spectrum.rst>`

----------------

For *input_file* see
:ref:`information about external inputs <input_file_from_external_source.rst>`


.. raw:: html

            </p>
        </div>
    </div>




