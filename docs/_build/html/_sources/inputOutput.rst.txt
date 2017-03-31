.. inputOutput-label:
Input/output files
===================

Both input and output files for Hazel are NetCDF files.

Input files
-----------

The input file constains the observations and information about the
observing position and boundary condition. The file consists of the
following variables:

-  lambda: vector of size *nlambda* containing the wavelength axis with
   respect to the center of the multiplet.

-  map: array of size *(npixel,8,nlambda)* containing the Stokes vector
   :math:`(I,Q,U,V)` and the associated standard deviation of the noise
   :math:`(\sigma_I,\sigma_Q,\sigma_U,\sigma_V)`.

-  boundary: array of size *(npixel,4)* containing the boundary
   condition for every inverted pixel.

-  height: vector of size *npixel* which contains the height of the
   slabs for every pixel.

-  obs\_theta: vector of size *npixel* which contains the observing
   angle :math:`\theta` for every pixel.

-  obs\_gamma: vector of size *npixel* which contains the observing
   angle :math:`\gamma` that defines the positive reference for Stokes
   :math:`Q` for every pixel.

-  mask: array of size *nx,ny* which tells whether this pixel will be
   inverted.

-  normalization: variable indicating whether the profiles are
   normalized to the peak amplitude or the continuum of Stokes
   :math:`I`.

-  pars: array of size *npixel,npars* which contains the initial value
   for the model parameters. These will be used to reinvert some pixels
   or, for instance, to refine the ambiguous solutions.

The routine ``gen_netcdf.pro`` on the directory ``IDL_routines`` and the
``genNetCDF.py`` on ``pyRoutines`` shows functions that generate such a
file by passing all the variables as parameters. The order of pars is
the following, depending on the number of slabs:

-  1-component (vector of size 8): :math:`B`, :math:`\theta_B`,
   :math:`\chi_B`, :math:`\tau`, :math:`v_\mathrm{dop}`, :math:`a`,
   :math:`v_\mathrm{mac}`, :math:`\beta`

-  2-component 1+1 with same field (vector of size 11): :math:`B`,
   :math:`\theta_B`, :math:`\chi_B`, :math:`\tau_1`, :math:`\tau_2`,
   :math:`v_\mathrm{dop}`, :math:`a`, :math:`v_\mathrm{mac1}`,
   :math:`v_\mathrm{mac2}`, :math:`\beta`, :math:`\beta_2`

-  2-component 1+1 with different field (vector of size 15):
   :math:`B_1`, :math:`\theta_{B1}`, :math:`\chi_{B1}`, :math:`B_2`,
   :math:`\theta_{B2}`, :math:`\chi_{B2}`, :math:`\tau_1`,
   :math:`\tau_2`, :math:`v_\mathrm{dop}`, :math:`v_\mathrm{dop2}`,
   :math:`a`, :math:`v_\mathrm{mac1}`, :math:`v_\mathrm{mac2}`,
   :math:`\beta`, :math:`\beta_2`

-  2-component 2 with different field with filling factor (vector of
   size 16): :math:`B_1`, :math:`\theta_{B1}`, :math:`\chi_{B1}`,
   :math:`B_2`, :math:`\theta_{B2}`, :math:`\chi_{B2}`, :math:`\tau_1`,
   :math:`\tau_2`, :math:`v_\mathrm{dop}`, :math:`v_\mathrm{dop2}`,
   :math:`a`, :math:`v_\mathrm{mac1}`, :math:`v_\mathrm{mac2}`,
   :math:`\mathrm{ff}`, :math:`\beta`, :math:`\beta_2`

Output files
------------

The results of the inversion are saved on two files defined on the
configuration file. The file with the inverted
profiles contains the following variables:

-  lambda: vector of size *nlambda* containing the wavelength axis with
   respect to the center of the multiplet.

-  map: array of size *(npixel,4,nlambda)* containing the synthetic
   Stokes vector :math:`(I,Q,U,V)` for every pixel.

The file with the inverted parameters contains the following variable:

-  map: array of size *(npixel,ncolumns)* containing the parameters of
   the inversion.

The number of columns depends on the selected model:

-  One-slab: nine columns with the vector
   :math:`(B,\theta_B,\chi_B,h,\tau,v_\mathrm{th},a,v_\mathrm{mac},\beta)`.

-  Two-slab with same magnetic field: eleven columns with the vector
   :math:`(B,\theta_B,\chi_B,h,[\tau]_1,[\tau]_2,v_\mathrm{th},a,[v_\mathrm{mac}]_1,
   [v_\mathrm{mac}]_2,\beta,\beta_2)`.

-  Two-slab with different magnetic field: fifteen columns with the
   vector
   :math:`([B]_1,[\theta_B]_1,[\chi_B]_1,[B]_2,[\theta_B]_2,[\chi_B]_2,
   h,[\tau]_1,[\tau]_2,[v_\mathrm{th}]_1,[v_\mathrm{th}]_2,a,[v_\mathrm{mac}]_1,[v_\mathrm{mac}]_2,\beta,\beta_2)`.

The file ``read_results.pro`` on the ``RunMPI`` directory shows how to
read the files from IDL.