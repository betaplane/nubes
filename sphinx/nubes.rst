*******************************
Cloud visualization with Mayavi
*******************************

This package contains an object-oriented, interactive version of a cloud visualization with Mayavi_, in the :mod:`.core` module, and a script version (file ``script.py``) which doesn't import anything from *this* module. On a computer with display and graphics drivers, the script version can be run (provided the `dependencies`_ are installed) as follows::

  wrf_file=<input_wrfout_file> mayavi2 -x script.py -o

where ``<input_wrfout_file>`` is to be replaced by the WRF output file on which to base the visualization. This will intermittently create a sequence of images files named ``scene_xxx.png`` where ``xxx`` is a 3-digit number indicating the order of each file in the sequence. This files will be deleted once a movie has been created from them via MoviePy_.

For portability and completely **offcreen** rendering, I have created a Singularity_ image which can be used in the following way - assuming, of course, that Singularity_ is installed::

  singularity run mayavi2.img script.py <input_wrfout_file>

with ``<input_wrfout_file>`` to be replaced as above.

.. Warning::
   The file specified as ``<input_wrfout_file>`` (as well as ``script.py``) has to be `reachable from within the singularity container <https://www.sylabs.io/guides/2.6/user-guide/quick_start.html#working-with-files>`_. In the case of ``script.py``, this is going to be the case if Singularity_ is run wihtout ``sudo`` and it lives somewhere in the user directory tree.

.. Note::

   I have currently two Singularity_ images, ``mayavi.img`` and ``mayavi2.img`` (both squashed), where the second one contains all the additional sections (runscript, environment, help) and the first the actual OS and installations.

   Inside the image, which is based on the intelpython/intelpython3_core docker image, there is a Conda_ environment called ``mayavi`` that contains all the necessary python modules.

Dependencies
------------

* Mayavi_
* numpy
* xarray
* wrf-python
* pyproj
* netCDF4_
* MoviePy_

Building Singularity
--------------------

If squashed images are desired, there is one dependency: `squashfs <http://squashfs.sourceforge.net/>`_, which at the time of this writing is simply built with the Makefile in ``squashfs-tools/``. (There was one other dependency I don't remember, maybe `libarchive <https://www.libarchive.org/>`_.) With the source from github::

  git clone https://github.com/sylabs/singularity.git

my configure line on ``diaguita`` was::

  CPPFLAGS="-I/HPC/arno/dia/gcc/include/" LDFLAGS="-L/HPC/arno/dia/gcc/lib" LD_LIBRARY_PATH=/HPC/arno/dia/gcc/lib ./configure --prefix=/HPC/arno/dia/gcc/ --sysconfdir=/HPC/arno/dia/gcc/etc/ --with-mksquashfs=/HPC/arno/dia/gcc/bin/mksquashfs


.. automodule:: core
   :members:

.. _netCDF4: http://unidata.github.io/netcdf4-python/

.. _SRTM: https://www2.jpl.nasa.gov/srtm/

.. _Blue Marble: https://visibleearth.nasa.gov/view_cat.php?categoryID=1484

.. _Mayavi: http://docs.enthought.com/mayavi/mayavi/index.html

.. _VTK: https://www.vtk.org/

.. _MoviePy: http://zulko.github.io/moviepy/

.. _Singularity: https://www.sylabs.io/singularity/

.. _Conda: https://conda.io/docs/

Links
-----

* `image textures <http://geoexamples.blogspot.com/2014/02/3d-terrain-visualization-with-python.html>`_

Issues and bug workarounds
^^^^^^^^^^^^^^^^^^^^^^^^^^

* `failure to illuminate rendered volume (i.e. clouds) before user interacts with the window <https://github.com/enthought/mayavi/issues/7>`_
  
  * a possible workaround in case this becomes relevant (it seems not to affect the 'scripted' animation::

      from tvtk.pyface.ui.qt4 import scene
      from tvtk.pyface.ui.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

      class tvtkBug(scene._VTKRenderWindowInteractor):
          def paintEvent(self, e):
          QVTKRenderWindowInteractor.paintEvent(self, e)
          scene._VTKRenderWindowInteractor = tvtkBug

      scene._VTKRenderWindowInteractor = tvtkBug

* main VTK_ bug resulting in the error message (and possibly a seg fault)::

    ERROR: In /root/VTK-8.1.1/Rendering/OpenGL2/vtkOpenGLRenderWindow.cxx, line 797
    vtkOSOpenGLRenderWindow (0x562b86863710): GL version 2.1 with the gpu_shader4 extension is not supported by your graphics driver but is required for the new OpenGL rendering backend. Please update your OpenGL driver. If you are using Mesa please make sure you have version 10.6.5 or later and make sure your driver in Mesa supports OpenGL 3.2.

  * https://www.vtk.org/pipermail/vtk-developers/2017-November/035592.html
  * http://vtk.1045678.n5.nabble.com/Offscreen-rendering-problems-on-headless-Ubuntu-td5746035.html

* `installer doesn't find source-compiled VTK <https://github.com/enthought/mayavi/issues/652>`_

  * solution is to delete 'vtk' from ``__requires__`` in the __init__.py file of the Mayavi_ source before installing with ``python setup.py install``
  * e.g.::

      git clone https://github.com/enthought/mayavi.git
      cd mayavi && git checkout 4.6.2 && sed -i '/vtk/d' mayavi/__init__.py && python setup.py install
