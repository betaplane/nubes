"""
Utility routines
================

.. NOTE::
    * Currently projection is hard-coded to 'lcc' (could also be read from netcdf file)

"""
#!/usr/bin/env python
import os
from importlib import import_module
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from warnings import catch_warnings, simplefilter

# hack for GDAL framework installation on MacOS
# import sys
# sys.path.insert(0, '/usr/local/lib/python3.6/site-packages/GDAL-2.1.3-py3.6-macosx-10.12-x86_64.egg/')
# import gdal


class Tools(object):
    """Class that holds a few methods to create the static files needed for the visualization:

        * :meth:`.dem` - to produce a DEM cropped (and optionally decimated) to the region corresponding to a WRF domain
        * :meth:`.image` - to produce an image (currently, a `Blue Marble`_ image) cropped to the WRF domain which is used as a texture for the DEM

    :param wrf_file: name of the WRF output netCDF4_ file (``wrfout_...``) from which to take the (cloud) data (will be accessible as :attr:`.nc`)
    :type wrf_file: :obj:`str`

    .. attribute:: nc

        input netCDF4_ WRF file containing the variables (e.g. ``CLDFRA``) to be visualized.

    """
    def __init__(self, wrf_file):
        self.nc = Dataset(wrf_file)
        X, Y = self.nc['XLONG'][0, : ,:], self.nc['XLAT'][0, :, :]
        self.affine = self.create_affine(X, Y)
        self.xbnds = X.min(), X.max()
        self.ybnds = Y.max(), Y.min()

    def __del__(self):
        try:
            self.nc.close()
        except: pass

    def proj(self):
        """Return projection parameters for the ``wrfout_...`` file (:attr:.nc).
        """
        pyproj = import_module('pyproj')
        return pyproj.Proj(
            lon_0 = self.nc.CEN_LON,
            lat_0 = self.nc.CEN_LAT,
            lat_1 = self.nc.TRUELAT1,
            lat_2 = self.nc.TRUELAT2,
            proj = 'lcc'
        )

    def create_affine(self, x, y):
        """Return a function that maps the geographical coordinates (latitude and longitude) of input data (e.g. the ``wrfout_...`` file or the input DEM) to a Cartesian grid spanned by integer index vectors (e.g. as returned by a call to :meth:`numpy.mgrid`). The returnd function is to be called as ``func(x, y)`` where ``x`` and ``y`` correspond to :class:`arrays <numpy.ndarray>` of longitudes and latitudes.

        :param x: array of longitudes
        :type x: :class:`~numpy.ndarray`
        :param y: array of latitudes
        :type y: :class:`~numpy.ndarray`
        :returns: function that performs an affine mapping from geographic to integer-gridded coordinates
        :rtype: :obj:`callable`
        """
        proj = self.proj()
        x, y = proj(x, y)
        m, n = x.shape
        j, i = np.mgrid[:m, :n]

        C = np.linalg.lstsq(
            np.r_[[x.flatten()], [y.flatten()], np.ones((1, m * n))].T,
            np.r_[[i.flatten()], [j.flatten()], np.ones((1, m * n))].T)[0].T
        b = C[:2, 2]
        A = C[:2, :2]

        def to_grid(lon, lat):
            lon, lat = proj(lon, lat)
            return A.dot(np.vstack((lon.flatten(), lat.flatten()))) + np.r_[[b]].T

        return to_grid

    def dem(self, file_path, decimate=12):
        """Return a netCDF4_ file containing a DEM cropped to the region corresponding to the input netCDF4_ WRF file (argument ``wrf_file`` to :class:`Tools` and :attr:`.nc`), and optionally downsampled.

        :param file_path: path to the input DEM file (in geotiff format, e.g. from the SRTM_)
        :param decimate: decimation factor for the resulting DEM (e.g. if ``12``, the resulting resolution is 1/12th of the original)
        :returns: cropped and optionally decimated DEM, with elevation as variable ``z`` and the affine-mapped coordinates (see :meth:`.create_affine` in variables ``x`` and ``y``, respectively)
        :rtype: :class:`xarray.DataArray`

        """
        gdal = import_module('gdal')
        ds = gdal.Open(file_path)
        g = ds.GetGeoTransform()
        b = ds.GetRasterBand(1)
        z = b.ReadAsArray()

        dx, dy, x = g[1], g[5], z
        if decimate > 1:
            sig = import_module('scipy.signal')
            x = sig.decimate(sig.decimate(z, decimate, axis=0), decimate, axis=1)
            dx = (g[1] * z.shape[1]) / x.shape[1]
            dy = (g[5] * z.shape[0]) / x.shape[0]

        Z = xr.DataArray(x, coords=[
            ('lat', g[3] + dy * (np.arange(x.shape[0]) + .5)),
            ('lon', g[0] + dx * (np.arange(x.shape[1]) + .5))
        ]).sel(lon=slice(*self.xbnds), lat=slice(*self.ybnds))

        return self.project_xarray(Z)

    def project_xarray(self, Z):
        x, y = self.affine(*np.meshgrid(Z.lon, Z.lat))
        y = y + self.nc.dimensions['south_north'].size
        return xr.Dataset({
            'z': Z,
            'x': (('lat', 'lon'), x.reshape(Z.shape)),
            'y': (('lat', 'lon'), y.reshape(Z.shape))
        })

    def image(self, file_path):
        """Return an image cropped to the region corresponding to the input netCDF4_ WRF file (argument ``wrf_file`` to :class:`Tools` and :attr:`.nc`). This function currently assumes that the image is a `Blue Marble`_ region ``B2`` image - this is important because these images are not geotiffs and I map the pixel coordinates based on the boundary of the ``B2`` region.

        :param file_path: file path to the original image
        :returns: cropped image
        :rtype: :class:`PIL.Image.Image`

        .. NOTE::

            The image needs, potentially, to be rotated and flipped (90 deg CCW and flipped horizontally) before applying it as texture.
        """
        Image = import_module('PIL.Image')
        Image.MAX_IMAGE_PIXELS = None
        im = Image.open(im_path)
        # blue marble next gen is 1/4th of an arc minute (1/dx = 240)
        i = (np.vstack((self.xbnds, self.ybnds)).T * [1, -1] + [90, 0]) * 240
        return im.crop(np.hstack((np.floor(i[0, :]), np.ceil(i[1, :]))))


class Visualize(Tools):
    """Class that sets up an interactive Mayavi_ visualization - at this points of clouds, using volume rendering.

    Usage::

        v = Visualize(<wrfout_file>)
        v.anim()
        v.show()

    :param wrf_file: path to a WRF output netCDF4_ file whose data is to be visualized
    :param var_name: name of the variable to be visualized as a volume density (e.g. ``CLDFRA``)
    :param dem_file: path to a DEM in netCDF4_ format which represents the terrain, as produced by :meth:`.dem`
    :param image_file: path to an image file which is used as texture on the rendered DEM.

    """
    def __init__(self, wrf_file, var_name='CLDFRA', dem_file='dem_d03_12.nc', image_file='marble_d03_rot.jpg'):
        super().__init__(wrf_file)
        tvtk = import_module('tvtk.api')
        self.mlab = import_module('mayavi.mlab')

        ds = xr.open_dataset(dem_file)

        # here I 'fill' the DEM westward with zeros (over the ocean)
        dx = ds.lon.diff('lon').mean().item()
        i = np.arange(ds.lon.min(), self.xbnds[0], -dx)[:0:-1]
        Z = xr.DataArray(np.pad(ds.z, [(0, 0), (len(i), 0)], 'constant'),
                         coords = [('lat', ds.lat), ('lon', np.r_[i, ds.lon])])
        self.z = self.project_xarray(Z)

        # get the texture from image (extensions either 'jpg' or 'png')
        reader = {'.jpg': tvtk.tvtk.JPEGReader, '.png': tvtk.tvtk.PNGReader}[os.path.splitext(image_file)[1]]
        im = reader(file_name = image_file)
        self.tex = tvtk.tvtk.Texture(input_connection = im.output_port, interpolate = 1)

        self.fig = self.mlab.figure(size=(800, 800))

        # surf = self.mlab.surf(self.z.x.T, self.z.y.T, self.z.z.values.T/1000, colormap='gist_earth')
        surf = self.mlab.surf(self.z.x.T, self.z.y.T, self.z.z.values.T/1000, color=(1, 1, 1))
        surf.actor.enable_texture = True
        surf.actor.tcoord_generator_mode = 'plane'
        surf.actor.actor.texture = self.tex

        # self.mlab.view(-110, 80, 115, [50, 60, 1])
        self.mlab.view(-120, 65, 200, [40, 66, -6])

        ly, lx, self.nt = [self.nc.dimensions[n].size for n in ['south_north', 'west_east', 'Time']]
        self.xyz = [c.transpose(2, 1, 0) for c in np.mgrid[0:10:0.05, :ly, :lx]][::-1]
        self.vol = None

        @self.mlab.animate
        def anim():
            for i in range(self.nt):
                self.anim_func(i, var_name)
                yield

        self.anim = anim
        self.show = self.mlab.show

    def anim_func(self, i, var_name):
        wrf = import_module('wrf')
        iz = np.arange(0, 10, .05)
        if self.vol is not None:
            self.vol.remove()
        cld = wrf.vinterp(self.nc, wrf.getvar(self.nc, var_name, timeidx=int(i)), 'ght_msl', iz)
        xyzc = self.xyz + [cld.values.transpose(2, 1, 0)]
        self.vol = self.mlab.pipeline.volume(self.mlab.pipeline.scalar_field(*xyzc), color=(1, 1, 1), figure=self.fig)

    def write_movie(self, var_name='CLDFRA', file_name='mov.mp4', fps=6):
        """Write a movie file to disk, based on the interactive animation set up by this class.

        :param var_name: name of the variable to be visualized as a volume density (e.g. ``CLDFRA``)
        :param file_name: path to the output (movie) file to be produced
        :param fps: frames per second

        """
        def anim(i):
            self.anim_func(i * fps, var_name)
            return self.mlab.screenshot(antialiased=True)
        mp = import_module('moviepy.editor')
        vc = mp.VideoClip(anim, duration=self.nt / fps)
        writer = {'.mp4': 'write_videofile', '.gif': 'write_gif'}[os.path.splitext(file_name)[1]]
        getattr(vc, writer)(file_name, fps=fps)
