from netCDF4 import Dataset
import xarray as xr
import numpy as np
from mayavi import mlab
from tvtk.api import tvtk
from pyproj import Proj
from traits.api import HasTraits
import ctypes

# VTK bug, see:
# https://www.vtk.org/pipermail/vtk-developers/2017-November/035592.html
# osm = ctypes.CDLL('/usr/local/lib/libOSMesa32.so', ctypes.RTLD_GLOBAL)

mlab.options.offscreen = True

class Nubes(HasTraits):
    @staticmethod
    def create_affine(nc):
        proj = Proj(lon_0=nc.CEN_LON, lat_0=nc.CEN_LAT, lat_1=nc.TRUELAT1, lat_2=nc.TRUELAT2, proj='lcc')
        x, y = proj(nc['XLONG'][0, : ,:], nc['XLAT'][0, :, :])
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

    def run(self):
        ds = xr.open_dataset('dem_d03_12.nc')
        nc = Dataset('wrfout_d03_2016-05-11_12:00:00')

        # here I 'fill' the DEM westward with zeros (over the ocean)
        dx = ds.lon.diff('lon').mean().item()
        i = np.arange(ds.lon.min(), nc['XLONG'][0, :, :].min(), -dx)[:0:-1]
        Z = xr.DataArray(np.pad(ds.z, [(0, 0), (len(i), 0)], 'constant'),
                         coords = [('lat', ds.lat), ('lon', np.r_[i, ds.lon])])

        affine = self.create_affine(nc)
        x, y = affine(*np.meshgrid(Z.lon, Z.lat))
        y = y + nc.dimensions['south_north'].size
        x, y = [v.reshape(Z.shape) for v in [x, y]]

        im = tvtk.JPEGReader(file_name='marble_d03_rot.jpg')
        tex = tvtk.Texture(input_connection=im.output_port, interpolate=1)

        surf = mlab.surf(x.T, y.T, Z.values.T/1000, color=(1, 1, 1))
        surf.actor.enable_texture = True
        surf.actor.tcoord_generator_mode = 'plane'
        surf.actor.actor.texture = tex

        mlab.view(-120, 60, 200, [55, 55, -9])

        mlab.savefig('test.png', size=(1200, 800))

if __name__ == '__main__':
    mlab.test_plot3d()
