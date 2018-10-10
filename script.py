from netCDF4 import Dataset
import xarray as xr
import numpy as np
import wrf
from mayavi import mlab
from tvtk.api import tvtk
from pyproj import Proj
from glob import glob
from os import remove, environ
from moviepy.editor import ImageSequenceClip
from ctypes import CDLL, RTLD_GLOBAL

fps = 6                                     # frames per second of the produced movie
movie_file = '/assets/mov.mp4'              # path to the file in which to save the produced movie
dem_file = '/assets/dem_d03_12.nc'          # path to static DEM file (in netcdf format)
image_file = '/assets/marble_d03_rot.jpg'   # path to texture image

intermediate_images = '/assets/scene'       # base name for the generated images (will be deleted after
                                            # movie has been generated)

# vertical resolution of volume rendering
top = 15                                    # upper boundary in km
vspace = 0.05                               # vertical resolution in km

# VTK bug affecting OFFSCREEN RENDERING, see:
# https://www.vtk.org/pipermail/vtk-developers/2017-November/035592.html
# and
# http://vtk.1045678.n5.nabble.com/Offscreen-rendering-problems-on-headless-Ubuntu-td5746035.html
# uncomment if used on-screen
osm = CDLL('/usr/local/lib/libOSMesa32.so', RTLD_GLOBAL)
mlab.options.offscreen = True

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

ds = xr.open_dataset(dem_file)
nc = Dataset(environ['wrf_file'])

# here I 'fill' the DEM westward with zeros (over the ocean)
dx = ds.lon.diff('lon').mean().item()
i = np.arange(ds.lon.min(), nc['XLONG'][0, :, :].min(), -dx)[:0:-1]
Z = xr.DataArray(np.pad(ds.z, [(0, 0), (len(i), 0)], 'constant'),
                 coords = [('lat', ds.lat), ('lon', np.r_[i, ds.lon])])

affine = create_affine(nc)
x, y = affine(*np.meshgrid(Z.lon, Z.lat))
y = y + nc.dimensions['south_north'].size
x, y = [v.reshape(Z.shape) for v in [x, y]]

im = tvtk.JPEGReader(file_name=image_file)
tex = tvtk.Texture(input_connection=im.output_port, interpolate=1)

surf = mlab.surf(x.T, y.T, Z.values.T/1000, color=(1, 1, 1))
surf.actor.enable_texture = True
surf.actor.tcoord_generator_mode = 'plane'
surf.actor.actor.texture = tex

mlab.view(-120, 60, 200, [55, 55, -9])

iz = np.arange(0, top, vspace)
ly, lx, nt = [nc.dimensions[n].size for n in ['south_north', 'west_east', 'Time']]
tr = lambda x: x.transpose(2, 1, 0)
z, y, x = np.mgrid[slice(0, top, vspace), :ly, :lx]

# NOTE: rendering the first timestep twice (and writing it to file!!!) seems to get around the issue with the lacking lighting in the first image
cld = wrf.vinterp(nc, wrf.getvar(nc, 'CLDFRA', timeidx=0), 'ght_msl', iz)
vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(tr(x), tr(y), tr(z), tr(cld.values)), color=(1, 1, 1), vmin=0, vmax=.7)
mlab.savefig('{}_{}.png'.format(intermediate_images, i), size=(1200, 800))
for i in range(nt):
    vol.remove()
    cld = wrf.vinterp(nc, wrf.getvar(nc, 'CLDFRA', timeidx=i), 'ght_msl', iz)
    vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(tr(x), tr(y), tr(z), tr(cld.values)), color=(1, 1, 1), vmin=0, vmax=.7)
    mlab.savefig('scene_{:03d}.png'.format(i), size=(1200, 800))

g = sorted(glob('scene_*png'))
mov = ImageSequenceClip(g, fps=fps)
mov.write_videofile(movie_file, fps=fps)
for f in g:
    remove(f)
