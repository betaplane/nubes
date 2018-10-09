from netCDF4 import Dataset
import xarray as xr
import numpy as np
import wrf
from mayavi import mlab
from tvtk.api import tvtk
from pyproj import Proj
from glob import glob
from os import remove
from moviepy.editor import ImageSequenceClip
from ctypes import CDLL, RTLD_GLOBAL

fps = 6
outfile = 'mov.mp4'

# see https://github.com/enthought/mayavi/issues/7
# from tvtk.pyface.ui.qt4 import scene
# from tvtk.pyface.ui.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

# class tvtkBug(scene._VTKRenderWindowInteractor):
#     def paintEvent(self, e):
#         QVTKRenderWindowInteractor.paintEvent(self, e)

# scene._VTKRenderWindowInteractor = tvtkBug

# VTK bug, see:
# https://www.vtk.org/pipermail/vtk-developers/2017-November/035592.html
# and
# http://vtk.1045678.n5.nabble.com/Offscreen-rendering-problems-on-headless-Ubuntu-td5746035.html
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

ds = xr.open_dataset('dem_d03_12.nc')
nc = Dataset('wrfout_d03_2016-05-11_12:00:00')

# here I 'fill' the DEM westward with zeros (over the ocean)
dx = ds.lon.diff('lon').mean().item()
i = np.arange(ds.lon.min(), nc['XLONG'][0, :, :].min(), -dx)[:0:-1]
Z = xr.DataArray(np.pad(ds.z, [(0, 0), (len(i), 0)], 'constant'),
                 coords = [('lat', ds.lat), ('lon', np.r_[i, ds.lon])])

affine = create_affine(nc)
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

top = 15
vspace = 0.05
iz = np.arange(0, top, vspace)
ly, lx, nt = [nc.dimensions[n].size for n in ['south_north', 'west_east', 'Time']]
tr = lambda x: x.transpose(2, 1, 0)
z, y, x = np.mgrid[slice(0, top, vspace), :ly, :lx]

vol = None
for i in range(nt):
    try:
        vol.remove()
    except: pass
    cld = wrf.vinterp(nc, wrf.getvar(nc, 'CLDFRA', timeidx=i), 'ght_msl', iz)
    vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(tr(x), tr(y), tr(z), tr(cld.values)), color=(1, 1, 1), vmin=0, vmax=.7)
    mlab.savefig('scene_{:03d}.png'.format(i), size=(1200, 800))

g = sorted(glob('scene_*png'))
mov = ImageSequenceClip(g, fps=fps)
mov.write_videofile(outfile, fps=fps)
for f in g:
    remove(f)

# make_frame(0)
# from tvtk.pyface import light_manager
# vol.scene.light_manager =
# mov = VideoClip(make_frame, duration=nt / fps)
# mov.write_videofile('test.mp4', fps=fps)
# mov.write_gif('test.gif', fps=6)
# mlab.savefig('test.png', size=(800, 800))
