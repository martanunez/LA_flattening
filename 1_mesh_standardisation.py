"""
    Copyright (c) - Marta Nunez Garcia
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General
    Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option)
    any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details. You should have received a copy of the GNU General Public License along with this
    program. If not, see <http://www.gnu.org/licenses/>.
"""
"""
    Standardise mesh: clip pulmonary veins (PVs), left atrial appendage (LAA), and mitral valve (MV)
    Launch GUI and ask the user to select 5 seeds close to the ending points of: LSPV, LIPV, RIPV, RSPV & LAA (IN THIS ORDER)

    Input: mesh and several parameters related to clipping
    Output: clipped mesh and auxiliary files containing info about seeds, clipping planes, etc.
    Example usage: python 1_mesh_standardisation.py --meshfile data/mesh.vtk --pv_dist 5 --laa_dist 5 --vis 1 
"""

from clip_aux_functions import *
import os
import sys
import xlsxwriter
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--meshfile', type=str, metavar='PATH', help='path to input mesh')
parser.add_argument('--pv_dist', type=int, default=3, help='PV clipping distance (mm)')
parser.add_argument('--laa_dist', type=int, default=3, help='LAA clipping distance (mm)')
parser.add_argument('--maxslope', type=int, default=5, help='Anything above this is ostium')
parser.add_argument('--clspacing', type=float, default=0.4, help='Resample the centerline with this spacing')
parser.add_argument('--skippointsfactor', type=float, default=0.1, help='Percentage of points to ignore at beginning of centerline')
parser.add_argument('--highslope', type=float, default=1.2, help='Above this slope we start counting')
parser.add_argument('--bumpcriterion', type=float, default=0.05, help='Ostium if slope higher than highslope and above bump criterion')
parser.add_argument('--pvends', type=int, default=1, help='Enforce the centerline to reach the end boundary of the surface.')
parser.add_argument('--vis', type=int, default=1, help='Set to 1 to visualise clipping results overlaid with original mesh')
parser.add_argument('--save', type=int, default=0, help='Set to 0 to remove intermediate results (centerlines, clippoints, etc.)')
args = parser.parse_args()


fileroot = os.path.dirname(args.meshfile)
filename = os.path.basename(args.meshfile)
filenameroot = os.path.splitext(filename)[0]

if os.path.isfile(args.meshfile):
    surface = readvtk(args.meshfile)
else:
    sys.exit('ERROR: input file does not exist. Please, specify a valid path.')

# compute seeds (save polydata as vtp)
seedsfile = os.path.join(fileroot, filenameroot + '_clipseeds.vtp')
outseedsfile = os.path.join(fileroot, filenameroot + '_clipseeds.csv')
if not os.path.exists(seedsfile):
    laa_seedon = 1
    print('Please, select exactly 5 seeds in this order: \n 1. RSPV \n 2. RIPV \n 3. LIPV \n 4. LSPV \n 5. LAA')
    select_seeds(surface, 'GTLabels', seedsfile, args.vis, laa_seedon)
else:
    print('Seeds already selected, using those ones. To compute new seeds, delete file with suffix clip_seeds.vtp')
seeds_to_csv(seedsfile, 'GTLabels', [77, 76, 78, 79, 36], outseedsfile)

# Compute centerlines
outfile = os.path.join(fileroot, filenameroot + '_')
pv_LAA_centerlines(args.meshfile, outseedsfile, outfile, args.pvends)

# label PVs automatically
outfile = os.path.join(fileroot, filenameroot + '_')
clip_veins_sections_and_LAA(args.meshfile, outfile, args.clspacing, args.maxslope, args.skippointsfactor, args.highslope, args.bumpcriterion)

# clip PV end points
sufixfile = os.path.join(fileroot, filenameroot + '_')
inputfile = os.path.join(fileroot, filenameroot + '_autolabels.vtp')
inputsurface = readvtp(inputfile)

# use special distance (dist_LAA) for the LAA (label=37). Save the info about the clipping planes
stdmesh, clip_planes = clip_vein_endpoint_and_LAA_save_planes(inputsurface, sufixfile, args.pv_dist, 37, args.laa_dist)

plane_file = os.path.join(fileroot, filenameroot + '_clip_planes.xlsx')
workbook = xlsxwriter.Workbook(plane_file)
worksheet = workbook.add_worksheet()
onoff = 1
for i in range(5):
    point = clip_planes[2*i, 0:3]
    normal = clip_planes[2*i+1, 0:3]
    for j in range(3):
        worksheet.write(2*i, j, point[j])
        worksheet.write(2*i+1, j, normal[j])
workbook.close()

# close small holes that can appear after clipping the veins & copy original scalar info in the filled holes
# in same cases the normals flip. Have in mind...
max_hole_size = 4  # empirical, be careful to do not close pv ostiums, check visually.
stdmesh_closed = fillholes(stdmesh, max_hole_size)
print('\n')
transfer_all_scalar_arrays(surface, stdmesh_closed)

if args.vis > 0:
    visualise_default(stdmesh_closed, surface, 'STD mesh', 'autolabels', 36, 79)
writevtk(stdmesh_closed, os.path.join(fileroot, filenameroot + '_clipped.vtk'))

# Apply crinkle clip to the veins. It does not cut cells.
array_labels = np.zeros(surface.GetNumberOfPoints())
locator = vtk.vtkPointLocator()
locator.SetDataSet(stdmesh_closed)
locator.BuildLocator()
for p in range(surface.GetNumberOfPoints()):
    point = surface.GetPoint(p)
    closestpoint_id = locator.FindClosestPoint(point)
    dist = euclideandistance(point, stdmesh_closed.GetPoint(closestpoint_id))
    if dist > 0.05:   # empirical distance
        array_labels[p] = 1
newarray = numpy_to_vtk(array_labels)
newarray.SetName("pv")
surface.GetPointData().AddArray(newarray)
m_ccliped = pointthreshold(surface, "pv", 0, 0)
transfer_array(stdmesh_closed, m_ccliped, 'autolabels', 'autolabels')
writevtk(cleanpolydata(m_ccliped), os.path.join(fileroot, filenameroot + '_crinkle_clipped.vtk'))

# MV clip, auto
w = [0.95, 0.05, 0.0]
o_file = os.path.join(fileroot, filenameroot + '_clipped_mitral')
surfaceclipped = find_mitral_cylinder_pvs(stdmesh, 'autolabels', o_file, 0.35, w, 0)
if args.vis > 0:
    visualise_color(surfaceclipped, surface, 'Mitral Valve clip')
writevtk(cleanpolydata(surfaceclipped), o_file + '.vtk')

if args.save == 0:  # remove intermediate results (centerlines, clippoints, etc.)
    if sys.platform == "linux" or sys.platform == "linux2":
        os.system("rm " + os.path.join(fileroot, filenameroot + '*clbranch*.vtp'))
        os.system("rm " + os.path.join(fileroot, filenameroot + '*clvein*.vtp'))
        os.system("rm " + os.path.join(fileroot, filenameroot + '*clraw*.vtp'))
        os.system("rm " + os.path.join(fileroot, filenameroot + '*clsection*.vtp'))
        os.system("rm " + os.path.join(fileroot, filenameroot + '*clippointid*.csv'))
    if sys.platform == "win32":
        os.system("del " + os.path.join(fileroot, filenameroot + '*clbranch*.vtp'))
        os.system("del " + os.path.join(fileroot, filenameroot + '*clvein*.vtp'))
        os.system("del " + os.path.join(fileroot, filenameroot + '*clraw*.vtp'))
        os.system("del " + os.path.join(fileroot, filenameroot + '*clsection*.vtp'))
        os.system("del " + os.path.join(fileroot, filenameroot + '*clippointid*.csv'))
