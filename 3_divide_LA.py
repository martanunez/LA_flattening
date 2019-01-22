# Parcellate mesh creating appropriate paths to divide the LA in the 5 pieces considered in our regional flattening.
# Launch GUI and ask the user to select 8 seeds in this order:
# 4 in RSPV, RIPV, LIPV, LSPV. Approx in the center of the filled holes.
# 4 in the MV in this order: going down from RSPV, RIPV, LIPV, and LAA.

# Input: mesh with clipped & filled holes corresponding to PVs and LAA; and clipped MV.
# Output: dividing paths saved as separated polydatas.
# Usage: python 3_divide_LA.py data/mesh_clipped_c.vtk

from clip_aux_functions import *
import sys
import os
import numpy as np

meshfile = sys.argv[1]
fileroot = os.path.dirname(meshfile)
filename = os.path.basename(meshfile)
filenameroot = os.path.splitext(filename)[0]

if not os.path.exists(meshfile):
    print('ERROR: Input file not found')
    exit()

outputfile = os.path.join(fileroot, filenameroot + '_seeds.vtk')   # selected by user
outputfile2 = os.path.join(fileroot, filenameroot + '_seeds_for_flat.vtk')  # modified seeds, used for flattening

surface = readvtk(meshfile)
nseeds = 9
labels = [0, 1, 2, 3, 4, 5, 6, 7, 8]

if not os.path.exists(outputfile):
    print('Select exactly 9 seeds in this order: \n 1. RSPV\n 2. RIPV\n 3. LIPV \n 4. LSPV \n 5. LAA and, \n 6. 4 seeds close to the MV contour (starting in the end of the line connecting RSPV with MV)')
    seeds = seed_interactor(surface)
    if not seeds.GetNumberOfIds() == nseeds:
        print('You should select exactly', nseeds, ' seeds. Try again!')
        seeds = seed_interactor(surface)

    # create the mesh with seeds
    newpoints = vtk.vtkPoints()
    newvertices = vtk.vtkCellArray()
    labels_array = vtk.vtkDoubleArray()
    labels_array.SetName('seed_label')

    for s in range(seeds.GetNumberOfIds()):
        label = labels[s]
        point = surface.GetPoint(seeds.GetId(s))
        pid = newpoints.InsertNextPoint(point)
        labels_array.InsertNextValue(label)
        # Create the topology of the point (a vertex)
        newvertices.InsertNextCell(1)
        newvertices.InsertCellPoint(pid)

    pointspd = vtk.vtkPolyData()
    pointspd.SetPoints(newpoints)
    pointspd.SetVerts(newvertices)
    pointspd.GetPointData().AddArray(labels_array)
    writevtk(pointspd, outputfile)
else:
    print('Using pre-computed seeds. To compute new ones delete file with suffix _seeds.vtk and run again')

# Find corresponding seeds in the complete mesh
seeds_poly = readvtk(outputfile)
locator = vtk.vtkPointLocator()
locator.SetDataSet(surface)
locator.BuildLocator()

id_1 = locator.FindClosestPoint(seeds_poly.GetPoint(0))
id_2 = locator.FindClosestPoint(seeds_poly.GetPoint(1))
id_3 = locator.FindClosestPoint(seeds_poly.GetPoint(2))
id_4 = locator.FindClosestPoint(seeds_poly.GetPoint(3))
id_laa = locator.FindClosestPoint(seeds_poly.GetPoint(4))
id_mv1 = locator.FindClosestPoint(seeds_poly.GetPoint(5))
id_mv2 = locator.FindClosestPoint(seeds_poly.GetPoint(6))
id_mv3 = locator.FindClosestPoint(seeds_poly.GetPoint(7))
id_mv4 = locator.FindClosestPoint(seeds_poly.GetPoint(8))

# detect MV contour
edges = extractboundaryedge(surface)
cont = extractlargestregion(edges)
edge_cont_ids = get_ordered_cont_ids_based_on_distance(cont).astype(int)

# Find closest point IN the contour corresponding to id_top, id_bottom, id_left, id_right
locator_cont = vtk.vtkPointLocator()
locator_cont.SetDataSet(cont)
locator_cont.BuildLocator()

id_mv1_cont = locator_cont.FindClosestPoint(surface.GetPoint(id_mv1))
id_mv2_cont = locator_cont.FindClosestPoint(surface.GetPoint(id_mv2))
id_mv3_cont = locator_cont.FindClosestPoint(surface.GetPoint(id_mv3))
id_mv4_cont = locator_cont.FindClosestPoint(surface.GetPoint(id_mv4))

# find corresponding ordered points in the COMPLETE mesh (before I was using only the contour)
locator_complete = vtk.vtkPointLocator()
locator_complete.SetDataSet(surface)
locator_complete.BuildLocator()
cont_mv_ids = np.zeros(edge_cont_ids.shape[0]) - 1
for i in range(cont_mv_ids.shape[0]):
    p = cont.GetPoint(edge_cont_ids[i])
    cont_mv_ids[i] = locator_complete.FindClosestPoint(p)

id_5 = locator_complete.FindClosestPoint(cont.GetPoint(id_mv1_cont))
id_6 = locator_complete.FindClosestPoint(cont.GetPoint(id_mv2_cont))
id_7 = locator_complete.FindClosestPoint(cont.GetPoint(id_mv3_cont))
id_8 = locator_complete.FindClosestPoint(cont.GetPoint(id_mv4_cont))

# save final seeds to be used in the flattening
poly_points = vtk.vtkPolyData()
points = vtk.vtkPoints()
points.SetNumberOfPoints(9)
points.SetPoint(0, surface.GetPoint(id_1))
points.SetPoint(1, surface.GetPoint(id_2))
points.SetPoint(2, surface.GetPoint(id_3))
points.SetPoint(3, surface.GetPoint(id_4))
points.SetPoint(4, surface.GetPoint(id_laa))
points.SetPoint(5, surface.GetPoint(id_5))
points.SetPoint(6, surface.GetPoint(id_6))
points.SetPoint(7, surface.GetPoint(id_7))
points.SetPoint(8, surface.GetPoint(id_8))
poly_points.SetPoints(points)
# poly_points.Update()
writevtk(poly_points, outputfile2)

# create and write inter-seeds paths
path1 = find_create_path(surface, id_1, id_2)
path2 = find_create_path(surface, id_2, id_3)
path3 = find_create_path(surface, id_3, id_4)
path4 = find_create_path(surface, id_4, id_1)
path5 = find_create_path(surface, id_1, id_5)
path6 = find_create_path(surface, id_2, id_6)
path7 = find_create_path(surface, id_3, id_7)
path_laa1 = find_create_path(surface, id_4, id_laa)
path_laa2 = find_create_path(surface, id_laa, id_8)
path_laa3 = find_create_path(surface, id_laa, id_1)

writevtk(path1, os.path.join(fileroot, filenameroot + 'path1.vtk'))
writevtk(path2, os.path.join(fileroot, filenameroot + 'path2.vtk'))
writevtk(path3, os.path.join(fileroot, filenameroot + 'path3.vtk'))
writevtk(path4, os.path.join(fileroot, filenameroot + 'path4.vtk'))
writevtk(path5, os.path.join(fileroot, filenameroot + 'path5.vtk'))
writevtk(path6, os.path.join(fileroot, filenameroot + 'path6.vtk'))
writevtk(path7, os.path.join(fileroot, filenameroot + 'path7.vtk'))
writevtk(path_laa1, os.path.join(fileroot, filenameroot + 'path_laa1.vtk'))
writevtk(path_laa2, os.path.join(fileroot, filenameroot + 'path_laa2.vtk'))
writevtk(path_laa3, os.path.join(fileroot, filenameroot + 'path_laa3.vtk'))
