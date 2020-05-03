import vtk
import math
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
import os
import sys
from scipy import sparse
import scipy.sparse.linalg as linalg_sp
from scipy.sparse import vstack, hstack, coo_matrix, csc_matrix
import csv

###     Input/Output    ###
def readvtk(filename):
    """Read VTK file"""
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def readvtp(filename):
    """Read VTP file"""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def writevtk(surface, filename, type='ascii'):
    """Write binary or ascii VTK file"""
    writer = vtk.vtkPolyDataWriter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        writer.SetInputData(surface)
    else:
        writer.SetInput(surface)
    writer.SetFileName(filename)
    if type == 'ascii':
        writer.SetFileTypeToASCII()
    elif type == 'binary':
        writer.SetFileTypeToBinary()
    writer.Write()

def writevtp(surface, filename):
    """Write VTP file"""
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        writer.SetInputData(surface)
    else:
        writer.SetInput(surface)
    writer.SetFileName(filename)
#    writer.SetDataModeToBinary()
    writer.Write()

###     Math    ###
def euclideandistance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2 + (point1[2] - point2[2])**2)

def normvector(v):
    return math.sqrt(dot(v, v))

def angle(v1, v2):
    return math.acos(dot(v1, v2) / (normvector(v1) * normvector(v2)))

def acumvectors(point1, point2):
    return [point1[0] + point2[0], point1[1] + point2[1], point1[2] + point2[2]]

def subtractvectors(point1, point2):
    return [point1[0] - point2[0], point1[1] - point2[1], point1[2] - point2[2]]

def dividevector(point, n):
    nr = float(n)
    return [point[0]/nr, point[1]/nr, point[2]/nr]

def multiplyvector(point, n):
    nr = float(n)
    return [nr*point[0], nr*point[1], nr*point[2]]

def sumvectors(vect1, scalar, vect2):
    return [vect1[0] + scalar*vect2[0], vect1[1] + scalar*vect2[1], vect1[2] + scalar*vect2[2]]

def cross(v1, v2):
    return [v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]]

def dot(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def normalizevector(v):
    norm = normvector(v)
    return [v[0] / norm, v[1] / norm, v[2] / norm]

###     Mesh processing     ###

def cleanpolydata(polydata):
    cleaner = vtk.vtkCleanPolyData()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        cleaner.SetInputData(polydata)
    else:
        cleaner.SetInput(polydata)
    cleaner.Update()
    return cleaner.GetOutput()

def fillholes(polydata, size):
    """Fill mesh holes smaller than 'size' """
    filler = vtk.vtkFillHolesFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        filler.SetInputData(polydata)
    else:
        filler.SetInput(polydata)
    filler.SetHoleSize(size)
    filler.Update()
    return filler.GetOutput()

def pointthreshold(polydata, arrayname, start=0, end=1, alloff=0):
    """ Clip polydata according to given thresholds in scalar array"""
    threshold = vtk.vtkThreshold()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        threshold.SetInputData(polydata)
    else:
        threshold.SetInput(polydata)
    threshold.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, arrayname)
    threshold.ThresholdBetween(start, end)
    if (alloff):
        threshold.AllScalarsOff()
    threshold.Update()
    surfer = vtk.vtkDataSetSurfaceFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        surfer.SetInputData(threshold.GetOutput())
    else:
        surfer.SetInput(threshold.GetOutput())
    surfer.Update()
    return surfer.GetOutput()

def cellthreshold(polydata, arrayname, start=0, end=1):
    threshold = vtk.vtkThreshold()
    threshold.SetInputData(polydata)
    threshold.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,arrayname)
    threshold.ThresholdBetween(start,end)
    threshold.Update()

    surfer = vtk.vtkDataSetSurfaceFilter()
    surfer.SetInputConnection(threshold.GetOutputPort())
    surfer.Update()
    return surfer.GetOutput()

def roundpointarray(polydata, name):
    """Round values in point array"""
    # get original array
    array = polydata.GetPointData().GetArray(name)
    # round labels
    for i in range(polydata.GetNumberOfPoints()):
        value = array.GetValue(i)
        array.SetValue(i, round(value))
    return polydata

def planeclip(surface, point, normal, insideout=1):
    """Clip surface using plane given by point and normal"""
    clipplane = vtk.vtkPlane()
    clipplane.SetOrigin(point)
    clipplane.SetNormal(normal)
    clipper = vtk.vtkClipPolyData()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        clipper.SetInputData(surface)
    else:
        clipper.SetInput(surface)
    clipper.SetClipFunction(clipplane)

    if insideout == 1:
        # print 'insideout ON'
        clipper.InsideOutOn()
    else:
        # print 'insideout OFF'
        clipper.InsideOutOff()
    clipper.Update()
    return clipper.GetOutput()

def cutdataset(dataset, point, normal):
    """Similar to planeclip but using vtkCutter instead of vtkClipPolyData"""
    cutplane = vtk.vtkPlane()
    cutplane.SetOrigin(point)
    cutplane.SetNormal(normal)
    cutter = vtk.vtkCutter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        cutter.SetInputData(dataset)
    else:
        cutter.SetInput(dataset)
    cutter.SetCutFunction(cutplane)
    cutter.Update()
    return cutter.GetOutput()

def pointset_centreofmass(polydata):
    centre = [0, 0, 0]
    for i in range(polydata.GetNumberOfPoints()):
        point = [polydata.GetPoints().GetPoint(i)[0],
          polydata.GetPoints().GetPoint(i)[1],
          polydata.GetPoints().GetPoint(i)[2]]
        centre = acumvectors(centre,point)
    return dividevector(centre, polydata.GetNumberOfPoints())

def seeds_to_csv(seedsfile, arrayname, labels, outfile):
    """Read seeds from VTP file, write coordinates in csv"""
    # f = open(outfile, 'wb')
    f = open(outfile, 'w')
    allseeds = readvtp(seedsfile)
    for l in labels:
        currentseeds = pointthreshold(allseeds, arrayname, l, l, 0)

        # currentpoint = currentseeds.GetPoint(0)   # only 1
        currentpoint = pointset_centreofmass(currentseeds)

        line = str(currentpoint[0]) + ',' + str(currentpoint[1]) + ',' + str(currentpoint[2]) + '\n'
        # line = currentpoint[0] + b',' + currentpoint[1] + b',' + currentpoint[2] + b'\n'
        f.write(line)
    f.close()

def point2vertexglyph(point):
    """Create glyph from points to visualise them"""
    points = vtk.vtkPoints()
    points.InsertNextPoint(point[0], point[1], point[2])
    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
    glyph = vtk.vtkVertexGlyphFilter()
    glyph.SetInputConnection(poly.GetProducerPort())
    glyph.Update()
    return glyph.GetOutput()

def generateglyph(polyIn, scalefactor=2):
    vertexGlyphFilter = vtk.vtkGlyph3D()
    sphereSource = vtk.vtkSphereSource()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        vertexGlyphFilter.SetSourceData(sphereSource.GetOutput())
        vertexGlyphFilter.SetInputData(polyIn)
    else:
        vertexGlyphFilter.SetSource(sphereSource.GetOutput())
        vertexGlyphFilter.SetInput(polyIn)
    vertexGlyphFilter.SetColorModeToColorByScalar()
    vertexGlyphFilter.SetSourceConnection(sphereSource.GetOutputPort())
    vertexGlyphFilter.ScalingOn()
    vertexGlyphFilter.SetScaleFactor(scalefactor)
    vertexGlyphFilter.Update()
    return vertexGlyphFilter.GetOutput()

def linesource(p1, p2):
    """Create vtkLine from coordinates of 2 points"""
    source = vtk.vtkLineSource()
    source.SetPoint1(p1[0], p1[1], p1[2])
    source.SetPoint2(p2[0], p2[1], p2[2])
    return source.GetOutput()

def append(polydata1, polydata2):
    """Define new polydata appending polydata1 and polydata2"""
    appender = vtk.vtkAppendPolyData()
    appender.AddInput(polydata1)
    appender.AddInput(polydata2)
    appender.Update()
    return appender.GetOutput()

def extractcells(polydata, idlist):
    """Extract cells from polydata whose cellid is in idlist."""
    cellids = vtk.vtkIdList()  # specify cellids
    cellids.Initialize()
    for i in idlist:
        cellids.InsertNextId(i)

    extract = vtk.vtkExtractCells()  # extract cells with specified cellids
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        extract.SetInputData(polydata)
    else:
        extract.SetInput(polydata)
    extract.AddCellList(cellids)
    extraction = extract.GetOutput()

    geometry = vtk.vtkGeometryFilter()  # unstructured grid to polydata
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        geometry.SetInputData(extraction)
    else:
        geometry.SetInput(extraction)
    geometry.Update()
    return geometry.GetOutput()

def extractboundaryedge(polydata):
    edge = vtk.vtkFeatureEdges()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        edge.SetInputData(polydata)
    else:
        edge.SetInput(polydata)
    edge.FeatureEdgesOff()
    edge.NonManifoldEdgesOff()
    edge.Update()
    return edge.GetOutput()

def extractlargestregion(polydata):
    """Keep only biggest region"""
    surfer = vtk.vtkDataSetSurfaceFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        surfer.SetInputData(polydata)
    else:
        surfer.SetInput(polydata)
    surfer.Update()

    cleaner = vtk.vtkCleanPolyData()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        cleaner.SetInputData(surfer.GetOutput())
    else:
        cleaner.SetInput(surfer.GetOutput())
    cleaner.Update()

    connect = vtk.vtkPolyDataConnectivityFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        connect.SetInputData(cleaner.GetOutput())
    else:
        connect.SetInput(cleaner.GetOutput())
    connect.SetExtractionModeToLargestRegion()
    connect.Update()

    cleaner = vtk.vtkCleanPolyData()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        cleaner.SetInputData(connect.GetOutput())
    else:
        cleaner.SetInput(connect.GetOutput())
    cleaner.Update()
    return cleaner.GetOutput()

def countregions(polydata):
    """Count number of connected components/regions"""
    # preventive measures: clean before connectivity filter to avoid artificial regionIds
    surfer = vtk.vtkDataSetSurfaceFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        surfer.SetInputData(polydata)
    else:
        surfer.SetInput(polydata)
    surfer.Update()

    cleaner = vtk.vtkCleanPolyData()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        cleaner.SetInputData(surfer.GetOutput())
    else:
        cleaner.SetInput(surfer.GetOutput())
    cleaner.Update()

    connect = vtk.vtkPolyDataConnectivityFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        connect.SetInputData(cleaner.GetOutput())
    else:
        connect.SetInput(cleaner.GetOutput())
    connect.Update()
    return connect.GetNumberOfExtractedRegions()

def extractclosestpointregion(polydata, point=[0, 0, 0]):
    # NOTE: preventive measures: clean before connectivity filter
    # to avoid artificial regionIds
    # It slices the surface down the middle
    surfer = vtk.vtkDataSetSurfaceFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        surfer.SetInputData(polydata)
    else:
        surfer.SetInput(polydata)
    surfer.Update()

    cleaner = vtk.vtkCleanPolyData()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        cleaner.SetInputData(surfer.GetOutput())
    else:
        cleaner.SetInput(surfer.GetOutput())
    cleaner.Update()

    connect = vtk.vtkPolyDataConnectivityFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        connect.SetInputData(cleaner.GetOutput())
    else:
        connect.SetInput(cleaner.GetOutput())
    connect.SetExtractionModeToClosestPointRegion()
    connect.SetClosestPoint(point)
    connect.Update()
    return connect.GetOutput()

def extractconnectedregion(polydata, regionid):
    """Extract connected region with label = regionid """
    surfer = vtk.vtkDataSetSurfaceFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        surfer.SetInputData(polydata)
    else:
        surfer.SetInput(polydata)
    surfer.Update()

    cleaner = vtk.vtkCleanPolyData()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        cleaner.SetInputData(surfer.GetOutput())
    else:
        cleaner.SetInput(surfer.GetOutput())
    cleaner.Update()

    connect = vtk.vtkPolyDataConnectivityFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        connect.SetInputData(cleaner.GetOutput())
    else:
        connect.SetInput(cleaner.GetOutput())

    connect.SetExtractionModeToAllRegions()
    connect.ColorRegionsOn()
    connect.Update()
    surface = pointthreshold(connect.GetOutput(), 'RegionId', float(regionid), float(regionid))
    return surface

def get_connected_edges(polydata):
    """Extract all connected regions"""
    connect = vtk.vtkPolyDataConnectivityFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        connect.SetInputData(polydata)
    else:
        connect.SetInput(polydata)
    connect.SetExtractionModeToAllRegions()
    connect.ColorRegionsOn()
    connect.Update()
    return connect

def find_create_path(mesh, p1, p2):
    """Get shortest path (using Dijkstra algorithm) between p1 and p2 on the mesh. Returns a polydata"""
    dijkstra = vtk.vtkDijkstraGraphGeodesicPath()
    if vtk.vtkVersion().GetVTKMajorVersion() > 5:
        dijkstra.SetInputData(mesh)
    else:
        dijkstra.SetInput(mesh)
    dijkstra.SetStartVertex(p1)
    dijkstra.SetEndVertex(p2)
    dijkstra.Update()
    return dijkstra.GetOutput()

def compute_geodesic_distance(mesh, id_p1, id_p2):
    """Compute geodesic distance from point id_p1 to id_p2 on surface 'mesh'
    It first computes the path across the edges and then the corresponding distance adding up point to point distances)"""
    path = find_create_path(mesh, id_p1, id_p2)
    total_dist = 0
    n = path.GetNumberOfPoints()
    for i in range(n-1):   # Ids are ordered in the new polydata, from 0 to npoints_in_path
        p0 = path.GetPoint(i)
        p1 = path.GetPoint(i+1)
        dist = math.sqrt(math.pow(p0[0]-p1[0], 2) + math.pow(p0[1]-p1[1], 2) + math.pow(p0[2]-p1[2], 2) )
        total_dist = total_dist + dist
    return total_dist, path

def transfer_array(ref, target, arrayname, targetarrayname):
    """Transfer scalar array using closest point approximation"""
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(ref)
    locator.BuildLocator()

    refarray = ref.GetPointData().GetArray(arrayname)  # get array from reference

    numberofpoints = target.GetNumberOfPoints()
    newarray = vtk.vtkDoubleArray()
    newarray.SetName(targetarrayname)
    newarray.SetNumberOfTuples(numberofpoints)
    target.GetPointData().AddArray(newarray)

    # go through each point of target surface, determine closest point on surface, copy value
    for i in range(target.GetNumberOfPoints()):
        point = target.GetPoint(i)
        closestpoint_id = locator.FindClosestPoint(point)
        value = refarray.GetValue(closestpoint_id)
        newarray.SetValue(i, value)
    return target

def transfer_all_scalar_arrays(m1, m2):
    """ Transfer all scalar arrays from m1 to m2"""
    for i in range(m1.GetPointData().GetNumberOfArrays()):
        print('Transferring scalar array: {}'.format(m1.GetPointData().GetArray(i).GetName()))
        transfer_array(m1, m2, m1.GetPointData().GetArray(i).GetName(), m1.GetPointData().GetArray(i).GetName())

def transfer_all_scalar_arrays_by_point_id(m1, m2):
    """ Transfer all scalar arrays from m1 to m2 by point id"""
    for i in range(m1.GetPointData().GetNumberOfArrays()):
        print('Transferring scalar array: {}'.format(m1.GetPointData().GetArray(i).GetName()))
        m2.GetPointData().AddArray(m1.GetPointData().GetArray(i))

def get_ordered_cont_ids_based_on_distance(mesh):
    """ Given a contour, get the ordered list of Ids (not ordered by default).
    Open the mesh duplicating the point with id = 0. Compute distance transform of point 0
    and get a ordered list of points (starting in 0) """
    m = vtk.vtkMath()
    m.RandomSeed(0)
    # copy the original mesh point by point
    points = vtk.vtkPoints()
    polys = vtk.vtkCellArray()
    cover = vtk.vtkPolyData()
    nver = mesh.GetNumberOfPoints()
    points.SetNumberOfPoints(nver+1)

    new_pid = nver  # id of the duplicated point
    added = False

    for j in range(mesh.GetNumberOfCells()):
        # get the 2 point ids
        ptids = mesh.GetCell(j).GetPointIds()
        cell = mesh.GetCell(j)
        if (ptids.GetNumberOfIds() != 2):
            # print "Non contour mesh (lines)"
            break

        # read the 2 involved points
        pid0 = ptids.GetId(0)
        pid1 = ptids.GetId(1)
        p0 = mesh.GetPoint(ptids.GetId(0))   # returns coordinates
        p1 = mesh.GetPoint(ptids.GetId(1))

        if pid0 == 0:
            if added == False:
                # Duplicate point 0. Add gaussian noise to the original point
                new_p = [p0[0] + m.Gaussian(0.0, 0.0005), p0[1] + m.Gaussian(0.0, 0.0005), p0[2] + m.Gaussian(0.0, 0.0005)]
                points.SetPoint(new_pid, new_p)
                points.SetPoint(pid1, p1)
                polys.InsertNextCell(2)
                polys.InsertCellPoint(pid1)
                polys.InsertCellPoint(new_pid)
                added = True
            else:  # act normal
                points.SetPoint(ptids.GetId(0), p0)
                points.SetPoint(ptids.GetId(1), p1)
                polys.InsertNextCell(2)
                polys.InsertCellPoint(cell.GetPointId(0))
                polys.InsertCellPoint(cell.GetPointId(1))
        elif pid1 == 0:
            if added == False:
                new_p = [p1[0] + m.Gaussian(0.0, 0.0005), p1[1] + m.Gaussian(0.0, 0.0005), p1[2] + m.Gaussian(0.0, 0.0005)]
                points.SetPoint(new_pid, new_p)
                points.SetPoint(pid0, p0)
                polys.InsertNextCell(2)
                polys.InsertCellPoint(pid0)
                polys.InsertCellPoint(new_pid)
                added = True
            else:  # act normal
                points.SetPoint(ptids.GetId(0), p0)
                points.SetPoint(ptids.GetId(1), p1)
                polys.InsertNextCell(2)
                polys.InsertCellPoint(cell.GetPointId(0))
                polys.InsertCellPoint(cell.GetPointId(1))

        else:
            points.SetPoint(ptids.GetId(0), p0)
            points.SetPoint(ptids.GetId(1), p1)
            polys.InsertNextCell(2)
            polys.InsertCellPoint(cell.GetPointId(0))
            polys.InsertCellPoint(cell.GetPointId(1))

    if added == False:
        print('Warning: I have not added any point, list of indexes may not be correct.')
    cover.SetPoints(points)
    cover.SetPolys(polys)
    if not vtk.vtkVersion.GetVTKMajorVersion() > 5:
        cover.Update()
    # compute distance from point with id 0 to all the rest
    npoints = cover.GetNumberOfPoints()
    dists = np.zeros(npoints)
    for i in range(npoints):
        [dists[i], poly] = compute_geodesic_distance(cover, int(0), i)
    list_ = np.argsort(dists).astype(int)
    return list_[0:len(list_)-1]    # skip last one, duplicated

def define_pv_segments_proportions(t_v5, t_v6, t_v7, alpha):
    """define number of points of each pv hole segment to ensure appropriate distribution"""
    props = np.zeros([4, 3])
    props[0, 0] = np.divide(1.0, 4.0)  # proportion of the total number of points of the pv contour according to the proportion of circle
    props[0, 1] = np.divide(1.0, 4.0) + t_v5 * np.divide(1.0, 2.0*np.pi)
    props[0, 2] = 1.0 - props[0, 0] - props[0, 1]
    # print('Proportions sum up:', props[0, 0]+props[0, 1]+props[0, 2])
    props[1, 0] = np.divide(t_v6, 2.0*np.pi) - np.divide(1.0, 2.0)
    props[1, 2] = np.divide(1.0, 4.0)  # s3
    props[1, 1] = 1.0 - props[1, 0] - props[1, 2]   # s2
    # print('Proportions sum up:', props[1, 0]+props[1, 1]+props[1, 2])
    props[2, 0] = np.divide(1.0, 4.0)
    props[2, 1] = np.divide(t_v7, 2.0*np.pi) - props[2, 0]
    props[2, 2] = 1.0 - props[2, 0] - props[2, 1]
    # print('Proportions sum up:', props[1, 0]+props[1, 1]+props[1, 2])
    props[3, 0] = np.divide(1.0, 4.0)   # a bit more if the LAA is displaced to the left
    props[3, 1] = np.divide(1.0, 2.0)   # a bit less if the LAA is displaced to the left
    # props[3, 0] = np.divide(1.0, 4.0) + alpha * np.divide(1.0, 2.0*np.pi)  # a bit more if the LAA is displaced to the left
    # props[3, 1] = np.divide(1.0, 2.0) - alpha * np.divide(1.0, 2.0*np.pi)  # a bit less if the LAA is displaced to the left
    props[3, 2] = np.divide(1.0, 4.0)
    return props

def define_disk_template(rdisk, rhole_rspv, rhole_ripv, rhole_lipv, rhole_lspv, rhole_laa, xhole_center, yhole_center,
                         laa_hole_center_x, laa_hole_center_y, t_v5, t_v6, t_v7, t_v8):
    """Define target positions in the disk template, return coordinates (x,y) corresponding to:
    v1r, v1d, v1l, v2u, v2r, v2l, v3u, v3r, v3l, v4r, v4u, v4d, vlaad, vlaau, p5, p6, p7, p8 """
    coordinates = np.zeros([2, 18])
    complete_circumf_t = np.linspace(0, 2 * np.pi, 1000, endpoint=False)
    rspv_hole_x = np.cos(complete_circumf_t) * rhole_rspv + xhole_center[0]
    rspv_hole_y = np.sin(complete_circumf_t) * rhole_rspv + yhole_center[0]
    ripv_hole_x = np.cos(complete_circumf_t) * rhole_ripv + xhole_center[1]
    ripv_hole_y = np.sin(complete_circumf_t) * rhole_ripv + yhole_center[1]
    lipv_hole_x = np.cos(complete_circumf_t) * rhole_lipv + xhole_center[2]
    lipv_hole_y = np.sin(complete_circumf_t) * rhole_lipv + yhole_center[2]
    lspv_hole_x = np.cos(complete_circumf_t) * rhole_lspv + xhole_center[3]
    lspv_hole_y = np.sin(complete_circumf_t) * rhole_lspv + yhole_center[3]
    laa_hole_x = np.cos(complete_circumf_t) * rhole_laa + laa_hole_center_x
    laa_hole_y = np.sin(complete_circumf_t) * rhole_laa + laa_hole_center_y
    # define (x,y) positions where I put v5, v6, v7 and v8
    coordinates[0, 14] = np.cos(t_v5) * rdisk  # p5_x
    coordinates[1, 14] = np.sin(t_v5) * rdisk  # p5_y
    coordinates[0, 15] = np.cos(t_v6) * rdisk  # p6_x
    coordinates[1, 15] = np.sin(t_v6) * rdisk  # p6_y
    coordinates[0, 16] = np.cos(t_v7) * rdisk  # p7_x
    coordinates[1, 16] = np.sin(t_v7) * rdisk  # p7_y
    coordinates[0, 17] = np.cos(t_v8) * rdisk  # p8_x
    coordinates[1, 17] = np.sin(t_v8) * rdisk  # p8_y

    # define target points corresponding to the pv holes
    # RSPV (right (in the line connecting to MV; left (horizontal line), down, vertical line))
    coordinates[0, 0] = rspv_hole_x[
        np.abs(complete_circumf_t - t_v5).argmin()]  # v1r_x, x in rspv circumf where angle is pi/4
    coordinates[1, 0] = rspv_hole_y[np.abs(complete_circumf_t - t_v5).argmin()]
    coordinates[0, 1] = rspv_hole_x[np.abs(complete_circumf_t - (3 * np.pi / 2)).argmin()]
    coordinates[1, 1] = rspv_hole_y[np.abs(complete_circumf_t - (3 * np.pi / 2)).argmin()]
    coordinates[0, 2] = rspv_hole_x[(np.abs(complete_circumf_t - np.pi)).argmin()]
    coordinates[1, 2] = rspv_hole_y[(np.abs(complete_circumf_t - np.pi)).argmin()]
    # RIPV
    coordinates[0, 3] = ripv_hole_x[np.abs(complete_circumf_t - (np.pi / 2)).argmin()]  # x in ripv circumf UP
    coordinates[1, 3] = ripv_hole_y[np.abs(complete_circumf_t - (np.pi / 2)).argmin()]
    coordinates[0, 4] = ripv_hole_x[np.abs(complete_circumf_t - t_v6).argmin()]
    coordinates[1, 4] = ripv_hole_y[np.abs(complete_circumf_t - t_v6).argmin()]
    coordinates[0, 5] = ripv_hole_x[np.abs(complete_circumf_t - (np.pi)).argmin()]
    coordinates[1, 5] = ripv_hole_y[np.abs(complete_circumf_t - (np.pi)).argmin()]
    # LIPV
    coordinates[0, 6] = lipv_hole_x[np.abs(complete_circumf_t - (np.pi / 2)).argmin()]
    coordinates[1, 6] = lipv_hole_y[np.abs(complete_circumf_t - (np.pi / 2)).argmin()]
    coordinates[0, 7] = lipv_hole_x[complete_circumf_t.argmin()]  # angle = 0
    coordinates[1, 7] = lipv_hole_y[complete_circumf_t.argmin()]
    coordinates[0, 8] = lipv_hole_x[np.abs(complete_circumf_t - t_v7).argmin()]
    coordinates[1, 8] = lipv_hole_y[np.abs(complete_circumf_t - t_v7).argmin()]
    # LSPV
    coordinates[0, 9] = lspv_hole_x[complete_circumf_t.argmin()]  # angle = 0
    coordinates[1, 9] = lspv_hole_y[complete_circumf_t.argmin()]
    coordinates[0, 10] = lspv_hole_x[np.abs(complete_circumf_t - (np.pi / 2)).argmin()]
    coordinates[1, 10] = lspv_hole_y[np.abs(complete_circumf_t - (np.pi / 2)).argmin()]
    coordinates[0, 11] = lspv_hole_x[np.abs(complete_circumf_t - (3 * np.pi / 2)).argmin()]
    coordinates[1, 11] = lspv_hole_y[np.abs(complete_circumf_t - (3 * np.pi / 2)).argmin()]
    # LAA
    coordinates[0, 12] = laa_hole_x[np.abs(complete_circumf_t - (3 * np.pi / 2)).argmin()]
    coordinates[1, 12] = laa_hole_y[np.abs(complete_circumf_t - (3 * np.pi / 2)).argmin()]
    coordinates[0, 13] = laa_hole_x[np.abs(complete_circumf_t - t_v8).argmin()]  # angle = pi/2 + pi/8
    coordinates[1, 13] = laa_hole_y[np.abs(complete_circumf_t - t_v8).argmin()]
    return coordinates

def get_coords(c):
    """Given all coordinates in a matrix, identify and return them separately"""
    return c[0,0], c[1,0], c[0,1], c[1,1], c[0,2], c[1,2], c[0,3], c[1,3], c[0,4], c[1,4], c[0,5], c[1,5], c[0,6], c[1,6], c[0,7], c[1,7], c[0,8], c[1,8], c[0,9], c[1,9], c[0,10], c[1,10], c[0,11], c[1,11], c[0,12], c[1,12], c[0,13], c[1,13], c[0,14], c[1,14], c[0,15], c[1,15], c[0,16], c[1,16], c[0,17], c[1,17]

def extract_LA_contours(m_open, filename, save=False):
    """Given LA with clipped PVs, LAA and MV identify and classify all 5 contours using 'autolabels' array.
    Save contours if save=True"""
    edges = extractboundaryedge(m_open)
    conn = get_connected_edges(edges)
    poly_edges = conn.GetOutput()
    if save==True:
        writevtk(poly_edges, filename[0:len(filename) - 4] + '_detected_edges.vtk')

    print('Detected {} regions'.format(conn.GetNumberOfExtractedRegions()))
    if conn.GetNumberOfExtractedRegions() != 6:
        print(
            'WARNING: the number of contours detected is not the expected. The classification of contours may be wrong')

    for i in range(6):
        print('Detecting region index: {}'.format(i))
        c = pointthreshold(poly_edges, 'RegionId', i, i)
        autolabels = vtk_to_numpy(c.GetPointData().GetArray('autolabels'))
        counts = np.bincount(autolabels.astype(int))
        mostcommon = np.argmax(counts)

        if mostcommon == 36:  # use the most repeated label since some of they are 36 (body). Can be 36 more common in the other regions?
            print('Detected MV')
            if save == True:
                writevtk(c, filename[0:len(filename) - 4] + '_cont_mv.vtk')
            cont_mv = c
        if mostcommon == 37:
            print('Detected LAA')
            if save == True:
                writevtk(c, filename[0:len(filename) - 4] + '_cont_laa.vtk')
            cont_laa = c
        if mostcommon == 76:
            print('Detected RSPV')
            if save == True:
                writevtk(c, filename[0:len(filename) - 4] + '_cont_rspv.vtk')
            cont_rspv = c
        if mostcommon == 77:
            print('Detected RIPV')
            if save == True:
                writevtk(c, filename[0:len(filename) - 4] + '_cont_ripv.vtk')
            cont_ripv = c
        if mostcommon == 78:
            print('Detected LSPV')
            if save == True:
                writevtk(c, filename[0:len(filename) - 4] + '_cont_lspv.vtk')
            cont_lspv = c
        if mostcommon == 79:
            print('Detected LIPV')
            if save == True:
                writevtk(c, filename[0:len(filename) - 4] + '_cont_lipv.vtk')
            cont_lipv = c
    return cont_rspv, cont_ripv, cont_lipv, cont_lspv, cont_mv, cont_laa

def build_locators(mesh, m_open, cont_rspv, cont_ripv, cont_lipv, cont_lspv, cont_laa):
    """Build different locators to find corresponding points between different meshes (open/closed, open/contours, etc)"""
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(mesh)  # clipped + CLOSED - where the seeds are marked
    locator.BuildLocator()

    locator_open = vtk.vtkPointLocator()
    locator_open.SetDataSet(m_open)
    locator_open.BuildLocator()

    locator_rspv = vtk.vtkPointLocator()
    locator_rspv.SetDataSet(cont_rspv)
    locator_rspv.BuildLocator()

    locator_ripv = vtk.vtkPointLocator()
    locator_ripv.SetDataSet(cont_ripv)
    locator_ripv.BuildLocator()

    locator_lipv = vtk.vtkPointLocator()
    locator_lipv.SetDataSet(cont_lipv)
    locator_lipv.BuildLocator()

    locator_lspv = vtk.vtkPointLocator()
    locator_lspv.SetDataSet(cont_lspv)
    locator_lspv.BuildLocator()

    locator_laa = vtk.vtkPointLocator()
    locator_laa.SetDataSet(cont_laa)
    locator_laa.BuildLocator()
    return locator, locator_open, locator_rspv, locator_ripv, locator_lipv, locator_lspv, locator_laa

def read_paths(filename, npaths):
    """read the paths (lines) defined in the 3D mesh using 3_divide_LA.py"""
    for i in range(npaths):
        if os.path.isfile(filename[0:len(filename)-4]+'path'+ str(i+1) +'.vtk')==False:
            sys.exit('ERROR: dividing line path' + str(i+1) + ' not found. Run 3_divide_LA.py')
        else:
            if i == 0:
                path1 = readvtk(filename[0:len(filename) - 4] + 'path' + str(i + 1) + '.vtk')
            elif i == 1:
                path2 = readvtk(filename[0:len(filename) - 4] + 'path' + str(i + 1) + '.vtk')
            elif i == 2:
                path3 = readvtk(filename[0:len(filename) - 4] + 'path' + str(i + 1) + '.vtk')
            elif i == 3:
                path4 = readvtk(filename[0:len(filename) - 4] + 'path' + str(i + 1) + '.vtk')
            elif i == 4:
                path5 = readvtk(filename[0:len(filename) - 4] + 'path' + str(i + 1) + '.vtk')
            elif i == 5:
                path6 = readvtk(filename[0:len(filename) - 4] + 'path' + str(i + 1) + '.vtk')
            elif i == 6:
                path7 = readvtk(filename[0:len(filename) - 4] + 'path' + str(i + 1) + '.vtk')
    # read laa related paths: line from lspv to laa and from laa to mv
    if os.path.isfile(filename[0:len(filename)-4] + 'path_laa1.vtk')==False:
        sys.exit('ERROR: dividing line path_laa1 not found. Run 3_divide_LA.py')
    else:
        path_laa1 = readvtk(filename[0:len(filename)-4] + 'path_laa1.vtk')

    if os.path.isfile(filename[0:len(filename)-4] + 'path_laa2.vtk')==False:
        sys.exit('ERROR: dividing line path_laa2 not found. Run 3_divide_LA.py')
    else:
        path_laa2 = readvtk(filename[0:len(filename)-4] + 'path_laa2.vtk')

    if os.path.isfile(filename[0:len(filename)-4] + 'path_laa3.vtk')==False:
        sys.exit('ERROR: dividing line path_laa3 not found. Run 3_divide_LA.py')
    else:
        path_laa3 = readvtk(filename[0:len(filename)-4] + 'path_laa3.vtk')
    return path1, path2, path3, path4, path5, path6, path7, path_laa1, path_laa2, path_laa3

def get_mv_contour_ids(cont_mv, locator_open):
    """Obtain ids of the MV contour"""
    edge_cont_ids = get_ordered_cont_ids_based_on_distance(cont_mv)
    mv_cont_ids = np.zeros(edge_cont_ids.size)
    for i in range(mv_cont_ids.shape[0]):
        p = cont_mv.GetPoint(edge_cont_ids[i])
        mv_cont_ids[i] = locator_open.FindClosestPoint(p)
    return mv_cont_ids

def identify_segments_extremes(path1, path2, path3, path4, path5, path6, path7, path_laa1, path_laa2, path_laa3,
                               locator_open, locator_rspv, locator_ripv, locator_lipv, locator_lspv, locator_laa,
                               cont_rspv, cont_ripv, cont_lipv, cont_lspv, cont_laa,
                               v5, v6, v7, v8):
    """Identify ids in the to_be_flat mesh corresponding to the segment extremes: v1d, v1r, ect."""
    # start with segments of PVs because they will modify the rest of segments (we try to have uniform number of points in the 3 segments of the veins)
    # first identify ALL pv segments extremes (v1d, v2u etc.)

    # s1 - Find ids corresponding to v1d and v2u as intersection of rspv (ripv) contour and path1
    dists1 = np.zeros(path1.GetNumberOfPoints())
    dists2 = np.zeros(path1.GetNumberOfPoints())
    for i in range(path1.GetNumberOfPoints()):
        p = path1.GetPoint(i)
        dists1[i] = euclideandistance(p, cont_rspv.GetPoint(locator_rspv.FindClosestPoint(p)))
        dists2[i] = euclideandistance(p, cont_ripv.GetPoint(locator_ripv.FindClosestPoint(p)))
    v1d_in_path1 = np.argmin(dists1)
    v2u_in_path1 = np.argmin(dists2)
    v1d = locator_open.FindClosestPoint(path1.GetPoint(v1d_in_path1))
    v2u = locator_open.FindClosestPoint(path1.GetPoint(v2u_in_path1))

    # s2 - Find 2l and v3r
    dists1 = np.zeros(path2.GetNumberOfPoints())
    dists2 = np.zeros(path2.GetNumberOfPoints())
    for i in range(path2.GetNumberOfPoints()):
        p = path2.GetPoint(i)
        dists1[i] = euclideandistance(p, cont_ripv.GetPoint(locator_ripv.FindClosestPoint(p)))
        dists2[i] = euclideandistance(p, cont_lipv.GetPoint(locator_lipv.FindClosestPoint(p)))
    v2l_in_path2 = np.argmin(dists1)
    v3r_in_path2 = np.argmin(dists2)
    v2l = locator_open.FindClosestPoint(path2.GetPoint(v2l_in_path2))
    v3r = locator_open.FindClosestPoint(path2.GetPoint(v3r_in_path2))

    # s3 - Find v3u and v4d
    dists1 = np.zeros(path3.GetNumberOfPoints())
    dists2 = np.zeros(path3.GetNumberOfPoints())
    for i in range(path3.GetNumberOfPoints()):
        p = path3.GetPoint(i)
        dists1[i] = euclideandistance(p, cont_lipv.GetPoint(locator_lipv.FindClosestPoint(p)))
        dists2[i] = euclideandistance(p, cont_lspv.GetPoint(locator_lspv.FindClosestPoint(p)))
    v3u_in_path3 = np.argmin(dists1)
    v4d_in_path3 = np.argmin(dists2)
    v3u = locator_open.FindClosestPoint(path3.GetPoint(v3u_in_path3))
    v4d = locator_open.FindClosestPoint(path3.GetPoint(v4d_in_path3))

    # s4 - Find v4r and v1l
    dists1 = np.zeros(path4.GetNumberOfPoints())
    dists2 = np.zeros(path4.GetNumberOfPoints())
    for i in range(path4.GetNumberOfPoints()):
        p = path4.GetPoint(i)
        dists1[i] = euclideandistance(p, cont_lspv.GetPoint(locator_lspv.FindClosestPoint(p)))
        dists2[i] = euclideandistance(p, cont_rspv.GetPoint(locator_rspv.FindClosestPoint(p)))
    v4r_in_path4 = np.argmin(dists1)
    v1l_in_path4 = np.argmin(dists2)
    v4r = locator_open.FindClosestPoint(path4.GetPoint(v4r_in_path4))
    v1l = locator_open.FindClosestPoint(path4.GetPoint(v1l_in_path4))

    # find ids in the MV
    id_v5 = locator_open.FindClosestPoint(v5)
    id_v6 = locator_open.FindClosestPoint(v6)
    id_v7 = locator_open.FindClosestPoint(v7)
    id_v8 = locator_open.FindClosestPoint(v8)

    # Next 4 segments: s5, s6, s7, s8 : FROM pvs (v1r,v2r,v3l,v4l) TO points in the MV
    dists1 = np.zeros(path5.GetNumberOfPoints())
    for i in range(path5.GetNumberOfPoints()):
        p = path5.GetPoint(i)
        dists1[i] = euclideandistance(p, cont_rspv.GetPoint(locator_rspv.FindClosestPoint(p)))
    v1r_in_path5 = np.argmin(dists1)
    v1r = locator_open.FindClosestPoint(path5.GetPoint(v1r_in_path5))

    # s6
    dists1 = np.zeros(path6.GetNumberOfPoints())
    for i in range(path6.GetNumberOfPoints()):
        p = path6.GetPoint(i)
        dists1[i] = euclideandistance(p, cont_ripv.GetPoint(locator_ripv.FindClosestPoint(p)))
    v2r_in_path6 = np.argmin(dists1)
    v2r = locator_open.FindClosestPoint(path6.GetPoint(v2r_in_path6))

    # s7
    dists1 = np.zeros(path7.GetNumberOfPoints())
    for i in range(path7.GetNumberOfPoints()):
        p = path7.GetPoint(i)
        dists1[i] = euclideandistance(p, cont_lipv.GetPoint(locator_lipv.FindClosestPoint(p)))
    v3l_in_path7 = np.argmin(dists1)
    v3l = locator_open.FindClosestPoint(path7.GetPoint(v3l_in_path7))

    # S8a -> segment from v4 (lspv) to LAA
    dists1 = np.zeros(path_laa1.GetNumberOfPoints())
    dists2 = np.zeros(path_laa1.GetNumberOfPoints())
    for i in range(path_laa1.GetNumberOfPoints()):
        p = path_laa1.GetPoint(i)
        dists1[i] = euclideandistance(p, cont_lspv.GetPoint(locator_lspv.FindClosestPoint(p)))
        dists2[i] = euclideandistance(p, cont_laa.GetPoint(locator_laa.FindClosestPoint(p)))
    v4u_in_pathlaa1 = np.argmin(dists1)
    vlaad_in_pathlaa1 = np.argmin(dists2)
    v4u = locator_open.FindClosestPoint(path_laa1.GetPoint(v4u_in_pathlaa1))
    vlaad = locator_open.FindClosestPoint(path_laa1.GetPoint(vlaad_in_pathlaa1))

    # S8b -> segment from LAA to V8 (MV)
    dists1 = np.zeros(path_laa2.GetNumberOfPoints())
    for i in range(path_laa2.GetNumberOfPoints()):
        p = path_laa2.GetPoint(i)
        dists1[i] = euclideandistance(p, cont_laa.GetPoint(locator_laa.FindClosestPoint(p)))
    vlaau_in_pathlaa2 = np.argmin(dists1)
    vlaau = locator_open.FindClosestPoint(path_laa2.GetPoint(vlaau_in_pathlaa2))

    # aux point vlaar (connecting laa and rspv - auxiliary to know laa contour direction)
    dists1 = np.zeros(path_laa3.GetNumberOfPoints())
    for i in range(path_laa3.GetNumberOfPoints()):
        p = path_laa3.GetPoint(i)
        dists1[i] = euclideandistance(p, cont_laa.GetPoint(locator_laa.FindClosestPoint(p)))
    vlaar_in_pathlaa3 = np.argmin(dists1)
    vlaar = locator_open.FindClosestPoint(path_laa3.GetPoint(vlaar_in_pathlaa3))
    return v1r, v1d, v1l, v2u, v2r, v2l, v3u, v3r, v3l, v4r, v4u, v4d, vlaad, vlaau, vlaar, id_v5, id_v6, id_v7, id_v8

def get_rspv_segments_ids(cont_rspv, locator_open, v1l, v1d, v1r, propn_rspv_s1, propn_rspv_s2, propn_rspv_s3):
    """ Return 3 arrays with ids of each of the 3 segments in rspv contour.
        Return also the modified (to have proportional number of points in the segments) extreme ids"""
    edge_cont_rspv = get_ordered_cont_ids_based_on_distance(cont_rspv)
    rspv_cont_ids = np.zeros(edge_cont_rspv.size)
    for i in range(rspv_cont_ids.shape[0]):
        p = cont_rspv.GetPoint(edge_cont_rspv[i])
        rspv_cont_ids[i] = locator_open.FindClosestPoint(p)
    pos_v1l = int(np.where(rspv_cont_ids == v1l)[0])
    rspv_ids = np.append(rspv_cont_ids[pos_v1l:rspv_cont_ids.size], rspv_cont_ids[0:pos_v1l])
    pos_v1d = int(np.where(rspv_ids == v1d)[0])
    pos_v1r = int(np.where(rspv_ids == v1r)[0])
    if pos_v1r < pos_v1d:   # flip
        aux = np.zeros(rspv_ids.size)
        for i in range(rspv_ids.size):
            aux[rspv_ids.size - 1 - i] = rspv_ids[i]
        # mantain the v1l as the first one (after the flip is the last one)
        flipped = np.append(aux[aux.size - 1], aux[0:aux.size - 1])
        rspv_ids = flipped.astype(int)
    rspv_s1 = rspv_ids[0:int(np.where(rspv_ids == v1d)[0])]
    rspv_s2 = rspv_ids[int(np.where(rspv_ids == v1d)[0]): int(np.where(rspv_ids == v1r)[0])]
    rspv_s3 = rspv_ids[int(np.where(rspv_ids == v1r)[0]): rspv_ids.size]

    # # correct to have proportional segments length
    # s1_prop_length = round(propn_rspv_s1*len(rspv_ids))
    # s2_prop_length = round(propn_rspv_s2*len(rspv_ids))
    # s3_prop_length = round(propn_rspv_s3*len(rspv_ids))
    # v1l_prop = v1l   # stays the same, reference
    # v1d_prop = rspv_ids[int(s1_prop_length)]
    # v1r_prop = rspv_ids[int(s1_prop_length + s2_prop_length)]
    # rspv_s1_prop = rspv_ids[0:int(s1_prop_length)]
    # rspv_s2_prop = rspv_ids[int(s1_prop_length): int(s1_prop_length + s2_prop_length)]
    # rspv_s3_prop = rspv_ids[int(s1_prop_length + s2_prop_length): rspv_ids.size]

    # INTERMEDIATE (final) solution. Offset
    s1_prop_length = round(propn_rspv_s1 * len(rspv_ids))
    s2_prop_length = round(propn_rspv_s2 * len(rspv_ids))
    s3_prop_length = round(propn_rspv_s3 * len(rspv_ids))
    v1l_prop = v1l   # stays the same, reference
    rspv_s1_offset = round((s1_prop_length - rspv_s1.size)/2)        # If negative, I'll shorten s1 in that case
    v1d_prop = rspv_ids[int(rspv_s1.size + rspv_s1_offset)]
    rspv_s1_prop = rspv_ids[0:int(rspv_s1.size + rspv_s1_offset)]
    new_s2_size = rspv_s2.size - rspv_s1_offset   # initial minus points now given to s1
    rspv_s2_offset = np.floor((s2_prop_length - new_s2_size)/2)    # I will add an offset of half the difference. Floor, otherwise s3 is always shorter since it is the remaining part
    v1r_prop = rspv_ids[int(rspv_s1_prop.size + new_s2_size + rspv_s2_offset)]
    rspv_s2_prop = rspv_ids[int(rspv_s1.size + rspv_s1_offset):int(rspv_s1.size + rspv_s1_offset + new_s2_size + rspv_s2_offset)]
    rspv_s3_prop = rspv_ids[int(rspv_s1.size + rspv_s1_offset + new_s2_size + rspv_s2_offset): rspv_ids.size]
    # # print('RSPV original lengths', rspv_s1.size, rspv_s2.size, rspv_s3.size)
    # # print('Proportional lengths', rspv_s1_prop.size, rspv_s2_prop.size, rspv_s3_prop.size)
    return rspv_ids, rspv_s1_prop, rspv_s2_prop, rspv_s3_prop, v1l_prop, v1d_prop, v1r_prop

def get_ripv_segments_ids(cont_ripv, locator_open, v2l, v2r, v2u, propn_ripv_s1, propn_ripv_s2, propn_ripv_s3):
    """ Return 3 arrays with ids of each of the 3 segments in ripv contour.
        Return also the modified (to have proportional number of points in the segments) extreme ids"""
    edge_cont_ripv = get_ordered_cont_ids_based_on_distance(cont_ripv)
    ripv_cont_ids = np.zeros(edge_cont_ripv.size)
    for i in range(ripv_cont_ids.shape[0]):
        p = cont_ripv.GetPoint(edge_cont_ripv[i])
        ripv_cont_ids[i] = locator_open.FindClosestPoint(p)
    pos_v2l = int(np.where(ripv_cont_ids == v2l)[0])
    ripv_ids = np.append(ripv_cont_ids[pos_v2l:ripv_cont_ids.size], ripv_cont_ids[0:pos_v2l])
    pos_v2r = int(np.where(ripv_ids == v2r)[0])
    pos_v2u = int(np.where(ripv_ids == v2u)[0])
    if pos_v2u < pos_v2r:  # flip
        aux = np.zeros(ripv_ids.size)
        for i in range(ripv_ids.size):
            aux[ripv_ids.size - 1 - i] = ripv_ids[i]
        flipped = np.append(aux[aux.size - 1], aux[0:aux.size - 1])
        ripv_ids = flipped.astype(int)
    ripv_s1 = ripv_ids[0:int(np.where(ripv_ids == v2r)[0])]
    ripv_s2 = ripv_ids[int(np.where(ripv_ids == v2r)[0]): int(np.where(ripv_ids == v2u)[0])]
    ripv_s3 = ripv_ids[int(np.where(ripv_ids == v2u)[0]): ripv_ids.size]

    # # # correct to have proportional segments length
    # s1_prop_length = round(propn_ripv_s1 * len(ripv_ids))
    # s2_prop_length = round(propn_ripv_s2 * len(ripv_ids))
    # s3_prop_length = round(propn_ripv_s3 * len(ripv_ids))
    # v2l_prop = v2l  # stays the same, reference
    # v2r_prop = ripv_ids[int(s1_prop_length)]
    # v2u_prop = ripv_ids[int(s1_prop_length + s2_prop_length)]
    # ripv_s1_prop = ripv_ids[0:int(s1_prop_length)]
    # ripv_s2_prop = ripv_ids[int(s1_prop_length): int(s1_prop_length + s2_prop_length)]
    # ripv_s3_prop = ripv_ids[int(s1_prop_length + s2_prop_length): ripv_ids.size]

    # INTERMEDIATE solution.
    s1_prop_length = round(propn_ripv_s1 * len(ripv_ids))
    s2_prop_length = round(propn_ripv_s2 * len(ripv_ids))
    s3_prop_length = round(propn_ripv_s3 * len(ripv_ids))
    v2l_prop = v2l   # stays the same, reference
    ripv_s1_offset = round((s1_prop_length - ripv_s1.size)/2)
    v2r_prop = ripv_ids[int(ripv_s1.size + ripv_s1_offset)]
    ripv_s1_prop = ripv_ids[0:int(ripv_s1.size + ripv_s1_offset)]
    new_s2_size = ripv_s2.size - ripv_s1_offset
    ripv_s2_offset = np.floor((s2_prop_length - new_s2_size)/2)
    v2u_prop = ripv_ids[int(ripv_s1_prop.size + new_s2_size + ripv_s2_offset)]
    ripv_s2_prop = ripv_ids[int(ripv_s1.size + ripv_s1_offset):int(ripv_s1.size + ripv_s1_offset + new_s2_size + ripv_s2_offset)]
    ripv_s3_prop = ripv_ids[int(ripv_s1.size + ripv_s1_offset + new_s2_size + ripv_s2_offset): ripv_ids.size]
    # print('RIPV original lengths', ripv_s1.size, ripv_s2.size, ripv_s3.size)
    # print('Proportional lengths', ripv_s1_prop.size, ripv_s2_prop.size, ripv_s3_prop.size)
    return ripv_ids, ripv_s1_prop, ripv_s2_prop, ripv_s3_prop, v2l_prop, v2r_prop, v2u_prop

def get_lipv_segments_ids(cont_lipv, locator_open, v3r, v3u, v3l, propn_lipv_s1, propn_lipv_s2, propn_lipv_s3):
    """ Return 3 arrays with ids of each of the 3 segments in lipv contour.
        Return also the modified (to have proportional number of points in the segments) extreme ids"""
    edge_cont_lipv = get_ordered_cont_ids_based_on_distance(cont_lipv)
    lipv_cont_ids = np.zeros(edge_cont_lipv.size)
    for i in range(lipv_cont_ids.shape[0]):
        p = cont_lipv.GetPoint(edge_cont_lipv[i])
        lipv_cont_ids[i] = locator_open.FindClosestPoint(p)
    pos_v3r = int(np.where(lipv_cont_ids == v3r)[0])
    lipv_ids = np.append(lipv_cont_ids[pos_v3r:lipv_cont_ids.size], lipv_cont_ids[0:pos_v3r])
    pos_v3u = int(np.where(lipv_ids == v3u)[0])
    pos_v3l = int(np.where(lipv_ids == v3l)[0])
    if pos_v3l < pos_v3u:  # flip
        aux = np.zeros(lipv_ids.size)
        for i in range(lipv_ids.size):
            aux[lipv_ids.size - 1 - i] = lipv_ids[i]
        flipped = np.append(aux[aux.size - 1], aux[0:aux.size - 1])
        lipv_ids = flipped.astype(int)
    lipv_s1 = lipv_ids[0:int(np.where(lipv_ids == v3u)[0])]
    lipv_s2 = lipv_ids[int(np.where(lipv_ids == v3u)[0]): int(np.where(lipv_ids == v3l)[0])]
    lipv_s3 = lipv_ids[int(np.where(lipv_ids == v3l)[0]): lipv_ids.size]

    # # # correct to have proportional segments length
    # s1_prop_length = round(propn_lipv_s1 * len(lipv_ids))
    # s2_prop_length = round(propn_lipv_s2 * len(lipv_ids))
    # s3_prop_length = round(propn_lipv_s3 * len(lipv_ids))
    # v3r_prop = v3r  # stays the same, reference
    # v3u_prop = lipv_ids[int(s1_prop_length)]
    # v3l_prop = lipv_ids[int(s1_prop_length + s2_prop_length)]
    # lipv_s1_prop = lipv_ids[0:int(s1_prop_length)]
    # lipv_s2_prop = lipv_ids[int(s1_prop_length): int(s1_prop_length + s2_prop_length)]
    # lipv_s3_prop = lipv_ids[int(s1_prop_length + s2_prop_length): lipv_ids.size]

    # INTERMEDIATE solution.
    s1_prop_length = round(propn_lipv_s1 * len(lipv_ids))
    s2_prop_length = round(propn_lipv_s2 * len(lipv_ids))
    s3_prop_length = round(propn_lipv_s3 * len(lipv_ids))
    v3r_prop = v3r   # stays the same, reference
    lipv_s1_offset = round((s1_prop_length - lipv_s1.size)/2)
    v3u_prop = lipv_ids[int(lipv_s1.size + lipv_s1_offset)]
    lipv_s1_prop = lipv_ids[0:int(lipv_s1.size + lipv_s1_offset)]
    new_s2_size = lipv_s2.size - lipv_s1_offset
    lipv_s2_offset = np.floor((s2_prop_length - new_s2_size)/2)
    v3l_prop = lipv_ids[int(lipv_s1_prop.size + new_s2_size + lipv_s2_offset)]
    lipv_s2_prop = lipv_ids[int(lipv_s1.size + lipv_s1_offset):int(lipv_s1.size + lipv_s1_offset + new_s2_size + lipv_s2_offset)]
    lipv_s3_prop = lipv_ids[int(lipv_s1.size + lipv_s1_offset + new_s2_size + lipv_s2_offset): lipv_ids.size]
    # print('LIPV original lengths', lipv_s1.size, lipv_s2.size, lipv_s3.size)
    # print('Proportional lengths', lipv_s1_prop.size, lipv_s2_prop.size, lipv_s3_prop.size)
    return lipv_ids, lipv_s1_prop, lipv_s2_prop, lipv_s3_prop, v3r_prop, v3u_prop, v3l_prop

def get_lspv_segments_ids(cont_lspv, locator_open, v4r, v4u, v4d, propn_lspv_s1, propn_lspv_s2, propn_lspv_s3):
    """ Return 3 arrays with ids of each of the 3 segments in lspv contour.
        Return also the modified (to have proportional number of points in the segments) extreme ids"""
    edge_cont_lspv = get_ordered_cont_ids_based_on_distance(cont_lspv)
    lspv_cont_ids = np.zeros(edge_cont_lspv.size)
    for i in range(lspv_cont_ids.shape[0]):
        p = cont_lspv.GetPoint(edge_cont_lspv[i])
        lspv_cont_ids[i] = locator_open.FindClosestPoint(p)
    pos_v4r = int(np.where(lspv_cont_ids == v4r)[0])
    lspv_ids = np.append(lspv_cont_ids[pos_v4r:lspv_cont_ids.size], lspv_cont_ids[0:pos_v4r])
    pos_v4u = int(np.where(lspv_ids == v4u)[0])
    pos_v4d = int(np.where(lspv_ids == v4d)[0])
    if pos_v4d < pos_v4u:   # flip
        aux = np.zeros(lspv_ids.size)
        for i in range(lspv_ids.size):
            aux[lspv_ids.size - 1 - i] = lspv_ids[i]
        flipped = np.append(aux[aux.size - 1], aux[0:aux.size - 1])
        lspv_ids = flipped.astype(int)
    lspv_s1 = lspv_ids[0:int(np.where(lspv_ids == v4u)[0])]
    lspv_s2 = lspv_ids[int(np.where(lspv_ids == v4u)[0]): int(np.where(lspv_ids == v4d)[0])]
    lspv_s3 = lspv_ids[int(np.where(lspv_ids == v4d)[0]): lspv_ids.size]

    ## correct to have proportional segments length
    # s1_prop_length = round(propn_lspv_s1*len(lspv_ids))
    # s2_prop_length = round(propn_lspv_s2*len(lspv_ids))
    # s3_prop_length = round(propn_lspv_s3*len(lspv_ids))
    # v4r_prop = v4r   # stays the same, reference
    # v4u_prop = lspv_ids[int(s1_prop_length)]
    # v4d_prop = lspv_ids[int(s1_prop_length + s2_prop_length)]
    # lspv_s1_prop = lspv_ids[0:int(s1_prop_length)]
    # lspv_s2_prop = lspv_ids[int(s1_prop_length): int(s1_prop_length + s2_prop_length)]
    # lspv_s3_prop = lspv_ids[int(s1_prop_length + s2_prop_length): lspv_ids.size]

    # INTERMEDIATE solution.
    s1_prop_length = round(propn_lspv_s1*len(lspv_ids))
    s2_prop_length = round(propn_lspv_s2*len(lspv_ids))
    s3_prop_length = round(propn_lspv_s3*len(lspv_ids))
    v4r_prop = v4r   # stays the same, reference
    lspv_s1_offset = round((s1_prop_length - lspv_s1.size)/2)
    v4u_prop = lspv_ids[int(lspv_s1.size + lspv_s1_offset)]
    lspv_s1_prop = lspv_ids[0:int(lspv_s1.size + lspv_s1_offset)]
    new_s2_size = lspv_s2.size - lspv_s1_offset
    lspv_s2_offset = np.floor((s2_prop_length - new_s2_size)/2)
    v4d_prop = lspv_ids[int(lspv_s1_prop.size + new_s2_size + lspv_s2_offset)]
    lspv_s2_prop = lspv_ids[int(lspv_s1.size + lspv_s1_offset):int(lspv_s1.size + lspv_s1_offset + new_s2_size + lspv_s2_offset)]
    lspv_s3_prop = lspv_ids[int(lspv_s1.size + lspv_s1_offset + new_s2_size + lspv_s2_offset): lspv_ids.size]
    # print('LSPV Original lengths', lspv_s1.size, lspv_s2.size, lspv_s3.size)
    # print('Proportional lengths', lspv_s1_prop.size, lspv_s2_prop.size, lspv_s3_prop.size)
    return lspv_ids, lspv_s1_prop, lspv_s2_prop, lspv_s3_prop, v4r_prop, v4u_prop, v4d_prop

def get_laa_segments_ids(cont_laa, locator_open, vlaau, vlaad, vlaar):
    """ Return 2 arrays with ids of each of the 2 segments in LAA contour."""
    edge_cont_laa = get_ordered_cont_ids_based_on_distance(cont_laa)
    laa_cont_ids = np.zeros(edge_cont_laa.size)
    for i in range(laa_cont_ids.shape[0]):
        p = cont_laa.GetPoint(edge_cont_laa[i])
        laa_cont_ids[i] = locator_open.FindClosestPoint(p)
    pos_vlaad = int(np.where(laa_cont_ids == vlaad)[0])  # intersection of laa contour and path 8a (from lspv to laa)
    laa_ids = np.append(laa_cont_ids[pos_vlaad:laa_cont_ids.size], laa_cont_ids[0:pos_vlaad])

    pos_vlaar = int(np.where(laa_ids == vlaar)[0])
    pos_vlaau = int(np.where(laa_ids == vlaau)[0])
    if pos_vlaau < pos_vlaar:  # flip
        aux = np.zeros(laa_ids.size)
        for i in range(laa_ids.size):
            aux[laa_ids.size - 1 - i] = laa_ids[i]
        flipped = np.append(aux[aux.size - 1], aux[0:aux.size - 1])
        laa_ids = flipped.astype(int)

    laa_s1 = laa_ids[0:int(np.where(laa_ids == vlaau)[0])]
    laa_s2 = laa_ids[int(np.where(laa_ids == vlaau)[0]): laa_ids.size]
    return laa_ids, laa_s1, laa_s2

def get_segment_ids_in_to_be_flat_mesh(path, locator, intersect_end, intersect_beginning):
    s = np.zeros(path.GetNumberOfPoints())
    for i in range(path.GetNumberOfPoints()):
        p = path.GetPoint(i)
        s[i] = int(locator.FindClosestPoint(p))
    intersect_wlast = np.intersect1d(s, intersect_end)   # find repeated values (s1 merges with rspv contour)
    nlasts_to_delete = len(intersect_wlast)
    index1 = np.arange(len(s) - nlasts_to_delete, len(s))
    final_s = np.delete(s, index1)

    intersect_wfirst = np.intersect1d(final_s, intersect_beginning)
    nfirst_to_delete = len(intersect_wfirst)
    index2 = np.arange(0, nfirst_to_delete)
    s = np.delete(final_s, index2)
    return s

def define_boundary_positions(rdisk, rhole_rspv, rhole_ripv, rhole_lipv, rhole_lspv, rhole_laa, xhole_center, yhole_center, laa_hole_center_x, laa_hole_center_y,
                              s9size, s10size, s11size, s12size, pv_laa_segment_lengths, t_v5, t_v6, t_v7, t_v8):
    """Define BOUNDARY target (x0,y0) coordinates given template parameters (hole radii and positions) and number of points of segments"""
    p_bound = s9size + s10size + s11size + s12size + np.sum(pv_laa_segment_lengths)
    x0_bound = np.zeros(int(p_bound))
    y0_bound = np.zeros(int(p_bound))
    # start with BOUNDARY (disk contour) 4 segments of the mv <-> contour of the disk
    # s9: left
    ind1 = 0
    ind2 = s9size
    t = np.linspace(-(2*np.pi - t_v6), t_v5, s9size+1, endpoint=True)   # +1 because later I will exclude the last point
    # flip to have clock wise direction in the angle
    aux = np.zeros(t.size)
    for i in range(t.size):
        aux[t.size-1-i] = t[i]
    t = aux
    final_t = t[0:len(t)-1]  # exclude extreme, only one, last
    x0_bound[ind1: ind2] = np.cos(final_t) * rdisk
    y0_bound[ind1: ind2] = np.sin(final_t) * rdisk

    # s10: bottom
    ind1 = ind2
    ind2 = ind2 + s10size
    t = np.linspace(t_v7, t_v6, s10size+1, endpoint=True)
    # flip to have clock wise direction in the angle
    aux = np.zeros(t.size)
    for i in range(t.size):
        aux[t.size-1-i] = t[i]
    t = aux
    final_t = t[0:len(t)-1]  # exclude extreme, only one, last
    x0_bound[ind1: ind2] = np.cos(final_t) * rdisk
    y0_bound[ind1: ind2] = np.sin(final_t) * rdisk

    # s11: left - from v7 to v8
    ind1 = ind2
    ind2 = ind2 + s11size
    t = np.linspace(t_v8, t_v7, s11size+1, endpoint=True)
    # flip to have clock wise direction in the angle
    aux = np.zeros(t.size)
    for i in range(t.size):
        aux[t.size-1-i] = t[i]
    t = aux
    final_t = t[0:len(t)-1]  # exclude extreme, only one, last
    x0_bound[ind1: ind2] = np.cos(final_t) * rdisk
    y0_bound[ind1: ind2] = np.sin(final_t) * rdisk

    # s12: top
    ind1 = ind2
    ind2 = ind2 + s12size
    t = np.linspace(t_v5, t_v8, s12size+1, endpoint=True)
    # flip to have clock wise direction in the angle
    aux = np.zeros(t.size)
    for i in range(t.size):
        aux[t.size-1-i] = t[i]
    t = aux
    final_t = t[0:len(t)-1]  # exclude extreme, only one, last
    x0_bound[ind1: ind2] = np.cos(final_t) * rdisk
    y0_bound[ind1: ind2] = np.sin(final_t) * rdisk

    # PV HOLES
    # RSPV, starts in pi
    # rspv_s1
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[0, 0]
    t = np.linspace(np.pi, 3*np.pi/2, pv_laa_segment_lengths[0, 0]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_rspv + xhole_center[0]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_rspv + yhole_center[0]
    # rspv_s2
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[0,1]
    t = np.linspace(3*np.pi/2, t_v5 + 2*np.pi, pv_laa_segment_lengths[0, 1]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_rspv + xhole_center[0]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_rspv + yhole_center[0]
    # rspv_s3
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[0,2]
    t = np.linspace(t_v5, np.pi, pv_laa_segment_lengths[0, 2]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_rspv + xhole_center[0]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_rspv + yhole_center[0]

    # RIPV, starts in pi
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[1,0]
    t = np.linspace(np.pi, t_v6, pv_laa_segment_lengths[1, 0]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_ripv + xhole_center[1]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_ripv + yhole_center[1]
    # ripv_s2
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[1,1]
    t = np.linspace(t_v6, 2*np.pi + np.pi/2, pv_laa_segment_lengths[1, 1]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_ripv + xhole_center[1]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_ripv + yhole_center[1]
    # ripv_s3
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[1,2]
    t = np.linspace(np.pi/2, np.pi, pv_laa_segment_lengths[1, 2]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_ripv + xhole_center[1]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_ripv + yhole_center[1]

    # LIPV, starts in 0
    # lipv_s1
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[2, 0]
    t = np.linspace(0, np.pi/2, pv_laa_segment_lengths[2, 0]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_lipv + xhole_center[2]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_lipv + yhole_center[2]
    # lipv_s2
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[2, 1]
    t = np.linspace(np.pi/2, t_v7, pv_laa_segment_lengths[2, 1]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_lipv + xhole_center[2]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_lipv + yhole_center[2]
    # lipv_s3
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[2, 2]
    t = np.linspace(t_v7, 2*np.pi, pv_laa_segment_lengths[2, 2]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_lipv + xhole_center[2]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_lipv + yhole_center[2]

    # LSPV, starts in 0
    # lspv_s1
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[3, 0]
    t = np.linspace(0, np.pi/2, pv_laa_segment_lengths[3, 0]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_lspv + xhole_center[3]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_lspv + yhole_center[3]
    # lspv_s2
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[3, 1]
    t = np.linspace(np.pi/2, 3*np.pi/2, pv_laa_segment_lengths[3, 1]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_lspv + xhole_center[3]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_lspv + yhole_center[3]
    # lspv_s3
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[3, 2]
    t = np.linspace(3*np.pi/2, 2*np.pi, pv_laa_segment_lengths[3, 2]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_lspv + xhole_center[3]
    y0_bound[ind1: ind2] = np.sin(t) * rhole_lspv + yhole_center[3]

    # LAA hole, circumf
    # laa s1, starts in 3*pi/2
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[4, 0]
    t = np.linspace(3*np.pi/2, t_v8 + 2*np.pi, pv_laa_segment_lengths[4, 0]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_laa + laa_hole_center_x
    y0_bound[ind1: ind2] = np.sin(t) * rhole_laa + laa_hole_center_y
    # laa s2
    ind1 = ind2
    ind2 = ind2 + pv_laa_segment_lengths[4, 1]
    t = np.linspace(t_v8, 3*np.pi/2, pv_laa_segment_lengths[4, 1]+1, endpoint=True)  # skip last one later
    t = t[0:len(t)-1]
    x0_bound[ind1: ind2] = np.cos(t) * rhole_laa + laa_hole_center_x
    y0_bound[ind1: ind2] = np.sin(t) * rhole_laa + laa_hole_center_y
    return x0_bound, y0_bound


def define_constraints_positions(s1, s2, s3, s4, s5, s6, s7, s8a, s8b, v1l_x, v1l_y, v1d_x, v1d_y, v1r_x, v1r_y, v2l_x,
                                 v2l_y, v2r_x, v2r_y, v2u_x, v2u_y, v3r_x, v3r_y, v3u_x, v3u_y, v3l_x, v3l_y,
                                 v4r_x, v4r_y, v4u_x, v4u_y, v4d_x, v4d_y, vlaad_x, vlaad_y, vlaau_x, vlaau_y, p5_x,
                                 p5_y, p6_x, p6_y, p7_x, p7_y, p8_x, p8_y):
    """Define target (x0,y0) coordinates of regional constraints given segments and template parameters (extreme coordinates of segments)"""
    p_const = s1.shape[0] + s2.shape[0] + s3.shape[0] + s4.shape[0] + s5.shape[0] + s6.shape[0] + s7.shape[0] + s8a.shape[0] + s8b.shape[0]
    x0_const = np.zeros(p_const)
    y0_const = np.zeros(p_const)
    # s1, vert line, right
    ind1 = 0
    ind2 = s1.shape[0]
    # vert line
    x0_const[ind1:ind2] = v1d_x
    aux = np.linspace(v1d_y, v2u_y, s1.shape[0] + 2, endpoint=True)
    y0_const[ind1:ind2] = aux[1:aux.size - 1]  # skip first and last

    # s2,  bottom line
    ind1 = ind2
    ind2 = ind2 + s2.shape[0]
    # crosswise lines (all with direction starting in the PV ending in the MV). General rule:
    # m = (y2-y1)/(x2-x1)
    # b = y - m*x
    # y = m*x + b   (any x and y in the line)
    aux = np.linspace(v2l_x, v3r_x, s2.size + 2, endpoint=True)
    x0_const[ind1: ind2] = aux[1:aux.size - 1]
    m = (v3r_y - v2l_y) / (v3r_x - v2l_x)
    b = v3r_y - m * v3r_x
    aux2 = m * aux + b
    y0_const[ind1: ind2] = aux2[1:aux2.size - 1]

    # s3, vert line left
    ind1 = ind2
    ind2 = ind2 + s3.shape[0]
    x0_const[ind1: ind2] = v3u_x
    aux = np.linspace(v3u_y, v4d_y, s3.shape[0] + 2, endpoint=True)
    y0_const[ind1: ind2] = aux[1:aux.size - 1]

    # s4, hori top line
    ind1 = ind2
    ind2 = ind2 + s4.shape[0]
    aux = np.linspace(v4r_x, v1l_x, s4.shape[0] + 2, endpoint=True)
    x0_const[ind1: ind2] = aux[1:aux.size - 1]
    m = (v1l_y - v4r_y) / (v1l_x - v4r_x)
    b = v4r_y - m * v4r_x
    aux2 = m * aux + b
    y0_const[ind1: ind2] = aux2[1:aux2.size - 1]

    # s5 - line crosswise line from v1r to v5
    ind1 = ind2
    ind2 = ind2 + s5.shape[0]
    m = (p5_y - v1r_y) / (p5_x - v1r_x)
    b = v1r_y - m * v1r_x
    aux = np.linspace(v1r_x, p5_x, s5.shape[0] + 2, endpoint=True)
    aux2 = m * aux + b
    x0_const[ind1: ind2] = aux[1:aux.size - 1]
    y0_const[ind1: ind2] = aux2[1:aux2.size - 1]

    # s6 - line crosswise line from v2r to v6
    ind1 = ind2
    ind2 = ind2 + s6.shape[0]
    m = (p6_y - v2r_y) / (p6_x - v2r_x)
    b = v2r_y - m * v2r_x
    aux = np.linspace(v2r_x, p6_x, s6.shape[0] + 2, endpoint=True)
    aux2 = m * aux + b
    x0_const[ind1: ind2] = aux[1:aux.size - 1]
    y0_const[ind1: ind2] = aux2[1:aux2.size - 1]

    # s7 - line crosswise line from v3l to v7
    ind1 = ind2
    ind2 = ind2 + s7.shape[0]
    m = (p7_y - v3l_y) / (p7_x - v3l_x)
    b = v3l_y - m * v3l_x
    aux = np.linspace(v3l_x, p7_x, s7.shape[0] + 2, endpoint=True)
    aux2 = m * aux + b
    x0_const[ind1: ind2] = aux[1:aux.size - 1]
    y0_const[ind1: ind2] = aux2[1:aux2.size - 1]

    # # s8a  - vertical line from lspv (v4u) to laa
    # ind1 = ind2
    # ind2 = ind2 + s8a.shape[0]    # vertical line
    # aux = np.linspace(v4u_y, vlaad_y, s8a.shape[0] + 2, endpoint=True)
    # x0_const[ind1: ind2] = xhole_center[3]
    # y0_const[ind1: ind2] = aux[1:aux.size-1]

    # s8a  - crosswise line from lspv (v4u) to laa
    ind1 = ind2
    ind2 = ind2 + s8a.shape[0]
    aux = np.linspace(v4u_x, vlaad_x, s8a.shape[0] + 2, endpoint=True)
    x0_const[ind1: ind2] = aux[1:aux.size - 1]
    m = (vlaad_y - v4u_y) / (vlaad_x - v4u_x)
    b = v4u_y - m * v4u_x
    aux2 = m * aux + b
    y0_const[ind1: ind2] = aux2[1:aux2.size - 1]

    # s8b- line crosswise line from vlaau to v8
    ind1 = ind2
    ind2 = ind2 + s8b.shape[0]
    m = (p8_y - vlaau_y) / (p8_x - vlaau_x)
    b = vlaau_y - m * vlaau_x
    if p8_x > vlaau_x:
        print('Warning: v8 is greater (in absolute value) than v_laa_up, consider select a different angle for point V8')
    aux = np.linspace(vlaau_x, p8_x, s8b.shape[0] + 2, endpoint=True)
    aux2 = m * aux + b
    x0_const[ind1: ind2] = aux[1:aux.size - 1]
    y0_const[ind1: ind2] = aux2[1:aux2.size - 1]
    return x0_const, y0_const


def ExtractVTKPoints(mesh):
    """Extract points from vtk structures. Return the Nx3 numpy.array of the vertices."""
    n = mesh.GetNumberOfPoints()
    vertex = np.zeros((n, 3))
    for i in range(n):
        mesh.GetPoint(i, vertex[i, :])
    return vertex


def ExtractVTKTriFaces(mesh):
    """Extract triangular faces from vtkPolyData. Return the Nx3 numpy.array of the faces (make sure there are only triangles)."""
    m = mesh.GetNumberOfCells()
    faces = np.zeros((m, 3), dtype=int)
    for i in range(m):
        ptIDs = vtk.vtkIdList()
        mesh.GetCellPoints(i, ptIDs)
        if ptIDs.GetNumberOfIds() != 3:
            raise Exception("Nontriangular cell!")
        faces[i, 0] = ptIDs.GetId(0)
        faces[i, 1] = ptIDs.GetId(1)
        faces[i, 2] = ptIDs.GetId(2)
    return faces


def ComputeLaplacian(vertex, faces):
    """Calculates the laplacian of a mesh
    vertex 3xN numpy.array: vertices
    faces 3xM numpy.array: faces"""
    n = vertex.shape[1]
    m = faces.shape[1]

    # compute mesh weight matrix
    W = sparse.coo_matrix((n, n))
    for i in np.arange(1, 4, 1):
        i1 = np.mod(i - 1, 3)
        i2 = np.mod(i, 3)
        i3 = np.mod(i + 1, 3)
        pp = vertex[:, faces[i2, :]] - vertex[:, faces[i1, :]]
        qq = vertex[:, faces[i3, :]] - vertex[:, faces[i1, :]]
        # normalize the vectors
        pp = pp / np.sqrt(np.sum(pp ** 2, axis=0))
        qq = qq / np.sqrt(np.sum(qq ** 2, axis=0))

        # compute angles
        ang = np.arccos(np.sum(pp * qq, axis=0))
        W = W + sparse.coo_matrix((1 / np.tan(ang), (faces[i2, :], faces[i3, :])), shape=(n, n))
        W = W + sparse.coo_matrix((1 / np.tan(ang), (faces[i3, :], faces[i2, :])), shape=(n, n))

    # compute laplacian
    d = W.sum(axis=0)
    D = sparse.dia_matrix((d, 0), shape=(n, n))
    L = D - W
    return L


def flat(m, boundary_ids, x0, y0):
    """Conformal flattening fitting boundary to (x0,y0) coordinate positions"""
    vertex = ExtractVTKPoints(m).T
    faces = ExtractVTKTriFaces(m).T
    n = vertex.shape[1]
    L = ComputeLaplacian(vertex, faces)

    L = L.tolil()
    L[boundary_ids, :] = 0
    for i in range(boundary_ids.shape[0]):
        L[boundary_ids[i], boundary_ids[i]] = 1

    Rx = np.zeros(n)
    Rx[boundary_ids] = x0
    Ry = np.zeros(n)
    Ry[boundary_ids] = y0
    L = L.tocsr()

    result = np.zeros((Rx.size, 2))
    result[:, 0] = linalg_sp.spsolve(L, Rx)  # x
    result[:, 1] = linalg_sp.spsolve(L, Ry)  # y

    pd = vtk.vtkPolyData()
    pts = vtk.vtkPoints()

    pts.SetNumberOfPoints(n)
    for i in range(n):
        pts.SetPoint(i, result[i, 0], result[i, 1], 0)

    pd.SetPoints(pts)
    pd.SetPolys(m.GetPolys())
    pd.Modified()
    return pd


def flat_w_constraints(m, boundary_ids, constraints_ids, x0_b, y0_b, x0_c, y0_c):
    """ Conformal flattening fitting boundary points to (x0_b,y0_b) coordinate positions
    and additional contraint points to (x0_c,y0_c).
    Solve minimization problem using quadratic programming: https://en.wikipedia.org/wiki/Quadratic_programming"""

    penalization = 1000
    vertex = ExtractVTKPoints(m).T    # 3 x n_vertices
    faces = ExtractVTKTriFaces(m).T
    n = vertex.shape[1]
    L = ComputeLaplacian(vertex, faces)
    L = L.tolil()
    L[boundary_ids, :] = 0.0     # Not conformal there
    for i in range(boundary_ids.shape[0]):
         L[boundary_ids[i], boundary_ids[i]] = 1

    L = L*penalization

    Rx = np.zeros(n)
    Ry = np.zeros(n)
    Rx[boundary_ids] = x0_b * penalization
    Ry[boundary_ids] = y0_b * penalization

    L = L.tocsr()
    # result = np.zeros((Rx.size, 2))

    nconstraints = constraints_ids.shape[0]
    M = np.zeros([nconstraints, n])   # M, zero rows except 1 in constraint point
    for i in range(nconstraints):
        M[i, constraints_ids[i]] = 1
    dx = x0_c
    dy = y0_c

    block1 = hstack([L.T.dot(L), M.T])

    zeros_m = coo_matrix(np.zeros([len(dx),len(dx)]))
    block2 = hstack([M, zeros_m])

    C = vstack([block1, block2])

    prodx = coo_matrix([L.T.dot(Rx)])
    dxx = coo_matrix([dx])
    cx = hstack([prodx, dxx])

    prody = coo_matrix([L.T.dot(Ry)])
    dyy = coo_matrix([dy])
    cy = hstack([prody, dyy])

    solx = linalg_sp.spsolve(C, cx.T)
    soly = linalg_sp.spsolve(C, cy.T)

    # print('There are: ', len(np.argwhere(np.isnan(solx))), ' nans')
    # print('There are: ', len(np.argwhere(np.isnan(soly))), ' nans')
    if len(np.argwhere(np.isnan(solx))) > 0:
        print('WARNING!!! matrix is singular. It is probably due to the convergence of 2 different division lines in the same point.')
        print('Trying to assign different 2D possition to same 3D point. Try to create new division lines or increase resolution of mesh.')

    pd = vtk.vtkPolyData()
    pts = vtk.vtkPoints()

    pts.SetNumberOfPoints(n)
    for i in range(n):
        pts.SetPoint(i, solx[i], soly[i], 0)

    pd.SetPoints(pts)
    pd.SetPolys(m.GetPolys())
    pd.Modified()
    return pd
