# -*- coding: utf-8 -*-
"""
Created on Mon Apr 09 17:15:00 2014

@author: Fabio Cesar Canesin
@email: Fabio.canesin@gmail.com
@desc: Conversor malha programa disciplina COC-MEF I para VTU
"""

import re
import numpy as np
import StringIO
import vtk
from vtk.util import numpy_support as npvtk

class Frozendict(dict):
    """ Immutable dictionary implementation
    """
    def _immutable(self, *args, **kws):
        """ Raise error if try to change value
        """
        raise TypeError('object is immutable')

    __setitem__ = _immutable
    __delitem__ = _immutable
    clear = _immutable
    update = _immutable
    setdefault = _immutable
    pop = _immutable
    popitem = _immutable

#Dictionary of VTK to VISAGE element number correpondence
MJC2VTK = {
# MJC : VTK
    #None: (5, 3, 1), # TRIANGLE - 3 NODES
    #None: 6, # TRIANGLE STRIP
    #None: (9, 4, 1), # QUAD - 4 NODES
    #None: 10, # TETRA - 4 NODES
    #None: 11, # VOXEL - 8 NODES 90 DEGREE AXIS PARALLEL BRICK
    #None: 12, # HEXAHEDRON - 8 NODES
    #None: 13, # WEDGE - 6 NODES
    #None: 14, # PYRAMID - 5 NODES
    #None: 15, # PENTAGONAL PRISM - 10 NODES
    #None: 16, # HEXAGONAL PRISM - 12 NODES
    #None: (21, 3, 2), # QUADRATIC EDGE - 3 NODES
    #None: (22, 6, 2), # QUADRATIC TRIANGLE - 6 NODES
    'Q8': (23, 8, 2), # QUADRATIC QUAD - 8 NODES
    #None: 24, # QUADRATIC TETRA - 10 NODES
    #None: 25, # QUADRATIC HEXAHEDRON - 20 NODES
    #None: 26, # QUADRATIC WEDGE - 15 NODES
    #None: 27, # QUADRATIC PYRAMID - 13 NODES
    #None: 28, # BIQUADRATIC QUAD - 9 NODES
    #None: 29, # TRIQUADRATIC HEXAHEDRON - 27 NODES
    #None: 30, # QUADRATIC LINEAR QUAD - 6 NODES
    #None: 31, # QUADRATIC LINEAR WEDGE - 12 NODES
    #None: 32, # BIQUADRATIC QUADRATIC WEDGE - 18 NODES
    #None: 33, # BIQUADRATIC QUADRATIC HEXAHEDRON - 24 NODES
    #None: 34, # BIQUADRATIC TRIANGLE - 7 NODES
    }

MJC2VTK = Frozendict(MJC2VTK)
#Elementos implementados

class Grid(object):
    """ Class that holds geometry data and grid for core test
    """
    def __init__(self, meshpath=None):
        """ Input mesh address to read
        """
        if meshpath is None:
            raise Exception("Missing needed geometric description")
        self.mesh = readfile(meshpath)

    def SaveVTK(self, outfile="Output"):
        """ Save grid in vtk structured format.

            Parameters
            ----------
            grid: vtkStructuredGrid : grid used in computation
            outfile: string : output file name in system
                            Defaults to "Output"
            Return
            ------
            pass
            save file directly on outfile location
        """
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(outfile + ".vtu")
        writer.SetInput(self.mesh)
        writer.Write()

def parsedata(lhs, rhs):
    """ Soma duas strings substituindo characteres merda por tabulação
    """
    if len(lhs) <= len(rhs):
        lhs = re.sub(r'\s+', ' ', lhs, flags=re.IGNORECASE)
    rhs = re.sub(r'\s+', ' ', rhs, flags=re.IGNORECASE)
    return '\n'.join((lhs, rhs))

def loadarray(dtp, datalist):
    """ Carrega um array de arrays segundo pad e lista de dados
    """
    nds = reduce(parsedata, datalist).strip()
    nds = StringIO.StringIO('\n'.join((dtp['names'], nds)))
    return np.genfromtxt(nds, names=True, dtype=dtp['type']).view(np.recarray)

def read_scalar_results(grid, data):
    """ Read results from scalr fields
    """
    nnodes = grid.GetNumberOfPoints()
    #ncells = grid.GetNumberOfCells()

    for idx, line in enumerate(data):
        mnodal_scalar_fields = re.match(r".*\.NODAL\.SCALAR$", line)
        mnodal_scalar_data = re.match(r".*\.NODAL\.SCALAR.DATA$", line)
        if mnodal_scalar_fields:
            nodal_scalar_fields = (int(data[idx+1].strip()), data[idx+2])
        if mnodal_scalar_data:
            dttype = {'names': 'id '+nodal_scalar_fields[1],
                      'type': tuple([int]+nodal_scalar_fields[0]*[float])}
            nodal_scalar_data = loadarray(dttype, data[idx+2:idx+nnodes+2])

    nodal_scalar_fields = re.findall(r"'(\w+\s*\w*)'", nodal_scalar_fields[1])
    for field in nodal_scalar_fields:
        vtkfield = npvtk.numpy_to_vtk(np.ascontiguousarray(nodal_scalar_data[field]), 1)
        vtkfield.SetName(field)
        grid.GetPointData().AddArray(vtkfield)
    pass

def readfile(meshpath):
    """ Generate half-cylindrical mesh for core test simulation.

        Parameters
        ----------
        meshpath : string : path to mesh file

        Return
        ------
        mesh : vtkUnstructuredGrid object.
    """
    arquivo = open(meshpath, 'r').readlines()
    data = []
    #Filtra as linhas em branco que só fodem a vida
    for idx, linha in enumerate(arquivo):
        if re.match(r'^\s*$', linha):
            pass
        else:
            data.append(linha.strip())
    del arquivo

    for idx, linha in enumerate(data):
        if linha == '%NODE':
            npts = int(data[idx+1])
        if linha == '%NODE.COORD':
            nodesidx = idx + 2
        if linha == '%ELEMENT':
            nelem = int(data[idx+1])
        if '%ELEMENT.' in linha:
            etype = linha.partition('.')[-1]
            elemsidx = idx + 2
        if linha == '%RESULT':
            residx = idx

    #Short cut se o tipo de elemetno nao puder ser lido
    if etype not in MJC2VTK:
        raise Exception(u"Tipo de elemento não implementado")

    #Alocando memória para a malha
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(npts)
    cells = vtk.vtkCellArray()

    def read_nodes():
        """ Read nodes from data list and provide data structure
        """
        dtype = {'names': 'id x y z', 'type': (int, float, float, float)}
        nds = loadarray(dtype, data[nodesidx:nodesidx+npts])
        nds.sort(axis=0, order=['id'])
        return nds

    def read_elements():
        """ Read elements from data list and provide data structure
        """
        nstring = 'id mat thick ordem conect'+(MJC2VTK[etype][1]-1)*' {}'
        dtype = {'names': nstring.format(*range(MJC2VTK[etype][1]-1)),
                 'type': (int, int, float, float, str(MJC2VTK[etype][1])+'int')}
        els = loadarray(dtype, data[elemsidx:elemsidx+nelem])
        els.sort(axis=0, order=['id'])
        return els
    #Read mesh structure
    nodes = read_nodes()
    elements = read_elements()
    #Free some memory by continuing from results only
    data = data[residx:]

    # Loop nos nos e preenche as coordenadas
    for i in xrange(npts):
        points.SetPoint(nodes.id[i]-1, nodes.x[i], nodes.y[i], nodes.z[i])
    #Loop nos elementos e depois loop na conectividade
    for i in xrange(nelem):
        if MJC2VTK[etype][0] == 23:
            cellobj = vtk.vtkQuadraticQuad()

        if MJC2VTK[etype][2] is 2:
            nnodestmp = MJC2VTK[etype][1]
            ordemconect = tuple(range(0, nnodestmp, 2)+range(1, nnodestmp, 2))
        else:
            ordemconect = tuple(range(MJC2VTK[etype][1]))

        for idx, j in enumerate(ordemconect):
            cellobj.GetPointIds().SetId(idx, elements.conect[i][j]-1)
        cells.InsertNextCell(cellobj)

    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.SetCells(MJC2VTK[etype][0], cells)
    read_scalar_results(grid, data)
    return grid

if __name__ == '__main__':
	""" Exemplo de uso
	"""
    GRID = Grid('n.pos')
    GRID.SaveVTK("noutf2")
