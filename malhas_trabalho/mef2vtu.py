# -*- coding: utf-8 -*-
"""
Created on Mon Apr 09 17:15:00 2014

@author: Fabio Cesar Canesin
@email: Fabio.canesin@gmail.com
@desc: Conversor malha programa disciplina COC-MEFI para VTU
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

#Dicionario de correspondencia tipo de elementos para MEF e VTK
MJC2VTK = {
# MJC : VTK
	#Numero no programa -- ( Numero no VTK, Numero de nós no elemento, 1-linear 2-quadratico)
    3: (5, 3, 1), # TRIANGLE - 3 NODES
    #None: 6, # TRIANGLE STRIP
    4: (9, 4, 1), # QUAD - 4 NODES
    #None: 10, # TETRA - 4 NODES
    #None: 11, # VOXEL - 8 NODES 90 DEGREE AXIS PARALLEL BRICK
    #None: 12, # HEXAHEDRON - 8 NODES
    #None: 13, # WEDGE - 6 NODES
    #None: 14, # PYRAMID - 5 NODES
    #None: 15, # PENTAGONAL PRISM - 10 NODES
    #None: 16, # HEXAGONAL PRISM - 12 NODES
    #None: (21, 3, 2), # QUADRATIC EDGE - 3 NODES
    #None: (22, 6, 2), # QUADRATIC TRIANGLE - 6 NODES
    #None: (23, 8, 2), # QUADRATIC QUAD - 8 NODES
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

def loadarray(dtp, datalist):
    """ Carrega um array de arrays segundo pad e lista de dados
    """
    datastring = ''
    for line in datalist:
        datastring += line+'\n'
    nds = StringIO.StringIO(datastring)
    return np.genfromtxt(nds, names=True, dtype=dtp['type']).view(np.recarray)

def read_scalar_results(grid, data):
    """ Read results from scalr fields
    """
    nnodes = grid.GetNumberOfPoints()
    # Vamos ler os campos 'id DX DY DZ'
    nodal_scalar_fields = "DX DY DZ"

    for idx, line in enumerate(data):
        # Caso estejamos no linha que tem a palavra 'nvec'
        if 'nvec' in line.split():
            # Descrever o tipo de dado que existe nesse resultado
            dttype = {'names': 'id '+nodal_scalar_fields, # Concatenacao de strings
                          'type': tuple([int]+3*[float])} # Definindo estrutura de dados
            # Agora carregamos os dados, apartir de idx+1 até o final dos numero de nós
            desloc = loadarray(dttype, [dttype['names']]+data[idx+1:idx+1+nnodes])
            #for field in nodal_scalar_fields.split():
            vecdata = np.array([desloc.DX, desloc.DY, desloc.DZ])
            vtkfield = npvtk.numpy_to_vtk(np.ascontiguousarray(vecdata.T), 1)
            vtkfield.SetName("Deslocamento")
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

        if 'coor' in linha.split():
            npts = int(linha.split()[1]) #Numero de coordenadas (nos/pontos)
            nodesidx = idx + 1 #Posicao em data que comecam os nos(pontos)
        if 'elem' in linha.split():
            nelem = int(linha.split()[1]) #Numero de elementos
            etype = int(data[idx+1].split()[1]) #Nome do elemento no MJC2VTK
            elemsidx = idx + 1 #Posicao em data que comeca a conectividade
        if 'nvec' in linha.split():
            residx = idx+1 #Posicao em data que comeca os resultados

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
        nds = loadarray(dtype, [dtype['names']]+data[nodesidx:nodesidx+npts])
        return nds

    def read_elements():
        """ Read elements from data list and provide data structure
        """
        nstring = 'id type conect'+(MJC2VTK[etype][1]-1)*' {}'
        dtype = {'names': nstring.format(*range(MJC2VTK[etype][1]-1)),
                 'type': (int, int, str(MJC2VTK[etype][1]) + 'int')}
        els = loadarray(dtype, [dtype['names']]+data[elemsidx:elemsidx+nelem])
        return els

    #Read mesh structure
    nodes = read_nodes()
    elements = read_elements()
    #Free some memory by continuing from results only
    data = data[residx-1:] #'no no no elem elem elem elem --> res res res ...'
    # Loop nos nos e preenche as coordenadas
    for i in xrange(npts):
        points.SetPoint(nodes[i][0]-1, nodes[i][1], nodes[i][2], nodes[i][3])
    #Loop nos elementos e depois loop na conectividade
    for i in xrange(nelem):
        #Verifica se e um quadrangular quadratico
        if MJC2VTK[etype][0] == 23:
            cellobj = vtk.vtkQuadraticQuad()
        #Verifica se e um triangulo
        if MJC2VTK[etype][0] == 5:
            cellobj = vtk.vtkTriangle()
        #Verifica se e um quadrilatero
        if MJC2VTK[etype][0] == 9:
            cellobj = vtk.vtkQuad()
		#Verifica se e quadratico
        if MJC2VTK[etype][2] is 2:
            nnodestmp = MJC2VTK[etype][1]
            ordemconect = tuple(range(0, nnodestmp, 2)+range(1, nnodestmp, 2))
        #Caso o elemento seja linear:
        else:
            # !! CHECAR ordem da conectividade
            ordemconect = tuple(range(MJC2VTK[etype][1])) # (0, 1, 2)

        for idx, j in enumerate(ordemconect): # (0[0], 1[1], 2[2])
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

    casos = [(5, 5), (9, 9), (17, 17), (33, 33), (65, 65)]
    tipos_malha = [True, False] # Quadrangular e triangular

    # quadrado_unitario(11, 11, 'testando', carga=-1.0, quad=tipo)
    for caso in casos:
        ni, nj = caso
        for tipo in tipos_malha:
            nome = '{}_{}x{}'.format('quad' if tipo else 'tria', ni-1, nj-1)
            GRID = Grid(nome+'.out')
            GRID.SaveVTK(nome)
