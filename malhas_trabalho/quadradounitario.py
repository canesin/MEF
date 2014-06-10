# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:09:23 2014

@author: Fabio Cesar Canesin
@email: Fabio.canesin@gmail.com
@desc: Gerador de malhas para trabalho COC-MEFI
"""

import numpy as np

def quadrado_unitario(ni, nj, outfile, carga=-1.0, quad=True, nu=0.3, E=1.0):
    """
        Gera e exporta quadrado unitario de ni x nj nós, quad ou tri
         /\ y(j)
          |
          |
          |-----ni------
          |            |
          |            |
          |            nj
          |            |
          |            |
         __________________________> x (i)
     """
    x = np.linspace(0, 1, ni)
    y = np.linspace(1, 0, nj)
     
    xi, yj = np.meshgrid(x, y, sparse=False, indexing='ij')
    with open(outfile + '.dat', 'w') as f:
        cabecalho = "{0:>5} {1:>4}    1    {2}    2    2\n".format(ni*nj,
                         (ni-1)*(nj-1) if quad else (2*(ni-1)*(nj-1)),
                         4 if quad else 3)
        f.write(cabecalho)                         
        f.write("  1    {}\n".format(3 if quad else 1))
        f.write("{:>5.1F} {:>5.2F}\n".format(E, nu))
        idx = 1
        for i in xrange(ni):
            for j in xrange(nj):
                node = '{0:>5} {1:>9.3F} {2:>9.3F} {3:>9.3F}\n'
                f.write(node.format(idx, xi[j][i], yj[j][i], 0.0))
                idx +=1
        if not quad:
            nelem = 1
            for j in xrange(nj-1):
                for i in xrange(ni-1):
                    conn = '{:>5} {:>4} {:>4} {:>4}    1\n'
                    pad = j*ni
                    f.write(conn.format(nelem, (i+1)+pad,
                                               (i+1)+ni+pad,
                                               (i+1)+ni+1+pad))
                    nelem += 1                    
                    f.write(conn.format(nelem, (i+1)+pad,
                                               (i+1)+ni+1+pad,
                                               (i+1)+1+pad))
                    nelem += 1
        if quad:
            nelem = 1
            for j in xrange(nj-1):
                for i in xrange(ni-1):
                    conn = '{:>5} {:>4} {:>4} {:>4} {:>4}    1\n'
                    pad = j*ni
                    f.write(conn.format(nelem, (i+1)+pad,
                                               (i+1)+ni+pad,
                                               (i+1)+ni+1+pad,
                                               (i+1)+1+pad))
                    nelem += 1
        for j in xrange(nj):
            pad = j*ni
            constraints = '{:>5}    1    1\n'
            f.write(constraints.format(pad+1))
        f.write('    0    0    0\n')
        for i in xrange(1, ni):
            force = '{:>5}     0.0 {: .5E}\n'
            if i != ni-1:
                nload = carga/(ni-1.0)
            else:
                nload = carga/(2.0*(ni-1.0))
            f.write(force.format(i+1, nload))
        f.write('    0     0.0      0.0\n')


casos = [(5, 5), (9, 9), (17, 17), (33, 33), (65, 65)]
tipos_malha = [True, False] # Quadrangular e triangular

# quadrado_unitario(11, 11, 'testando', carga=-1.0, quad=tipo)
for caso in casos:
    ni, nj = caso
    for tipo in tipos_malha:
        nome = '{}_{}x{}'.format('quad' if tipo else 'tria', ni-1, nj-1)
        quadrado_unitario(ni, nj, nome, carga=-1.0, quad=tipo)
