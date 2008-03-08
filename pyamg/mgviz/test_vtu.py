""" Using the examples in the examples directory"""

#import os
#from scipy import rand, array
#from scipy.io import loadmat
#from pyamg.gallery import load_example
from scipy.testing import *
import xml.parsers.expat
from numpy import array, uint32
from mgviz import write_vtu

class TestWriteVtu(TestCase):
    def setUp(self):
        #cdata = ({'5':random.random((E2V.shape[0],1))}, {'5':2*random.random((E2V.shape[0],1))})
        #data = zeros((Vert.shape[0],1))
        #data[5:10]=1
        #pdata=concatenate((random.random((Vert.shape[0],1)),data),1)
        cases = []
        class mesh:
            file='test.vtu'
            Vert=None
            E2V=None
            pdata=None
            cdata=None
        mesh=mesh()
        mesh.file = 'test.vtu'

        # 1 triangle
        mesh.Vert = array([[0.0,0.0],
                           [0.0,1.0],
                           [1.0,1.0]])
        E2V = array([[0,2,1]],uint32)
        mesh.Cells = {'5':E2V}
        mesh.pdata = None
        mesh.cdata = None
        cases.append(mesh)
        
        # 2 triangles
        mesh.Vert = array([[0.0,0.0],
                           [1.0,0.0],
                           [0.0,1.0],
                           [1.0,1.0]])
        E2V = array([[0,3,2],
                     [0,1,3]],uint32)
        mesh.Cells = {'5':E2V}
        mesh.pdata = None
        mesh.cdata = None
        cases.append(mesh)

        # 8 triangles
        mesh.Vert = array([[0.0,0.0],
                           [1.0,0.0],
                           [2.0,0.0],
                           [0.0,1.0],
                           [1.0,1.0],
                           [2.0,1.0],
                           [0.0,2.0],
                           [1.0,2.0],
                           [2.0,2.0]])
        E2V = array([[0,4,3],
                     [0,1,4],
                     [1,5,4],
                     [1,2,5],
                     [3,7,6],
                     [3,4,7],
                     [4,8,7],
                     [4,5,8]],uint32)
        mesh.Cells = {'5':E2V}
        mesh.pdata = None
        mesh.cdata = None
        cases.append(mesh)

        self.cases=cases

    def test_xml(self):
        for mesh in self.cases:
            try:
                write_vtu(mesh.Vert,mesh.Cells,mesh.file,mesh.pdata,mesh.cdata)
            except:
                assert False, 'cannot write test.vtu'
            try:
                parser = xml.parsers.expat.ParserCreate()
                parser.ParseFile(open(mesh.file, 'r'))
            except Exception, ex:
                assert False, 'problem: %s' % (ex)

if __name__ == '__main__':
    nose.run(argv=['', __file__])
