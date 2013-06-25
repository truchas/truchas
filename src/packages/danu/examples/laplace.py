#!/usr/bin/env python

'''
LaPlace Solver on a fixed grid using Danu for output
Solver based on example laplace.py file I found on
the SciPy web site.

Need numpy and scipy 
'''

import sys, os
try:
  import numpy
except ImportError:
  print 'This script requires numpy'
  raise ImportError


try:
  import Danu
except ImportError:
  raise ImportError, 'This script requires the Danu Python module. Append the install path to PYTHONPATH'

class Grid:
    
    """A simple grid class that stores the details and solution of the
    computational grid."""
    
    def __init__(self, nx=10, ny=10, xmin=0.0, xmax=1.0,
                 ymin=0.0, ymax=1.0):
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.nx, self.ny = nx, ny
        self.dx = float(xmax-xmin)/(nx-1)
        self.dy = float(ymax-ymin)/(ny-1)
        self.u = numpy.zeros((nx, ny), 'd')
        # used to compute the change in solution in some of the methods.
        self.old_u = self.u.copy()        

    def setBC(self, l, r, b, t):        
        """Sets the boundary condition given the left, right, bottom
        and top values (or arrays)"""        
        self.u[0, :] = l
        self.u[-1, :] = r
        self.u[:, 0] = b
        self.u[:,-1] = t
        self.old_u = self.u.copy()

    def setBCFunc(self, func):
        """Sets the BC given a function of two variables."""
        xmin, ymin = self.xmin, self.ymin
        xmax, ymax = self.xmax, self.ymax
        x = numpy.arange(xmin, xmax + self.dx*0.5, self.dx)
        y = numpy.arange(ymin, ymax + self.dy*0.5, self.dy)
        self.u[0 ,:] = func(xmin,y)
        self.u[-1,:] = func(xmax,y)
        self.u[:, 0] = func(x,ymin)
        self.u[:,-1] = func(x,ymax)

    def computeError(self):        
        """Computes absolute error using an L2 norm for the solution.
        This requires that self.u and self.old_u must be appropriately
        setup."""        
        v = (self.u - self.old_u).flat
        return numpy.sqrt(numpy.dot(v,v))

    def setupOutput(self,h5file='laplace-test.h5'):
        
        import os
        if os.path.exists(h5file):
          os.remove(h5file)
        output=Danu.Output(h5file,'w')
        sim=output.add_simulation('Heat Equation')

        return
        
    def dumpMesh(self,h5file='laplace-test.h5'):
        
        output=Danu.Output(h5file,'a')
        mesh=Danu.Mesh(output,'Heat Equation Cell Center Mesh', Danu.STRUCTURED_MESH, Danu.QUAD_ELEM)
        x=numpy.linspace(self.xmin,self.xmax,self.nx)
        y=numpy.linspace(self.ymin,self.ymax,self.ny)
        mesh.write_coordinates(x,y)
        mesh.set_attribute('X Coordinate Column',0)
        mesh.set_attribute('Y Coordinate Column',1)

        return

    def dumpResult(self,h5file='laplace-test.h5', iteration=0, t=0.0):

        output=Danu.Output(h5file,'w')
        sim=output.get_simulation('Heat Equation')
        seq=sim.get_nextSequence(iteration,t)
        seq.data_write('u solution', self.u)
        
        return
         

        


class LaplaceSolver:

    """A simple Laplacian solver"""
    
    def __init__(self, grid):
        self.grid = grid

    def timeStep(self, dt=0.0):

        """Takes a time step using a numeric expressions."""
        g = self.grid
        dx2, dy2 = g.dx**2, g.dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        u = g.u
        g.old_u = u.copy()

        # The actual iteration
        u[1:-1, 1:-1] = ((u[0:-2, 1:-1] + u[2:, 1:-1])*dy2 + 
                         (u[1:-1,0:-2] + u[1:-1, 2:])*dx2)*dnr_inv
        
        return g.computeError()

    def solve(self, n_iter=0, eps=1.0e-16):        
        """Solves the equation given an error precision -- eps.  If
        n_iter=0 the solving is stopped only on the eps condition.  If
        n_iter is finite then solution stops in that many iterations
        or when the error is less than eps whichever is earlier.
        Returns the error if the loop breaks on the n_iter condition
        and returns the iterations if the loop breaks on the error
        condition."""        
        err = self.timeStep()
        count = 1

        while err > eps:
            if n_iter and count >= n_iter:
                self.grid.dumpResult(iteration=count)
                return err
            err = self.timeStep()
            if count%10 == 0:
              self.grid.dumpResult(iteration=count)
            count = count + 1

        return count


if __name__ == '__main__':

  # Define the grid
  grid=Grid()

  # Set BC
  grid.setBC(l=-1.0,r=1.0,t=1.0,b=-1.0)

  # Setup the output file
  grid.setupOutput()

  # Dump out the Mesh
  grid.dumpMesh()

  solver=LaplaceSolver(grid)
  print solver.solve()
  print grid.u



