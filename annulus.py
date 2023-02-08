#! /usr/bin/env python3
# 
# This is an example of using gmsh (with boundary tags) in nutils.
# In this script we solve the Laplace equation :math:`- k ΔT = 0` on a 
# annulus ring domain :math:`Ω` with Dirichlet boundary :math:`Γ`,
# subject to boundary conditions:
#
# .. math::   T &= T_1    && Γ_{\rm inner},
#
#             T &= T_2    && Γ_{\rm outer},
#
# where k denotes the thermal conductivity with unit W/(m·K).
# On the Dirichlet boundaries, the temperature is set to T1 
# and T2, respectively. The exact solution is:
#
# .. math:: T_{\rm exact} = ( ( ln( 0.1 ) - 0.5 ln( x_0^2 + x_1^2 ) ) / ( ln( 0.1 ) - ln( 0.15 ) ) ) ( T_2 - T_1 ) + T_1.
#
# We start by importing the necessary modules.

import numpy, treelog
from nutils import mesh, function, solver, export, cli
from nutils.expression_v2 import Namespace
from matplotlib import collections

# main
def main(fname: str, degree: int, k: float, T1: float, T2: float):
    '''
    Laplace problem on a ring.

    .. arguments::

       fname [meshfiles/annulus0.msh]
         Name of gmsh file.
       degree [1]
         Polynomial degree.
       k [0.25]
         Diffusivity constant in (w/m-K).
       T1 [70.0]
         Inner ring temperature (K).
       T2 [20.0]
         Outer ring temperature (K).
    '''

    # construct mesh
    domain, geom = mesh.gmsh(fname)

    # create namespace
    ps = Namespace()
    ps.x  = geom
    ps.define_for('x', gradient='∇', normal='n', jacobians=('dV', 'dS'))
    ps.k  = k
    ps.T1 = T1
    ps.T2 = T2

    # construct basis function
    ps.basis = domain.basis('std', degree=degree)

    # discretized solution
    ps.T = function.dotarg('lhs', ps.basis)

    # define weakform
    res = domain.integral('k ∇_i(basis_n) ∇_i(T) dV' @ ps, degree=degree*2)

    # define dirichlet constraints
    sqr  = domain.boundary['inner'].integral('(T - T1)^2 dS' @ ps, degree=degree*2)
    sqr += domain.boundary['outer'].integral('(T - T2)^2 dS' @ ps, degree=degree*2)
    cons = solver.optimize('lhs', sqr, droptol=1e-15)

    # assembled and solved
    lhs = solver.solve_linear('lhs', res, constrain=cons)

    # post-process solution 
    bezier = domain.sample('bezier', 9)
    x, T   = bezier.eval(['x_i', 'T'] @ ps, lhs=lhs)
    triplot(bezier, x, T, 'solution')

    # error
    ps.Texact = '( ( ln( 0.1 ) - 0.5 ln( x_0^2 + x_1^2 ) ) / ( ln( 0.1 ) - ln( 0.15 ) ) ) ( T2 - T1 ) + T1 '
    ps.err    = 'T - Texact'
    
    # post-process absolute error
    bezier = domain.sample('bezier', 9)
    x, e   = bezier.eval(['x_i', 'err'] @ ps, lhs=lhs)
    triplot(bezier, x, e, 'error')

    # norm error
    L2err = domain.integral('( err )^2 dV' @ ps, degree=degree*2).eval(lhs=lhs)**.5
    H1err = domain.integral('( ( err )^2 + ∇_i(err) ∇_i(err)) dV' @ ps, degree=degree*2).eval(lhs=lhs)**.5
    treelog.user('L2 error: {:.2e}'.format(L2err))
    treelog.user('H1 error: {:.2e}'.format(H1err))

    return L2err, H1err

def convergence(ncases=4, degree=1, k=0.25, T1=70, T2=20):
    
    # define mesh sizes 
    allh  = [0.01, 0.005, 0.0025, 0.00125]
    h     = numpy.array(allh[:ncases])

    # define empty arrays
    L2err = numpy.empty(len(h))
    H1err = numpy.empty(len(h))

    # loop over number of cases
    with treelog.iter.fraction('mesh', range(ncases)) as counter:
        for icase in counter:
            fname = 'meshfiles/annulus'+str(icase)+'.msh'
            # compute solution and errors
            L2err[icase], H1err[icase] = main(fname=fname, degree=degree ,k=k, T1=T1, T2=T2)

    # plot convergence
    with export.mplfigure('convergence.png',dpi=300) as fig:
        ax = fig.add_subplot(111)
        ax.set_xlabel(r'1/h')
        ax.set_ylabel(r'error')
        # plot computed errors
        ax.loglog(1/h, L2err, 'bo-', label='L²-norm')
        ax.loglog(1/h, H1err, 'ro-', label='H¹-norm')
        slope_triangle(fig, ax, 1/h, L2err)
        slope_triangle(fig, ax, 1/h, H1err)
        ax.legend(loc='lower left')
        ax.grid(which='both', axis='both', linestyle=":") 
        fig.set_tight_layout(True)

def triplot(bezier, points, value, name):

    with export.mplfigure(name+'.png', dpi=300) as fig:
        ax  = fig.add_subplot(111, title=name)
        im  = ax.tripcolor(points[:,0], points[:,1], bezier.tri, value, cmap='jet')
        fig.colorbar(im)
        ax.add_collection(collections.LineCollection(points[bezier.hull,:2], colors='k', linewidths=0.5, alpha=0.5))
        ax.autoscale(enable=True, axis='both', tight=True)
        ax.set_aspect('equal')

# plot slope triangle
def slope_triangle(fig, ax, x, y):

    # choose the points (from last) to show slope
    i, j = (-2, -1) if x[-1] < x[-2] else (-1, -2)
    if not all(numpy.isfinite(x[-2:])) or not all(numpy.isfinite(y[-2:])):
      treelog.warning('Not plotting slope triangle for +/-inf or nan values')
      return

    from matplotlib import transforms
    shifttrans = ax.transData + transforms.ScaledTranslation(0, 0.1, fig.dpi_scale_trans)
    xscale, yscale = ax.get_xscale(), ax.get_yscale()

    # delta() checks if either axis is log or lin scaled
    delta = lambda a, b, scale: numpy.log10(float(a) / b) if scale=='log' else float(a - b) if scale=='linear' else None
    slope = delta(y[-2], y[-1], yscale) / delta(x[-2], x[-1], xscale)
    if slope in (numpy.nan, numpy.inf, -numpy.inf):
      treelog.warning(f'Cannot draw slope triangle with slope: {str(slope)}, drawing nothing')
      return slope

    # handle positive and negative slopes correctly
    xtup, ytup = ((x[i], x[j], x[i]), (y[j], y[j], y[i]))
    a, b = (2/3., 1/3.)
    xval = a*x[i] + b*x[j]
    yval = b*y[i] + a*y[j]

    ax.fill(xtup, ytup, color='0.9', edgecolor='k', transform=shifttrans)
    slopefmt='{0:.1f}'
    ax.text(xval, yval, slopefmt.format(slope), horizontalalignment='center', verticalalignment='center', transform=shifttrans, fontsize=12)

cli.choose(main,convergence)
