import pyvista
import numpy as np
from dolfinx import geometry, plot


def evaluate_on_points(function,points):
    mesh = function.function_space.mesh
    comm = mesh.comm
    # comm
    bb_tree = geometry.BoundingBoxTree(mesh, mesh.topology.dim)
    cells = []
    points_on_proc = []
    # Find cells whose bounding-box collide with the the points
    cell_candidates = geometry.compute_collisions(bb_tree, points.T)
    # Choose one of the cells that contains the point
    colliding_cells = geometry.compute_colliding_cells(mesh, cell_candidates, points.T)
    for i, point in enumerate(points.T):
        if len(colliding_cells.links(i)) > 0:
            points_on_proc.append(point)
            cells.append(colliding_cells.links(i)[0])

    points_on_proc = np.array(points_on_proc)
    if len(points_on_proc) > 0:
        values_on_proc = function.eval(points_on_proc, cells)
        point_data_proc = [points_on_proc.T, values_on_proc]
    else:
        point_data_proc = None

    point_data = comm.gather(point_data_proc, root=0)

    if comm.rank == 0:
        point_data = list(filter(None, point_data))
        points = np.concatenate([data_proc[0].T for data_proc in point_data])
        values = np.concatenate([data_proc[1] for data_proc in point_data])
    else:
        points = None
        values = None

    return values

def warp_plot_2d(u,cell_field=None,field_name="Field",factor=1.,backend="none",**kwargs):
    #"ipyvtklink", "panel", "ipygany", "static", "pythreejs", "none"
    msh = u.function_space.mesh
    
    # Create plotter and pyvista grid
    plotter = pyvista.Plotter()

    topology, cell_types, geometry = plot.create_vtk_mesh(msh)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    # Attach vector values to grid and warp grid by vector
    values = np.zeros((geometry.shape[0], 3), dtype=np.float64)
    values[:, :len(u)] = u.x.array.real.reshape((geometry.shape[0], len(u)))
    grid["u"] = values
    warped_grid = grid.warp_by_vector("u", factor=factor)
    if cell_field is not None:
        warped_grid.cell_data[field_name] = cell_field.vector.array
        warped_grid.set_active_scalars(field_name)
    plotter.add_mesh(warped_grid,**kwargs)
    #plotter.show_axes()
    plotter.camera_position = 'xy'
    
    return plotter


if __name__ == "__main__":

    from dolfinx import fem, mesh
    from mpi4py import MPI


    # Create a mesh and a first function

    mesh = mesh.create_unit_square(
        MPI.COMM_WORLD, 8, 8, dolfinx.mesh.CellType.triangle
    )
    V_u = fem.VectorFunctionSpace(mesh, ("CG", 1))
    V_f = fem.FunctionSpace(mesh, ("CG", 1))
    u = fem.Function(V_u)
    u.interpolate(lambda x: [1 + x[0] ** 2, 2 * x[1] ** 2])

    #pyvista.set_jupyter_backend(backend)
    pyvista.OFF_SCREEN = True
    plotter = warp_plot_2d(u,factor=1,show_edges=True)
    if pyvista.OFF_SCREEN:
        pyvista.start_xvfb(wait=0.1)
        plotter.screenshot("deflection.png",transparent_background=True)
    else:
        plotter.show()
    
    x_s = np.linspace(0,1,100)
    y_s = 0.2 * np.ones_like(x_s)
    z_s = 0.0 * np.ones_like(x_s)
    points = np.array([x_s,y_s,z_s])
    values = evaluate_on_points(u,points)
    print(values)

