# How to install FEniCS

## On Linux and Macosx using Anaconda

Instructions on how to install FEniCS-dolfinx are available https://fenicsproject.org/download/.

We suggest and support the installation method using anaconda. Anaconda is a useful package manager for python.

The following install instruction work only on Linux on Macos. If you have Windows, see below.

1. Install Anaconda from https://www.anaconda.com/products/distribution.
    - If you are on Macosx: choose the install form your platform:
        - Intel processors: https://repo.anaconda.com/archive/Anaconda3-2022.05-MacOSX-x86_64.pkg
        - M1 processors: https://repo.anaconda.com/archive/Anaconda3-2022.05-MacOSX-arm64.pkg

2. Open a new terminal and go to the directory containing the file `fenicsx-0.5.1.yaml`. You will find this file in the present git repository.

3. You should be now in the 'base' environment and your command prompt should show '(base)'.
To be sure to use updated version of the package and avoid further conflicts, let us update the `base` environment with the following command:
```
conda update -n base -c defaults conda
```

4. Create a new conda environment from the file `fenicsx-0.5.1.yaml`
    ```
    conda env create --file fenicsx-0.5.1.yml
    ```

5. You have now installed fenics in the conda environment `fenicsx-0.5.1`. To use it you must activate the environment with the following command
    ```
    conda activate fenicsx-0.5.1
    ```

After the first installation, you need only step 5 above (`conda activate fenicsx-0.5.1`) to use fenicsx on a terminal.

You can find further help on conda [here](https://docs.conda.io/projects/conda/en/latest/_downloads/843d9e0198f2a193a3484886fa28163c/conda-cheatsheet.pdf)

### Troubleshooting
- If you get `CondaValueError: prefix already exists: /Users/maurini/opt/anaconda3/envs/fenicsx-0.5.1` try with
    ```
    conda env create --file fenicsx-0.5.1.yml --force
    ```

- If you do have issues with pyvista try setting `pyvista.set_jupyter_backend("none")`. If this does not work, just do not use pyvista: comment the related lines of code, save to xdmf format and visualize results in paraview.

## On Windows

FEniCS is not distributed for Windows boxes. For Windows 10, the preferred option is the [Windows subsystem for linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install-win10).
Install the Ubuntu distribution as WSL, then refer to the section above inside the Ubuntu WSL

If you cannot install the WSL, you will need to create Ubuntu virtual machine. Get in touch with your instructors in that case.
Then inside the ubuntu virtual machine follows the instruction for conda installation


## Cloud-based  Google colab installations
You can run python programs jupyter notebooks and FEniCS on online servers. The basic service is free and can be a solution if all other installation systems fail.

* You need a google account
* You can use this example to start https://colab.research.google.com/github/fem-on-colab/fem-on-colab/blob/main/fenicsx/test-dolfinx.ipynb#scrollTo=infectious-train
* You can save the notebooks and your working environment on your google drive

We suggest to use this only as emergency solution
* *advantages:* You do not need to install anything on your machine and you do note use the resources of your machine

* *disadvantages:* It can be slow. You can use only jupyter notebooks, you share your data with google or microsoft, you do not a full control of the system, you need to be online with a good network connection.

## How to test the installation

Run following commands in your terminal:

```
wget https://fenicsproject.org/docs/dolfin/latest/python/_downloads/598330b504d63e359baad030e1010987/demo_poisson.py
python3 demo_poisson.py
```
or copy and paste this code in the notebook and exectute the cell

```
import numpy as np
import ufl
from dolfinx import fem, io, mesh, plot
from ufl import ds, dx, grad, inner
from mpi4py import MPI
from petsc4py.PETSc import ScalarType
msh = mesh.create_rectangle(comm=MPI.COMM_WORLD,
                            points=((0.0, 0.0), (2.0, 1.0)), n=(32, 16),
                            cell_type=mesh.CellType.triangle,)
V = fem.FunctionSpace(msh, ("Lagrange", 1))
facets = mesh.locate_entities_boundary(msh, dim=1,
                                       marker=lambda x: np.logical_or(np.isclose(x[0], 0.0),
                                                                      np.isclose(x[0], 2.0)))
dofs = fem.locate_dofs_topological(V=V, entity_dim=1, entities=facets)
bc = fem.dirichletbc(value=ScalarType(0), dofs=dofs, V=V)
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
x = ufl.SpatialCoordinate(msh)
f = 10 * ufl.exp(-((x[0] - 0.5) ** 2 + (x[1] - 0.5) ** 2) / 0.02)
g = ufl.sin(5 * x[0])
a = inner(grad(u), grad(v)) * dx
L = inner(f, v) * dx + inner(g, v) * ds
problem = fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()
with io.XDMFFile(msh.comm, "out_poisson/poisson.xdmf", "w") as file:
    file.write_mesh(msh)
    file.write_function(uh)
```
