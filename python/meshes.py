import gmsh
import numpy as np
from mpi4py import MPI
from dolfinx.io import gmshio, XDMFFile


def generate_mesh_with_crack(
    Lx=1, Ly=1, Lcrack=0.3, lc=0.015, refinement_ratio=10, dist_min=0.05, dist_max=0.2
):
    # For further documentation see
    # - gmsh tutorials, e.g. see https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorial/python/t10.py
    # - dolfinx-gmsh interface https://github.com/FEniCS/dolfinx/blob/master/python/demo/gmsh/demo_gmsh.py
    #
    gmsh.initialize()
    gdim = 2
    mesh_comm = MPI.COMM_WORLD
    model_rank = 0
    
    if mesh_comm.rank == model_rank:
        model = gmsh.model()
        model.add("Rectangle")
        model.setCurrent("Rectangle")
        p1 = model.geo.addPoint(0.0, 0.0, 0, lc)
        p2 = model.geo.addPoint(Lcrack, 0.0, 0, lc)
        p3 = model.geo.addPoint(Lx, 0, 0, lc)
        p4 = model.geo.addPoint(Lx, Ly, 0, lc)
        p5 = model.geo.addPoint(0, Ly, 0, lc)
        l1 = model.geo.addLine(p1, p2)
        l2 = model.geo.addLine(p2, p3)
        l3 = model.geo.addLine(p3, p4)
        l4 = model.geo.addLine(p4, p5)
        l5 = model.geo.addLine(p5, p1)
        cloop1 = model.geo.addCurveLoop([l1, l2, l3, l4, l5])
        surface_1 = model.geo.addPlaneSurface([cloop1])

        model.mesh.field.add("Distance", 1)
        model.mesh.field.setNumbers(1, "NodesList", [p2])
        # model.mesh.field.setNumber(1, "NNodesByEdge", 100)
        # model.mesh.field.setNumbers(1, "EdgesList", [2])
        #
        # SizeMax -                     /------------------
        #                              /
        #                             /
        #                            /
        # SizeMin -o----------------/
        #          |                |    |
        #        Point         DistMin  DistMax

        model.mesh.field.add("Threshold", 2)
        model.mesh.field.setNumber(2, "IField", 1)
        model.mesh.field.setNumber(2, "LcMin", lc / refinement_ratio)
        model.mesh.field.setNumber(2, "LcMax", lc)
        model.mesh.field.setNumber(2, "DistMin", dist_min)
        model.mesh.field.setNumber(2, "DistMax", dist_max)
        model.mesh.field.setAsBackgroundMesh(2)

        model.geo.synchronize()
        surface_entities = [model[1] for model in model.getEntities(2)]
        model.addPhysicalGroup(2, surface_entities, tag=5)
        model.setPhysicalName(2, 2, "Rectangle surface")
        model.mesh.generate(gdim)
    
    msh, mt, ft = gmshio.model_to_mesh(model, mesh_comm, model_rank,gdim=gdim)
    msh.name = "crack"
    mt.name = f"{msh.name}_cells"
    ft.name = f"{msh.name}_facets"
    return msh, mt, ft


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    from pathlib import Path

    Lcrack = 0.5
    dist_min = 0.1
    dist_max = 0.3
    msh, mt, ft = generate_mesh_with_crack(
        Lcrack=Lcrack,
        Ly=0.5,
        lc=0.1,  # caracteristic length of the mesh
        refinement_ratio=10,  # how much it is refined at the tip zone
        dist_min=dist_min,  # radius of tip zone
        dist_max=dist_max,  # radius of the transition zone
    )

    with XDMFFile(msh.comm, "out_gmsh/mesh.xdmf", "w") as file:
        file.write_mesh(msh)
        msh.topology.create_connectivity(1, 2)
        file.write_meshtags(mt, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry")
        file.write_meshtags(ft, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry")
