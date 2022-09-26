# Homework (Devoir Maison)

## General information

- We will assign task to complete progressively during the course. 
- We will ask you to submit your work only at the end of the class.
- We will evaluate your work only at the end of the class.
- Submitting the homework is optional
- The final grade for the class will be calculated as follows

```
max(0.2 * homework_grade + 0.8 * final_exam_grade, final_exam_grade)
```

## First assigment (20/09/2022)

- Install FEniCSx on your computer following `INSTALL.md`
- Check that you can run without errors `00-Mesh/Mesh.ipynb` and `01-LinearElasticity/LinearElasticity.ipynb`
- Write the Boundary Value Problem corresponding to `01-LinearElasticity/LinearElasticity.ipynb`, including:
   - strong formulation with Dirichlet and Neumann boundary conditions
   - weak formulation 

## Second assigment (27/09/2022)

<img src="plateBCs.png"
     alt="Clamped plate with a crack"
     style="float: center; width:200px;" />

- Modify the boundary conditions to simulate the case of the plate in the figure above. We consider the case where the horizontal displacement is null at the top. You can use the symmetry as done in the previous example to model only half of the plate. This configuration corresponds to an experimental test that will consider later in this DM. 
 