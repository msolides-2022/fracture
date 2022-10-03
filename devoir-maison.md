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
- Check that you can run without errors `00-Mesh.ipynb` and `01-LinearElasticity.ipynb`
- Write the Boundary Value Problem corresponding to `01-LinearElasticity.ipynb`, including:
   - strong formulation with Dirichlet and Neumann boundary conditions
   - weak formulation 

## Second assigment (27/09/2022)

<img src="figures/plateBCs.png"
     alt="Clamped plate with a crack"
     style="float: center; width:200px;" />

- Modify the boundary conditions to simulate the case of the plate in the figure above. We consider the case where the horizontal displacement is null at the top. You can use the symmetry as done in the previous example to model only half of the plate. This configuration corresponds to an experimental test that will consider later in this DM. 
 
 ## Third assigment (04/10/2022)

- Using dimensional analysis, show that $K_{I}=\frac{E \Delta}{\sqrt{l}}K^{*}(a/l)$
- Calculate $K_I$ for a PMMA plate (E = 3000 MPa, $\nu$ = 0.4) of sizes l = 10 cm,  2h = 45 cm, thickness e =2 mm, one value of $a \in [0.5-3]$ cm, loading $\Delta = 1$ mm

