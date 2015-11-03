2D BEUT Matlab program

The 2D Boundary Element Unstructured Transission-line (BEUT) method is a hybrid electromagnetic simulation scheme which combines the Boundary Element Method (BEM) and Unstructured Transmission-Line Modelling (UTLM) method, to give an unconditionally stable time domain solver which has perfect radiating boundaries.

Please note, this implementation in Matlab is intended for research/educational use, to show the concept, but at the expense of computational cost! I have created a more efficient, parallelised C++ program to compute the BEM operators which will run faster (speed depends on number of cores), useful if you decide to do any complex structures.



The general program flow is as follows:

1. Create/load mesh

2. Convert mesh to "UTLMClass" and "MeshBoundary" objects

3. Calculate BEM operators (using "MeshBoundary" object in the Matlab or C++ solver)

4. Set material parameters

5. Create excitation

5. Define probe locations

6. Timestepping loop (Marchin-on-in-Time)

7. View results



To install, place the entire contents of +BEUT inside a folder which can be seen in the Matlab search path.

There are 5 subfolders:
+BEM
+Excitation
+Main
+Meshing
+UTLM

The +Main folder contains examples on how to use BEUT. There are also +Main folders in the other folders to demonstrate the use of UTLM and BEM as indivdual solvers. The +Meshing/+Main folder is the first point of call for creating a mesh or loading a custom 2D mesh. For demonstration and testing of individual classes and functions, refer to the +Demo folders.

Generally, you will only want to open and run scripts (in +Main or +Demo folders), for which the names begin with a capital letter. Classes also start with a capital letter, but Matlab helpfully changes the icon for these files. Functions begin with a small letter.

If using the C++ program to compute BEM operators, you will need to change the path which stores the resulting .mat files. This can done in +BEUT/CFolder.


Some functions have not been written by me, but I appreciate the original authors for their help! Please inspect these files/folders for affiliations and copyrights:
+distmesh
read_wobj
lgquad (originally lgwt)
