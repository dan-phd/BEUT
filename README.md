# 2D BEUT Matlab program

## Introduction
The 2D Boundary Element Unstructured Transmission-line (BEUT) method is a hybrid electromagnetic simulation scheme which combines the Boundary Element Method (BEM) and Unstructured Transmission-Line Modelling (UTLM) method, to give an unconditionally stable time domain solver which has perfectly radiating boundaries.

Please note, this implementation in Matlab is intended for research/educational use, to show the concept, but at the expense of computational cost! I have created a more efficient, parallelised C++ program to compute the BEM operators which will run faster (speed depends on number of cores), useful if you decide to do any complex structures. This can be found at my [2DTDBEM] repository.

## Installation
To install, place the entire contents of *+BEUT* inside a Matlab search path.

## Contents
The *BEUT* project repository has the following directory structure:
* *BEUT* - license and readme files
 * *+BEUT* - project folder
  * *+BEM* - container for the BEM solver
   * *+Analytical* - analytical solutions for comparison
   * *+Demo* - examples demonstrating BEM components
   * *+Main* - main code for BEM test cases
  * *+Excitation* - container for the excitation classes
   * *+Demo* - examples demonstrating the excitation classes
   * *+Main* - main code for the BEUT test cases
  * *+Meshing* - container for the meshing classes
   * *+distmesh* - a third party mesher
   * *+Main* - main code for creating or loading meshes
   * *meshes* - contains .mat meshes ready to be used
    * *unconverted* - contins meshes to be converted
  * *+UTLM* - container for the UTLM solver
   * *@UTLMClass* - the UTLM class definition and its methods
   * *+Analytical* - analytical solutions for comparison
   * *+Main* - main code for UTLM test cases

The *+BEUT/+Main* folder contains examples on how to use BEUT. There are *+Main* folders elsewhere in the project to demonstrate the use of UTLM and BEM as individual solvers. The *+Meshing/+Main* folder is the first point of call for creating a mesh or loading a custom 2D mesh. For demonstration and testing of individual classes and functions, refer to the *+Demo* folders.
Generally, you will only want to open and run scripts in the *+Main* and *+Demo* folders.

If using the [2DTDBEM] C++ program to compute BEM operators, you will need to change the path which stores the resulting *.mat* files; This can done by modifying the path string in *+BEUT/CFolder.m*.

## Program flow
Once you are all set up, the general flow of a typical program is as follows:

1. Create/load mesh and save as a *UTLMClass* object (for UTLM) and a *MeshBoundary* object (for BEM)
3. Calculate BEM operators (using either the Matlab or C++ solver)
4. Set material parameters
5. Create excitation
5. Define probe locations
6. Roll through the timestepping loop (Marching-on-in-Time)
7. View results

## Credits
Some functions have not been written by me, but I would like to give thanks to the original authors. Please inspect these files/folders for affiliations and copyrights:
* *+distmesh*
* *read_wobj.m*
* *lgquad.m* (originally *lgwt.m*)


[2DTDBEM]: https://github.com/dan-phd/2DTDBEM