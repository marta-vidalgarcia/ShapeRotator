# ShapeRotator 
## An R tool for standardised rigid rotations of articulated Three-Dimensional structures with application for geometric morphometrics
### Welcome to the ShapeRotator Wiki!
#### The quantification of complex morphological patterns typically involves comprehensive shape and size analyses, usually obtained by gathering morphological data from all the structures that capture the phenotypic diversity of an organism or object. Articulated structures are a critical component of overall phenotypic diversity, but data gathered from these structures are difficult to incorporate in to modern analyses because of the complexities associated with jointly quantifying 3D shape in multiple structures. While there are existing methods for analysing shape variation in articulated structures in Two-Dimensional (2D) space, these methods do not work in 3D, a rapidly growing area of capability and research. Here we describe a simple geometric rigid rotation approach that removes the effect of random translation and rotation, enabling the morphological analysis of 3D articulated structures. Our method is based on Cartesian coordinates in 3D space so it can be applied to any morphometrics problem that also uses 3D coordinates (e.g. spherical harmonics). We demonstrate the method by applying it to a landmark-based data set for analysing shape variation using geometric morphometrics. We have developed an R tool (ShapeRotator) so that the method can be easily implemented in the commonly used R package geomorph and MorphoJ software. This method will be a valuable tool for 3D morphological analyses in articulated structures by allowing an exhaustive examination of shape and size diversity.

`library(devtools)`

`install_github("marta-vidalgarcia/ShapeRotator")`
`library(ShapeRotator)`

***

Please check out the Wiki to know more about the functions available within the ShapeRotator R tool and the basic steps required in order to successfully implement the rotation on a dataset of 3D coordinates. ShapeRotator allows the rigid rotation of sets of both landmarks and semilandmarks used in geometric morphometric analyses, enabling morphometric analyses of complex objects, articulated structures, or multiple parts within an object or specimen.



***

If you have any questions you can contact Marta Vidal-García at **marta.vidalga@gmail.com**



***

Citation:
M. Vidal-García, L. Bandara and J.S. Keogh. (2018) ShapeRotator: An R tool for standardized rigid rotations of articulated three-dimensional structures with application for geometric morphometrics. Ecology and Evolution. 8:4669–4675.
DOI: 10.1002/ece3.4018
