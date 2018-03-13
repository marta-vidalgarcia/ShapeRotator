# ShapeRotator: 
### An R tool for standardised rigid rotations of articulated Three-Dimensional structures with application for geometric morphometrics


The quantification of complex morphological patterns typically involves comprehensive shape and size analyses, usually obtained by gathering morphological data from all the structures that capture the phenotypic diversity of an organism or object. Articulated structures are a critical component of overall phenotypic diversity, but data gathered from these structures are difficult to incorporate in to modern analyses because of the complexities associated with jointly quantifying 3D shape in multiple structures. 
While there are existing methods for analysing shape variation in articulated structures in Two-Dimensional (2D) space, these methods do not work in 3D, a rapidly growing area of capability and research.  
Here we describe a simple geometric rigid rotation approach that removes the effect of random translation and rotation, enabling the morphological analysis of 3D articulated structures. Our method is based on Cartesian coordinates in 3D space so it can be applied to any morphometric problem that also uses 3D coordinates (e.g. spherical harmonics). We demonstrate the method by applying it to a landmark-based data set for analysing shape variation using geometric morphometrics. 
We have developed an R tool (ShapeRotator) so that the method can be easily implemented in the commonly used R package geomorph and MorphoJ software.  This method will be a valuable tool for 3D morphological analyses in articulated structures by allowing an exhaustive examination of shape and size diversity. 

### To install the ShapeRotator R-tool from Github using devtools:
install.packages("devtools")

devtools::install_github("marta-vidalgarcia/ShapeRotator")

### Translate
data.1_t  <- translate (T=data_1, landmark=landmarkA)
