### SHAPEROTATOR R package

library(devtools)
install_github("marta-vidalgarcia/ShapeRotator")

### INTERNAL FUNCTIONS


# Obtain the length of a vector u
norm_3D <- function (u) 
{
  return( sqrt(u[1]^2 + u[2]^2 + u[3]^2));
}

# Return the unit renormalisation of the vector from v
unit_3D <- function (v)
{ 
  return ( v/norm_3D(v) );
}

# Compute the dot product between vectors u and v
dot_3D <- function (u,v)
{
  d = u[1] * v[1] + u[2] * v[2] + u[3] * v[3];
  return (d); 
}

# Take the cross product between the vectors u and v
cross_3D <- function (u,v)
{
  w = c( u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1] );
  return (w);
}

# Compute the angle between two vectors
angle_3D <- function (u,v)
{
  a = acos( dot_3D(u,v)/ (norm_3D(u) * norm_3D(v)) );
  return (a);
}

# Return M * v, the action of the matrix M on a vector v.
mmult_3D <- function (M, v)
{
  w = c( M[1,1] * v[1] + M [1, 2] * v[2] + M [1,3] * v[3], 
         M [2,1] * v[1] + M [2,2] * v[2] + M[2,3] * v[3],
         M [3,1] * v[1] + M [3,2] * v[2] + M[3,3] * v[3]);
  return (w); 
}


# Return a rotation matrix given an axis v and an angle t:
rotmat_3D <- function (v, t) 
{ 
  # Create a 3 x 3 matrix with 0's.
  R = matrix (0, 3, 3);
  
  # normalise v to u, because the rotation matrix needs a unit vector
  u = unit_3D (v); 
  
  # The first entry of the rotation matrix is: 
  R[1,] = c (cos(t) + u[1]^2 *(1 - cos(t)), u[1]* u[2] * (1 - cos(t)) - u[3] * sin(t), u[1]* u[3]* (1 - cos(t)) + u[2] *sin(t) );
  R[2,] = c (u[2]* u[1] * (1- cos(t)) + u[3] * sin(t), cos(t) + u[2]^2 *(1 - cos(t)), u[2]*u[3] * (1 - cos(t)) - u[1] *sin(t) );
  R[3,] = c ( u[3]* u[1] * (1 - cos(t)) -  u[2] * sin(t), u[3]* u[2]* (1 - cos(t)) + u[1] * sin(t), cos(t) + u[3]^2 * (1 - cos(t)) );
  
  return (R);
}

# Given a rotation matrix and a list of vectors, return the accompanying
# list of rotated vectors.

rotveclist_3D <- function (R, vlist)
{
  # Create an empty vector list for the rotation
  rvlist <-matrix(NA, nrow = dim(vlist)[1], ncol = 3) # Define size of the matrix to fill in the loop
  
  # Run through and perform the rotatio 
  for( i in 1: dim(vlist) [1]  ){
    rvlist [i, ] = mmult_3D (R, vlist [i, ]);
  }
  
  return (rvlist);
}



X_to_Z_rot_angle_3D = function (v) 
{ 
  return (atan (-v[3]/v[2]) ); 
}

Y_to_Z_rot_angle_3D = function (v) 
{ 
  return (atan (v[3]/v[1]) ); 
}


# function vector.angle (vector from angle)
vector.angle <- function (angle)
{
  tangent <- tan(angle*pi/180);
  if (angle == 90) {
    vector <- c(0, 1, 0)
  } else if (angle == 0) {
    vector <- c(1, 0, 0)
  } else if (angle == 180) {
    vector <- c(-1, 0, 0)
  } else if ((angle < 90) && (angle > 0)) {
    vector <- c(1, tangent, 0)
  } else if ((90 < angle) && (angle < 180)) {
    vector <- c(-1, -tangent, 0)
  } else if ((180 < angle) && (angle < 270)) {
    vector <- c(-1, -tangent, 0)
  } else if ((270 < angle) && (angle < 360)) {
    vector <- c(1, tangent, 0)
  } else if (angle == 360) {
    vector <- c(1, 0, 0)
  } else
    vector <- NaN
  return (vector);
}


######
