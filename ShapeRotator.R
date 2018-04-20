
######################################
### SHAPEROTATOR R TOOL ###
######################################


## DEBUG
immediate_on = TRUE

shaperotator3D_warning_on <- function ()
{
	immediate_on = TRUE
}


shaperotator3D_warning_on <- function ()
{
	immediate_on = FALSE
}



######################################
### INTERNAL FUNCTIONS ###
######################################

# Obtain the length of a vector u
norm_3D <- function (u)
{
  return( sqrt(u[1]^2 + u[2]^2 + u[3]^2));}

# Return the unit renormalisation of the vector from v
unit_3D <- function (v)
{
  return ( v/norm_3D(v) );
}

# Compute the dot product between vectors u and v
dot_3D <- function (u,v){
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
  R[3,] = c ( u[3]* u[1] * (1 - cos(t)) - u[2] * sin(t), u[3]* u[2]* (1 - cos(t)) + u[1] * sin(t), cos(t) + u[3]^2 * (1 - cos(t)) );
   return (R);
}

# Given a rotation matrix and a list of vectors, return the accompanying
# list of rotated vectors.

rotveclist_3D <- function (R, vlist)
{
  # Create an empty vector list for the rotation
  rvlist <-matrix(NA, nrow = dim(vlist)[1], ncol = 3) # Define size of the matrix to fill in the loop
   # Run through and perform the rotatio 
  for( i in 1: dim(vlist) [1] ){
    rvlist [i, ] = mmult_3D (R, vlist [i, ]);
  }
   return (rvlist);
}

# Compute the angle between [-\pi, pi) given ajacent and opposite sides
angle_tan <- function(adj, opp)
{
	# When the vector is at pi/2 or -pi/2, can't feed a vector that's zero though.
	if ( isTRUE (adj == 0) ) {
		if ( isTRUE (opp == 0) ) {
			angle <- NaN		# Zero vector, no angle.
		} else if ( isTRUE (opp < 0 ) ) {
			angle <- - pi/2		#
		} else {
			angle <- pi/2
		}

		return (angle)
	}

	# Outside of the pi/2 region.
	angle <- atan (opp/adj)

	return (angle)
}


# This returns the angle of the vector in the Y-Z plane. A rotation around the X axis of this angle
# brings the Z component to zero.
X_to_Z_rot_angle_3D = function (v)
{
	return (angle_tan (v[2],v[3]) )
}

# The angle in the X-Z plane. A rotation around the Y axis of this angle brings the Y component to zero.
Y_to_Z_rot_angle_3D = function (v)
{
	return (angle_tan (v[1],v[3]) )
}


# function vector.angle (vector from angle)
vector.angle <- function (angle)
{
	theta <- angle* pi/180;
	return( c(cos(theta), sin(theta), 0) ); # This simply produces a *unit* vector to the specified angle theta (in radians)
}


##############################
### FUNCTIONS SHAPEROTATOR ###
##############################

# function translate
# Takes an array T, an index in the array landmark, and translates it to the origin
# The default is to set the origin to (0,0,0), but this can be specified to be something else
translate <- function (T, landmark, origin = c(0,0,0))
{
  translist <- array(data = NA, dim = c(dim(T)[1], dim(T)[2], dim(T)[3]),dimnames = list(NULL, NULL, c(dimnames(T)[[3]])));
  for( i in 1:dim(T)[3] ){
	for (j in 1: dim(T)[1]) {
	    translist [j,,i] <- T [j,,i] - T[landmark,,i] + origin
	}
  }
  return (translist);
}


# Function that rotates an array, given a landmark and a vector (to which the rotation should take place).
# The resulting array will be such that the landmark arr[land,] will actually be on the same line, including orientation
# as the vector vec.
rotate.orientation <- function(arr, land, vec)
{
	# Create an array of the right dimensions
	rot <- matrix(NA, nrow = dim(arr)[1], ncol = 3)
	# Lets normalise and use unit vectors (vec  can be non-unital)
	uvec <- unit_3D(vec)
	# Rotation axis
	axis <- cross_3D(arr[land,], uvec)
	# ROtation matrix
	rotmat <- rotmat_3D(axis, angle_3D(arr[land,], uvec))
	# Rotate and return
	rot <- rotveclist_3D(rotmat, arr)

	# Check if rotmat[land,] is in the same orentation as vec.
	# Logic: svec which is the rotated vector should be the same as uvec if it is in the right orientation.
	# Unfortunately, due to numerical error, we can't just equality, which is equivalent to asking that ||svec - uvec|| is close to zero.
	# However, the numerically stable solution is to check that || svec - uvec || > || svec + uvec ||, which means that you have the wrong orientation.
	# In that case, we rotate again on the axis given by the variable "axis", but this time, we go around \pi radians.
	svec <- unit_3D(rot[land,])
	if ( isTRUE (norm_3D(svec - uvec) > norm_3D(svec + uvec)) ) {
		rotmat <- rotmat_3D(axis, pi)
		rot <- rotveclist_3D(rotmat, rot)
	}

	return(rot)
}



### SIMPLE.ROTATION
# function rotation for a single-point articulation (rotate translated datasets)
# Returns a list where $rotated1 gives the first object and $rotated2 gives the second
simple.rotation <- function(data.1, data.2, land.a, land.b, land.c, land.d, land.e, land.f, angle)
{

  # The first data set
  rot.array.1 <- array(data = NA, dim = c(dim(data.1)[1], dim(data.1)[2], dim(data.1)[3]),dimnames = list(NULL, NULL, c(dimnames(data.1)[[3]])));

  # This is the angle which land.b will face to land.e
  angle.v <-vector.angle(angle)


   #Rotating for data.1
  for( i in 1:(dim(data.1)[3]) ){
	# Untranslated specimen?
	if ( ! isTRUE(  all.equal((data.1[,,i])[land.a,] , c(0,0,0)) ) ) {
		warning (sprintf("Landmark A is not at the origin (data.1[,,%i])[%i,] = (%f,%f,%f)?", i, land.a, (data.1[,,i])[land.a,1], (data.1[,,i])[land.a,2], (data.1[,,i])[land.a,3]), immediate. = immediate_on)
	}

	# First, we rotate so that land.b rests on the axis (0,1,0).
	rot.array.1[,,i] <- rotate.orientation(data.1[,,i], land.b, c(0,1,0))

	#Sanity check: we need to have rot.array.1[,,i][land.b,] = (0, y, 0)
	if ( ! isTRUE(  all.equal( (rot.array.1[,,i])[land.b,1] , 0)) || ! isTRUE(  all.equal( (rot.array.1[,,i])[land.b,3] , 0)) ) {
		warning (sprintf("First rotation gone wrong in (rot.array.1[,,%i])[%i,] = (%f,%f,%f)?", i, land.b, (rot.array.1[,,i])[land.b,1], (rot.array.1[,,i])[land.b,2], (rot.array.1[,,i])[land.b,3]), immediate. = immediate_on)
	}


	# Rotate landmark C to the X-Y plane - we do so by computing the cross product between the unit.landmark.c.pos.proj vector and (-1,0,0),
	# which in the X-Y plane forms an angle given by the angle between these two vectors:
	unit.landmark.c.pos.proj = unit_3D( c( (rot.array.1[,,i])[land.c,1], 0, (rot.array.1[,,i])[land.c,3]))	# Unit vector with landmark.c.pos projected to the X-Z plane.
	rotmat <- rotmat_3D( cross_3D(unit.landmark.c.pos.proj, c(-1,0,0)), angle_3D(unit.landmark.c.pos.proj, c(-1,0,0)) )
	rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])


	# Let us now check that the landmark C is in the right spot. That is, we want C_x < B_x, if not, we rotate by angle pi.
	if ( isTRUE((rot.array.1[,,i])[land.b,1] < (rot.array.1[,,i])[land.c,1]) ) {
		rotmat <- rotmat_3D( c(0,1,0), pi)
		rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
	}

	# Sanity check:
	if ( ! isTRUE(  all.equal( (rot.array.1[,,i])[land.c,3] , 0)) || isTRUE((rot.array.1[,,i])[land.b,1] < (rot.array.1[,,i])[land.c,1])) {
		warning (sprintf("Second rotation gone wrong in (rot.array.1[,,%i])[%i,] = (%f,%f,%f)?", i, land.c, (rot.array.1[,,i])[land.c,1], (rot.array.1[,,i])[land.c,2], (rot.array.1[,,i])[land.c,3]), immediate. = immediate_on)
	}
# 	else  {
#		warning (sprintf("OK %i", i), immediate. = immediate_on)
#	}


	# Now, we rigidly move this object to sit at an angle given by "angle" to the X-axis.
	rot.array.1[,,i] <- rotate.orientation(rot.array.1[,,i], land.b, angle.v)
   }


  # The second data set
  rot.array.2 <- array(data = NA, dim = c(dim(data.2)[1], dim(data.2)[2], dim(data.2)[3]),dimnames = list(NULL, NULL, c(dimnames(data.2)[[3]])));

  for( i in 1:(dim(data.2)[3]) ){
	# Untranslated specimen?
	if ( ! isTRUE(  all.equal( (data.2[,,i])[land.d,] , c(0,0,0)) ) ) {
		warning (sprintf("Landmark D is not at the origin (data.2[,,%i])[%i,] = (%f,%f,%f)?", i, land.d, (data.2[,,i])[land.d,1], (data.2[,,i])[land.d,2], (data.2[,,i])[land.d,3]), immediate. = immediate_on)
	}

	# First, we rotate so that land.e rests on the axis (0,1,0).
	rot.array.2[,,i] <- rotate.orientation(data.2[,,i], land.e, c(0,1,0))


	#Sanity check: we need to have rot.array.1[,,i][land.e,] = (0, y, 0)
	if ( ! isTRUE(  all.equal( (rot.array.2[,,i])[land.e,1] , 0)) || ! isTRUE(  all.equal( (rot.array.2[,,i])[land.e,3] , 0)) ) {
		warning (sprintf("First rotation gone wrong in (rot.array.2[,,%i])[%i,] = (%f,%f,%f)?", i, land.e, (rot.array.2[,,i])[land.e,1], (rot.array.2[,,i])[land.e,2], (rot.array.2[,,i])[land.e,3]), immediate. = immediate_on)
	}


	# Rotate landmark F to the X-Y plane:
	unit.landmark.f.pos.proj = unit_3D( c( (rot.array.2[,,i])[land.f,1], 0, (rot.array.2[,,i])[land.f,3]))	# Unit vector with landmark.f.pos projected to the X-Z plane.
	rotmat <- rotmat_3D( cross_3D(unit.landmark.f.pos.proj, c(1,0,0)), angle_3D(unit.landmark.f.pos.proj, c(1,0,0)) )
	rot.array.2[,,i] <- rotveclist_3D ( rotmat, rot.array.2[,,i])


	# Let us now check that the landmark E is in the right spot. That is, we want E_x < F_x, if not, we rotate so that landmark F sits correctly on the X-Y plane.
	if ( isTRUE ((rot.array.2[,,i])[land.e,1] > (rot.array.2[,,i])[land.f,1])) {
		rotmat <- rotmat_3D( c(0,1,0), pi)
		rot.array.2[,,i] <- rotveclist_3D ( rotmat, rot.array.2[,,i])
	}

	# This is just a sanity check:
	if ( ! isTRUE(  all.equal( (rot.array.2[,,i])[land.f,3] , 0)) || isTRUE ((rot.array.2[,,i])[land.e,1] > (rot.array.2[,,i])[land.f,1]) ) {
		warning (sprintf("Second rotation gone wrong in (rot.array.2[,,%i])[%i,] = (%f,%f,%f)?", i, land.f, (rot.array.2[,,i])[land.f,1], (rot.array.2[,,i])[land.f,2], (rot.array.2[,,i])[land.f,3]), immediate. = immediate_on)
	}

	# Now, we rigidly move this object so that landmark E sits on the (1,0,0) axis.
	rot.array.2[,,i] <- rotate.orientation(rot.array.2[,,i], land.e, c(1,0,0))
   }

  # Join the two rotated datasets
  return ( list("rotated1" = rot.array.1, "rotated2" = rot.array.2))

}


### DOUBLE.ROTATION
# function rotation for a two-points articulation (2 different rotations for each translated dataset)
# Returns a list where $rotated1 gives the first object and $rotated2 gives the second
double.rotation <- function(data.1, data.2, land.a, land.b, land.c, land.d, land.e, land.f, land.g, land.h, angle)
{
	# First rotation.
	rot.array.1 <- array(data = NA, dim = c(dim(data.1)[1], dim(data.1)[2], dim(data.1)[3]),dimnames = list(NULL, NULL, c(dimnames(data.1)[[3]])));
	angle.v <-vector.angle(angle)


	for( i in 1:(dim(data.1)[3]) ){
		# Landmark A should be sitting on (0,0,0)

		# Untranslated specimen?
		if ( ! isTRUE(  all.equal((data.1[,,i])[land.a,] , c(0,0,0)) ) ) {
			warning (sprintf("Landmark A is not at the origin (data.1[,,%i])[%i,] = (%f,%f,%f)?", i, land.a, (data.1[,,i])[land.a,1], (data.1[,,i])[land.a,2], (data.1[,,i])[land.a,3]), immediate. = immediate_on)
		}

		# Rotate so that landmark C is sitting on the Z axis
		rot.array.1[,,i] <- rotate.orientation(data.1[,,i], land.c, c(0,0,1))


		#Sanity check: we need to have rot.array.1[,,i][land.c,] = (0, 0, y)
		if ( ! isTRUE(  all.equal( (rot.array.1[,,i])[land.c,1] , 0)) || ! isTRUE(  all.equal( (rot.array.1[,,i])[land.c,2] , 0)) ) {
			warning (sprintf("First rotation gone wrong in (rot.array.1[,,%i])[%i,] = (%f,%f,%f)?", i, land.c, (rot.array.1[,,i])[land.c,1], (rot.array.1[,,i])[land.c,2], (rot.array.1[,,i])[land.c,3]), immediate. = immediate_on)
		}

		# Now, we want to move landmark B to the X-Z axis.
		unit.landmark.b.pos.proj = unit_3D( c( (rot.array.1[,,i])[land.b,1], (rot.array.1[,,i])[land.b,2],0) )	# Unit vector with landmark.b.pos projected to the X-Y plane.
		rotmat <- rotmat_3D( cross_3D(unit.landmark.b.pos.proj, c(1,0,0)), angle_3D(unit.landmark.b.pos.proj, c(1,0,0)) )
		rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])

		# We want B_x > 0.
		if ( isTRUE((rot.array.1[,,i])[land.b,1] < 0) ) {
			rotmat <- rotmat_3D( c(0,0,1), pi)
			rot.array.1[,,1] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
		}


		# Let us now check that the landmark C is in the right spot. That is, we want B_y < D_y, if not, we rotate by angle pi in the X axis.
		if ( isTRUE((rot.array.1[,,i])[land.d,2] < (rot.array.1[,,i])[land.b,2]) ) {
			rotmat <- rotmat_3D( c(1,0,0), pi)
			rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
		}

		# Landmark B has to be 0 on the Y axis, and we want to make sure that the orientation is right
		if ( ! isTRUE(  all.equal( (rot.array.1[,,i])[land.b,2] , 0)) || isTRUE((rot.array.1[,,i])[land.d,2] < (rot.array.1[,,i])[land.b,2]) || isTRUE(rot.array.1[,,i][land.b,1] < 0) ) {
			warning (sprintf("Second rotation gone wrong in (rot.array.1[,,%i])[%i,] = (%f,%f,%f)?",
			i, land.b, (rot.array.1[,,i])[land.b,1], (rot.array.1[,,i])[land.b,2], (rot.array.1[,,i])[land.b,3]), immediate. = immediate_on)
		}

		# Rotate body to be the correct angle in the X-Y plane (so around the Z axis)
		rotmat <- rotmat_3D( c(0,0,1), angle * pi/180)
		rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
  	}


	# Establishing the rotation axes for data.2
	rot.array.2 <- array(data = NA, dim = c(dim(data.2)[1], dim(data.2)[2], dim(data.2)[3]),dimnames = list(NULL, NULL, c(dimnames(data.2)[[3]])));


	#Rotating for data.2 (ROTATION number 1)
	for( i in 1:(dim(data.2)[3]) ){
		# Untranslated specimen?
		if ( ! isTRUE(  all.equal((data.2[,,i])[land.e,] , c(0,0,0)) ) ) {
			warning (sprintf("Landmark E is not at the origin (data.2[,,%i])[%i,] = (%f,%f,%f)?", i, land.e, (data.2[,,i])[land.e,1], (data.2[,,i])[land.e,2], (data.2[,,i])[land.e,3]), immediate. = immediate_on)
		}

		# Rotate so that landmark F is sitting on the Z axis in the same direction as landmark C
		rot.array.2[,,i] <- rotate.orientation(data.2[,,i], land.f, unit_3D(rot.array.1[,,i][land.c,]) )

		#Sanity check: we need to have rot.array.1[,,i][land.c,] in the same direction as rot.array.2[,,i][land.f,]
		if ( ! isTRUE(  all.equal(unit_3D(rot.array.1[,,i][land.c,]), unit_3D(rot.array.2[,,i][land.f,]))) ) {
			warning (sprintf("Rotation of landmark F gone wrong?  (rot.array.1[,,%i])[%i,] = (%f,%f,%f), (rot.array.2[,,%i])[%i,] = (%f,%f,%f)",
				i, land.c, (rot.array.1[,,i])[land.c,1], (rot.array.1[,,i])[land.c,2], (rot.array.1[,,i])[land.c,3],
				i, land.f, (rot.array.2[,,i])[land.f,1], (rot.array.2[,,i])[land.f,2], (rot.array.2[,,i])[land.f,3]),
				immediate. = immediate_on)
		}


		# Now, we want to move landmark G to the X-Z axis.
		unit.landmark.g.pos.proj = unit_3D( c( (rot.array.2[,,i])[land.g,1], (rot.array.2[,,i])[land.g,2],0) )	# Unit vector with landmark.g.pos projected to the X-Y plane.
		rotmat <- rotmat_3D( cross_3D(unit.landmark.g.pos.proj, c(1,0,0)), angle_3D(unit.landmark.g.pos.proj, c(1,0,0)) )
		rot.array.2[,,i] <- rotveclist_3D ( rotmat, rot.array.2[,,i])

		# We want G_x > 0, because it should be in the same direction as B_x.
		if ( isTRUE((rot.array.2[,,i])[land.g,1] < 0) ) {
			rotmat <- rotmat_3D( c(0,0,1), pi)
			rot.array.2[,,1] <- rotveclist_3D ( rotmat, rot.array.2[,,i])
		}

		# This completely constrains the problem, we should now have G_y > H_y (because correcting for this would mean that landmark C ~ landmark F).
		# Nevertheless, we should check and issue a warning.

		if (  isTRUE(  rot.array.2[,,i][land.g,2] < rot.array.2[,,i][land.h,2]) ){
			warning (sprintf("Unsatisfied constraint in second rotation? (rot.array.2[,,%i])[%i,] = (%f,%f,%f), (rot.array.2[,,%i])[%i,] = (%f,%f,%f)",
				i, land.g, (rot.array.2[,,i])[land.g,1], (rot.array.2[,,i])[land.g,2], (rot.array.2[,,i])[land.g,3],
				i, land.h, (rot.array.2[,,i])[land.h,1], (rot.array.2[,,i])[land.h,2], (rot.array.2[,,i])[land.h,3]),
				immediate. = immediate_on)
		}
  	}

	return ( list("rotated1" = rot.array.1, "rotated2" = rot.array.2))

}



#####################
### PLOTTING CODE ###
#####################


#FUNCTION plot.rotation

plot.rotation.3D <- function (joined.data, data.1, data.2, specimen.num)
{ specimen <- joined.data[,,specimen.num]
colour <- c((rep.int(x=c(1), times=max(dim(data.1[,,1])))), (rep.int(x=c(2), times=max(dim(data.2[,,1])))))
lim=max(abs(specimen))
plot3d(specimen[,1], specimen[,2], specimen[,3], col = colour, xlim=c(-lim*1.1,lim*1.1), ylim=c(-lim*1.1, lim*1.1), zlim=c(-lim*1.1,lim*1.1));
}


# PLOT each dataset separately
plot.3D <- function (data, specimen.num, colour)
{ specimen <- data[,,specimen.num]
lim=max(abs(specimen))
plot3d(specimen[,1], specimen[,2], specimen[,3], col = colour, xlim=c(-lim*1.1,lim*1.1), ylim=c(-lim*1.1, lim*1.1), zlim=c(-lim*1.1,lim*1.1));
}


################################
### JOINING AND SORTING CODE ###
################################

# Join two arrays in a simple way: simply concatenate data.2[,,i] to data.1[,,i] and return a bigger array.
# Typically, this should be used after issuing a match.datasets
join.arrays <- function(data.1, data.2)
{
	if ( isTRUE( dim(data.1)[3] != dim(data.1)[3] ) ) {
		warning ( sprintf("Datasets have different number of samples: data.1 - %i and data.2 - %i, merging only the minimum of the two.", dim(data.1)[3], dim(data.2)[3]) )
	}

	new_array <- array(data = NA, dim = c(dim(data.1)[1]+dim(data.2)[1], dim(data.1)[2], min ( dim(data.1)[3], dim(data.2)[3]) ), dimnames = list(NULL, NULL, c(dimnames(data.1)[[3]])));

	for (i in 1:dim(new_array)[3]) {
		new_array[,,i] <- rbind(data.1[,,i], data.2[,,i])
	}

	return (new_array)
}

# Return two datasets, with data.1[i] matched with data.2[j] where j is given by the same species name as data.1[i].
# Returns a list with $matched1 and $matched2 as the parameters.
match.datasets <- function(data.1, data.2)
{

	if ( isTRUE( dim(data.1)[3] != dim(data.1)[3] ) ) {
		warning ( sprintf("Datasets have different number of samples: data.1 - %i and data.2 - %i, merging only the minimum of the two.", dim(data.1)[3], dim(data.2)[3]) )
	}

	# In case we have wrong lenghts, we still continue
	min.length = min( dim(data.1)[3], dim(data.2)[3] )

	dim.list.1 <-dimnames(data.1)[[3]]
	dim.list.2 <-dimnames(data.2)[[3]]

	m <- as.vector(match(dim.list.1, dim.list.2))

	new.data.1  <- array(data = NA, dim = c(dim(data.1)[1], dim(data.1)[2], min.length),dimnames = list(NULL, NULL, c(dimnames(data.1)[[3]])))
	new.data.2  <- array(data = NA, dim = c(dim(data.2)[1], dim(data.2)[2], min.length),dimnames = list(NULL, NULL, c(dimnames(data.2)[[3]])))

	for (i in 1:min.length) {
		if ( !isTRUE(is.na(m[i])) ) {
			# If we don't hit an NA in m[i], then we add the data, otherwise, we don't!
			new.data.1[,,i] <- data.1[,,i]
		 	new.data.2[,,i] <- data.2[,,m[i]]
    		}
	}


	# Delete the unmatched
	delete.NA <- c( which(is.na(m) == TRUE) )

	# The way we're assigning doesn't copy the dimnames over. We forcibly copy the dimnames of new.data.1 to new.data.2,
	# Because match.datasets should be an idempotent function.
	dimnames(new.data.2)[[3]] <- dimnames (new.data.1)[[3]]


	# If delete.NA is empty  then we return everything.
	if ( isTRUE (length(delete.NA) == 0 ) ) {
	 	return (list("matched1" = new.data.1[,,], "matched2" = new.data.2[,,]))
	} else {
	# delete.NA has some values in it, we have to delete those rows.
		return (list("matched1" = new.data.1[,,-delete.NA], "matched2" = new.data.2[,,-delete.NA]))
	}
}


## THese should be obselete.
# Join two datasets
join.datasets <- function(data.1, data.2)
{

  new_array <- join.data.noname (data.1, data.2)


  dim.list.1 <-dimnames(data.1)[[3]]
  dim.list.2 <-dimnames(data.2)[[3]]

  m <- as.vector(match(dim.list.1, dim.list.2))
  n <- c(1:length(dim.list.1))

  delete.NA <- c(which(is.na(m)==TRUE))
  return (new_array[,,-delete.NA])
}

# Join data
join.data.noname <- function(data.1, data.2)
{
  dim.list.1 <-dimnames(data.1)[[3]]
  dim.list.2 <-dimnames(data.2)[[3]]

  m <- as.vector(match(dim.list.1, dim.list.2))
  n <- c(1:length(dim.list.1))

  new_array <- array(data = NA, dim = c(dim(data.1)[1]+dim(data.2)[1], dim(data.1)[2], dim(data.1)[3]), dimnames = list(NULL, NULL, c(dimnames(data.1)[[3]])));

  for (i in 1:length(dim.list.1)){
    if ( !isTRUE(is.na(m[i])) ) {
      new_array[,,i] <- rbind(data.1[,,i], data.2[,,m[i]])
    }
  }

  return (new_array)
}


######################################
### DISTORTION AND SANITY CHECKING ###
######################################

# Expects translated data tdata and a rotated array rdata
# First entry: length of tdata points from the origin
# Second entry: length of rdata points from the origin
# Third entry: difference of the first and second
# Forth entry: testing for whehther this is close to zero (upto precision error)
compute.distortion <- function(tdata, rdata)
{

	dist <- array(data = NA, dim = c(dim(tdata)[1], 4, dim(tdata)[3]),dimnames = list(NULL, NULL, c(dimnames(tdata)[[3]])));

	# If the dimensions aren't the same, we just return the NA's.
	if ( isTRUE (dim(tdata) != dim(rdata) ) ) {
		return (dist)
	}


  	for( i in 1:dim(tdata)[3] ){
		for (j in 1: dim(tdata)[1]) {
		    	dist[j,1,i] <- norm_3D(tdata [j,,i])
		   	dist[j,2,i] <- norm_3D(rdata[j,,i])
		   	dist[j,3,i] <- abs(norm_3D(tdata [j,,i])  - norm_3D(rdata[j,,i]))
			dist[j,4,i] <- isTRUE(all.equal( dist[j,,i][3] , 0))

			# Nonzero distortion?
			if ( ! dist[j,,i][4] ) {
				warning (sprintf("Possible distortion in sample %i at landmark %i.", i, j), immediate. = immediate_on)
			}
		}
  	}

	return (dist)
}

