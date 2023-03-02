# Tropical-Balls
This repository provides functions and vignettes that may be used to construct maximum inscribed and minimum enclosing tropical balls with applications.  The measurements that the user obtains are based on a novel algorithm called the tropical Markov Chain Monte Carlo Hit-and-Run (HAR) sampler.  This sampler samples a tropical polytopes full-dimensional element uniformly.  It may also be modified to sample in a similar way as a HAR sampler which samples according to a normal distribution.

-Tropical Polytopes
Tropical polytopes, in all cases are represented as a matrix where a tropical point is a row in the matrix.  It is assumed that the first coordinate of any point is zero.

-Normalizing points and polytopes
There are several helper functions in the 'THAR_Funcs.R' script that will allow the user to normalize a tropical point or polytope.  Normalizing means subtracting the scalar value in the first element of any point from every element of that point.  This is because a point in the tropical projective torus that is modified by adding a scalar to each element of a point is the same point.

-Minimum Enclosing Tropical Ball
The minimum enclosing tropical ball is the smallest ball based on the tropical metric that encloses a polytope.

-Maximum Inscribed Tropical Ball
The maximum inscribed tropical ball is the tropical ball of largest radius defined by the tropica metric, fully enclosed inside a tropical polytope.

-Volume of a Tropical Polytope
The volume of a tropical polytope is calculated by sampling from an enclosing tropical polytope of known volume.  In this case the enclosing tropical polytope is a minimum enclosing ball.  By sampling from the minimum enclosing ball we can determine the proportion of points that also fall inside the tropical polytope of interest.  Multiplying this proportion by the volume of the minimum enclosing ball will give us the volume estimate.  Note that our volume estimates depend significantly on the volume of the minimum enclosing ball as well as the sample size.

-Uniform Sampling of a Tropical Polytope
A tropical simplex is a tropical polytope defined on e vertices where the tropical polytope is in R^e/R1, the tropical projective torus.  The tropical HAR sampler, samples the full-dimensional part, or e-trunk, of any tropical simplex uniformly.  If however, we apply it to a tropical polytope that has more than e vertices we cannot apply the HAR sampler directly.  Instead we sample each simplex defined by e vertices of the original tropical polytope.  However, the first step is to sample from the polytope's minimum enclosing ball and determine the number of points that fall inside of the tropical polytope as well as well as the points that fall inside of each identifiable tropical simplex.  We then identify the subset of tropical simplices that account for all sample points.  These simplices are then sampled according to the proportion of points accounted for by each tropical simplex.

-Rounding a Tropical Poltyope
By rounding, we mean removing all elements, or tentacles, of a tropical polytope that are of dimension e-2 or less.  The reason for this is to be able to sample from only the e-1 trunk of the tropical polytope since only the e-1 trunk contributes to the volume.  This is accomplished through use of Fukuda's double description algorithm, cddlib, and its R implementation, rcdd.

If any bugs are noted in the code, please contact Dave Barnhill and david.barnhill@nps.edu.
