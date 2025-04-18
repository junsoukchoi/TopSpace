% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_mcmc_TopSpace.R
\name{init.mcmc.TopSpace}
\alias{init.mcmc.TopSpace}
\title{Initialize the MCMC algorithm for TopSpace}
\usage{
init.mcmc.TopSpace(
  X,
  nbd.sizes,
  H,
  K,
  coord,
  kernel.params = c(0.01, 1),
  poly.degree = 20L
)
}
\arguments{
\item{X}{A list of vectors, each representing a local neighborhood and specifying the cell types observed within that neighborhood.}

\item{nbd.sizes}{A vector indicating the number of cells present in each local neighborhood.}

\item{H}{An integer indicating the total number of cell types observed in the data.}

\item{K}{An integer indicating the number of topics.}

\item{coord}{A matrix whose rows correspond to the spatial coordinates of each local neighborhood.}

\item{kernel.params}{A two-element vector specifying the hyperparameters of the modified exponential squared covariance kernel for the GP over spatially varying Dirichlet hyperparameters. The first element is the concentration parameter, where higher values yield greater concentration around the center. The second element is the smoothness parameter, where lower values result in a smoother process. The default values are (0.01, 1).}

\item{poly.degree}{An integer specifying the highest degree of Hermite polynomials that determines the dimension (length) of the truncated basis expansion for the GP over spatially varying Dirichlet hyperparameters. The default value is 10L.}
}
\value{
A list of initial values for the MCMC algorithm for TopSpace. Valid tags include 'Z', 'beta', 'alpha', 'B', 'n.eigen', 'lambda', 'Bases', and 'theta'. The value corresponding to each tag contains either the initial MCMC values or the eigendecomposition specification for the GP in TopSpace.
}
\description{
Initialize the MCMC algorithm for TopSpace
}
\examples{
# set random seed
set.seed(2024)

# load the pre-processed non-small cell lung cancer (NSCLC) dataset
data("NSCLC")

# select the image of the second patient
img = 2
X = NSCLC[[img]]$X
nbd_sizes = NSCLC[[img]]$nbd.sizes
coord = NSCLC[[img]]$coord

# specify the total number of cell types observed in the data: 
# CD19+ B-cells, CD14+ cells, CD8+ T-cells, CD4+ T-cells, CK+ cancer cells, and Others
H = 6

# specify the number of topics to be considered
K = 4
   
# obtain initial values for the MCMC algorithm for TopSpace
init = init.mcmc.TopSpace(X, nbd_sizes, H, K, coord, kernel.params =  c(0.01, 1))
print(init)
}
