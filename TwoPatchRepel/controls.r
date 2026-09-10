# Set time #and space discretisation
# Duration of run max_t
# Number of runs N
# Random seed (may be NULL - then not reset)
# Initial x as x0 or (if NULL) uniform on xmin,xmax
# Initial s as s0 or (if NULL) uniform on states
controls <- list(
  dt=0.1,
  max_t=10^4,
  xmin=-50,
  xmax=50,
  dx=0.01,
  delta=10^-3,
  x0=NULL,
  s0=NULL,
  seed=0,
  N=10
)
