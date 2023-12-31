% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rast_to_df.R
\name{rast_to_df}
\alias{rast_to_df}
\title{rast_to_df}
\usage{
rast_to_df(x, time_dim)
}
\arguments{
\item{x}{Rast stack resulted from detect3event and re-organized until 'nc_no_NA' (see detect3 function)}

\item{time_dim}{The time dimension. It should be a single value to use as the origin for the date indexes}
}
\value{
This function will convert a raster object to a dataframe.
}
\description{
Function for converting a raster to a dataframe.
}
\keyword{internal}
