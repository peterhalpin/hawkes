#-----------------------------------------------------------------------------
# Last update 05/09/2017
# Functions for estimation of Hawkes process

# To Do:
	# 1) Make set_kernel return an object with attributes governing the kernel, its arguments, its box contstaints.

	# 2) Sort out one place to put "shit lists." I want one function that checks for all for all bad args, and if that function has been called already, then it isnt called again.

	# 3) Time series data object? Instead of passing various combinations of out inp, time, and events to each function, just pass the TS data object and pull the necessary args from it.

	#4) Sort out non-montonicity of the EM. Lapham reported a similar phenom when using n_k to approximate \sum_k F(T-t_k) in the compensator. So it appears that the problem can be caused by this kind of approximation error. One possible sources of this problem is the rounding error imposed in compute_I.

	# Another plausible source is that the posterior probabilties are wrong in Halpin & De Boeck 2013. Here are three versions of the likelihood in eqation A(6).

	# As stated on .p23 of Daley & Vera-Jones, the likelihood of singe point from a Poisson process on the interval (0, T] is just \lambda(t_k) exp{- \Lambda(T)}. This is th probability of a single point at t_k and no points in the rest of the interval.

	# If the interval under consideration is just t_k, then the exp term is equal to 1 and we get just the CIF. This is probably the closest thing to the probability of the point t_k.

  # The conditional probabilties, \lambda(t_k) / \Lambda(T) the probability of a point at t_k, conditional on the number of points in the interval [0, T].

	#Each of these leads to a different expresion for the posterior distribution

	#But I used \lambda(T_k) / \Lambda(T) for some goddam reason. So the posterior probability of a point belong to a given parent is

	#\lambda_m(t_k) exp{- \Lambda_m(T)} * \Lambda_m(T) / sum_r \Lambda_r(T)

	# all over the sum over m of that, or

	# \lambda_m(t_k) exp{- \Lambda_m(T)} * \Lambda_m(T)
  # --------------------------------------------------------
	# \sum_m \lambda_m(t_k) exp{- \Lambda_m(T)} * \Lambda_m(T)

	# Where the sum over r drops out because it is a constant in all terms. Noting that \Lambda_m(\infty) = 1, the term  exp{- \Lambda_m(T)} * \Lambda_m(T) might be approximated as 1/e, which would leave

	# \lambda_m(t_k)
	# ---------------------
	# \sum_m \lambda_m(t_k)

	# This last is the quantity used in Halpin & De Boeck 2013. Since the rounding error reported Lapham was due to replacing \Lambda(T) with \Lamda(\infty) it seems that this "approximation" may be the source of the current problems. Yet Veen & Schoenberg use the last expression, and so does Lapham!
#-----------------------------------------------------------------------------

require(dplyr)
# devtools::document("R")

#-----------------------------------------------------------------------------
#' Creates a point process object (S3).
#'
#'  A \code{pp} object is just a list whose components are the event times of each margin of a point process, plus some additional attributes. Role is analgous to \code{stats::ts}.
#'
#' @param list_times a list, each component containing the event times for a margin of the point process.
#' @param n_margins if \code{list_times} is missing, used to create an empty \code{pp} object.
#' @param end_time if \code{list_times}  is missing, used to create an empty \code{pp} object.
#' @return A \code{pp} object (S3).
#' @export

pp <- function(list_times, n_margins = length(list_times), end_time = max(unlist(list_times))) {

  if (missing(end_time)) end_time <- 0
	if (missing(list_times)) {
	if (missing(n_margins)) n_margins <- 1
		list_times <- vector("list", n_margins)
	}

  if (is.null(names(list_times))) {
    names(list_times) <- paste0("X", 1:length(list_times))
  }

	events <- structure(list_times, class = "pp")
	attr(events, "n_margins") <- length(list_times)
	attr(events, "n_events") <- unlist(lapply(list_times, length))
	attr(events, "end_time") <- lapply(list_times, function(x) max(x, 0)) %>% unlist %>% {max(c(.,end_time))}
	events
}


#-----------------------------------------------------------------------------
#' Print function for \code{pp} objects.
#'
#' @param pp_obj a \code{pp} object.
#' @export

print.pp <- function(pp_obj) {
	cat("Point Process: \n")
	cat("n_margins =", attr(pp_obj, "n_margins"), "\n")
	cat("n_events per margin =", attr(pp_obj, "n_events"), "\n")
	cat("Time interval = [0, ", attr(pp_obj, "end_time"), "] \n", sep = "")
  for (i in 1:attr(pp_obj, "n_margins")) {
    print(pp_obj[i])
    #cat(names(pp)[i], "=", round(pp[[i]][1:4], 3),
    #  "...", round(pp[[i]][attr(pp, "n_events")[i]], 3), "\n" )
  }
}


#-----------------------------------------------------------------------------
#' Default plot method for \code{pp} objects.
#'
#'Cumulative event count plot function. To do: add waiting times plots.
#'
#' @param pp_obj a \code{pp} object.
#' @return A plot.
#' @export

plot.pp <-function(pp_obj) {
  m <- attr(pp_obj, "n_margins")
	n <- attr(pp_obj, "n_events")

  plot(pp_obj[[1]], 1:n[1], type = "n", main = "Number of Events by Time",
    ylab = "Events", xlab = "Time", ylim = c(0, max(n)),
    xlim = c(0, max(unlist(pp_obj))), las = 1)

  axis(side = 4, las = 1)
  par(mar = c(5,4,4,4) + 1)
  #mtext("Events", side = 4, line = 3, cex = par("cex.lab"))
  legend(x = "topleft", legend = attr(pp_obj, "names"),
    col = c(3:(m+2)), lty = 1, lwd = 2, border = "n")

  for (i in 1:m) {
		points(pp_obj[[i]], 1:n[i], type = "l", col = i+2, lwd = 2)
	}

  # Under construction: Event times
  # plot(pp[[1]], rep(1, length(pp[[1]])), type = "n", main = "Event Times", ylab = "", xlab = "Time", pch = 124, yaxt = "n", ylim = c(-.5, -(m+.5)))
  #
  # for (i in 1:m) {
  #   points(pp[[i]], rep(-i, length(pp[[i]])), pch = 124, col = i+2)
  # }
}


#-----------------------------------------------------------------------------
#' Computes residual waiting times of a Hawkes process.
#'
#' Residuals are computed using the time change theorem (see e.g., Daley and Verrra-Jones, 2013, sec. 7.4). If \code{parms} is missing, computes residuals for homogeneous Poisson process. To do: Re-write this to work as a default \code{resid} method on \code{EM} objects.

#' @param pp_obj a \code{pp} object.
#' @param parms parameters of the Hawkes process formatted as described for \code{set_parms}. If missing, output is for homogeneous Poisson process.
#' @param kernel_type currently only implemented for \code{"dgamma"}.
#' @return A list of length \code{attr(pp_obj, "n_margins")}, each component a sorted vector of residual interevent times for the corresponding margin.
#' @export

res <- function(pp_obj, parms, kernel_type = "dgamma") {

  k <- set_kernel(kernel_type)
  m <- attr(pp_obj, "n_margins")
  n <- attr(pp_obj, "n_events")

  if (missing(parms) || is.null(parms)) {
    parms <- set_parms(pp_obj)
    parms$mu <- n / attr(pp_obj, "end_time")
    parms$alpha <- parms$alpha*0
  }

  temp <- lapply(pp_obj, function(x) x*0)

  # Ugly loop: need to compute compensator for each time point (j) in each margin (i)
  for (i in 1:m) {
    for (j in 2:n[i]) {
      t <- pp_obj[[i]][j]
      temp[[i]][j] <- parms$mu[i] * t + sum(compensator(i, 1:m, pp_obj, parms, k, t))
    }
  }
  lapply(temp, function(x) sort(diff(x)))
}


#-----------------------------------------------------------------------------
#'  Goodness of fit plots and KS tests for a Hawkes process.
#'
#' If \code{parms} is missing, fit is reported for homogeneous Poisson process.

#' @param pp_obj a \code{pp} object.
#' @param parms parameters of the Hawkes process formatted as described for \code{set_parms}. If missing, output is for homogeneous Poisson process.
#' @param kernel_type currently only implemented for \code{"dgamma"}.

#' @return qq plot of residual waiting times with KS test(s) against \code{exp(1)}.
#' @export

gof <- function(pp_obj, parms, kernel_type = "dgamma") {

  resid <- res(pp_obj, parms, kernel_type = "dgamma")
  m <- length(resid)
  quant <- lapply(resid, function(x) qexp(ecdf(x)(x)))
  ks <- mapply(ks.test, resid, quant, SIMPLIFY = F)
  ks_p <- paste("; p =",  lapply(ks, function(x) round(x$p.value, 4)))
  ks_stat <-  paste("KS =", lapply(ks, function(x) round(x$statistic, 4)))
  legend_text <- paste(names(pp_obj), ": ", ks_stat, ks_p, sep = "")
  #lim <- c(0, max(unlist(resid)))
  lim <- c(0, 8)
  plot(resid[[1]], quant[[1]], type = "n",
    main = "QQ Plot of Residual Waiting Times", xlab = "Residual",
    ylab = "Exponential", xlim = lim, ylim = lim, las = 1)

  legend(x = 0, y = lim[2], legend = legend_text, col = c(3:(m+2)),
    lty = 1, lwd = 2, border = "n")

  abline(a = 0, b = 1)
  for (i in 1:m) {
    points(resid[[i]], quant[[i]], type = "p", pch = 1, col = (i+2))
  }
  #ks
}


#-----------------------------------------------------------------------------
#' Under construction: Define response kernel functions
#'
#' For user-chosen \code{kernel_type} in \code{c("dgamma", "dexp", "dlnorm")}, returns a list of repsponse kernel functions. Currently only implemented for \code{"dgamma"}. The goal is to set up various response kernel choices so that they have unifrom args and returns, and catch inadmissable args passed during M-step optimization.
#'
#' @param kernel_type currently only implemented for \code{"dgamma"}.
#' @return A named list of functions with components \code{c("Density", "Dist", "Quantile")}.
#' @export

set_kernel <- function(kernel_type = "dgamma") {
	warn <- "\'kernel_type\' must be one of c(\'dgamma\', \'dexp\', \'dlnorm\')."

	if (!is.character(kernel_type)) stop(warn)
	if (!kernel_type %in% c("dgamma", "dexp", "dlnorm")) stop(warn)

	if (kernel_type == "dgamma") {
		Density <- dgamma
		Dist <- pgamma
		Quantile <- qgamma
	}

	if (kernel_type == "dlnorm") {
		Density <- dlnorm
		Dist <- plnorm
		Quantile <- qlnorm
	}

	temp <- c(Density, Dist, Quantile)
	names(temp) <- c("Density", "Dist", "Quantile")
	structure(temp, kernel_type = kernel_type)
	temp
}

#-----------------------------------------------------------------------------
#' Under construction: Generate random starting values for a Hawkes process.
#'
#'  For a given \code{pp} object with a given \code{kernel_type}, generate random starting values for EM algorithm. Currently only implemented for \code{"dgamma"}. To Do: implement other kernels, generate parameters that follow equality constraints.
#'
#' @param pp_obj a \code{pp} object.
#' @param kernel_type currently only implemented for \code{"dgamma"}.
#' @param constraints currently not implemented.
#'
#' @return A named list of numeric arrays with components \code{c("mu", "alpha", "k_parm1", "k_parm2")} defined as follows. Let n_margins denote the \code{attr(pp_obj, "n_margins")}. Then \code{"mu"} is the n_margins-vector of baseline parameters; \code{"alpha"} is the n_margins-by-n_margins matrix of intensity parameters; \code{"k_parm1"} is the n_margins-by-n_margins matrix containing the first arg to the response kernels; \code{"k_parm2"} is the n_margins-by-n_margins matrix containing the second arg to the response kernels.
#' @export

set_parms <- function(pp_obj, kernel_type = "dgamma", constraints = NULL) {
	m <- attr(pp_obj, "n_margins")
	n <- attr(pp_obj, "n_events")
	end_time <- attr(pp_obj, "end_time")

  mu <- rep(.5, m)
  temp_ind <- which(n > 0)

  if (length(temp_ind)>0) {
    mu[temp_ind] <- lapply(mu[temp_ind], function(x) x<- runif (1, min = .001, max = n/end_time)) %>% unlist
  }

	alpha <- array(runif (m^2, min = .1, max = .9), dim = c(m,m))
	k_parm2 <- array(runif (m^2, min = .5, max = 2), dim = c(m,m))

	if (kernel_type == "dlnorm") {
    k_parm1 <- array(runif (m^2, min = -2, max = 2), dim = c(m,m))
	} else {
    #k_parm1 <- array(1, dim = c(m,m))
    k_parm1 <- array(runif (m^2, min = 1, max = 5), dim = c(m,m))
	}

	# Apply constraints
	if (!is.null(constraints)) {
		alpha[constraints == 0] <- 0
		k_parm1[constraints == 0] <- 0
		k_parm2[constraints == 0] <- 0
		constraints[constraints < 1] <- NA
		equal_ind <- which(duplicated(constraints, MARGIN = 0) & is.na(constraints) == F)

		while (length(equal_ind) > 0) {
			temp_ind <- which(constraints == constraints[equal_ind[1]])
			k_parm1[temp_ind] <- k_parm1[temp_ind[1]]
			k_parm2[temp_ind] <- k_parm2[temp_ind[1]]
			equal_ind <- equal_ind[!equal_ind %in% temp_ind]
		}
	}
	list("mu" = mu, "alpha" = alpha, "k_parm1" = k_parm1, "k_parm2" = k_parm2)
}

#-----------------------------------------------------------------------------
#' Formats output of \code{EM}
#'
#' Reformatrs \code{EM} output to be compatable with other functions that take the \code{parms} argument. Currently only implemented for \code{"dgamma"}.
#'
#' @param em_obj the output of a call to \code{EM}
#' @param solution a scalar indicating which row of \code{EM$Solutions} to read.
#'
#' @return A named list of numeric arrays with components \code{c("mu", "alpha", "k_parm1", "k_parm2")} defined as follows. Let n_margins denote the \code{attr(pp_obj, "n_margins")}. Then \code{"mu"} is the n_margins-vector of baseline parameters; \code{"alpha"} is the n_margins-by-n_margins matrix of intensity parameters; \code{"k_parm1"} is the n_margins-by-n_margins matrix containing the first arg to the response kernels; \code{"k_parm2"} is the n_margins-by-n_margins matrix containing the second arg to the response kernels.
#' @export

get_parms <- function(em_obj, solution = 1) {
  temp <- em_obj$Solutions
  if (is.null(dim(temp))) {
    temp <- t(as.matrix(temp))
  }
	temp <- temp[solution,-1]
  m <- length(temp[grep("mu", names(temp))])
  mu <- temp[1:m]
  grep("^k_parm1", names(temp))
  alpha <- matrix(temp[grep("alpha", names(temp))], nrow = m, ncol = m)
  k_parm1 <- matrix(temp[grep("^k_parm1", names(temp))], nrow = m, ncol = m)
  k_parm2 <- matrix(temp[grep("^k_parm2", names(temp))], nrow = m, ncol = m)
	list("mu" = mu, "alpha" = alpha, "k_parm1" = k_parm1, "k_parm2" = k_parm2)
}

#-----------------------------------------------------------------------------
#' Compute the intensity functions for a Hawkes process.
#'
#' For each event in an output process, compute its intensity under each event in an input process. This is the main workhorse of both the E and M step. Note that the return fromat would more intuitively be represented as n_output_events by n_input_events matrix, but this gets very large and sparse for longer time series. So instead of a matrix, the output is fromated two lists, each being itself a list of length n_output_events -- one list contains intensities, the other the indices, for input events that have non-negligible intensity. This is less intuitive in terms of storage but substantially faster (linear versus quadratic computational complexity; see Halpin 2013).
#'
#' @param out index for which component of \code{pp} object is the output process.
#' @param inp index for which component of \code{pp} object is the input process.
#' @param pp_obj a \code{pp} object.
#' @param parms parameters of the Hawkes process formatted as described for \code{set_parms}.
#' @param k optional: output from \code{set_kernel(kernel_type)}. Currently only implemented for \code{kernel_type = "dgamma"}.
#' @param Log logical: return \code{alpha*kernel} or \code{log(alpha) + log(kernel)} ?
#'
#' @return A named list with components \code{c("Intensity", "Compensator", "Parents")}, each of which are lists with \code{length(pp_obj[[out]])} components. \code{"Intensity"} contains the non-zero intensities for each output event, and \code{"Parents"} contains the parent index (i.e., from \code{pp_obj[[inp]])} for each value in \code{"Intensity"}. Similarly for \code{"Compensator"}.
#' @export

compute_I <- function(out, inp, pp_obj, parms, k, Log = F) {
	o_events <- pp_obj[[out]]
	i_events <- pp_obj[[inp]]
	end_time <- attr(pp_obj, "end_time")
	alpha <- parms$alpha[out, inp]
	k_parm1 <- parms$k_parm1[out, inp]
	k_parm2 <- parms$k_parm2[out, inp]

	if (missing(k)) k <- set_kernel()
	k_dist <- function(x)  alpha * k$Dist(x, k_parm1, k_parm2)

	if (Log == F) {
		k_density <- function(x) alpha * k$Density(x, k_parm1, k_parm2)
	} else {
		k_density <- function(x) log(alpha) + k$Density(x, k_parm1, k_parm2, log = T)
	}

	# Controls rounding error; see Halpin 2013.
	error <- 1e-10
	delta <- k$Quantile(1 - error, k_parm1, k_parm2)

	# For each output, get indices of inputs that are possible parents
	parent_ind <- lapply(o_events, function(x) which(i_events < x & i_events > x - delta))

	# For each output, get time lag from each possible parent
	parent_lag <- mapply(function(x, y) o_events[x] - i_events[y], 1:length(o_events), parent_ind)

	# For each output, get intensity under each possible parent
	temp_I <- lapply(parent_lag, k_density)
	temp_I <- lapply(temp_I, function(x) if (length(x) == 0) { x <- 0 } else { x <- x })

	# For each parent of each output, get compensator (this seems inefficient)
	temp_C <- lapply(parent_ind, function(x) k_dist(end_time - i_events[x]))
	temp_C <- lapply(temp_C, function(x) if (length(x) == 0) { x <- 0 } else { x <- x })

	temp <- list(temp_I, temp_C, parent_ind)
	names(temp) <- c("Intensity", "Compensator", "Parents")
	temp
}

#-----------------------------------------------------------------------------
#' Wrapper for \code{compute_I$Intensity}.
#'
#' Takes an m-dimensional multivariate Hawkes process and computes the intensities for indicated output and input processes. Main difference with \code{compute_I} is that \code{inp} can be a vector.
#'
#' @param out index for which component of \code{pp} object is the output process.
#' @param inp a vector of indices for which components of \code{pp} object are the input process.
#' @param pp_obj a \code{pp} object.
#' @param parms parameters of the Hawkes process formatted as described for \code{set_parms}.
#' @param k output from \code{set_kernel(kernel_type)}. Currently only implemented for \code{kernel_type = "dgamma"}.
#' @param Log logical: return \code{alpha*kernel} or \code{log(alpha) + log(kernel)} ?
#'
#' @return A list of \code{length(inp)} components, where the j-th element is \code{compute_I} for \code{pp_obj[[out]]} and \code{pp_obj[[inp[j]]]}.
#' @export

intensity <- function(out, inp, pp_obj, parms, k, Log = F) {
	n <- length(inp)
	temp <- vector("list", n)

	for (j in 1:n) {
		temp[[j]] <- compute_I(out, j, pp_obj, parms, k, Log)
	}
	names(temp) <- paste(out, inp, sep = ":")
	temp
}

#-----------------------------------------------------------------------------
#' Wrapper for \code{compute_I$Compensator}.
#'
#' Similar to \code{intensity}. Takes an m-dimensional multivariate process and computes the compensator for indicated output and input processes.
#'
#' @param out index for which component of \code{pp} object is the output process.
#' @param inp a vector of indices for which components of \code{pp} object are the input process.
#' @param pp_obj a \code{pp} object.
#' @param parms parameters of the Hawkes process formatted as described for \code{set_parms}.
#' @param k output from \code{set_kernel(kernel_type)}. Currently only implemented for \code{kernel_type = "dgamma"}.
#' @param # t: override default time interval \code{[0, attr(events, "end_time"]} with a different upper bound.
#'
#' @return vector of \code{length(inp)}, where the j-th element is the (scalar) value of the compensator for \code{events[[out]]} and \code{events[[inp[j]]]}.
#' @export

compensator <- function(out, inp, pp_obj, parms, k, t = NULL) {
	if (!is.null(t)) {attr(pp_obj, "end_time") <- t}

	n_inp <- length(inp)
	i_events <- pp_obj[inp]
	end_time <- attr(pp_obj, "end_time")

	temp_ind <- lapply(pp_obj, function(x) max(c(which(x <= end_time), 0)))
	temp <- rep(0, times = n_inp)

	for (j in 1:n_inp) {
		if (temp_ind[[j]] > 0) {
			temp[j] <- with(parms, alpha[out, j] * sum(k$Dist(end_time - i_events[[j]], k_parm1[out, j], k_parm2[out, j])))
		}
	}
	names(temp) <- paste(out, inp, sep = ":")
	temp
}

#-----------------------------------------------------------------------------
#' Normalizes the output of \code{intensity} to be probabilities.
#'
#' Similar to \code{intensity}, except the intensities are rescaled to probabilities.
#'
#' @param out index for which component of \code{pp} object is the output process.
#' @param inp a vector of indices for which components of \code{pp} object are the input process.
#' @param pp_obj a \code{pp} object.
#' @param parms parameters of the Hawkes process formatted as described for \code{set_parms}.
#' @param k output from \code{set_kernel(kernel_type)}. Currently only implemented for \code{kernel_type = "dgamma"}.
#'
#' @return A list of \code{length(inp)} components, where the j-th element is \code{compute_I} for \code{pp_obj[[out]]} and \code{pp_obj[[inp[j]]]}, normalized over the input processes.
#' @export

probs <-function(out, inp, pp_obj, parms, k) {
  list_I <- intensity(out, inp, pp_obj, parms, k)
	n <- attr(pp_obj, "n_events")[out]
	m <- length(inp)
	mu <- parms$mu[out]
	end_time <- attr(pp_obj, "end_time")
	temp_I <- lapply(list_I, function(x) x$Intensity)
	temp_C <- lapply(list_I, function(x) x$Compensator)

	# Allows for different methods of computing posteriors
	p_sum <- function(list_1, list_2)
	{
		mapply(function(x, y) x, list_1, list_2)
	}

	temp_P <- mapply(p_sum, temp_I, temp_C, SIMPLIFY = F)

	total <- lapply(temp_P, function(x) lapply(x, sum)) %>% unlist %>% array(dim = c(n, m)) %>% cbind(mu) %>% apply(1, sum)

	temp_P <- lapply(temp_P, function(r) mapply(function(x, y) x / y, r, total, SIMPLIFY = F))

	# Storage: How to do this with mapply? list_I is list of lists...
	temp <- vector("list", m)

	for (i in 1:m) {
	  temp[[i]] <- list(temp_P[[i]],list_I[[i]]$Parents)
		names(temp[[i]]) <- c("Prob", "Parents")
	}
	names(temp) <- names(list_I)
	temp
}


#-----------------------------------------------------------------------------
#' Computes the Q function the EM algorithm.
#'
#' Computes obejctive function of the M step of the EM algorithm, for a single input / output process combination
#'
#' @param out index for which component of \code{pp} object is the output process.
#' @param inp index for which component of \code{pp} object is the input process.
#' @param pp_obj a \code{pp} object.
#' @param parms parameters of the Hawkes process formatted as described for \code{set_parms}.
#' @param k output from \code{set_kernel(kernel_type)}. Currently only implemented for \code{kernel_type = "dgamma"}.
#' @param list_P output from \code{probs}.
#' @return  The (scalar) value of the Q function.
#' @export

compute_Q <- function(out, inp, pp_obj, parms, k, list_P) {
	m <- attr(pp_obj, "n_margins")
	n_events <- attr(pp_obj, "n_events")[out]
	parms$alpha <- array(compute_alpha(out, inp, pp_obj, parms, k, list_P), dim = c(m,m))
	temp_I <- compute_I(out, inp, pp_obj, parms, k, Log = T)
	temp_C <- compensator(out, inp, pp_obj, parms, k)
  P <- list_P[[paste(out, inp, sep = ":")]]

	# Loop to weight temp_I$Intensity[[i]] by P$Prob[[i]] -- slippery indexing.
	for (i in 1:n_events) {
		m <- length(temp_I$Parent[[i]])

		# Do nothing if temp_I$Parent[[i]] is empty
		if (m > 0) {
			temp_P <- rep(0, times = m)
			n <- temp_I$Parent[[i]] %in% P$Parent[[i]]

			# Get non-zero P$Prob[[i]] for each temp_I$Parent[[i]]
			if (length(n) > 0) {
				r <- P$Parent[[i]] %in% temp_I$Parent[[i]]
				temp_P[n] <- P$Prob[[i]][r]  # slippery
			}

			# Finally, apply the weights
			temp_I$Intensity[[i]] <- temp_I$Intensity[[i]] * temp_P
		}
	}
	sum(unlist(temp_I$Intensity)) - temp_C
}


#-----------------------------------------------------------------------------
#' Under construction: Wrapper for \code{compute_Q}.
#'
#' Returns sum of \code{compute_Q} over different the output/input combinations given by \code{out_inp}. Currently only implemented for constraints over a single output process (i.e., only 1 row of \code{out_inp} can have non-zero entries). Used to apply model constraints on response kernels. Yuck.

#' @param out_inp n_margins-by-n_margins binary matrix in which \code{out_inp_ij = 1} denotes the output/input combinations on which to call compute_Q. When only one element of \code{out_inp} is equal is 1, this is equivalent to calling compute_Q directly on the output/input pair denoted by that element.
#' @param pp_obj a \code{pp} object.
#' @param parms parameters of the Hawkes process formatted as described for \code{set_parms}.
#' @param k output from \code{set_kernel(kernel_type)}. Currently only implemented for \code{kernel_type = "dgamma"}.
#' @param list_P output from \code{probs}.
#' @return  The (scalar) sum of the values of the Q functions.
#' @export

Q_fun <- function(out_inp, pp_obj, parms, k, list_P) {
  ind <- which(out_inp == 1, arr.ind = T)
	total <- 0

	for (i in 1:nrow(ind)) {
		out <- ind[i,1]
		inp <- ind[i,2]
		total <- total + compute_Q(out, inp, pp_obj, parms, k, list_P)
	}
	total
}


#-----------------------------------------------------------------------------
#' Computes the baserate parameter of the Hawkes process.
#'
#'  Computes the (constant) baserate for a single output process.
#'
#' @param out index for which component of \code{pp} object is the output process.
#' @param inp a vector of indices for which components of \code{pp} object are the input process.
#' @param pp_obj a \code{pp} object.
#' @param list_P output from \code{probs}.
#' @return The (scalar) value of the (constant) baserate.
#' @export

compute_mu <- function(out, inp, pp_obj, list_P) {
  n_events <- attr(pp_obj, "n_events")[out]
	end_time <- attr(pp_obj, "end_time")
	p_ind <- paste(out, inp, sep = ":")
	sum_P <- lapply(list_P[p_ind], function(x) x[[1]]) %>% unlist  %>% sum
	(n_events - sum_P) / end_time
}


#-----------------------------------------------------------------------------
#' Computes the intensity parameter of the Hawkes process.
#'
#' Compute the intensity parameters for single intput/output process combination.
#'
#' @param out index for which component of \code{pp} object is the output process.
#' @param inp index for which component of \code{pp} object is the input process.
#' @param pp_obj a \code{pp} object.
#' @param k output from \code{set_kernel(kernel_type)}. Currently only implemented for \code{kernel_type = "dgamma"}.
#' @param list_P output from \code{probs}.
#' @return The (scalar) value of the intensity parameter.
#' @export

compute_alpha <- function(out, inp, pp_obj, parms, k, list_P) {
	end_time <- attr(pp_obj, "end_time")
	p_ind <- paste(out, inp, sep = ":")
	sum_P <- list_P[[p_ind]]$Prob %>% unlist %>% sum
	sum_P / with(parms, sum(k$Dist(end_time - pp_obj[[inp]], k_parm1[out, inp], k_parm2[out, inp])))
}


#-----------------------------------------------------------------------------
#' Computes the incomplete data log likelihood of the Hawkes process.
#'
#' Computes the contribution of a single margin to the incomplete data log likelihood of the Hawkes process.
#'
#' @param out index for which component of \code{pp} object is the output process.
#' @param inp a vector of indices for which components of \code{pp} object are the input process.
#' @param pp_obj a \code{pp} object.
#' @param parms parameters of the Hawkes process formatted as described for \code{set_parms}.
#' @param k output from \code{set_kernel(kernel_type)}. Currently only implemented for \code{kernel_type = "dgamma"}.
#' @return The (scalar) value of incomplete data log likelihood for the output margin.
#' @export

logL <- function(out, inp, pp_obj, parms, k) {
	end_time <- attr(pp_obj, "end_time")
	temp_C <- compensator(out, inp, pp_obj, parms, k)
	temp_I <- intensity(out, inp, pp_obj, parms, k) %>% lapply(function(x) (x$Intensity))

	# Store as df. Simplify this??
	temp_df <- lapply(temp_I, function(x) lapply(x, function(x) sum(unlist(x)) )) %>% lapply(function(x) unlist(x)) %>% data.frame()

	temp_df$mu <- parms$mu[out]
	temp_out <- apply(temp_df, 1, sum) %>% log %>% sum
	temp_out - sum(temp_C) - parms$mu[out] * end_time
}


#----------------------------------------------------------------------------
#' Runs EM algorithm for Multivariate Hawkes process.
#'
#' The only required argument is \code{pp_obj}. Repeated calls are made to \code{optim} with \code{method = "Nelder-Mead"} and its default arguments. Nelder-Mead may attempt non-admissilbe parameter values for the kernel densities of the Hawkes process (e.g., a negative scale/shape parameter for the gamma distribution), in which case \code{optim} will throw a warning message -- but this is still faster than optim's \code{method = "L-BFGS-B"}. To Do: 1) catch \code{optim} warnings and build box constraints into kernels. 2) Implement equality constraints over response functions. While the source code has many (too many) lines of code directed at equality constraints over response functions, this is currently not implemented robustly so the \code{constraints} argument is omitted from the function call. Also need to improve \code{Trace = T} output.
#'
#' @param pp_obj a \code{pp} object.
#' @param kernel_type currently only implemented for \code{kernel_type "dgamma"}  (optional).
#' @param starts starting values for parameters of the Hawkes process formatted as described for \code{set_parms}. If omitted randomly generated starting values are used.
#' @param nstarts How many times to run the EM algorithm (to address non-convex likelihood). Default = 1. Note that if \code{nstarts > 1} then \code{starts} should be omitted, or else each different run of the algorithm will converge to the same value.
#' @param conv convergenve criterion of EM algorithm (absolute change in complete data log likelihood on successive steps M steps). Default = 1e-4.
#' @param maxit maximum number of iterations of the EM algorithm. Default = 500.
#' @param Trace (logical) if \code{True}, some graphical convergence diagnostics are displayed during optimization.
#' @return  A named list with components \code{c("Solutions", "History")}. \code{"Solutions"} is a matrix with \code{nstarts} rows, the first column of which cotains the convergence value of the incomplete data loglikelihood for those starting values, with the remaining rows containing the parameters estimates. \code{"History"} is an \code{nstarts}-by-\code{maxit} matrix that contains the convergence history for each set of starting values. Both matrices are ordered from the best to worst solutions.
#' @export

EM <- function(pp_obj, kernel_type, starts, nstarts, conv, maxit, Trace) {

	# Need to set up proper "catch list" for missing and bad arguments...
	warn <- " The arg \'pp_obj\' must be an object of class 'pp'. See function 'pp'."

	if (missing(pp_obj)) stop(cat(warn))
	if (class(pp_obj) != "pp") stop(cat(warn))

	# Set up constants
	if (missing(conv)) conv <- 1e-4
	if (missing(nstarts)) nstarts <- 1
	if (missing(maxit)) maxit <- 500
	if (missing(kernel_type)) kernel_type <- "dgamma"
	if (missing(Trace)) Trace <- F

	end_time <- attr(pp_obj, "end_time")
	m <- attr(pp_obj, "n_margins")
	n <- attr(pp_obj, "n_events")
	k <- set_kernel(kernel_type)  # Don't use k as an index!


  #Helper functions; require args from EM env, not for external use
	logL_fun <- function(parms) {
		temp <- 0
		for (i in 1:m) {
			inp <- which(constraints[i,] != 0)
			temp <- temp + logL(i, inp, pp_obj, parms, k)
		}
		temp
	}

	mu_fun <- function(list_P) {
		temp <- rep(0, m)
		for (i in 1:m) {
			inp <- which(constraints[i, ] != 0)
		 	temp[i] <- compute_mu(i, inp, pp_obj, list_P)
	 	}
		temp
	}

	alpha_fun <- function(i, parms, list_P) {
		ind <- which(constraints != 0, arr.ind = T)
		out <- ind[i,1]
		inp <- ind[i,2]
		compute_alpha(out, inp, pp_obj, parms, k, list_P)
	}

	apply_constraints <- function() {
		out_list <- vector("list", n_Q)
		names(out_list) <- paste0("set", 1:n_Q)
		temp <- array(0, dim = c(m,m))

		for (i in 1:n_Q) {
			temp <- temp*0
			temp_ind <- which(constraints == constraints[unique_ind[i]])
			temp[temp_ind] <- 1
			out_list[[i]] <- temp
		}
		out_list
	}

	P_fun <- function(i, parms) {
		inp <- which(constraints[i,] != 0)
		probs(i, inp, pp_obj, parms, k)
	}

  # Need to change this up to apply different constraints :(
	Q_optim <- function(par, equal, list_P) {

		#k_parm1 <- array(1, dim = c(m,m))
		#k_parm2 <- array(par, dim = c(m,m))
		k_parm1 <- array(par[1], dim = c(m,m))
		k_parm2 <- array(par[2], dim = c(m,m))
		o_parms <- list("k_parm1" = k_parm1, "k_parm2" = k_parm2)

		Q_fun(equal, pp_obj, o_parms, k, list_P)
	}

	# Set up constraints; unique_ind and n_Q index the unique Q-functions
	#if (missing(constraints))
  constraints <- array(1:m^2, dim = c(m,m))
	unique_ind <- which(!duplicated(constraints, MARGIN = 0) & constraints > 0)
	n_Q <- length(unique_ind)
  equal <- apply_constraints()

 	# Set up storage and control
	History <- matrix(0, nrow = nstarts, ncol = maxit+1)
	Solutions <- matrix(0, nrow = nstarts, ncol = (3 * m^2 + m + 1))
	temp_P <- vector("list", m)

	# Set up Trace
	if (Trace == T) {
		best <- -1e7
		par(mfrow = c(m,m))
	}

	# nstarts loop
	for (i in 1:nstarts) {
		iter <- 1
		epsilon <- 10
		parms <- set_parms(pp_obj, kernel_type, constraints)
		History[i,(iter)] <- logL_fun(parms)

		# EM loop
		while (epsilon > conv & iter <= maxit) {

			# E step:
			for (j in 1:m) {
				temp_P[[j]] <- P_fun(j, parms)
			}
		  list_P <- unlist(temp_P, recursive = F)

			# M step (a): maximumize Q; skip on first iter to help avoid bad starts? # j = 1
			#if (iter > 1) {
				for (j in 1:n_Q) {
					#par <- with(parms, k_parm2[unique_ind[j]])

					par <- with(parms, c(k_parm1[unique_ind[j]], k_parm2[unique_ind[j]]))

					S <- optim(par = par,
	      	 				 fn = Q_optim,
									 equal = equal[[j]],
              		 list_P = list_P,
									 method  = "Nelder-Mead",
         					 control = list(fnscale = -1, trace = 1))

			    #parms$k_parm2[equal[[j]] == 1] <- S$par
					parms$k_parm1[equal[[j]] == 1] <- S$par[1]
					parms$k_parm2[equal[[j]] == 1] <- S$par[2]

			  }
			#}

			# M step (b): update mu and alpha
			parms$mu <- mu_fun(list_P)

			for (j in which(constraints != 0)) {
				 parms$alpha[j] <- alpha_fun(j, parms, list_P)
			}

			# Trace
			current <-  logL_fun(parms)
			History[i,(iter+1)] <- current

			if (Trace == T) {
				if (current > best) {best <- current}

				main <- c( paste("Best = ", round(best,5 ), "; Current = ", round(current, 5), sep = ""), rep("", times = m^2 - 1) )

				for (j in 1:m^2) {
					with(parms, curve(alpha[j] * dgamma(x, shape = k_parm1[j], rate = k_parm2[j]), to = 30, ylab = "Response Kernel", main = main[j] ))
				}
			}

			# Convergence
			epsilon <- (History[i,(iter+1)] - History[i,iter])
			iter <-iter + 1

		} # end EM loop

	  Solutions[i, 1] <- History[i, iter]
		Solutions[i, -1] <- unlist(parms)

	} # end nstarts loop

	ind <- order(Solutions[,1], decreasing = T)
	colnames(Solutions) <- c("logL", names(unlist(parms)))
	colnames(History) <- paste0("iter", 0:maxit)
	temp <- list(Solutions[ind,], History[ind, ])
	names(temp) <- c("Solutions", "History")
	temp
}



# logH <- function(parms, pp_obj, time, constraints, k = set_kernel()) {
#
#---------------------------------------------------------------------------
# #  Wrapper for logL to be used when computing Hessian with numDeriv. Only supports full constraints.
#
# # Arguments
# # parms is c(mu, alpha, k_parm_1, k_parm_2); k_parm1 and k_parm2 scalars
#
# # Returns
# # The value of incomplete data log likelihood for the output margin
#
#---------------------------------------------------------------------------
#
# m <- length(pp_obj)
# mu <- parms[1:2]
# alpha <- array(parms[3:6], dim = c(m,m))
#
# k_parm1 <- array(parms[7], dim = c(m,m))
# k_parm2 <- array(parms[8], dim = c(m,m))
# temp <- 0
#
# for (i in 1:m) {
# 	inp <- which(constraints[i,] != 0)
# 	temp <- temp + logL(i, inp, events, time, mu, alpha, k_parm1, k_parm2, k)
# 	}
# temp
# }
