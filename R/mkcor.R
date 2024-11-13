#' Fit the correlated Mk model by expectation-maximization
#' 
#' Returns a function that can be used to fit the correlated Mk model
#' to a set of observed character state data on a phylogeny using the 
#' expectation-maximization algorithm.
#'
#' @param x A data.frame with character state data. One column should be
#' named 'tip.label' and another should be named 'state.id'. Character state
#' id's should be stored as factors. If a terminal node state is
#' ambiguous just list it in multiple rows with one row for each possible
#' state. If a terminal node state is unknown it can simply be omitted from the
#' data.frame.
#' @param y A data.frame for the second character formatted like \code{x}.
#' @param phy An \code{ephylo} object.
#' @return A function that can be used to fit a model of correlated evolution
#' of \code{x} and \code{y} by expectation-maximization. The parameters to the
#' function are:
#' \itemize{
#' \item \code{par} A vector of initial parameter estimates. The first two values are
#'   rates of independent change in x and y, respectively. If given, the third 
#'   value is the rate of simultaneous change in x and y.
#' \item \code{maxiter} The maximum number of EM iterations to perform (default 100).
#' \item \code{tol} Stop the EM algorithm when the log likelihood improves by less
#'   than this amount between consecutive iterations (default 0.001).
#'}
#' The function returns a list with the following components:
#' \itemize{
#' \item \code{par} The maximum likelihood estimates of the overall rates of
#'   character evolution and character correlation.
#' \item \code{value} The log likelihood of the fitted model.
#' \item \code{event.counts} The expected number of different types of character state
#'   changes summed over all branches.
#' \item \code{state.probs} Ancestral state probabilities for each node in \code{phy}.
#'   E.g., \code{state.probs[i,j]} is the probability that node \eqn{i} is in state 
#'   \eqn{j}.
#' \item \code{branch.counts} An array with the expected number of different types of 
#'   character state changes on each branch. E.g., \code{branch.counts[i,j,u]} returns
#'   the expected number of evolutionary transitions from state \eqn{i} to state \eqn{j}
#'   on the branch leading to node \eqn{u}.
#' \item \code{logL} The log likelihood for each iteration of the EM algorithm. The
#'   final value corresponds to the parameter estimates stored in the other
#'   list components.
#' \item \code{states} A two column matrix enumerating the states in the model as
#'   combinations of the different states in \code{x} and \code{y}.
#'}
#' @seealso \code{\link{mkcor_fit_em}}
mkcor_em = function(x, y, phy)
{
    stopifnot(inherits(phy, "ephylo"))
    stopifnot(is.data.frame(x))
    stopifnot(is.data.frame(y))
    stopifnot(all(names(x) %in% c("tip.label", "state.id")))
    stopifnot(all(names(y) %in% c("tip.label", "state.id")))
    stopifnot(is.character(x$tip.label))
    stopifnot(is.factor(x$state.id))
    stopifnot(is.character(y$tip.label))
    stopifnot(is.factor(y$state.id))
    num_states_x = nlevels(x$state.id)
    num_states_y = nlevels(y$state.id)
    num_states = num_states_x * num_states_y
    num_tips = phy$num.tips
    num_nodes = phy$num.nodes
    left_child = phy$left.child
    right_sib = phy$right.sib
    postorder = phy$postorder.nodes
    preorder = phy$preorder.nodes
    parent = phy$parent
    brlen = phy$brlen
    brlen[phy$root] = 0
    node_loglik = matrix(0, num_states, num_nodes)

    # we can recover the joint state id as
    #
    #   ([y-1] + [x-1]*num_states_y) + 1
    #
    # because we are numbering the states in column-major
    # order of a num_states_y by num_states_x matrix (i.e.,
    # column numbers correspond to states in x and row numbers
    # correspond to states in y).
    #
    # e.g.,
    #     x 1 2 3 4  
    # y 1 [ 1 4 7 10 ]
    #   2 [ 2 5 8 11 ]
    #   3 [ 3 6 9 12 ]
    state_pairs = setNames(
        expand.grid(levels(y$state.id), levels(x$state.id))[,2:1], c("X", "Y"))

    x = split(as.integer(x$state.id), x$tip.label)[phy$tip.label]
    y = split(as.integer(y$state.id), y$tip.label)[phy$tip.label]

    for (tip in union(names(x), names(y)))
    {
        u = match(tip, phy$tip.label)
        node_loglik[, u] = -Inf
        for (xl in x[[tip]])
        {
            for (yl in y[[tip]])
            {
                i = ((yl - 1) + (xl - 1)*num_states_y) + 1
                node_loglik[i, u] = 0
            }
        }
    }
    for (tip in setdiff(names(x), names(y))) # missing y
    {
        u = match(tip, phy$tip.label)
        node_loglik[, u] = -Inf
        for (xl in x[[tip]])
        {
            for (yl in 1:num_states_y)
            {
                i = ((yl - 1) + (xl - 1)*num_states_y) + 1
                node_loglik[i, u] = 0
            }
        }
    }
    for (tip in setdiff(names(y), names(x))) # missing x
    {
        u = match(tip, phy$tip.label)
        node_loglik[, u] = -Inf
        for (xl in 1:num_states_x)
        {
            for (yl in 1:y[[tip]])
            {
                i = ((yl - 1) + (xl - 1)*num_states_y) + 1
                node_loglik[i, u] = 0
            }
        }
    }

    function(par, tol=0.001, maxiter=100L)
    {
        num_pars = length(par)
        stopifnot(num_pars == 2 || num_pars == 3)
        stopifnot(all(par >= 0))
        ans = .Call(
            C_mkcor_em,
            num_states,
            num_states_x,
            num_states_y,
            as.numeric(par),
            as.numeric(tol),
            as.integer(maxiter),
            node_loglik,
            num_tips,
            num_nodes,
            postorder,
            preorder,
            parent,
            left_child,
            right_sib,
            brlen
        )

        r = ans[[1]][3] / sqrt(sum(ans[[1]][c(1,3)])*sum(ans[[1]][c(2,3)]))
        ans[[1]] = c(ans[[1]], r)
        names(ans[[1]]) = c("rate.x", "rate.y", "rate.xy", "corr.xy")
        names(ans[[2]]) = c("N.x", "N.y", "N.xy")
        ans[[3]] = t(ans[[3]])
        ans[[6]] = state_pairs
        ret = vector("list", 7L)
        ret[[1]] = ans[[1]]
        ret[[2]] = tail(ans[[5]], 1L)
        ret[3:7] = ans[2:6]
        setNames(ret, c("par", "value", "event.counts", "state.probs", 
            "branch.counts", "logL", "states"))
    }
}


#' Fit the correlated Mk model by expectation-maximization
#' 
#' Return parameter estimates of the correlated Mk model given a set of 
#' observed character state data using expectation-maximization.
#'
#' @param x A data.frame with character state data. One column should be
#' named 'tip.label' and another should be named 'state.id'. Character state
#' id's should be stored as factors. If a terminal node state is
#' ambiguous just list it in multiple rows with one row for each possible
#' state. If a terminal node state is unknown it can simply be omitted from the
#' data.frame.
#' @param y A data.frame for the second character formatted like \code{x}.
#' @param phy An \code{ephylo} object.
#' @param correlated A boolean indicating if the correlated model (TRUE) or
#' the independent model (FALSE) should be fit.
#' @param par_init Optional initial parameter estimates to seed the algorithm.
#' If not \code{NULL}, any value for \code{correlated} and \code{num_fits} is
#' ignored and a single EM chain is run until termination.
#' @param maxiter The maximum number of EM iterations to perform (default 100).
#' @param tol Stop the EM algorithm when the log likelihood improves by less
#'   than this amount between consecutive iterations (default 0.001).
#' @param lower_init A lower bound for initial rate parameters (>= 0) used to
#' seed the fitting algorithm (default 0).
#' @param upper_init An upper bound for initial rate parameters used to seed
#' the fitting algorithm (default 0.001).
#' @param ndeps A scalar that multiplies the parameter estimates. The resulting
#' vector is then passed as the \code{ndeps} control parameter when calculating
#' the Hessian matrix. See \code{optim} for details. The Hessian is
#' used for calculating standard errors of parameter estimates.
#' @param num_fits Return the best (highest log likelihood) of \code{num_fits} 
#' models.
#' @return A list with the following components:
#' \itemize{
#' \item \code{par} The maximum likelihood estimates of the overall rates of
#'   character evolution and character correlation along with their standard
#'   errors.
#' \item \code{value} The log likelihood of the best-fitting model.
#' \item \code{event.counts} The expected number of different types of character state
#'   changes summed over all branches.
#' \item \code{state.probs} Ancestral state probabilities for each node in \code{phy}.
#'   E.g., \code{state.probs[i,j]} is the probability that node \eqn{i} is in state 
#'   \eqn{j}.
#' \item \code{branch.counts} An array with the expected number of different types of 
#'   character state changes on each branch. E.g., \code{branch.counts[i,j,u]} returns
#'   the expected number of evolutionary transitions from state \eqn{i} to state \eqn{j}
#'   on the branch leading to node \eqn{u}.
#' \item \code{logL} The log likelihood for each iteration of the EM algorithm. The
#'   final value corresponds to the parameter estimates stored in the other
#'   list components.
#' \item \code{states} A two column matrix enumerating the states in the model as
#'   combinations of the different states in \code{x} and \code{y}.
#'}
#' @seealso The \code{vignette("CorrelatedMk")}
mkcor_fit_em = function(x, y, phy, correlated=TRUE, par_init=NULL, maxiter=100L, 
    tol=0.001, lower_init=0, upper_init=0.001, ndeps=0.3, num_fits=1L)
{
    EM = mkcor_em(x, y, phy)
    LIK = mkcor_loglik(x, y, phy)
    if (!is.null(par_init))
    {
        stopifnot(is.numeric(par_init))
        stopifnot(length(par_init) == 2L || length(par_init) == 3L)
        stopifnot(all(par_init >= 0))
        num_pars = length(par_init)
        fit = EM(par_init, maxiter=maxiter, tol=tol)
    }
    else
    {
        num_pars = if (correlated) 3L else 2L
        best_score = -Inf
        for (i in 1:num_fits)
        {
            tmp = EM(runif(num_pars, lower_init, upper_init), maxiter=maxiter,
                tol=tol)
            if (tmp$value > best_score)
            {
                best_score = tmp$value
                fit = tmp
            }
        }
    }
    hess = try(optimHess(fit$par[1:num_pars], LIK, 
        control=list(fnscale=-1, ndeps=ndeps*fit$par[1:num_pars])), silent=TRUE)
    V = try(solve(-hess), silent=TRUE)
    if (!inherits(V, "try-error"))
        cV = try(chol(V), silent=TRUE)
    else
        cV = V
    # asymptotic standard error estimate
    if (!inherits(V, "try-error") && !inherits(cV, "try-error"))
    {
        se = sqrt(diag(V))
    }
    else
    {
        se = rep(NA_real_, num_pars)   
    }
    
    if (num_pars == 3L)
    {
        rx = fit$par[1]
        ry = fit$par[2]
        rxy = fit$par[3]
        r = fit$par[4]
        # asymptotic standard error estimate using delta method
        # https://en.wikipedia.org/wiki/Delta_method#Multivariate_delta_method
        # these are the derivatives of the correlation coefficient w.r.t.
        # \lambda_x (rx), \lambda_y (ry), and \lambda_{xy} (rxy)
        grad = c(
            -(rxy*(ry+rxy))/(2*((rx+rxy)*(ry+rxy))^(3/2)),
            -(rxy*(rx+rxy))/(2*((rx+rxy)*(ry+rxy))^(3/2)),
            (rxy*(rx+ry)+2*rx*ry)/(2*((rx+rxy)*(ry+rxy))^(3/2))
        )
        if (!inherits(V, "try-error") && !inherits(cV, "try-error"))
        {
            se = c(se, c(sqrt(t(grad) %*% V %*% grad)))
        }
        else
        {
            se = c(se, NA_real_)
        }

        fit$par = data.frame(
            par=unname(fit$par), 
            se=se,
            row.names=c("rate.x", "rate.y", "rate.xy", "corr.xy")
        )

    }
    else
    {

        fit$par = data.frame(
            par=unname(fit$par),
            se=c(se, 0, 0),
            row.names=c("rate.x", "rate.y", "rate.xy", "corr.xy")
        )

    }
    
    return (fit)
}


#' Calculate the log likelihood of the correlated Mk model
#' 
#' Returns a function that calculates the log likelihood of parameters of
#' the correlated Mk model given a set of observed character state data.
#'
#' @param x A data.frame with character state data. One column should be
#' named 'tip.label' and another should be named 'state.id'. Character state
#' id's should be stored as factors. If a terminal node state is
#' ambiguous just list it in multiple rows with one row for each possible
#' state. If a terminal node state is unknown it can simply be omitted from the
#' data.frame.
#' @param y A data.frame for the second character formatted like \code{x}.
#' @param phy An \code{ephylo} object.
#' @return A function that calculates the log likelihood of parameters of
#' the correlated Mk model given the observed data. The function takes a
#' single argument:
#' \itemize{
#' \item \code{par} 
#'   A vector of transition rate parameters. The first two values are
#'   rates of independent change in x and y, respectively. If given, the third 
#'   value is the rate of simultaneous change in x and y.
#'}
#' @seealso \code{\link{mkcor_fit}}
mkcor_loglik = function(x, y, phy)
{
    stopifnot(inherits(phy, "ephylo"))
    stopifnot(is.data.frame(x))
    stopifnot(is.data.frame(y))
    stopifnot(all(names(x) %in% c("tip.label", "state.id")))
    stopifnot(all(names(y) %in% c("tip.label", "state.id")))
    stopifnot(is.character(x$tip.label))
    stopifnot(is.factor(x$state.id))
    stopifnot(is.character(y$tip.label))
    stopifnot(is.factor(y$state.id))
    num_states_x = nlevels(x$state.id)
    num_states_y = nlevels(y$state.id)
    num_states = num_states_x * num_states_y
    num_tips = phy$num.tips
    num_nodes = phy$num.nodes
    left_child = phy$left.child
    right_sib = phy$right.sib
    postorder = phy$postorder.nodes
    preorder = phy$preorder.nodes
    parent = phy$parent
    brlen = phy$brlen
    brlen[phy$root] = 0
    node_loglik = matrix(0, num_states, num_nodes)

    x = split(as.integer(x$state.id), x$tip.label)[phy$tip.label]
    y = split(as.integer(y$state.id), y$tip.label)[phy$tip.label]

    for (tip in union(names(x), names(y)))
    {
        u = match(tip, phy$tip.label)
        node_loglik[, u] = -Inf
        for (xl in x[[tip]])
        {
            for (yl in y[[tip]])
            {
                i = ((yl - 1) + (xl - 1)*num_states_y) + 1
                node_loglik[i, u] = 0
            }
        }
    }
    for (tip in setdiff(names(x), names(y))) # missing y
    {
        u = match(tip, phy$tip.label)
        node_loglik[, u] = -Inf
        for (xl in x[[tip]])
        {
            for (yl in 1:num_states_y)
            {
                i = ((yl - 1) + (xl - 1)*num_states_y) + 1
                node_loglik[i, u] = 0
            }
        }
    }
    for (tip in setdiff(names(y), names(x))) # missing x
    {
        u = match(tip, phy$tip.label)
        node_loglik[, u] = -Inf
        for (xl in 1:num_states_x)
        {
            for (yl in 1:y[[tip]])
            {
                i = ((yl - 1) + (xl - 1)*num_states_y) + 1
                node_loglik[i, u] = 0
            }
        }
    }

    function(par)
    {
        stopifnot(all(par >= 0))
        stopifnot(length(par) == 2 || length(par) == 3)
        .Call(
            C_mkcor_loglik,
            num_states,
            num_states_x,
            num_states_y,
            as.numeric(par),
            node_loglik,
            num_tips,
            num_nodes,
            postorder,
            left_child,
            right_sib,
            brlen
        )
    }
}


#' Fit the correlated Mk model with L-BFGS-B 
#' 
#' Return parameter estimates of the correlated Mk model given a set of 
#' observed character state data using the L-BFGS-B numerical optimization
#' algorithm.
#'
#' @param x A data.frame with character state data. One column should be
#' named 'tip.label' and another should be named 'state.id'. Character state
#' id's should be stored as factors. If a terminal node state is
#' ambiguous just list it in multiple rows with one row for each possible
#' state. If a terminal node state is unknown it can simply be omitted from the
#' data.frame.
#' @param y A data.frame for the second character formatted like \code{x}.
#' @param phy An \code{ephylo} object.
#' @param correlated A boolean indicating if the correlated model (TRUE) or
#' the independent model (FALSE) should be fit.
#' @param par_init Optional initial parameter estimates to seed the algorithm.
#' If not \code{NULL}, any value for \code{correlated} and \code{num_fits} is
#' ignored and a single optimization is run until termination.
#' @param lower A lower bound for rate parameters (>= 0, default 1e-5).
#' @param upper An upper bound for rate parameters (default 1).
#' @param lower_init A lower bound for initial rate parameters (>= 0) used to
#' seed the fitting algorithm (default 0).
#' @param upper_init An upper bound for initial rate parameters used to seed
#' the fitting algorithm (default 0.001).
#' @param ndeps A scalar that multiplies the parameter estimates. The resulting
#' vector is then passed as the \code{ndeps} control parameter when calculating
#' the Hessian matrix. See \code{optim} for details.
#' @param num_fits Return the best (highest log likelihood) of \code{num_fits} 
#' models.
#' @return A list with the following components:
#' \itemize{
#' \item \code{par} The maximum likelihood estimates of the overall rates of
#'   character evolution and character correlation along with their standard
#'   errors.
#' \item \code{value} The log likelihood of the best-fitting model.
#' \item \code{counts} The number of likelihood evaluations made by \code{optim}.
#' \item \code{convergence} The convergence code returned by \code{optim}.
#' \item \code{message} Additional information returned by \code{optim}.
#'}
#' @seealso The \code{vignette("CorrelatedMk")}
mkcor_fit = function(x, y, phy, correlated=TRUE, par_init=NULL, lower=1e-5,
     upper=1, lower_init=0, upper_init=0.001, ndeps=0.3, num_fits=1L)
{
    LIK = mkcor_loglik(x, y, phy)
    if (!is.null(par_init))
    {
        stopifnot(is.numeric(par_init))
        stopifnot(length(par_init) == 2L || length(par_init) == 3L)
        stopifnot(all(par_init >= 0))
        num_pars = length(par_init)
        fit = optim(par_init, LIK, method="L-BFGS-B", lower=lower, upper=upper, 
            control=list(fnscale=-1))
    }
    else
    {
        num_pars = if (correlated) 3L else 2L
        best_score = -Inf
        for (i in 1:num_fits)
        {
            tmp = optim(runif(num_pars, lower_init, upper_init), LIK, 
                method="L-BFGS-B", lower=lower, upper=upper, 
                control=list(fnscale=-1))
            if (tmp$value > best_score)
            {
                best_score = tmp$value
                fit = tmp
            }
        }
    }
    hess = try(optimHess(fit$par, LIK, 
        control=list(fnscale=-1, ndeps=ndeps*fit$par)), silent=TRUE)
    V = try(solve(-hess), silent=TRUE)
    if (!inherits(V, "try-error"))
        cV = try(chol(V), silent=TRUE)
    else
        cV = V
    # asymptotic standard error estimate
    if (!inherits(V, "try-error") && !inherits(cV, "try-error"))
    {
        se = sqrt(diag(V))
    }
    else
    {
        se = rep(NA_real_, num_pars)   
    }
    
    if (num_pars == 3L)
    {
        rx = fit$par[1]
        ry = fit$par[2]
        rxy = fit$par[3]
        r = rxy / sqrt((rx+rxy)*(ry+rxy))
        # asymptotic standard error estimate using delta method
        # https://en.wikipedia.org/wiki/Delta_method#Multivariate_delta_method
        # these are the derivatives of the correlation coefficient w.r.t.
        # \lambda_x (rx), \lambda_y (ry), and \lambda_{xy} (rxy)
        grad = c(
            -(rxy*(ry+rxy))/(2*((rx+rxy)*(ry+rxy))^(3/2)),
            -(rxy*(rx+rxy))/(2*((rx+rxy)*(ry+rxy))^(3/2)),
            (rxy*(rx+ry)+2*rx*ry)/(2*((rx+rxy)*(ry+rxy))^(3/2))
        )
        if (!inherits(V, "try-error") && !inherits(cV, "try-error"))
        {
            se = c(se, c(sqrt(t(grad) %*% V %*% grad)))
        }
        else
        {
            se = c(se, NA_real_)
        }

        fit$par = data.frame(
            par=c(fit$par, r), 
            se=se,
            row.names=c("rate.x", "rate.y", "rate.xy", "corr.xy")
        )
    }
    else
    {
        fit$par = data.frame(
            par=c(fit$par, 0, 0),
            se=c(se, 0, 0),
            row.names=c("rate.x", "rate.y", "rate.xy", "corr.xy")
        )
    }

    return (fit)
}


#' Simulate data under the correlated Mk model
#' 
#' @param par A vector of transition rate parameters. The first two values are
#' rates of independent change in X and Y, respectively. The third 
#' value, if given, is the rate of simultaneous change in X and Y.
#' @param num_states_x The number of states in X.
#' @param num_states_y The number of states in Y.
#' @param phy An \code{ephylo} object.
#' @return A list with named components `x' and `y', which are data.frame
#' objects containing the simulated tip states. Each data.frame has a
#' `node.state' attribute containing the simulated states at internal
#' nodes.
simulate_mkcor = function(par, num_states_x, num_states_y, phy)
{
    stopifnot(inherits(phy, "ephylo"))
    stopifnot(all(par > 0))
    stopifnot(num_states_x > 1)
    stopifnot(num_states_y > 1)
    num_states_x = as.integer(num_states_x)
    num_states_y = as.integer(num_states_y)
    num_states = as.integer(num_states_x * num_states_y)
    num_tips = phy$num.tips
    num_nodes = phy$num.nodes
    preorder = phy$preorder.nodes
    parent = phy$parent
    brlen = phy$brlen
    ans = .Call(
        C_mkcor_simulate,
        num_states,
        num_states_x,
        num_states_y,
        par,
        num_tips,
        num_nodes,
        preorder,
        parent,
        brlen
    )
    x = structure(data.frame(
        tip.label=phy$tip.label, 
        state.id=factor(ans[[1]][1:num_tips])), 
        node.state=ans[[1]][(num_tips+1):num_nodes])
    y = structure(data.frame(
        tip.label=phy$tip.label, 
        state.id=factor(ans[[2]][1:num_tips])),
        node.state=ans[[2]][(num_tips+1):num_nodes])
    list(x=x, y=y)
}
