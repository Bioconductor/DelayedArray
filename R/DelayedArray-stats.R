### =========================================================================
### Statistical methods for DelayedArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The Normal Distribution
###
### All these methods return a DelayedArray object of the same dimensions
### as their first argument.
###

setMethod("dnorm", "DelayedArray",
    function(x, mean=0, sd=1, log=FALSE)
        stash_DelayedUnaryIsoOp(x, "dnorm",
            Rargs=list(mean=mean, sd=sd, log=log))
)

setMethod("pnorm", "DelayedArray",
    function(q, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE)
        stash_DelayedUnaryIsoOp(q, "pnorm",
            Rargs=list(mean=mean, sd=sd, lower.tail=lower.tail, log.p=log.p))
)

setMethod("qnorm", "DelayedArray",
    function(p, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE)
        stash_DelayedUnaryIsoOp(p, "qnorm",
            Rargs=list(mean=mean, sd=sd, lower.tail=lower.tail, log.p=log.p))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The Binomial Distribution
###
### All these methods return a DelayedArray object of the same dimensions
### as their first argument.
###

setMethod("dbinom", "DelayedArray",
    function(x, size, prob, log=FALSE)
        stash_DelayedUnaryIsoOp(x, "dbinom",
            Rargs=list(size=size, prob=prob, log=log))
)

setMethod("pbinom", "DelayedArray",
    function(q, size, prob, lower.tail=TRUE, log.p=FALSE)
        stash_DelayedUnaryIsoOp(q, "pbinom",
            Rargs=list(size=size, prob=prob,
                       lower.tail=lower.tail, log.p=log.p))
)

setMethod("qbinom", "DelayedArray",
    function(p, size, prob, lower.tail=TRUE, log.p=FALSE)
        stash_DelayedUnaryIsoOp(p, "qbinom",
            Rargs=list(size=size, prob=prob,
                       lower.tail=lower.tail, log.p=log.p))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The Poisson Distribution
###
### All these methods return a DelayedArray object of the same dimensions
### as their first argument.
###

setMethod("dpois", "DelayedArray",
    function(x, lambda, log=FALSE)
        stash_DelayedUnaryIsoOp(x, "dpois",
            Rargs=list(lambda=lambda, log=log))
)

setMethod("ppois", "DelayedArray",
    function(q, lambda, lower.tail=TRUE, log.p=FALSE)
        stash_DelayedUnaryIsoOp(q, "ppois",
            Rargs=list(lambda=lambda, lower.tail=lower.tail, log.p=log.p))
)

setMethod("qpois", "DelayedArray",
    function(p, lambda, lower.tail=TRUE, log.p=FALSE)
        stash_DelayedUnaryIsoOp(p, "qpois",
            Rargs=list(lambda=lambda, lower.tail=lower.tail, log.p=log.p))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The Logistic Distribution
###
### All these methods return a DelayedArray object of the same dimensions
### as their first argument.
###

setMethod("dlogis", "DelayedArray",
    function(x, location=0, scale=1, log=FALSE)
        stash_DelayedUnaryIsoOp(x, "dlogis",
            Rargs=list(location=location, scale=scale, log=log))
)

setMethod("plogis", "DelayedArray",
    function(q, location=0, scale=1, lower.tail=TRUE, log.p=FALSE)
        stash_DelayedUnaryIsoOp(q, "plogis",
            Rargs=list(location=location, scale=scale,
                       lower.tail=lower.tail, log.p=log.p))
)

setMethod("qlogis", "DelayedArray",
    function(p, location=0, scale=1, lower.tail=TRUE, log.p=FALSE)
        stash_DelayedUnaryIsoOp(p, "qlogis",
            Rargs=list(location=location, scale=scale,
                       lower.tail=lower.tail, log.p=log.p))
)

