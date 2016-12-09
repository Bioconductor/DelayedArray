### =========================================================================
### Statistical methods for DelayedArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Logistic
###
### All these methods return a DelayedArray object of the same dimensions
### as their first argument.
###

setMethod("dlogis", "DelayedArray",
    function(x, location=0, scale=1, log=FALSE)
        register_delayed_op(x, "dlogis",
            Rargs=list(location=location, scale=scale, log=log))
)

setMethod("plogis", "DelayedArray",
    function(q, location=0, scale=1, lower.tail=TRUE, log.p=FALSE)
        register_delayed_op(q, "plogis",
            Rargs=list(location=location, scale=scale,
                       lower.tail=lower.tail, log.p=log.p))
)

setMethod("qlogis", "DelayedArray",
    function(p, location=0, scale=1, lower.tail=TRUE, log.p=FALSE)
        register_delayed_op(p, "qlogis",
            Rargs=list(location=location, scale=scale,
                       lower.tail=lower.tail, log.p=log.p))
)

