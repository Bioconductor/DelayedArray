setClass("ConstantArray",
    contains="DelayedArray",
    representation(seed="ConstantArraySeed")
)

setClass("ConstantMatrix", contains=c("ConstantArray", "DelayedMatrix"))

setMethod("matrixClass", "ConstantArray", function(x) "ConstantMatrix")

setMethod("DelayedArray", "ConstantArraySeed", function(seed) new_DelayedArray(seed, Class="ConstantArray"))

ConstantArray <- function(dim, value) {
    DelayedArray(ConstantArraySeed(dim, value))
}

# Copied from the corresponding methods for RleArray/RleMatrix.
setAs("ConstantArray", "ConstantMatrix", function(from) new2("ConstantMatrix", from))

setAs("ConstantMatrix", "ConstantArray", function(from) from) 
