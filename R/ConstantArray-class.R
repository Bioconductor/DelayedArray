setClass("ConstantArray",
    contains="DelayedArray",
    representation(seed="ConstantArraySeed")
)

setClass("ConstantMatrix",
    contains="DelayedMatrix",
    representation(seed="ConstantArraySeed")
)

setMethod("matrixClass", "ConstantArray", function(x) "ConstantMatrix")

setMethod("DelayedArray", "ConstantArraySeed", function(seed) new_DelayedArray(seed, Class="ConstantMatrix"))

ConstantArray <- function(dim, value) {
    DelayedArray(ConstantArraySeed(dim, value))
}
