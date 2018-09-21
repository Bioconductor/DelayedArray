aperm2 <- DelayedArray:::aperm2

test_aperm2 <- function()
{
    a <- array(1:6000, c(1, 40, 1, 50, 3))

    checkException(aperm2(a, 0))
    checkException(aperm2(a, 6))
    checkException(aperm2(a, c(4:1, 4)))
    checkException(aperm2(a, 1:4))

    for (perm in list(1:5, 5:1, c(4, 5, 2, 1, 3)))
        checkIdentical(base::aperm(a, perm), aperm2(a, perm))

    checkIdentical(`dim<-`(a, dim(a)[-1]), aperm2(a, 2:5))
    checkIdentical(`dim<-`(a, dim(a)[-3]), aperm2(a, c(1:2, 4:5)))
    checkIdentical(`dim<-`(a, dim(a)[-c(1, 3)]), aperm2(a, c(2, 4:5)))

    target <- base::aperm(`dim<-`(a, dim(a)[-1]), 4:1)
    checkIdentical(target, aperm2(a, 5:2))

    target <- base::aperm(`dim<-`(a, dim(a)[-3]), c(4, 2, 1, 3))
    checkIdentical(target, aperm2(a, c(5, 2, 1, 4)))

    current <- aperm2(a, c(NA, NA, 5, 2, NA, 1, 4, NA))
    target <- `dim<-`(target, c(1, 1, 3, 40, 1, 1, 50, 1))
    checkIdentical(target, current)
    checkIdentical(drop(target), drop(current))
}

