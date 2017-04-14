#setRealizationBackend("RleArray")
#setRealizationBackend("HDF5Array")

DEFAULT_BLOCK_SIZE <- DelayedArray:::DEFAULT_BLOCK_SIZE

Arith_members <- c("+", "-", "*", "/", "^", "%%", "%/%")
Compare_members <- c("==", "!=", "<=", ">=", "<", ">")
Logic_members <- c("&", "|")  # currently untested

a1 <- array(sample(5L, 150, replace=TRUE), c(5, 10, 3))  # integer array
a2 <- a1 + runif(150) - 0.5                              # numeric array
m2 <- matrix(runif(60), ncol=6)                          # numeric matrix

block_sizes1 <- c(12L, 20L, 50L, 15000L)
block_sizes2 <- 2L * block_sizes1

test_DelayedMatrix_Ops <- function()
{
    test_delayed_Ops_on_matrix <- function(.Generic, m, M) {
        GENERIC <- match.fun(.Generic)

        target_current <- list(
            list(GENERIC(m, m[ , 1]), GENERIC(M, M[ , 1])),
            list(GENERIC(m[ , 2], m), GENERIC(M[ , 2], M))
        )
        for (i in seq_along(target_current)) {
            target <- target_current[[i]][[1L]]
            current <- target_current[[i]][[2L]]
            checkIdentical(target, as.matrix(current))
            checkIdentical(t(target), as.matrix(t(current)))
            checkIdentical(target[-2, 8:5], as.matrix(current[-2, 8:5]))
            checkIdentical(t(target[-2, 8:5]), as.matrix(t(current[-2, 8:5])))
            checkIdentical(target[-2, 0], as.matrix(current[-2, 0]))
            checkIdentical(t(target[-2, 0]), as.matrix(t(current[-2, 0])))
            checkIdentical(target[0, ], as.matrix(current[0, ]))
            checkIdentical(t(target[0, ]), as.matrix(t(current[0, ])))
        }

        target_current <- list(
            list(GENERIC(t(m), 8:-1), GENERIC(t(M), 8:-1)),
            list(GENERIC(8:-1, t(m)), GENERIC(8:-1, t(M))),

            list(GENERIC(t(m), m[1 , ]), GENERIC(t(M), M[1 , ])),
            list(GENERIC(m[2 , ], t(m)), GENERIC(M[2 , ], t(M))),

            list(GENERIC(t(m), m[1 , 6:10]), GENERIC(t(M), M[1 , 6:10])),
            list(GENERIC(m[2 , 8:7], t(m)), GENERIC(M[2 , 8:7], t(M)))
        )
        for (i in seq_along(target_current)) {
            target <- target_current[[i]][[1L]]
            current <- target_current[[i]][[2L]]
            checkIdentical(target, as.matrix(current))
            checkIdentical(target[1:3 , ], as.matrix(current[1:3 , ]))
            checkIdentical(target[ , 1:3], as.matrix(current[ , 1:3]))
            checkIdentical(t(target), as.matrix(t(current)))
            checkIdentical(t(target)[1:3 , ], as.matrix(t(current)[1:3 , ]))
            checkIdentical(t(target)[ , 1:3], as.matrix(t(current)[ , 1:3]))
            checkIdentical(target[8:5, -2], as.matrix(current[8:5, -2]))
            checkIdentical(t(target[8:5, -2]), as.matrix(t(current[8:5, -2])))
            checkIdentical(target[0, -2], as.matrix(current[0, -2]))
            checkIdentical(t(target[0, -2]), as.matrix(t(current[0, -2])))
            checkIdentical(target[ , 0], as.matrix(current[ , 0]))
            checkIdentical(t(target[ , 0]), as.matrix(t(current[ , 0])))
        }
    }

    a <- a2
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA

    toto <- function(x) t((5 * x[ , 1:2] ^ 3 + 1L) * log(x)[, 10:9])[ , -1]

    m <- a[ , , 2]
    M <- realize(m)
    checkIdentical(toto(m), as.array(toto(M)))
    ## "Logic" members currently untested.
    for (.Generic in c(Arith_members, Compare_members))
        test_delayed_Ops_on_matrix(.Generic, m, M)

    A <- realize(a)[ , , 2]
    M <- drop(A)
    checkIdentical(toto(m), as.array(toto(M)))
    for (.Generic in c(Arith_members, Compare_members))
        test_delayed_Ops_on_matrix(.Generic, m, M)
}

test_DelayedMatrix_mult <- function()
{
    m <- m2
    m[2, 4] <- NA
    m[5, 4] <- Inf
    m[6, 3] <- -Inf
    M <- realize(m)

    Lm <- rbind(rep(1L, 10), rep(c(1L, 0L), 5), rep(-100L, 10))
    Rm <- rbind(Lm + 7.05, 0.1 * Lm)

    on.exit(options(DelayedArray.block.size=DEFAULT_BLOCK_SIZE))
    for (block_size in block_sizes2) {
        options(DelayedArray.block.size=block_size)
        P <- Lm %*% M
        checkEquals(Lm %*% m, as.matrix(P))
        P <- M %*% Rm
        checkEquals(m %*% Rm, as.matrix(P))
    }
}

