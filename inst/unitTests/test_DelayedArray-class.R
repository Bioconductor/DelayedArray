#setAutoRealizationBackend("RleArray")
#setAutoRealizationBackend("HDF5Array")

test_DelayedArray_constructor <- function()
{
    check_DelayedMatrix <- function(m, M) {
        checkTrue(is(M, "DelayedMatrix"))
        checkTrue(validObject(M, complete=TRUE))
        checkIdentical(typeof(m), type(M))
        checkIdentical(typeof(m), type(M[ , -2]))
        checkIdentical(typeof(m), type(M[ , 0]))
        checkIdentical(typeof(m), type(M[0 , ]))
        checkIdentical(dim(m), dim(M))
        checkIdentical(rownames(m), rownames(M))
        checkIdentical(colnames(m), colnames(M))
    }

    ## 6-col seed

    DF1 <- DataFrame(aa=1:5,
                     bb=letters[1:5],
                     cc=Rle(c(-0.25, 7.1), 3:2),
                     dd=factor(c("E", "B", "A", "A", "B"), levels=LETTERS),
                     ee=Rle(factor(c("d", "b"), levels=letters[1:5]), 2:3),
                     ff=c(FALSE, TRUE, TRUE, FALSE, FALSE))

    df1 <- as.data.frame(DF1)
    m1 <- as.matrix(df1)      # character matrix

    M1a <- DelayedArray(m1)   # matrix seed
    check_DelayedMatrix(m1, M1a)
    checkIdentical(m1, as.matrix(M1a))
    checkIdentical(m1, as.array(M1a))

    M1b <- DelayedArray(df1)  # data.frame seed
    rownames(M1b) <- NULL
    check_DelayedMatrix(m1, M1b)

    M1c <- DelayedArray(DF1)  # DataFrame seed
    check_DelayedMatrix(m1, M1c)

    ## 5-col seed

    DF2 <- DF1[ , -2, drop=FALSE]
    df2 <- as.data.frame(DF2)
    m2 <- as.matrix(df2)      # character matrix

    M2a <- DelayedArray(m2)   # matrix seed
    check_DelayedMatrix(m2, M2a)
    checkIdentical(m2, as.matrix(M2a))
    checkIdentical(m2, as.array(M2a))

    M2b <- DelayedArray(df2)  # data.frame seed
    rownames(M2b) <- NULL
    check_DelayedMatrix(m2, M2b)

    M2c <- DelayedArray(DF2)  # DataFrame seed
    check_DelayedMatrix(m2, M2c)
    
    ## 4-col seed

    DF3 <- DF1[ , -c(2, 4), drop=FALSE]
    df3 <- as.data.frame(DF3)
    m3 <- as.matrix(df3)      # character matrix

    M3a <- DelayedArray(m3)   # matrix seed
    check_DelayedMatrix(m3, M3a)
    checkIdentical(m3, as.matrix(M3a))
    checkIdentical(m3, as.array(M3a))

    M3b <- DelayedArray(df3)  # data.frame seed
    rownames(M3b) <- NULL
    check_DelayedMatrix(m3, M3b)

    M3c <- DelayedArray(DF3)  # DataFrame seed
    check_DelayedMatrix(m3, M3c)

    ## 3-col seed

    DF4 <- DF1[ , -c(2, 4, 5), drop=FALSE]
    df4 <- as.data.frame(DF4)
    m4 <- as.matrix(df4)      # double matrix

    M4a <- DelayedArray(m4)   # matrix seed
    check_DelayedMatrix(m4, M4a)
    checkIdentical(m4, as.matrix(M4a))
    checkIdentical(m4, as.array(M4a))

    M4b <- DelayedArray(df4)  # data.frame seed
    rownames(M4b) <- NULL
    check_DelayedMatrix(m4, M4b)
    checkIdentical(m4, as.matrix(M4b))
    checkIdentical(m4, as.array(M4b))

    M4c <- DelayedArray(DF4)  # DataFrame seed
    check_DelayedMatrix(m4, M4c)
    checkIdentical(m4, as.matrix(M4c))
    checkIdentical(m4, as.array(M4c))

    ## 2-col seed

    DF5 <- DF1[ , c(1, 6), drop=FALSE]
    df5 <- as.data.frame(DF5)
    m5 <- as.matrix(df5)      # integer matrix

    M5a <- DelayedArray(m5)   # matrix seed
    check_DelayedMatrix(m5, M5a)
    checkIdentical(m5, as.matrix(M5a))
    checkIdentical(m5, as.array(M5a))

    M5b <- DelayedArray(df5)  # data.frame seed
    rownames(M5b) <- NULL
    check_DelayedMatrix(m5, M5b)
    checkIdentical(m5, as.matrix(M5b))
    checkIdentical(m5, as.array(M5b))

    M5c <- DelayedArray(DF5)  # DataFrame seed
    check_DelayedMatrix(m5, M5c)
    checkIdentical(m5, as.matrix(M5c))
    checkIdentical(m5, as.array(M5c))

    ## 1-col seed

    DF6 <- DF1[ , 6, drop=FALSE]
    df6 <- as.data.frame(DF6)
    m6 <- as.matrix(df6)      # logical matrix

    M6a <- DelayedArray(m6)   # matrix seed
    check_DelayedMatrix(m6, M6a)
    checkIdentical(m6, as.matrix(M6a))
    checkIdentical(m6, as.array(M6a))

    M6b <- DelayedArray(df6)  # data.frame seed
    rownames(M6b) <- NULL
    check_DelayedMatrix(m6, M6b)
    checkIdentical(m6, as.matrix(M6b))
    checkIdentical(m6, as.array(M6b))

    M6c <- DelayedArray(DF6)  # DataFrame seed
    check_DelayedMatrix(m6, M6c)
    checkIdentical(m6, as.matrix(M6c))
    checkIdentical(m6, as.array(M6c))
}

test_DelayedArray_subsetting <- function()
{
    a <- array(runif(78000), dim=c(600, 26, 5),
                             dimnames=list(NULL, LETTERS, letters[1:5]))
    A <- realize(a)
    checkTrue(is(A, "DelayedArray"))
    checkTrue(validObject(A, complete=TRUE))

    checkIdentical(a, as.array(A[ , , ]))
    checkException(A[ , 27, ], silent=TRUE)
    checkException(A[ , , "f"], silent=TRUE)

    target <- a[0, , ]
    checkIdentical(target, as.array(A[0, , ]))
    checkIdentical(target, as.array(A[integer(0), , ]))
    checkIdentical(target, as.array(A[NULL, , ]))

    target <- a[ , 0, ]
    checkIdentical(target, as.array(A[ , 0, ]))
    checkIdentical(target, as.array(A[ , integer(0), ]))
    checkIdentical(target, as.array(A[ , NULL, ]))

    target <- a[ , , 0]
    checkIdentical(target, as.array(A[ , , 0]))
    checkIdentical(target, as.array(A[ , , integer(0)]))
    checkIdentical(target, as.array(A[ , , NULL]))

    target <- a[ , 0, 0]
    checkIdentical(target, as.array(A[ , 0, 0]))
    checkIdentical(target, as.array(A[ , integer(0), integer(0)]))
    checkIdentical(target, as.array(A[ , NULL, NULL]))

    target <- a[0, , 0]
    checkIdentical(target, as.array(A[0, , 0]))
    checkIdentical(target, as.array(A[integer(0), , integer(0)]))
    checkIdentical(target, as.array(A[NULL, , NULL]))

    target <- a[0, 0, ]
    checkIdentical(target, as.array(A[0, 0, ]))
    checkIdentical(target, as.array(A[integer(0), integer(0), ]))
    checkIdentical(target, as.array(A[NULL, NULL, ]))

    target <- a[0, 0, 0]
    checkIdentical(target, as.array(A[0, 0, 0]))
    checkIdentical(target, as.array(A[integer(0), integer(0), integer(0)]))
    checkIdentical(target, as.array(A[NULL, NULL, NULL]))

    i <- c(FALSE, TRUE)
    target <- a[i, , ]
    checkIdentical(target, as.array(A[i, , ]))
    target <- a[i, 0, ]
    checkIdentical(target, as.array(A[i, 0, ]))
    checkIdentical(target, as.array(A[i, integer(0), ]))
    checkIdentical(target, as.array(A[i, NULL, ]))
    i <- c(FALSE, TRUE, TRUE, FALSE, FALSE)
    target <- a[i, , ]
    checkIdentical(target, as.array(A[i, , ]))
    target <- a[i, 0, ]
    checkIdentical(target, as.array(A[i, 0, ]))
    checkIdentical(target, as.array(A[i, integer(0), ]))
    checkIdentical(target, as.array(A[i, NULL, ]))

    k <- c(TRUE, FALSE)
    target <- a[ , , k]
    checkIdentical(target, as.array(A[ , , k]))
    target <- a[0, , k]
    checkIdentical(target, as.array(A[0, , k]))
    checkIdentical(target, as.array(A[integer(0), , k]))
    checkIdentical(target, as.array(A[NULL, , k]))
    target <- a[i, , k]
    checkIdentical(target, as.array(A[i, , k]))
    target <- a[i, 0, k]
    checkIdentical(target, as.array(A[i, 0, k]))
    checkIdentical(target, as.array(A[i, integer(0), k]))
    checkIdentical(target, as.array(A[i, NULL, k]))

    j <- c(20:5, 11:22)
    target <- a[0, j, ]
    checkIdentical(target, as.array(A[0, j, ]))
    checkIdentical(target, as.array(A[integer(0), j, ]))
    checkIdentical(target, as.array(A[NULL, j, ]))
    target <- a[0, j, 0]
    checkIdentical(target, as.array(A[0, j, 0]))
    checkIdentical(target, as.array(A[integer(0), j, integer(0)]))
    checkIdentical(target, as.array(A[NULL, j, NULL]))
    target <- a[99:9, j, ]
    checkIdentical(target, as.array(A[99:9, j, ]))

    k <- c("e", "b", "b", "d", "e", "c")
    target <- a[-3, j, k]
    checkIdentical(target, as.array(A[-3, j, k]))
    target <- a[-3, -j, k]
    checkIdentical(target, as.array(A[-3, -j, k]))

    target <- a[-3, -j, ]
    checkIdentical(target, as.array(A[-3, -j, ]))
    target <- a[-3, -j, -555555]
    checkIdentical(target, as.array(A[-3, -j, -555555]))

    target <- a[99:9, j, k]
    checkIdentical(target, as.array(A[99:9, j, k]))
    i <- c(150:111, 88:90, 75:90)
    target <- a[i, j, k]
    checkIdentical(target, as.array(A[i, j, k]))

    target1 <- a[99, j, k]
    checkIdentical(target1, as.matrix(A[99, j, k]))
    B1 <- A[99, j, k]
    checkIdentical(target1, as.array(B1))
    target <- target1[15:8, c("c", "b")]
    checkIdentical(target, as.array(B1[15:8, c("c", "b")]))
    target <- target1[0, ]
    checkIdentical(target, as.array(B1[0, ]))
    target <- target1[ , 0]
    checkIdentical(target, as.array(B1[ , 0]))
    target <- target1[0, 0]
    checkIdentical(target, as.array(B1[0, 0]))

    target2 <- a[i, 22, k]
    checkIdentical(target2, as.matrix(A[i, 22, k]))
    B2 <- A[i, 22, k]
    checkIdentical(target2, as.array(B2))
    target <- target2[15:8, c("c", "b")]
    checkIdentical(target, as.array(B2[15:8, c("c", "b")]))
    target <- target2[0, ]
    checkIdentical(target, as.array(B2[0, ]))
    target <- target2[ , 0]
    checkIdentical(target, as.array(B2[ , 0]))
    target <- target2[0, 0]
    checkIdentical(target, as.array(B2[0, 0]))

    target3 <- a[i, j, 5]
    checkIdentical(target3, as.matrix(A[i, j, 5]))
    B3 <- A[i, j, 5]
    checkIdentical(target3, as.array(B3))
    target <- target3[15:8, LETTERS[18:12]]
    checkIdentical(target, as.array(B3[15:8, LETTERS[18:12]]))
    target <- target3[0, ]
    checkIdentical(target, as.array(B3[0, ]))
    target <- target3[ , 0] 
    checkIdentical(target, as.array(B3[ , 0]))
    target <- target3[0, 0] 
    checkIdentical(target, as.array(B3[0, 0]))
}

