setMethod("t", "DelayedArray", S4Arrays::t.Array)
colsum <- S4Arrays::colsum
makeNindexFromArrayViewport <- S4Arrays::makeNindexFromArrayViewport
DummyArrayGrid <- S4Arrays::DummyArrayGrid
RegularArrayGrid <- S4Arrays::RegularArrayGrid
ArbitraryArrayGrid <- S4Arrays::ArbitraryArrayGrid
extract_array <- S4Arrays::extract_array
is_sparse <- S4Arrays::is_sparse
read_block <- S4Arrays::read_block
write_block <- S4Arrays::write_block

### Non-exported stuff (but used in some DelayedArray rev deps).
get_Nindex_lengths <- S4Arrays:::get_Nindex_lengths
set_dim <- S4Arrays:::set_dim
set_dimnames <- S4Arrays:::set_dimnames
subset_by_Nindex <- S4Arrays:::subset_by_Nindex
to_linear_index <- S4Arrays:::to_linear_index
bplapply2 <- S4Arrays:::bplapply2
