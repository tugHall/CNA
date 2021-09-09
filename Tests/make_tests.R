
### TESTs for all functions including main function - 'model'


# References --------------------------------------------------------------

##  the library is testthat
##  https://github.com/r-lib/testthat  

# Library -----------------------------------------------------------------


library(testthat)




# Examples ----------------------------------------------------------------


### Vector
test_that("Distinct roots", {
    
    roots <- sqrt( c( 1, 7, 12) )
    
    expect_that( roots, is_a("numeric") )
    expect_that( length(roots), equals(3) )
    expect_that( roots[1] < roots[2], is_true() )
})

### Function
test_that("trigonometric functions match identities", {
    expect_equal(sin(pi / 4), 1 / sqrt(2))
    expect_equal(cos(pi / 4), 1 / sqrt(2))
    expect_equal(tan(pi / 4), 1)
})




### function with Data.Frame as output
test_that("Data.Frame", {
    # manually created data
    dat <- iris[1:5, c("Species", "Sepal.Length")]
    
    # function
    myfun <- function(row, col, data) {
        data[row, col]
    }
    
    # result of applying function
    outdat <- myfun(1:5, c("Species", "Sepal.Length"), iris)
    
    # two versions of the same test
    expect_true(identical(dat, outdat))
    expect_identical(dat, outdat)
    expect_true(identical(dat, outdat))

})



