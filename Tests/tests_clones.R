
### TESTs for all functions including main function - 'model'
### Part II. Tests for clones

# References --------------------------------------------------------------

##  the library is testthat
##  https://github.com/r-lib/testthat  


### LIBRARIES and SOURCES ------------------------------------------------

# library for testing
library(testthat)

###  The simulation uses the functions and classes in the "Code/tugHall_2.1_functions.R" 

library(stringr)   # to use string data in input files
library(actuar)    # to use BIG NUMBERS in N_cell variable

source(file = "Code/tugHall_3.0_functions.R")
source(file = "Code/read_maps.R")


# Generate objects of simulation ------------------------------------------

### Parameters
if ( TRUE ){
    # Probabilities of processes
    
    E0 <<-  1E-4       # parameter in the division probability  
    F0 <<-  10         # parameter in the division probability  
    m0 <<-  1E-6       # mutation probability  
    uo <<-  0.5        # oncogene mutation probability  
    us <<-  0.5        # suppressor mutation probability  
    s0 <<-  10         # parameter in the sigmoid function  
    k0 <<-  0.2        # Environmental death probability  
    d0 <<-  0.35       # Initial probability to divide cells
    ### Additional parameters of simulation
    censore_n <<- 10^5       # Max cell number where the program forcibly stops
    censore_t <<- 100         # Max time where the program forcibly stops
    ### New parameters for CNA:
    m_dup  <<- 8E-8 # mutation probability for duplication
    m_del  <<- 1E-9 # mutation probability for deletion 
    lambda_dup  <<- 5000  # CNA duplication average length (of the geometrical distribution for the length)
    lambda_del  <<- 7000  # CNA deletion average length
    uo_dup  <<- 0.8 # Gene malfunction probability by CNA duplication for oncogene
    us_dup  <<- 0   # Gene malfunction probability by CNA duplication for suppressor
    uo_del  <<- 0   # Gene malfunction probability by CNA deletion    for oncogene
    us_del  <<- 0.8 # Gene malfunction probability by CNA deletion    for suppressor
    
}

gene_map  =  make_map( f_out = './Tests/GENE_MAP/current_test_map.txt', 
                       f_in = './Input/CCDS.current.txt')
### Objects
if ( TRUE ){
    genefile = 'Input/gene_hallmarks.txt'
    clonefile = 'Input/cloneinit.txt' 
    onco = oncogene$new()        # make the vector onco about the hallmarks
    onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
    hall = hallmark$new()        # make a vector hall with hallmarks parameters
    hall$read(genefile, onco$name)     # read from the genefile - 'gene_hallmarks.txt'
    env = environ$new(F0)               # new vector for average values of cells
    pnt = Point_Mutations$new()
    pnt_clones = NULL
    cna = CNA_Mutations$new()
    cna_clones = NULL
    mut_order  <-  0            #  mutation order to reproduce gene map
    assign("mut_order", mut_order, env=.GlobalEnv)
    assign("onco", onco, env=.GlobalEnv)
    assign("env", env, env=.GlobalEnv)
    assign("pnt", pnt, env=.GlobalEnv)
    assign("pnt_clones", pnt_clones, env=.GlobalEnv)
    assign("cna", cna, env=.GlobalEnv)
    assign("cna_clones", cna_clones, env=.GlobalEnv)
    assign("hall", hall, env=.GlobalEnv)
    assign("uo", uo, env=.GlobalEnv)
    assign("us", us, env=.GlobalEnv)
    

    clone1 = clone$new(gene_size=length(onco$cds_1),
                       m=m0, s=s0, k=k0, E=E0)          # clone1  -  empty object of clone
    clones = init_clones(clonefile, clone1)           # clones - the clones with hallmarks from cellfile - cellinit.txt - initial cells  
    onco_clones = init_onco_clones( onco, clones )    # onco_clones - the onco related to each clone in clones
    
}


# Tests for clones --------------------------------------------------------

### The function to get fields of object S4 and 
get_fields_S4  <-  function( objct, nm_objct){
  
  fl = paste0('./Tests/Clones/', nm_objct, '_fields.txt' )
  nms = sort(names( objct ) )
  nms = nms[ (which( nms == '.self') + 1) : length( nms ) ]
  nms = nms[ sapply( X = nms, FUN = function( x ) length( names( as.list( objct$field( x ) ) ) ) == 0 ) ]
  if ( length( which( nms == 'show') ) > 0 )  nms = nms[ -which( nms == 'show') ]
  flds  =  sapply( X = nms, FUN = function(x) objct$field( x ) )
  return( flds )
}

msg  =  'Check the function to generate and initialize an OncoGene object. 
        The functions: new()  and  read() of class oncogene'
cat( paste0(msg, '. \n') )
test_that( paste0( msg, ': \n'), {
    
    objct  =  onco 
    nm_objct  =  'onco'
    flds  =  get_fields_S4( objct, nm_objct )
    fl = paste0('./Tests/Clones/', nm_objct, '_fields.txt' )
    chk_data  =  readRDS( file = fl ) 
    
    # saveRDS( flds, file = fl ) 
    print('The fields of oncoGene object: ')
    print( flds )
    sapply( X = names(flds), FUN = function(x) expect_equal( chk_data[x], flds[x] )  )
})


msg  =  'Check the function to generate and initialize an clone object. 
        The function initialize of class clone'
cat( paste0(msg, '. \n') )
test_that( paste0( msg, ': \n'), {
  
  objct  =  clone1 
  nm_objct  =  'clone'
  flds  =  get_fields_S4( objct, nm_objct )
  fl = paste0('./Tests/Clones/', nm_objct, '_fields.txt' )

  chk_data  =  readRDS( file = fl ) 
  
  # saveRDS( flds, file = fl ) 
  print('The fields of clone object: ')
  print( flds )
  sapply( X = names( flds ), FUN = function(x) expect_equal( chk_data[x], flds[x] )  )
})


msg  =  'Check the function to generate and initialize an environ object. 
        The function initialize of class environ'
cat( paste0(msg, '. \n') )
test_that( paste0( msg, ': \n'), {
  
  objct  =  env 
  nm_objct  =  'env'
  flds  =  get_fields_S4( objct, nm_objct )
  fl = paste0('./Tests/Clones/', nm_objct, '_fields.txt' )
  chk_data  =  readRDS( file = fl ) 
  
  # saveRDS( flds, file = fl ) 

  print('The fields of environ object: ')
  print( flds )
  sapply( X = names( flds ), FUN = function(x) expect_equal( chk_data[x], flds[x] )  )
})



msg  =  'Check the function to generate and initialize an hallmark object. 
        The function initialize of class hallmark'
cat( paste0(msg, '. \n') )
test_that( paste0( msg, ': \n'), {
  
  objct  =  hall 
  nm_objct  =  'hall'
  flds  =  get_fields_S4( objct, nm_objct )
  fl = paste0('./Tests/Clones/', nm_objct, '_fields.txt' )

  chk_data  =  readRDS( file = fl ) 
  
  # saveRDS( flds, file = fl ) 
  print('The fields of hallmark object: ')
  print( flds )
  sapply( X = names( flds ), FUN = function(x) expect_equal( chk_data[x], flds[x] )  )
})





# Examples ----------------------------------------------------------------

if (FALSE){
  source('./Tests/tests_clones.R')
}

if (FALSE){
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
    
}

