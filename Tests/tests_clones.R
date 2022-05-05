
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

library(binaryLogic)  # to convert integer to bunary

source(file = "Code/tugHall_3.0_functions.R")
source(file = "Code/read_maps.R")

if ( !exists('Print_output_data') )  {
    pk = readline(prompt="Print output data of the tests? Press key: 1 - Yes, any others - No. ")
    if ( pk == '1' ) Print_output_data  =  TRUE else Print_output_data  =  FALSE
}

# Generate objects of simulation ------------------------------------------
define_files_names( sbdr_Input = '/Tests/Input', sbdr_Output = '/Tests/Output' )   
define_gene_location()
define_paramaters( read_fl = TRUE , file_name = './Tests/Input/parameters.txt' )
define_compaction_factor( read_fl = TRUE , file_name = './Tests/Input/CF.txt' )
print( 'Tested parameters are obtained from /Tests/Input/ folder.' )

# print_parameters()

# Define trial() function: trial_complex or trial_simple
if ( model_name != 'simplified' ){
  trial  =  trial_complex
} else {
  trial  =  trial_simple
}

### Objects
if ( TRUE ){
    genefile = 'Tests/Input/gene_hallmarks.txt'
    clonefile = 'Tests/Input/cloneinit.txt' 
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
    if ( Print_output_data ){
        print('The fields of oncoGene object: ')
        print( flds )
    }
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
    if ( Print_output_data ){
        print('The fields of clone object: ')
        print( flds )
    }
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
    if ( Print_output_data ){
        print('The fields of environ object: ')
        print( flds )
    }
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
    if ( Print_output_data ){
        print('The fields of hallmark object: ')
        print( flds )
    }
    sapply( X = names( flds ), FUN = function(x) expect_equal( chk_data[x], flds[x] )  )
})


# Functions of calculation for clones -------------------------------------

### UPDATE functions

### I)  Update Hallmarks
# it's used in the code as:
# lapply(clones,update_Hallmarks)

for ( i in 0:15 ){
    clone1 = clone$new(gene_size=length(onco$cds_1),
                        m=m0, s=s0, k=k0, E=E0) 
    b = as.binary( i , n = 4 )
    for( l in 1:4 ) if ( b[ l ] ) clone1$gene[ l ] = 1 else clone1$gene[ l ] = 0
    clones[[ i + 2 ]]  =  clone1
}

write_clones_hallmarks <- function( outfile = './Tests/Clones/Update_hallmarks.txt', 
                                    append = FALSE, clones, mark = 'Before_update' ) {
    data <- c( mark, 'Nmax', 'Ha', 'Him', 'Hi', 'Hd', 'Hb', 'mutden',
                      'driver_genes', 'pass_genes' )
    write(data, outfile, append = append, ncolumn=length(data), sep="\t")
    
    if (length(clones) > 0 ) {
        for ( i in 1:length( clones ) ) {
            clone1 = clones[[i]]
            data <- c(' - ', clone1$Nmax, clone1$Ha, clone1$Him, clone1$Hi, 
                      clone1$Hd, clone1$Hb, clone1$mutden, 
                      paste(clone1$gene, collapse =  ' '), paste(clone1$pasgene, collapse =  ' ') )
            write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
        }
    }
}

msg  =  'Check the function to calculate hallmarks values for driver genes.'
cat( paste0(msg, '. \n') )
test_that( paste0( msg, ': \n'), {
    
    fl = './Tests/Clones/Update_hallmarks.txt'
    write_clones_hallmarks( outfile = fl, append = FALSE,
                                        clones, mark = 'Before_update' )
    lapply(clones,update_Hallmarks)
    write_clones_hallmarks( outfile = fl, append = TRUE,
                            clones, mark = 'After_update' )
    
    chk_data = read.csv( file = './Tests/Clones/Update_hallmarks_check.txt', 
                         header = TRUE, sep = '\t' )
    
    data = read.csv( file = fl, header = TRUE, sep = '\t' )
    
    if ( Print_output_data ){
        print('The fields of clones related to hallmark: ')
        print( data )
    }
    expect_equal( chk_data, data )
})

clones = vector(mode = "list", length = 16)
for ( i in 0:15 ){
  clone1 = clone$new(gene_size=length(onco$cds_1),
                     m=m0, s=s0, k=k0, E=E0) 
  b = as.binary( i , n = 4 )
  for( l in 1:4 ) if ( b[ l ] ) clone1$pasgene[ l ] = 1 else clone1$pasgene[ l ] = 0
  clones[[ i+1 ]]  =  clone1
}

msg  =  'Check the function to calculate hallmarks values for passenger genes.'
cat( paste0(msg, '. \n') )
test_that( paste0( msg, ': \n'), {
  
  fl = './Tests/Clones/Update_hallmarks_passenger.txt'
  write_clones_hallmarks( outfile = fl, append = FALSE,
                          clones, mark = 'Before_update' )
  lapply(clones,update_Hallmarks)
  write_clones_hallmarks( outfile = fl, append = TRUE,
                          clones, mark = 'After_update' )
  
  
  chk_data = read.csv( file = './Tests/Clones/Update_hallmarks_passenger_check.txt', 
                       header = TRUE, sep = '\t' )
  
  data = read.csv( file = fl, header = TRUE, sep = '\t' )
  
  
  
  if ( Print_output_data ){
    print('The fields of clones related to hallmark: ')
    print( data )
  }
  expect_equal( chk_data, data )
})

# hall$updateEnviron(env, clones) 


### II) TRIAL() function
msg  =  'Check the function to calculate trial() function for passenger and driver genes'
cat( paste0(msg, '. \n') )

write_clones_trial <- function( outfile = './Tests/Clones/trial_clones.txt', N_clones_new, 
                                    append = TRUE, clones, onco_clones, mark = 'Before_trial' ) {
  data <- c( mark, 'N_new_clones', 'mutation rate', 'p0', 'N_cells', 'driver_genes', 'pass_genes', 'd', 'Hd' )
  write(data, outfile, append = append, ncolumn=length(data), sep="\t")
  
  if (length(clones) > 0 ) {
    for ( i in 1:length( clones ) ) {
      clone1 = clones[[ i ]]
      onco1  = onco_clones[[ i ]]
      data <- c(' - ', N_clones_new[[ i ]], clone1$m, onco1$p0_1, clone1$N_cells,  paste(clone1$gene, collapse =  ' '), 
                paste(clone1$pasgene, collapse =  ' ') , clone1$d, clone1$Hd )
      write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
    }
  }
}
fl  =  './Tests/Clones/trial_clones.txt'
if (file.exists(fl)) {
  #Delete file if it exists
  file.remove(fl)
}

for( rdr in -10:-2 ){
    m0 = 10^rdr
    
    # To fix a random choose:
    set.seed(123456)
    
    onco = oncogene$new()        # make the vector onco about the hallmarks
    onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
    
    clones = vector(mode = "list", length = 16)
    for ( i in 0:15 ){
      clone1 = clone$new(gene_size=length(onco$cds_1),
                         m=m0, s=s0, k=k0, E=E0) 
      b = as.binary( i , n = 4 )
      for( l in 1:4 ) if ( b[ l ] ) clone1$pasgene[ l ] = 1 else clone1$pasgene[ l ] = 0
      clones[[ i+1 ]]  =  clone1
    }
    # generate the same onco for all clones:
    onco_clones = init_onco_clones( onco, clones )
    # sapply(X = 1:16, FUN = function(x) ( onco_clones[[x]]$p0_1 ) )
    ### Trial without mutated genes (only passenger genes)
    
    cells_number <- sum_N_P_M(env, clones)                 # to calculate cells numbers - N,M 
    lapply(clones,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    hall$updateEnviron(env, clones) 
    # sapply(X = 1:16, FUN = function(x) ( clones[[x]]$N_cells ) )
    write_clones_trial( outfile = fl,  N_clones_new = rep(0,length(clones)), append = TRUE, 
                        clones, onco_clones, mark = 'Before_trial' )
    N_clones_new = unlist( mapply( trial, clones, onco_clones ) ) 
    # sapply(X = 1:16, FUN = function(x) ( clones[[x]]$N_cells ) )
    write_clones_trial( outfile = fl, N_clones_new, append = TRUE, 
                        clones, onco_clones, mark = 'After_trial' )
    
    for ( i in 0:15 ){
      clone1 = clone$new(gene_size=length(onco$cds_1),
                         m=m0, s=s0, k=k0, E=E0) 
      b = as.binary( i , n = 4 )
      for( l in 1:4 ) if ( b[ l ] ) clone1$gene[ l ] = 1 else clone1$gene[ l ] = 0
      clones[[ i + 1 ]]  =  clone1
    }
    cells_number <- sum_N_P_M(env, clones)             # to calculate cells numbers - N,M 
    lapply(clones,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    hall$updateEnviron(env, clones) 
    #sapply(X = 1:16, FUN = function(x) ( clones[[x]]$N_cells ) )
    write_clones_trial( outfile = fl,  N_clones_new = rep(0,length(clones)), append = TRUE, 
                        clones, onco_clones, mark = 'Before_trial' )
    N_clones_new = unlist( mapply( trial, clones, onco_clones ) ) 
    #sapply(X = 1:16, FUN = function(x) ( clones[[x]]$N_cells ) )
    write_clones_trial( outfile = fl, N_clones_new, append = TRUE, 
                        clones, onco_clones, mark = 'After_trial' )
    
}

test_that( paste0( msg, ': \n'), {
  
  chk_data = read.csv( file = './Tests/Clones/trial_clones_check.txt', stringsAsFactors = FALSE,
                       header = TRUE, sep = '\t' )
  
  data     = read.csv( file = './Tests/Clones/trial_clones.txt', stringsAsFactors = FALSE,
                      header = TRUE, sep = '\t' )

  if ( Print_output_data ){
    print('The fields of clones related to trial() function: ')
    print( data )
  }
  expect_equal( chk_data, data )
})


### III) trial_mutagenesis()  function
# The number of mutations for each NEW clone is num_mut

msg  =  'Check the function to calculate trial_mutagenesis() function for clones'
cat( paste0(msg, '. \n') )

fl  =  './Tests/Clones/trial_mutagenesis.txt'
if (file.exists(fl)) {
  #Delete file if it exists
  file.remove(fl)
}

onco = oncogene$new()        # make the vector onco about the hallmarks
onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
env = environ$new(F0)  
write_header( fl, env, onco ) 

for( rdr in -10:-2 ){
    m0 = 10^rdr
    
    # To fix a random choose:
    set.seed(123456)
    
    env = environ$new( F0 )   ### Again
    onco = oncogene$new( )        # make the vector onco about the hallmarks
    onco$read( genefile )          # read the input info to the onco from genefile - 'gene_cds2.txt'
    
    pnt = Point_Mutations$new()
    pnt_clones = NULL
    cna = CNA_Mutations$new()
    cna_clones = NULL
    
    clones = vector(mode = "list", length = 16)
    for ( i in 0:15 ){
        clone1 = clone$new(gene_size=length(onco$cds_1),
                           m=m0, s=s0, k=k0, E=E0) 
        b = as.binary( i , n = 4 )
        for( l in 1:4 ) if ( b[ l ] ) clone1$pasgene[ l ] = 1 else clone1$pasgene[ l ] = 0
        clones[[ i+1 ]]  =  clone1
    }
    # generate the same onco for all clones:
    onco_clones = init_onco_clones( onco, clones )
    
    cells_number <- sum_N_P_M(env, clones)             # to calculate cells numbers - N,M 
    lapply(clones,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    hall$updateEnviron(env, clones) 
    
    write_cloneout( fl, env, clones, isFirst = TRUE, onco_clones )
    
    for ( nn in 1:length( clones ) )  {
        trial_mutagenesis( clones[[nn]], num_mut = 5, onco_clones[[nn]]  )
    }
    env$T = env$T  +  1
    lapply(clones,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    hall$updateEnviron(env, clones) 
    write_cloneout( fl, env, clones, isFirst = TRUE, onco_clones )
    
    for ( i in 0:15 ){
        clone1 = clone$new(gene_size=length(onco$cds_1),
                           m=m0, s=s0, k=k0, E=E0) 
        b = as.binary( i , n = 4 )
        for( l in 1:4 ) if ( b[ l ] ) clone1$gene[ l ] = 1 else clone1$gene[ l ] = 0
        clones[[ i + 1 ]]  =  clone1
    }
    cells_number <- sum_N_P_M(env, clones)              # to calculate cells numbers - N,M 
    lapply(clones,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    hall$updateEnviron(env, clones) 
    
    write_cloneout( fl, env, clones, isFirst = TRUE, onco_clones )
    
    for ( nn in 1:length( clones ) )  {
      trial_mutagenesis( clones[[nn]], num_mut = 5, onco_clones[[nn]]  )
    }
    env$T = env$T  +  1
    lapply(clones,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    hall$updateEnviron(env, clones) 
    
    write_cloneout( fl, env, clones, isFirst = TRUE, onco_clones )
}

test_that( paste0( msg, ': \n'), {
    
    chk_data = read.csv( file = './Tests/Clones/trial_mutagenesis_check.txt', stringsAsFactors = FALSE,
                         header = TRUE, sep = '\t' )
    
    data     = read.csv( file = './Tests/Clones/trial_mutagenesis.txt', stringsAsFactors = FALSE,
                         header = TRUE, sep = '\t' )
    
    if ( Print_output_data ){
        print('All the fields of clones related to trial_mutagenesis() function: ')
        print( data )
    }
    expect_equal( chk_data, data, tolerance = 0.00001 )
})




