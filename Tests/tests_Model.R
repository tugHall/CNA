
### TESTs for function 'model'
### Part III. Tests for only function 'model'

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

Print_output_data  =  FALSE 



# Test for model() --------------------------------------------------------

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
    censore_t <<- 10         # Max time where the program forcibly stops
    ### New parameters for CNA:
    m_dup  <<- 8E-8 # mutation probability for duplication
    m_del  <<- 1E-9 # mutation probability for deletion 
    lambda_dup  <<- 5000  # CNA duplication average length (of the geometrical distribution for the length)
    lambda_del  <<- 7000  # CNA deletion average length
    uo_dup  <<- 0.8 # Gene malfunction probability by CNA duplication for oncogene
    us_dup  <<- 0   # Gene malfunction probability by CNA duplication for suppressor
    uo_del  <<- 0   # Gene malfunction probability by CNA deletion    for oncogene
    us_del  <<- 0.8 # Gene malfunction probability by CNA deletion    for suppressor
    
    genefile = 'Tests/Input/gene_hallmarks.txt'
    clonefile = 'Tests/Input/cloneinit.txt' 
    geneoutfile <- 'Output/geneout.txt'  # Gene Out file with Hallmarks 
    logoutfile <-  'Output/log.txt' 
}

gene_map  =  make_map( f_out = './Tests/GENE_MAP/current_test_map.txt', 
                       f_in = './Tests/Input/CCDS.current.txt')

msg  =  'Check the main function model() for different mutation rates with 10 time-steps for each simulation'
cat( paste0(msg, '. \n') )

rdr  =  -11:-5
for( nos in  1:7 ){
    m0 = 10^rdr[nos]
    
    subDir  =  paste0('./Tests/Model/', nos )
    if (! file.exists(subDir))  dir.create( subDir ) 
    
    # To fix a random choose:
    set.seed(123456)
    cloneoutfile <- paste0( subDir, '/cloneout.txt' )  # output information of simulation
    
    ### Simulation of the cancer cell/clone evolution:
    smlt = model(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0)
    
    # clones      =  smlt[[1]]
    # onco_clones =  smlt[[2]] 
    
    write_pnt_clones( pnt_clones, file_out  =  paste0( subDir, '/point_mutations.txt' )  )
    write_pnt_clones( cna_clones, file_out  =  paste0( subDir, '/CNA_mutations.txt' )  )
    
}


test_that( paste0( msg, ': \n'), {
    
    for( nos in  1:7 ){
        subDir        =  paste0('./Tests/Model/', nos )
        subDir_check  =  paste0('./Tests/Model/Check/', nos )
        for( nm in c( '/point_mutations.txt', '/CNA_mutations.txt', '/cloneout.txt' ) ){ 
            file_check  =  paste0( subDir_check, nm )
            file_data   =  paste0( subDir,       nm )
            
            if ( file.size( file_check ) > 10 ){
                chk_data = read.csv( file = file_check, stringsAsFactors = FALSE,
                                 header = TRUE, sep = '\t' )
            }  else  chk_data  =  NULL
            
            if ( file.size( file_data )  > 10 ){
                data     = read.csv( file = file_data, stringsAsFactors = FALSE,
                                 header = TRUE, sep = '\t' )
            }  else data = NULL
            
            if ( Print_output_data ){
                print('All the fields of clones related to model() function: ')
                print( data )
            }
            expect_equal( chk_data, data, tolerance = 0.00001, 
                          info = paste('Check simulation number ', nos, ' for file: ', nm ) )
        }
    }
})

