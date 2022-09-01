
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

if ( !exists('Print_output_data') )  {
    pk = readline(prompt="Print output data of the tests? Press key: 1 - Yes, any others - No. ")
    if ( pk == '1' ) Print_output_data  =  TRUE else Print_output_data  =  FALSE
}



# Test for model() --------------------------------------------------------
define_files_names( sbdr_Input = '/Tests/Input', sbdr_Output = '/Tests/Output' )   
define_gene_location()
define_paramaters( read_fl = TRUE , file_name = './Tests/Input/parameters.txt' )
define_compaction_factor( read_fl = TRUE , file_name = './Tests/Input/CF.txt' )
print( 'Tested parameters are obtained from /Tests/Input/ folder.' )
# print_parameters()

# To accelerate the calculations:
censor_time_step  <<- 10 

# Define trial() function: trial_complex or trial_simple
if ( model_name != 'simplified' ){
    trial  =  trial_complex
} else {
    trial  =  trial_simple
}


msg  =  'Check the main function model() for different mutation rates with 10 time-steps for each simulation'
cat( paste0(msg, '. \n') )
msg  =  'as well as VAF calculations and order of genes dysfunctions for each simulation'
cat( paste0(msg, '. \n') )

rdr  =  -11:-5
for( nos in  1:7 ){
    m0 = 10^rdr[nos]
    
    subDir  =  paste0('Tests/Model/', nos )
    if (! file.exists(subDir))  dir.create( subDir ) 
    
    # To fix a random choose:
    set.seed(123456)
    cloneoutfile <- paste0( subDir, '/cloneout.txt' )  # output information of simulation
    
    ### Simulation of the cancer cell/clone evolution:
    smlt = model(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m0, uo, us, s0, k0, censor_cells_number, censor_time_step, d0)
    
    clones      =  smlt[[1]]
    onco_clones =  smlt[[2]] 
    
    write_pnt_clones( pnt_clones, file_out  =  paste0( subDir, '/point_mutations.txt' )  )
    write_pnt_clones( cna_clones, file_out  =  paste0( subDir, '/CNA_mutations.txt' )  )
    
    source("Code/Functions_clones.R")
    get_flow_data(cloneoutfile, genefile, sbdr_Output = subDir )
    vf = get_VAF( file_name = paste0( subDir, '/VAF_data.txt' ) )
    VAF  =  get_rho_VAF( vf = vf, rho = c( 0.0, 0.1, 0.2, 0.5, 0.7, 0.9 ) , file_name = paste0( subDir, '/VAF.txt' ) )
    
    rdr_dysf  =  get_order_of_genes_dysfunction( pnt_mut = pnt_mut_B, file_name = paste0( subDir, '/order_genes_dysfunction.txt' ) )
    
}


test_that( paste0( msg, ': \n'), {
    
    for( nos in  1:7 ){
        subDir        =  paste0('./Tests/Model/', nos )
        subDir_check  =  paste0('./Tests/Model/Check/', nos )
        for( nm in c( '/cloneout.txt', '/point_mutations.txt', '/CNA_mutations.txt',  
                      '/VAF_data.txt', '/VAF.txt', '/order_genes_dysfunction.txt' ) ){ 
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

