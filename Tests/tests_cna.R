
### TESTs for all functions including main function - 'model'
### Part I. Tests for CNA and point mutation functions

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

if ( !exists('Print_output_data') )  {
  pk = readline(prompt="Print output data of the tests? Press key: 1 - Yes, any others - No. ")
  if ( pk == '1' ) Print_output_data  =  TRUE else Print_output_data  =  FALSE
}

# Generate objects of simulation ------------------------------------------

### Parameters
define_files_names( sbdr_Input = '/Tests/Input', sbdr_Output = '/Tests/Output' )   
define_gene_location()
define_paramaters( read_fl = TRUE , file_name = './Tests/Input/parameters.txt' )
define_compaction_factor( read_fl = TRUE , file_name = './Tests/Input/CF.txt' )
print( 'Tested parameters are obtained from /Tests/Input/ folder.' )

# Define trial() function: trial_complex or trial_simple
if ( model_name != 'simplified' ){
  trial  =  trial_complex
} else {
  trial  =  trial_simple
}



### Objects
if ( TRUE ){
    onco = oncogene$new()        # make the vector onco about the hallmarks
    onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
    hall = hallmark$new()        # make a vector hall with hallmarks parameters
    hall$read(genefile, onco$name)     # read from the genefile - 'gene_hallmarks.txt'
    env = environ$new(F0)               # new vector for average values of cells
    pnt = Point_Mutations$new()
    pnt_clones = NULL
    cna = CNA_Mutations$new()
    cna_clones = NULL
    mut_order  <<-  0            #  mutation order to reproduce gene map
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

# Test for gene maps reading ----------------------------------------------



cat('Check the reading of chromosomal location from CCDS.current.txt 
            which was getting from CCDS database. ')
test_that("Check the reading of chromosomal 
          location from CCDS.current.txt which was getting from CCDS database: ", {
    
    chk_data  =  readRDS( file = './Tests/GENE_MAP/test_map.txt' ) 

    # saveRDS( gene_map, file = './Tests/GENE_MAP/test_map.txt' ) 
    # print( 'The chromosomal locations for genes: ')
    if ( Print_output_data ){
        print('Output: ')
        print( gene_map )
    }
    expect_identical(chk_data, gene_map )
    # expect_true(identical(chk_data, gene_map ))
    
})


# saveRDS( list_len_cds, file = './Tests/GENE_MAP/length_CDS.txt')

cat('Check the calculation of lengths of genes and their CDS lengths. ')
test_that("Check the calculation of lengths of genes and their CDS length: ", {
              list_len_cds  =  get_len_cds_rna( gene_map )
              chk_data = readRDS( file = './Tests/GENE_MAP/length_CDS.txt')
              # print( 'The chromosomal locations for genes: ')
              if ( Print_output_data ){
                  print('Output: ')
                  print( list_len_cds )
              }
              expect_identical(chk_data, list_len_cds )
              expect_true(identical(chk_data, list_len_cds ))
              
          })



# Tests for CNA  ----------------------------------------------------------

# list_prob_len = get_cds_rna( gm = gene_map )
# saveRDS( list_prob_len, file = './Tests/CNA/prob_length.txt')
gene_map$pnts  = ''

cat('Check the calculation of probabilities based on lengths of genes and 
    their CDS lengths. The function get_cds_rna( gm ). ')
test_that("Check the calculation of probabilities based on lengths of genes and 
          their CDS lengths. The function get_cds_rna( gm ): ", {
    list_prob_len  =  get_cds_rna( gm = gene_map ) 
    chk_data = readRDS( file = './Tests/CNA/prob_length.txt')
    
    # saveRDS( list_prob_len, file = './Tests/CNA/prob_length.txt')
    
    if ( Print_output_data ){
        print( 'Given probabilities for each site: ')
        print( paste0('m0 = ', as.character( m0 ) ) )
        print( paste0('m_dup = ', as.character( m_dup ) ) )
        print( paste0('m_del = ', as.character( m_del ) ) )
        
        print('Output: ')
        print( list_prob_len )
    }
    expect_identical(chk_data, list_prob_len )
    # expect_true(identical(chk_data, list_prob_len ))
    
})


### Tests for add_duplication() function:

### Make data.frame for testing of deletion and duplication
if ( TRUE){
gm1  =  gene_map$Start[2]
gm2  =  gene_map$End[2]
gm12 =  gene_map$End[12]

St_En_Chr   =   data.frame( Start = c( gm1-50, gm1,    gm1+50, gm1+20, gm1, gm1+50,  gm1-50,  gm1,  gm1+50,  gm1+50 ), 
                            End   = c( gm2-50, gm2-50, gm2,    gm2-20, gm2, gm12-50, gm12-50, gm12, gm12-50, gm12  ), 
                            Chr = as.character( rep('5', 10 ) ) )
}

if ( FALSE ){
    lst_dupl  =  lapply( X = 1:10, FUN = function( x ) {
        add_duplication( gm = gene_map, Ref_start = St_En_Chr$Start[ x ], 
                         Ref_end = St_En_Chr$End[ x ], Chr = St_En_Chr$Chr[ x ] ) })
    
    saveRDS( lst_dupl, file = './Tests/CNA/lst_dupl.txt' ) 
}

cat('Check the modification of chromosomal locations due to duplication. 
    The function add_duplication( ). ')
test_that("Check the modification of chromosomal locations due to duplication. 
            The function add_duplication( ). : ", {
      lst_dupl  =  lapply( X = 1:10, FUN = function( x ) {
            add_duplication( gm = gene_map, Ref_start = St_En_Chr$Start[ x ], 
                             Ref_end = St_En_Chr$End[ x ], Chr = St_En_Chr$Chr[ x ] ) })
        
      chk_data = readRDS( file = './Tests/CNA/lst_dupl.txt')
      if ( Print_output_data ){
          print( 'Given the duplication for 10 sites: ')
          print( St_En_Chr )
      }
      expect_identical(chk_data, lst_dupl )
      # expect_true(identical(chk_data, list_prob_len ))
      
})


if ( FALSE ){
    lst_del  =  lapply( X = 1:10, FUN = function( x ) {
        add_deletion( gm = gene_map, Ref_start = St_En_Chr$Start[ x ], 
                         Ref_end = St_En_Chr$End[ x ], Chr = St_En_Chr$Chr[ x ] ) })
    
    saveRDS( lst_del, file = './Tests/CNA/lst_deletion.txt' ) 
}


cat('Check the modification of chromosomal locations due to deletion 
    The function add_deletion( ). ')
test_that("Check the modification of chromosomal locations due to deletion. 
            The function add_deletion( ). : ", {
        lst_del  =  lapply( X = 1:10, FUN = function( x ) {
            add_deletion( gm = gene_map, Ref_start = St_En_Chr$Start[ x ], 
                             Ref_end = St_En_Chr$End[ x ], Chr = St_En_Chr$Chr[ x ] ) })
        
        chk_data = readRDS( file = './Tests/CNA/lst_deletion.txt')
        if ( Print_output_data ){
            print( 'Given the deletions for 10 sites: ')
            print( St_En_Chr )
        }
        expect_identical(chk_data, lst_del )
})




if ( FALSE ){
    pos_pnts  =  1:10 * 5  +  gene_map$Start[2:11]
    lst_pnts  =  lapply( X = 1:10, FUN = function( x ) {
        add_pnt_mutation( gm = gene_map, pos_pnt = pos_pnts[x], Chr = '5' ) })
    saveRDS( lst_pnts, file = './Tests/CNA/lst_pnts.txt' ) 
}


cat( 'Check the modification of chromosomal locations due to point mutation 
    The function add_pnt_mutation( ). ' )
test_that( "Check the modification of chromosomal locations due to point mutation 
    The function add_pnt_mutation( ) : ", {
        pos_pnts  =  1:10 * 5  +  gene_map$Start[2:11]
        lst_pnts  =  lapply( X = 1:10, FUN = function( x ) {
            add_pnt_mutation( gm = gene_map, pos_pnt = pos_pnts[ x ], Chr = '5' ) })
        chk_data = readRDS( file = './Tests/CNA/lst_pnts.txt' )
        if ( Print_output_data ){
            print( 'Given the point mutations at 10 sites: ' )
            print( pos_pnts )
        }
        expect_identical( chk_data, lst_pnts )
} )



cat( 'Check the modification of chromosomal locations due to duplication after 
        the point mutations. The function add_duplication( ). ' )
test_that( "Check the modification of chromosomal locations due to duplication after 
        the point mutations. The function add_duplication( ): ", {
        pos_pnts  =  gene_map$End[2:11] - 1:10 * 5  
        gm = gene_map
        for( i in 1:10 ){
            gm = add_pnt_mutation( gm = gm, pos_pnt = pos_pnts[ i ], Chr = '5' ) 
        }
        
        lst_dup_and_pnts  =  lapply( X = 1:10, FUN = function( x ) {
            add_duplication( gm = gm, Ref_start = St_En_Chr$Start[ x ], 
                             Ref_end = St_En_Chr$End[ x ], Chr = St_En_Chr$Chr[ x ] ) })
        
        # saveRDS( lst_dup_and_pnts, file = './Tests/CNA/lst_dup_and_pnts.txt' ) 
        chk_data = readRDS( file = './Tests/CNA/lst_dup_and_pnts.txt')
        if ( Print_output_data ){
            print( 'Given the point mutations at 10 sites: ' )
            print( pos_pnts )
        
            print( 'And after that the different duplications: ')
            print( St_En_Chr )
        }
        expect_identical( chk_data, lst_dup_and_pnts )
    } )



cat( 'Check the modification of chromosomal locations due to deletion after 
        the point mutations. The function add_deletion( ). ' )
test_that( "Check the modification of chromosomal locations due to deletion after 
        the point mutations. The function add_deletion( ): ", {
        pos_pnts  =  gene_map$End[2:11] - 1:10 * 5  
        gm = gene_map
        for( i in 1:10 ){
            gm = add_pnt_mutation( gm = gm, pos_pnt = pos_pnts[ i ], Chr = '5' ) 
        }
        
        lst_del_and_pnts  =  lapply( X = 1:10, FUN = function( x ) {
            add_deletion( gm = gm, Ref_start = St_En_Chr$Start[ x ], 
                             Ref_end = St_En_Chr$End[ x ], Chr = St_En_Chr$Chr[ x ] ) })
        
        # saveRDS( lst_del_and_pnts, file = './Tests/CNA/lst_del_and_pnts.txt' ) 
        chk_data = readRDS( file = './Tests/CNA/lst_del_and_pnts.txt')
        if ( Print_output_data ){
            print( 'Given the point mutations at 10 sites: ' )
            print( pos_pnts )
        
            print( 'And after that the different deletions: ')
            print( St_En_Chr )
        }
        expect_identical( chk_data, lst_del_and_pnts )
} )

