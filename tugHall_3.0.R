
### LIBRARIES and SOURCES ------------------------------------------------


    ###  The simulation uses the functions and classes in the "Code/tugHall_2.1_functions.R" 
    
    library(stringr)   # to use string data in input files
    library(actuar)    # to use BIG NUMBERS in N_cell variable
    
    source(file = "Code/tugHall_3.0_functions.R")
    source(file = "Code/read_maps.R")

### SIMULATION -----------------------------------------------------------
### Simulation of the cancer cell/clone evolution:
  
    # source('./tugHall_3.0.R')
    define_files_names()    
    define_gene_location()
    define_paramaters( read_fl = TRUE , file_name = './Input/parameters.txt' )
    define_compaction_factor( read_fl = TRUE , file_name = './Input/CF.txt' )
    print_parameters()
    
    # Define trial() function: trial_complex or trial_simple
    if ( model_name != 'simplified' ){
        trial  =  trial_complex
    } else {
        trial  =  trial_simple
    }
    
    n_c  =  0 
    repeat{
        
        n_c  =  n_c + 1
        
        smlt = model(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0)
        
        if ( file.exists( cloneoutfile ) ) break
        if ( n_c  >=  n_repeat )           break
    }
    
    
    clones        <-  smlt[[ 1 ]]
    onco_clones   <-  smlt[[ 2 ]] 
    
    write_pnt_clones( pnt_clones, file = 'Output/point_mutations.txt' )
    write_pnt_clones( cna_clones, file = 'Output/CNA_mutations.txt' )
    # rm( list = ls() )


# Get VAF data  -----------------------------------------------------------
    source("Code/Functions_clones.R")
    get_flow_data(cloneoutfile, genefile )
    
    ### Also VAF data in the file 'Output/VAF_data.txt'
    vf = get_VAF()
    
    VAF  =  get_rho_VAF( vf = vf, rho = c( 0.0, 0.1, 0.2, 0.5, 0.7, 0.9 ) , file_name = './Output/VAF.txt' )
    
    

# Plot data ---------------------------------------------------------------

    source( './Code/my_plots.R' )
    
    rdr_dysf  =  get_order_of_genes_dysfunction()
    
    plot_order_dysfunction( rdr_dysf , pos = c(29,450), logscale = 'y', cex = 0.5 )
    
    

    plot_average_simulation_data() 
    
    
    
    