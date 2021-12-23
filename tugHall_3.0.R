
### LIBRARIES and SOURCES ------------------------------------------------


    ###  The simulation uses the functions and classes in the "Code/tugHall_2.1_functions.R" 
    
    library(stringr)   # to use string data in input files
    library(actuar)    # to use BIG NUMBERS in N_cell variable
    
    source(file = "Code/tugHall_3.0_functions.R")
    source(file = "Code/read_maps.R")
    


### MAKE INPUT CLONES ----------------------------------------------------

    # if you have a new format of gene file, please, use change of columns function like: 
    # genefile <- changeCol(genefile)
    
    ### Making of the input file for initial clones
    
    if (FALSE){
        
        
        x <- 1
        #xz <- data.frame(V1=x,V2=as.character("PIK3CA,APC,KRAS,TP53"), V3=rep.int(1000,length(x)))   # the id of clone, the genes, the number of cells in the clone
        xz <- data.frame(V1=x,V2=as.character(""), V3=rep.int(10^3,length(x)))   # the id of clone, the genes, the number of cells in the clone
        #xz[1,2] <- "PIK3CA"
        xz[2,] <- c(2,"APC",10^3)
        #xz[3,2] <- "KRAS"
        #xz[4,2] <- "TP53"
        xz$V2 <- as.character(xz$V2)
        
        write.table(xz,file = clonefile, col.names = FALSE,sep = "\t",row.names = FALSE)
    }



### SIMULATION -----------------------------------------------------------
### Simulation of the cancer cell/clone evolution:
  
    # source('./tugHall_3.0.R')
    define_files_names()    
    define_gene_location()
    define_paramaters()
    
    smlt = model(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0)
    
    clones      <- smlt[[1]]
    onco_clones <- smlt[[2]] 
    
    write_pnt_clones( pnt_clones, file = 'Output/point_mutations.txt' )
    write_pnt_clones( cna_clones, file = 'Output/CNA_mutations.txt' )
    # rm(list = ls())

    
### Get VAF data    
    source("Code/Functions_clones.R")
    get_flow_data(cloneoutfile, genefile )
    
    ### Also VAF data in the file 'Output/VAF_data.txt'
    vf = get_VAF()
    
    VAF  =  get_rho_VAF( vf = vf, rho = c( 0.1, 0.2, 0.5, 0.7, 0.9, 1 ) , file_name = './Output/VAF.txt' )
    
    
    
    
    
    
    
    