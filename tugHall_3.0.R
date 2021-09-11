
### LIBRARIES and SOURCES ------------------------------------------------


    ###  The simulation uses the functions and classes in the "Code/tugHall_2.1_functions.R" 
    
    library(stringr)   # to use string data in input files
    library(actuar)    # to use BIG NUMBERS in N_cell variable
    
    source(file = "Code/tugHall_3.0_functions.R")
    source(file = "Code/read_maps.R")
    

### Define the FOLDERS and files' names ---------------------------------------------------
    ## Create folders:  /Input, /Output and /Figures 
    
    mainDir <- getwd()
    subDir <- "Output"
    if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }
    
    subDir <- "Input"
    if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }
    
    subDir <- "Figures"
    if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }
    
    
    ### Files to output and input data
    genefile <- 'Input/gene_hallmarks.txt'    # gene file 
    clonefile <- 'Input/cloneinit.txt'     # initial Cells 
    
    ### Output files
    geneoutfile <- 'Output/geneout.txt'  # Gene Out file with Hallmarks 
    cloneoutfile <- 'Output/cloneout.txt'  # output information of simulation
    logoutfile <-  'Output/log.txt'      # log file to save the input information of simulation - "log.txt"
    ### Output/Weights.txt               # file with gene weights for hallmarks
    
### Define the gene map - chromosomal locations --------------------------

    
    ### Make a map of genes with sorting of start position for each chromosome:
    gene_map  <-   make_map(f_out    =  'Input/gene_map.txt', 
                             ls   =  c( 'CCDS4107.1', 'CCDS8702.1', 
                                        'CCDS43171.1', 'CCDS11118.1' ), 
                             f_in =  'Input/CCDS.current.txt' )
    gene_map  <-  order_gene_map( gene_map )  ### We have to be sure in the sorting of start position for each chromosome
    write.table(gene_map, file = 'Output/gene_MAP.txt', col.names = TRUE, 
                            sep = "\t", row.names = FALSE)                
    

### Define the PARAMETERS ------------------------------------------------

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
    censore_t <<- 30         # Max time where the program forcibly stops
    ### New parameters for CNA:
    m_dup  <<- 1E-8 # mutation probability for duplication
    m_del  <<- 1E-9 # mutation probability for deletion 
    lambda_dup  <<- 5000  # CNA duplication average length (of the geometrical distribution for the length)
    lambda_del  <<- 7000  # CNA deletion average length
    uo_dup  <<- 0.8 # Gene malfunction probability by CNA duplication for oncogene
    us_dup  <<- 0   # Gene malfunction probability by CNA duplication for suppressor
    uo_del  <<- 0   # Gene malfunction probability by CNA deletion    for oncogene
    us_del  <<- 0.8 # Gene malfunction probability by CNA deletion    for suppressor
    
    

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
smlt = model(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0)

### GET RESULTS ----------------------------------------------------------

if ( FALSE ){
    clones      <- smlt[[1]]
    onco_clones <- smlt[[2]] 
    
    write_pnt_clones( pnt_clones, file = 'Output/point_mutations.txt' )
    write_pnt_clones( cna_clones, file = 'Output/CNA_mutations.txt' )
    
    cn <- read.csv(file = 'Output/CNA_mutations.txt', sep = '\t')
    View(cn)
    pn <- read.csv(file = 'Output/point_mutations.txt', sep = '\t')
    View(pn)
    
    sapply(clones, FUN = function(x) x$field(name = 'id') ) 
    sapply(onco_clones, FUN = function(x) x$field(name = 'id') )
    
    unlist( sapply(clones, FUN = function(x) x$field(name = 'PointMut_ID') )  )
    unlist( sapply(clones, FUN = function(x) x$field(name = 'CNA_ID') )  )
    
    # check the CNA and point mutations
    unlist( sapply( clones, FUN = function(x) print( c(x$id, x$CNA_ID, x$PointMut_ID) ) ) )
    
}


### ANALYZE the RESULTS --------------------------------------------------
#### Analysis of the output data:

# Note: if output files have no data to plot Code/Analysis.R produces errors during plotting

# source("Code/Analysis_clones.R")

# In order to make report, please, use USER-GUIDE.Rmd to show results of simulation

# In order to improve the output plot, please, use Code/Functions.R and Code/Analysis.R scripts

