
### Define the FOLDERS and files' names ---------------------------------------------------
## Create folders:  /Input, /Output and /Figures 
define_files_names  <-  function(){    
  mainDir <- getwd()
  subDir <- "Output"
  if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }
  
  subDir <- "Input"
  if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }
  
  subDir <- "Figures"
  if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }
  
  
  ### Files to output and input data
  genefile <<- 'Input/gene_hallmarks.txt'    # gene file 
  clonefile <<- 'Input/cloneinit.txt'     # initial Cells 
  
  ### Output files
  geneoutfile <<- 'Output/geneout.txt'  # Gene Out file with Hallmarks 
  cloneoutfile <<- 'Output/cloneout.txt'  # output information of simulation
  logoutfile <<-  'Output/log.txt'      # log file to save the input information of simulation - "log.txt"
  ### Output/Weights.txt               # file with gene weights for hallmarks
}    
### Define the gene map - chromosomal locations --------------------------

define_gene_location  <-  function( file_input  =  'Input/CCDS.current.txt',
                                    genes_list  =  c( 'CCDS4107.1', 'CCDS8702.1', 
                                                      'CCDS43171.1', 'CCDS11118.1' ) ){    
  ### Make a map of genes with sorting of start position for each chromosome:
  gene_map  <<-   make_map(f_out    =  'Input/gene_map.txt', 
                           ls   =  genes_list, 
                           f_in =  file_input )
  gene_map  <<-  order_gene_map( gene_map )  ### We have to be sure in the sorting of start position for each chromosome
  
  write.table(gene_map, file = 'Output/gene_MAP.txt', col.names = TRUE, 
              sep = "\t", row.names = FALSE)                
}

### Define the PARAMETERS ------------------------------------------------

# Probabilities of processes
define_paramaters  <-  function( E0 =  1E-4, F0 =  10, m0 =  1E-7, uo =  0.9, us =  0.9, 
                                 s0 =  10, k0 =  0.12, d0 =  0.4, censore_n = 10^5,
                                 censore_t = 50, m_dup  = 1E-8, m_del  = 1E-8,
                                 lambda_dup  = 5000, lambda_del  = 7000, 
                                 uo_dup  = 0.8, us_dup  = 0.5, uo_del  = 0, us_del  = 0.8,
                                 CF  =  TRUE,  model  =  'proportional_metastatic', time_stop = 120, 
                                 read_fl = FALSE, file_name ='./Input/parameters.txt', 
                                 n_repeat = 1000 ){  
    if ( read_fl ){
        data_log  =  read.table( file = file_name, sep = '\t', stringsAsFactors = FALSE )
        names( data_log )  =  c( 'var', 'value' )
        # Model definition
        Compaction_factor  <<-  as.logical( data_log[ which( data_log$var == 'Compaction_factor' ), 2 ] )
        model_name         <<-  data_log[ which( data_log$var == 'model_name' ), 2 ]  
        time_stop          <<-  as.numeric( data_log[ which( data_log$var == 'time_stop' ), 2 ] )  # max time in seconds
        n_repeat           <<-  as.numeric( data_log[ which( data_log$var == 'n_repeat' ), 2 ] )  # max number of repetitions
        # Parameters:
        E0 <<-  as.numeric( data_log[ which( data_log$var == 'E' ), 2 ] )       # parameter in the division probability  
        F0 <<-  as.numeric( data_log[ which( data_log$var == 'F' ), 2 ] )       # parameter in the division probability  
        m0 <<-  as.numeric( data_log[ which( data_log$var == 'm0' ), 2 ] )     # mutation probability  
        uo <<-  as.numeric( data_log[ which( data_log$var == 'uo' ), 2 ] )        # oncogene mutation probability  
        us <<-  as.numeric( data_log[ which( data_log$var == 'us' ), 2 ] )        # suppressor mutation probability  
        s0 <<-  as.numeric( data_log[ which( data_log$var == 's' ), 2 ] )         # parameter in the sigmoid function  
        k0 <<-  as.numeric( data_log[ which( data_log$var == 'k' ), 2 ] )        # Environmental death probability  
        d0 <<-  as.numeric( data_log[ which( data_log$var == 'd0' ), 2 ] )      # Initial probability to divide cells
        ### Additional parameters of simulation
        censore_n <<- as.numeric( data_log[ which( data_log$var == 'censore_n' ), 2 ] )       # Max cell number where the program forcibly stops
        censore_t <<- as.numeric( data_log[ which( data_log$var == 'censore_t' ), 2 ] )       # Max time where the program forcibly stops
        ### New parameters for CNA:
        m_dup  <<- as.numeric( data_log[ which( data_log$var == 'm_dup' ), 2 ] ) # mutation probability for duplication
        m_del  <<- as.numeric( data_log[ which( data_log$var == 'm_del' ), 2 ] ) # mutation probability for deletion 
        lambda_dup  <<- as.numeric( data_log[ which( data_log$var == 'lambda_dup' ), 2 ] )  # CNA duplication average length (of the geometrical distribution for the length)
        lambda_del  <<- as.numeric( data_log[ which( data_log$var == 'lambda_del' ), 2 ] )  # CNA deletion average length
        uo_dup  <<- as.numeric( data_log[ which( data_log$var == 'uo_dup' ), 2 ] ) # Gene malfunction probability by CNA duplication for oncogene
        us_dup  <<- as.numeric( data_log[ which( data_log$var == 'us_dup' ), 2 ] )   # Gene malfunction probability by CNA duplication for suppressor
        uo_del  <<- as.numeric( data_log[ which( data_log$var == 'uo_del' ), 2 ] )   # Gene malfunction probability by CNA deletion    for oncogene
        us_del  <<- as.numeric( data_log[ which( data_log$var == 'us_del' ), 2 ] ) # Gene malfunction probability by CNA deletion    for suppressor
        
    } else {
        
        # Model definition:
        Compaction_factor  <<-  CF 
        model_name         <<-  model 
        # Parameters:
        E0 <<-  E0       # parameter in the division probability  
        F0 <<-  F0         # parameter in the division probability  
        m0 <<-  m0      # mutation probability  
        uo <<-  uo        # oncogene mutation probability  
        us <<-  us        # suppressor mutation probability  
        s0 <<-  s0         # parameter in the sigmoid function  
        k0 <<-  k0        # Environmental death probability  
        d0 <<-  d0       # Initial probability to divide cells
        ### Additional parameters of simulation
        censore_n <<- censore_n       # Max cell number where the program forcibly stops
        censore_t <<- censore_t         # Max time where the program forcibly stops
        time_stop  <<-  time_stop     # Max time in seconds of running after that the program forcibly stops
        n_repeat     <<-   n_repeat     # Max number of repetition of the program until the NON-ZERO output will be getting
        ### New parameters for CNA:
        m_dup  <<- m_dup # mutation probability for duplication
        m_del  <<- m_del # mutation probability for deletion 
        lambda_dup  <<- lambda_dup  # CNA duplication average length (of the geometrical distribution for the length)
        lambda_del  <<- lambda_del  # CNA deletion average length
        uo_dup  <<- uo_dup # Gene malfunction probability by CNA duplication for oncogene
        us_dup  <<- us_dup   # Gene malfunction probability by CNA duplication for suppressor
        uo_del  <<- uo_del   # Gene malfunction probability by CNA deletion    for oncogene
        us_del  <<- us_del # Gene malfunction probability by CNA deletion    for suppressor
    }
}    

# function to print GLOBAL parameters:
print_parameters  <-  function(){

    msg  =  c(
        'Model definition:  \n ' ,
        'Compaction_factor = ', Compaction_factor, '\n',
        'model_name  =  ', model_name, '\n', 
        'Parameters:  \n', 
        'parameter of the division probability E0 =  ', E0, '\n', 
        'another parameter of the division probability F0  = ',  F0, '\n',
        'mutation probability m0 =  ', m0, '\n', 
        'oncogene mutation probability uo = ', uo, '\n',   
        'suppressor mutation probability  us  =  ', us, '\n',  
        'parameter in the sigmoid function  s0  =  ', s0, '\n',   
        'Environmental death probability  k0 =  ',  k0, '\n', 
        'Initial probability to divide cells  d0  =  ',  d0, '\n', 
        'Additional parameters of simulation  \n ',
        'Max cell number where the program forcibly stops  censore_n  = ',  censore_n,  '\n',  
        'Max time steps where the program forcibly stops  censore_t  = ',  censore_t,  '\n',
        'Max time (in seconds) where the program forcibly stops time_stop  =  ',  time_stop,  '\n',
        'Max number of repetition of the program until the NON-ZERO output will be getting, n_repeat  =  ', n_repeat ,   '\n', 
        'New parameters for CNA:  \n', 
        'mutation probability for duplication  m_dup  =  ', m_dup ,  '\n',  
        'mutation probability for deletion',  m_del,  '\n',  
        'CNA duplication average length (of the geometrical distribution for the length)  lambda_dup  =  ', lambda_dup ,  '\n',  
        'CNA deletion average length  lambda_del  = ', lambda_del ,  '\n',  
        'Gene malfunction probability by CNA duplication for oncogene  uo_dup  =  ', uo_dup ,  '\n',  
        'Gene malfunction probability by CNA duplication for suppressor  us_dup  =  ', us_dup ,  '\n',  
        'Gene malfunction probability by CNA deletion for oncogene  uo_del  = ', uo_del ,  '\n',  
        'Gene malfunction probability by CNA deletion for suppressor  us_del  = ', us_del,  '\n',
        'Compaction factor is applied if variable Compaction_factor ==  TRUE \n',
        'Compaction factor for apoptosis hallmark CF$Ha = ', CF$Ha, ' \n', 
        'Compaction factor for angiogenesis hallmark CF$Hb = ', CF$Hb, ' \n', 
        'Compaction factor for growth/antigrowth hallmark CF$Hd = ', CF$Hd, ' \n', 
        'Compaction factor for immortalization hallmark CF$Hi = ', CF$Hi, ' \n', 
        'Compaction factor for invasion/metastasis hallmark CF$Him = ', CF$Him 
    )
    
    cat( paste0( msg, collapse = ' ' ) )
}

define_compaction_factor  <-  function( cf = data.frame( Ha = 1, Hb = 1, Hd = 1,
                                                         Hi = 1, Him = 1 ), 
                                        read_fl = TRUE , file_name = './Input/CF.txt' ){
    if ( read_fl ){
        data_log  =  read.table( file = file_name, sep = '\t', stringsAsFactors = FALSE )
        names( data_log )  =  c( 'var', 'value' )
        
        cf$Ha   =  data_log$value[ data_log$var == 'apoptosis' ]
        cf$Hb   =  data_log$value[ data_log$var == 'angiogenesis' ]
        cf$Hd   =  data_log$value[ data_log$var == 'growth' ]
        cf$Hi   =  data_log$value[ data_log$var == 'immortalization' ]
        cf$Him  =  data_log$value[ data_log$var == 'invasion' ]
    }
    
    CF  <<-  cf
}

#### The code:

#### I)  Define CLONE'S CLASSES ----------------------------------------------------------


# Clone class
clone <- setRefClass(
    # name of the class
    Class = "Clone",
    # Field
    fields = list(
        id = "numeric",          # identificator
        parent = "numeric",      # parent ID (for first - 0)
        N_cells = "numeric",     # Number of cells in clone 
        c = "numeric",           # split counter
        d = "numeric",           # probability of division
        i = "numeric",           # probability of hayflick limit
        m = "numeric",           # probability that gene normal function is destroyed due to epigenome abnormality
        a = "numeric",           # probability of apoptosis
        s = "numeric",           # coefficient in the of apoptosis
        k = "numeric",           # probability of cell death by environment
        E = "numeric",           # coefficient of friction term against to the split probability,
                                 # coefficient for determination the max number of cells that can exist in the primary tumor (Nmax = 1/E)
        Nmax = "numeric",        # the max number of cells that can exist in the primary tumor (Nmax = 1/E)
        im = "numeric",          # invasion/ metastasis probability
        Ha = "numeric",          # apoptosis probability difference (apoptosis gene x weight) 
        Him = "numeric",         # invasion/ metastasis probability difference (invasion/ metastasis gene x weight) 
        Hi = "numeric",          # mitotic restriction probability difference (immortalized gene x weight)
        Hd = "numeric",          # divide probability difference (cancer gene x weight)
        Hb = "numeric",          # friction coefficient (angiogenic gene x weight)
        gene = "numeric",        # flag for cancer gene function deficit (for drivers)
        pasgene = "numeric",     # flag for cancer gene as passenger dysfunction 
        # deleted: posdriver = "character", # position of cancer gene damage (function deficit)
        # deleted: pospasngr = "character", # position of cancer gene damage (maintenance of function)
        PointMut_ID  =  "numeric", # ID of point mutation in object PointMut / pnt 
        CNA_ID       =  "numeric", # ID of CNA mutation in object CNA
        mutden = "numeric",      # gene mutation density
        invasion = "logical",    # Wetting/Displacement flag:    TRUE: Wetting/Displacement      FALSE: Limited pattern
        birthday = "numeric"    # time step of birth of cell
        # lenCDS      = "numeric"     # length of all oncogenes of interest
        # lenRNA      = "numeric"      # length of all oncogenes of interest including introns
    ),

    # Method
    methods = list(
        # Initialize
        initialize = function(gene_size, id=1, parent=0, c=0, d=d0, i=1, m=m0,  N_cells = 1, 
                              mutden=0, a=0, k=k0, E=E0, Nmax=0, gene=NULL, pasgene=NULL,
                              PointMut_ID = 0, CNA_ID = 0,
                              # deleted: posdriver=NULL, pospasngr=NULL, 
                              invasion=FALSE, s=s0, birthday=0) {
            id <<- id
            parent <<- parent
            N_cells <<- N_cells
            c <<- c
            d <<- d
            i <<- i
            m <<- m
            s <<- s
            birthday <<- birthday
            mutden <<- mutden
            if (is.null(a)) {
                a <<- 1/(1+exp(-s*(mutden - 0.5)))
            } else {
                a <<- a
            }
            k <<- k
            E <<- E
            Nmax <<- 1.0 / E
            im <<- 0
            Ha <<- 0
            Him <<- 0
            Hi <<- 0
            Hd <<- 0
            Hb <<- 0

            if (is.null(gene)) {
                gene <<- rep(0, gene_size)
            } else {
                gene <<- gene
            }
            if (is.null(pasgene)) {
                pasgene <<- rep(0, gene_size)
            } else {
                pasgene <<- pasgene
            }
            # if (is.null(posdriver)) {
            #    posdriver <<- rep("", gene_size)
            # } else {
            #    posdriver <<- posdriver
            # }
            # if (is.null(pospasngr)) {
            #     pospasngr <<- rep("", gene_size)
            # } else {
            #     pospasngr <<- pospasngr
            # }
            PointMut_ID <<- PointMut_ID
            CNA_ID     <<- CNA_ID
            invasion   <<- invasion
            # lenCDS     <<-  sum( onco$cds )
            # lenRNA     <<-  integer(0)
        },
        # Apoptosis
        calcApoptosis = function() {
#            if (mutden <= sum(gene)/length(gene)) {
#                a1 = 1/(1+exp(s*(mutden - 0.5)))
#                mutden <<- sum(gene)/length(gene)
#                a2 = 1/(1+exp(s*(mutden - 0.5)))
#                a <<- a - (a1 - a2)

          mutden <<- sum(gene) / length(gene)
          a <<- 1 / ( 1 + exp( -1 * s * ( mutden - 0.5 ) ) )
          if ( a < 0 ) {
              a <<- 0
          }

        },
        # Aggregate
        calcMutden = function() {
            mutden <<- sum( gene ) / length( gene )
        }
    )
)

# Environ class
environ <- setRefClass(
    # the class name
    Class = "Environ",

    # Fields
    fields = list(
        T = "numeric",           # time counter
        N = "numeric",           # localized clones number
        M = "numeric",           # number of infiltrting / metastatic clones
        F = "numeric",           # a coeffitient (Nmax = F/E) that determines the maximal number of cells 
                                 # that can exist in the primary tumor when the hallmark is engraved  
        c = "numeric",           # average number of divisions 
        d = "numeric",           # mean value of spliting probability
        i = "numeric",           # average value of spliting probability
        a = "numeric",           # average value of apoptosis probability
        k = "numeric",           # average probability of cell death
        E = "numeric",           # average value of coefficients of friction term proportional to N, for splitting probability
        Nmax = "numeric",        # Maximal number of cells that can exist 
        im = "numeric",          # average value of invasion / metastasis probability
        Ha = "numeric",
        Him = "numeric",
        Hi = "numeric",
        Hd = "numeric",
        Hb = "numeric",
        type = "numeric",        # invasion / metastatic ratio
        gene = "numeric",        # cancer gene damage rate
        # posdriver = "character", # cancer gene damage position (function deficit)
        # pospasngr = "character", # cancer gene damage position (maintaince of function)
        mutden = "numeric",      # average mutation rate
        last_id = "numeric"
    ),

    # Methods
    methods = list(
        # Initialize
        initialize = function(F0) {
            T <<- 0
            N <<- 0
            M <<- 0
            F <<- F0
        }
    )
)


# Cancer gene
oncogene <- setRefClass(
    # name of the class
    Class = "OncoGene",

    # Fields
    fields = list(
        id = "numeric",       # identificator is same as in clone (key for clones)
        name  = "character",   # Cancer gene name list
        onsp  = "character",   # oncogene/suppressor indicator
        len   = "numeric",      # lengths of cancer genes
        cds_1 = "numeric",      # cancer gene CDS base length for parental chr 1
        cds_2 = "numeric",      # cancer gene CDS base length for parental chr 2
        rna_1 = "numeric",      # cancer gene RNA base number length for parental chr 1
        rna_2 = "numeric",      # cancer gene RNA base number length for parental chr 2
        p0_1  = "numeric",      # the probability of absent of mutations for parental chr 1 
        p0_2  = "numeric",      # the probability of absent of mutations for parental chr 2
        prob_1 = "numeric",     # prob = c(m0*sumCDS,  m_del*sumRNA,  m_dup*sumRNA) / sum(m0*sumCDS,  m_del*sumRNA,  m_dup*sumRNA) 
        prob_2 = "numeric",     # same for parental chr 2 
        sum_prob_1 = "numeric",  # sum_prob = sum(m0*sumCDS,  m_del*sumRNA,  m_dup*sumRNA) 
        sum_prob_2 = "numeric"  # same for parental chr 2
    ),

    # Methods
    methods = list(
        # read the configuration file
        read = function(file) {
            data = read.table(file, sep="\t")
            name0 = NULL
            onsp0 = NULL
            cds0 = NULL
            rna0 = NULL 
            for (i in 1:nrow(data)) {
                name <<- as.character(data[i, 1])
                if (!is.element(name, name0)) {
                    name0 = c(name0, name)
                    type = as.character(data[i, 3])
                    if (type == "?") {
                        if (runif(1) > 0.5) {
                            type = "o"
                        } else {
                            type = "s"
                        }
                    }
                    onsp0 = c(onsp0, type)
                    # changed the length CDS, now get from gene_map file
                    # cds0 = c(cds0, as.numeric(as.character(data[i, 2])))
                }
            }
            
            # Get length of CDS and gene's length from gene_map:
            for ( i in 1:length(name0) ) {
                w     =  which( gene_map$Gene == name0[ i ] )
                cds0  =  c( cds0, sum( gene_map$End[w]  -  gene_map$Start[w] + 1 ) )
                rna0  =  c( rna0, max( gene_map$End[w] ) - min( gene_map$Start[w]) + 1  )
            }
            
            id   <<- 1
            name <<- name0
            onsp <<- onsp0
            cds_1  <<- cds0
            cds_2  <<- cds0
            len  <<- length(name0)
            rna_1  <<- rna0 
            rna_2  <<- rna0 
            sum_prob_1 <<- sum( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) )
            sum_prob_2 <<- sum( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) )
            prob_1     <<- c( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) ) / sum_prob_1
            prob_2     <<- c( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) ) / sum_prob_2
            p0_1   <<- (1-m0)**sum(cds0) * (1-m_dup)**sum(rna0) * (1-m_del)**sum(rna0)
            p0_2   <<- (1-m0)**sum(cds0) * (1-m_dup)**sum(rna0) * (1-m_del)**sum(rna0)
        }
    )
)


hallmark <- setRefClass(
    # 
    Class = "HallMark",

    # 
    fields = list(
        Ha = "numeric",       # (Evading apoptosis)
        Hi = "numeric",       # (immortalization limit)
        Hd = "numeric",       # (Insensitivity to anti-growth signals || Self-sufficiency in growth signals)
        Hb = "numeric",       # (Sustained angiogenesis)
        Him = "numeric",      # (Tissue invasion & metastasis)
        Ha_w = "numeric",     # 
        Hi_w = "numeric",     # 
        Hd_w = "numeric",     # 
        Hb_w = "numeric",     # 
        Him_w = "numeric",    # 
        notHa = "numeric"
    ),

    # 
    methods = list(
        # 
        read = function(file, names, normalization  =  TRUE ) {
            # normalization is an indicator to normalize all Hallmarks values
            data <- read.table(file, sep="\t")
            Ha0 = NULL
            Hi0 = NULL
            Hd0 = NULL
            Hb0 = NULL
            Him0 = NULL
            Ha0_w = NULL
            Hi0_w = NULL
            Hd0_w = NULL
            Hb0_w = NULL
            Him0_w = NULL
            w_flg = FALSE
            if (ncol(data) >= 4) {
                w_flg = TRUE
            }
            # Acquire gene name and Hallmark coefficient by function from definition file
            for (i in 1:nrow(data)) {
                if (data[i, 2] == "apoptosis") {
                    Ha0 = c(Ha0, as.character(data[i, 1]))
                    if (w_flg) {
                        Ha0_w = c(Ha0_w, as.numeric(as.character(data[i, 4])))
                    }
                } else if (data[i, 2] == "immortalization") {
                    Hi0 = c(Hi0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hi0_w = c(Hi0_w, as.numeric(as.character(data[i, 4])))
                    }
                } else if (data[i, 2] == "anti-growth" | data[i, 2] == "growth") {
                    Hd0 = c(Hd0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hd0_w = c(Hd0_w, as.numeric(as.character(data[i, 4])))
                    }
                } else if (data[i, 2] == "angiogenesis") {
                    Hb0 = c(Hb0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hb0_w = c(Hb0_w, as.numeric(as.character(data[i, 4])))
                    }
                } else if (data[i, 2] == "invasion") {
                    Him0 = c(Him0, as.character(data[i, 1]))
                    if (w_flg) {
                        Him0_w = c(Him0_w, as.numeric(as.character(data[i, 4])))
                    }
                }
            }
            Ha <<- match(Ha0, names)
            notHa <<- setdiff(seq(1,length(names)),Ha)
            Hi <<- match(Hi0, names)
            Hd <<- match(Hd0, names)
            Hb <<- match(Hb0, names)
            Him <<- match(Him0, names)

            # if there is no Hallmark coefficient then generate a Hallmark coefficient as a random number - beta distribution
            if (!w_flg) {
                if (length(Ha) > 0) {
                    Ha_rnd = 1:length(Ha)
                } else {
                    Ha_rnd = c()
                }
                total0 = length(Ha)
                if (length(Hi) > 0) {
                    Hi_rnd = (total0 + 1):(total0 + length(Hi))
                } else {
                    Hi_rnd = c()
                }
                total0 = total0 + length(Hi)
                if (length(Hd) > 0) {
                    Hd_rnd = (total0 + 1):(total0 + length(Hd))
                } else {
                    Hd_rnd = c()
                }
                total0 = total0 + length(Hd)
                if (length(Hb) > 0) {
                    Hb_rnd = (total0 + 1):(total0 + length(Hb))
                } else {
                    Hb_rnd = c()
                }
                total0 = total0 + length(Hb)
                if (length(Him) > 0) {
                    Him_rnd = (total0 + 1):(total0 + length(Him))
                } else {
                    Him_rnd = c()
                }
                total = total0 + length(Him)
                # random = runif(total)
                random = rbeta(total, 0.01, 1)
                Ha0_w = random[Ha_rnd]
                Hi0_w = random[Hi_rnd]
                Hd0_w = random[Hd_rnd]
                Hb0_w = random[Hb_rnd]
                Him0_w = random[Him_rnd]
            }
            # Total by genetic mode 
            if ( normalization ){
                Ha_sum = sum(Ha0_w)
                Hi_sum = sum(Hi0_w)
                Hd_sum = sum(Hd0_w)
                Hb_sum = sum(Hb0_w)
                Him_sum = sum(Him0_w)
            } else {
                Ha_sum = 1
                Hi_sum = 1
                Hd_sum = 1
                Hb_sum = 1
                Him_sum = 1
            }
            
            Ha_w <<- Ha0_w / Ha_sum
            Hi_w <<- Hi0_w / Hi_sum
            Hd_w <<- Hd0_w / Hd_sum
            Hb_w <<- Hb0_w / Hb_sum
            Him_w <<- Him0_w / Him_sum
            
            if ( Compaction_factor ){
                Ha_w  <<-  CF$Ha  * Ha_w
                Hi_w  <<-  CF$Hi  * Hi_w
                Hd_w  <<-  CF$Hd  * Hd_w
                Hb_w  <<-  CF$Hb  * Hb_w
                Him_w <<-  CF$Him * Him_w
                
            }
            
        },
        # Change the cell variables
        # mode = 2 Corresponding (Hallmark) Gene Mode
        updateClone = function(clone1, F) {
          # Apoptosis
          clone1$calcApoptosis()
          clone1$Ha  =  sum( clone1$gene[Ha] * Ha_w )
          clone1$a   =  clone1$a - clone1$Ha
          if ( clone1$a < 0 ) {
              clone1$a = 0
          }
          # Not dead - Immortalized
          clone1$Hi = sum( clone1$gene[Hi] * Hi_w )
          clone1$i = 1 - clone1$Hi
          if ( clone1$i < 0 ) {
              clone1$i = 0
          }
          # Angiogenesis
          clone1$Hb = sum( clone1$gene[Hb] * Hb_w )
          
          clone1$E     =  E0  / (1 + F * clone1$Hb)
          clone1$Nmax  =  1.0 / clone1$E
          
          # Cancer gene, tumor suppressor gene
          clone1$Hd  = sum( clone1$gene[Hd] * Hd_w )
          clone1$Him = sum( clone1$gene[Him] * Him_w )
          
          clone1$d = clone1$Hd + d0     # d0 is defined in parameters 
          if ( clone1$d > 1 ) { clone1$d = 1 }
          if ( !clone1$invasion ) {
              clone1$d = clone1$d - clone1$E * env$N
              if ( clone1$d < 0 ) { clone1$d = 0 }
          }

                # Invasion metastasis
                if ( !clone1$invasion ) {
                  clone1$im = clone1$Him 
            } else {
              clone1$im = 1
            }
          },
        # Change the environment variables
        updateEnviron = function(env, clones) {
            sum_cell(env, clones)
            # env$M = sum()  # ceiling(length(clones) * env$type)
            # env$N = sum()  # length(clones)  - env$M
        }
    )
)


# Function to update Hallmark and variable after division or under initialization
update_Hallmarks <- function(clone1) {
    # Hallmark
    hall$updateClone(clone1, env$F)
}


#### II) CNA and Point Mutations: CLASSES -------------------------------------------------

# Class and functions related to point mutations of genes of interest
Point_Mutations <- setRefClass(
    # 
    Class = "Point_Mutations",
    # 
    fields = list(
        PointMut_ID	  = "numeric",         # ID of point mutation
        Allele        = "character",        # A or B allele
        Parental_1or2	= "numeric",         # 1 or 2
        Chr	          = "character",       # Chromosome name
        Ref_pos	      = "numeric",         # Reference position
        Phys_pos	    = "vector",          # Physical positions 
        Delta	        = "vector",          # Delta of positions 
        Copy_number	  = "numeric",         # Copy number of allele 
        Gene_name	    = "character",       # Gene's name
        MalfunctionedByPointMut  = "logical",       # True or False
        mut_order     = "numeric"          # order of mutation to reproduce the map 
    ),
    
    # 
    methods = list(
        # Initialization with default values
        initialize = function( ID = 1, ... ){
            PointMut_ID	  <<- ID       # ID of point mutation 
            Allele        <<- ""       # A or B allele
            Parental_1or2	<<-  0       # 1 or 2
            Chr	          <<-  ""       # Chromosome name
            Ref_pos	      <<-  0       # Reference position
            Phys_pos	    <<-  0       # Physical positions [from:to]
            Delta	        <<-  0       # Delta of positions [from:to]
            Copy_number	  <<-  0       # Copy number of allele 
            Gene_name	    <<-  ""       # Gene's name
            MalfunctionedByPointMut  <<- NA       # True of False 
            mut_order     <<-  0       # order of mutation
        },
  
        # Function to safe data to data frame df1
        safe = function() {
          df1 = data.frame( PointMut_ID = PointMut_ID,
                            Parental_1or2 = Parental_1or2,
                            Chr = Chr,
                            Ref_pos = Ref_pos,
                            Phys_pos = paste0( '[', paste(Phys_pos, collapse = ', ' ),  ']'),
                            Delta = paste0( '[', paste(Delta, collapse = ', ' ),  ']'),
                            Copy_number = Copy_number,
                            Gene_name = Gene_name,
                            MalfunctionedByPointMut = MalfunctionedByPointMut,
                            mut_order  =  mut_order,
                            stringsAsFactors = FALSE)
          # df = rbind( df, df1 )
          return( as.data.frame( df1 ) )
      }
        
    )
)


# Class and functions related to CNA mutations of genes of interest
CNA_Mutations <- setRefClass(
  # 
  Class = "CNA_Mutations",
  # 
  fields = list(
    CNA_ID	  = "numeric",         # ID of point mutation
    Parental_1or2	= "numeric",         # 1 or 2
    dupOrdel      = "character",         # duplication OR deletion indicator
    Chr	          = "character",       # Chromosome name
    Ref_start	    = "numeric",         # Reference start position
    Ref_end	      = "numeric",         # Reference end position
    # Phys_pos	    = "vector",          # Physical positions [from:to]
    # Delta	        = "vector",          # Delta of positions [from:to]
    # Copy_number	  = "numeric",         # Copy number of allele 
    Gene_names	    = "character",       # Genes' name
    MalfunctionedByCNA  = "logical",     # True of False 
    mut_order       = "numeric"          # Order of mutation
  ),
  
  # 
  methods = list(
    # Initialization with default values
    initialize = function( ID = 1, ... ){
      CNA_ID	      <<-  ID       # ID of point mutation 
      Parental_1or2	<<-  0       # 1 or 2
      dupOrdel      <<-  ""
      Chr	          <<-  ""       # Chromosome name
      Ref_start	    <<-  0       # Reference position
      Ref_end 	    <<-  0       # 
      Gene_names	  <<- ""       # Gene's name
      MalfunctionedByCNA  <<- NA       # True of False 
      mut_order     <<-  0
    },
    
    # Function to safe data to data frame df1
    safe = function() {
      df1 = data.frame( CNA_ID = CNA_ID,
                        Parental_1or2 = Parental_1or2,
                        dupOrdel = dupOrdel,
                        Chr = Chr,
                        Ref_start = Ref_start,
                        Ref_end = Ref_end,
                        Gene_names = Gene_names,
                        MalfunctionedByCNA = MalfunctionedByCNA,
                        mut_order  =  mut_order,
                        stringsAsFactors = FALSE)
      # df = rbind( df, df1 )
      return( as.data.frame( df1 ) )
    }
    
  )
)



### III) CNA and Point Mutations: Functions --------------------------------------

# Safe function for 1 point mutation
safe_pnt_mut  <-  function(pnt){
  return( pnt$safe() )
}

copy_pnt_no_mutation  <-  function( pnt1 ){
  ### The function to copy of pnt1 without mutation info for allele A
  pnt2  =  Point_Mutations$new()
  pnt2$PointMut_ID  =  pnt1$PointMut_ID
  pnt2$Allele = 'A'
  pnt2$Parental_1or2 = ifelse( pnt1$Parental_1or2 == 2, as.integer(1), as.integer(2) )
  pnt2$Chr  =  pnt1$Chr
  pnt2$Ref_pos  =  pnt1$Ref_pos
  pnt2$Phys_pos  =  NA
  pnt2$Delta     =  NA
  pnt2$Copy_number  =  1
  pnt2$Gene_name    =  pnt1$Gene_name
  pnt2$MalfunctionedByPointMut  =  NA
  pnt2$mut_order    =  pnt1$mut_order
  
  return( pnt2 )
}

copy_pnt  <-  function( pnt1 ){
  ### The function to copy of pnt1
  pnt2  =  Point_Mutations$new()
  pnt2$PointMut_ID  =  pnt1$PointMut_ID
  pnt2$Allele       = pnt1$Allele
  pnt2$Parental_1or2  =  pnt1$Parental_1or2
  pnt2$Chr  =  pnt1$Chr
  pnt2$Ref_pos  =  pnt1$Ref_pos
  pnt2$Phys_pos  =  pnt1$Phys_pos
  pnt2$Delta     =  pnt1$Delta
  pnt2$Copy_number  =  pnt1$Copy_number
  pnt2$Gene_name    =  pnt1$Gene_name
  pnt2$MalfunctionedByPointMut  =  pnt1$MalfunctionedByPointMut
  pnt2$mut_order    =  pnt1$mut_order
  
  return( pnt2 )
}

# Copy function for CNA
copy_CNA  <-  function( CNA1 ){
    ### The function to copy of CNA1
    CNA2  =  CNA_Mutations$new()
    CNA2$CNA_ID  =  CNA1$CNA_ID
    CNA2$Parental_1or2  =  CNA1$Parental_1or2
    CNA2$dupOrdel = CNA1$dupOrdel
    CNA2$Chr  =  CNA1$Chr
    CNA2$Ref_start  =  CNA1$Ref_start
    CNA2$Ref_end    =  CNA1$Ref_end
    CNA2$Gene_names =  CNA1$Gene_names
    CNA2$MalfunctionedByCNA  =  CNA1$MalfunctionedByCNA
    
    CNA2$mut_order  =  CNA1$mut_order
    
    return( CNA2 )
}

### Generation point mutation info
get_point_mutation <- function( onco1, gm_1_2 ){
    ### The function to get position of point mutation related to gene_map info for each chr
    prntl  =   sample( c(1,2), size = 1, replace = TRUE, prob = c( sum(onco1$cds_1), sum(onco1$cds_2) ))
    gene   =   ifelse( prntl == 1, 
      sample( onco1$name, size = 1, replace = TRUE, prob = onco1$cds_1/sum(onco1$cds_1) ),
      sample( onco1$name, size = 1, replace = TRUE, prob = onco1$cds_2/sum(onco1$cds_2) ))
    
    gm     =   gm_1_2[[ prntl ]]
    
    ### get position for  gene !!! + CNA later!!!
    sp  =  which( gm$Gene == gene)   ### get the name of related gene 
    cds =  ifelse( prntl == 1, 
                   onco1$cds_1[ which(onco1$name == gene ) ],
                   onco1$cds_2[ which(onco1$name == gene ) ] ) ### get CDS for the related gene
    p   =  sample( sp, size = 1, replace = TRUE, 
                  prob = gm[sp,'Len'] / sum( gm[sp,'Len'] ) )   ### get the block in gene_map
    pos =  sample( gm$Start[p]:gm$End[p], 1, replace=FALSE) 
    Chr =  gm$Chr[p]
    
    return( list( prntl, gene, pos, Chr ) )
}

### Generation point mutation info for the particular gene
get_point_mutation_for_gene <- function( onco1, gm_1_2, gene ){
  ### The function to get position of point mutation related to gene_map info for each chr
  prntl  =   sample( c(1,2), size = 1, replace = TRUE, prob = c( sum(onco1$cds_1), sum(onco1$cds_2) ))
  
  gm     =   gm_1_2[[ prntl ]]
  
  ### get position for  gene !!! + CNA later!!!
  sp  =  which( gm$Gene == gene)   ### get the name of related gene 
  cds =  ifelse( prntl == 1, 
                 onco1$cds_1[ which(onco1$name == gene ) ],
                 onco1$cds_2[ which(onco1$name == gene ) ] ) ### get CDS for the related gene
  p   =  sample( sp, size = 1, replace = TRUE, 
                 prob = gm[sp,'Len'] / sum( gm[sp,'Len'] ) )   ### get the block in gene_map
  pos =  sample( gm$Start[p]:gm$End[p], 1, replace=FALSE) 
  Chr =  gm$Chr[p]
  
  return( list( prntl, gene, pos, Chr ) )
}

generate_pnt  <-  function( prntl, gene, pos, onco1, Chr, mutation = NA ) {
    ### The function to generate object of point mutation pnt
    pnt0 = Point_Mutations$new()
    pnt0$Gene_name = gene
    pnt0$PointMut_ID =  ifelse( length(pnt_clones) == 0, 1, 
                                pnt_clones[[ length(pnt_clones) ]]$PointMut_ID + 2 ) 
    pnt0$Allele = 'B'  # Mutation occurs on allele B
    pnt0$Parental_1or2  =  as.integer(prntl)   
    pnt0$Chr = Chr
    pnt0$Ref_pos  = pos
    pnt0$Phys_pos = pos
    pnt0$Delta = 0
    pnt0$Copy_number = 1
    u = ifelse( onco1$onsp[ which(onco1$name == gene) ] == 'o', uo, us)
    if ( is.na( mutation ) ) {
        pnt0$MalfunctionedByPointMut  =  ifelse( (runif(1) < u), TRUE, FALSE )
    } else {
        pnt0$MalfunctionedByPointMut  =  TRUE 
    }
    mut_order  <<-  mut_order  +  1
    pnt0$mut_order  =  mut_order
    
    pnt_clones <<- c( pnt_clones, pnt0 )
    
    pnt1  =  copy_pnt_no_mutation( pnt0 )
    pnt_clones <<- c( pnt_clones, pnt1 )

    return( pnt0 )
}

# the function to generate pnt with coping all information from input pnt
generate_to_copy_pnt  <-  function( pnt ) {
    ### The function to generate object of point mutation pnt with data from input pnt
    pnt0 = Point_Mutations$new()
    pnt0$Gene_name = pnt$Gene_name
    pnt0$PointMut_ID =  ifelse( length(pnt_clones) == 0, 1, 
                                pnt_clones[[ length(pnt_clones) ]]$PointMut_ID + 2 ) 
    pnt0$Allele = pnt$Allele
    pnt0$Parental_1or2  =  pnt$Parental_1or2
    pnt0$Chr = pnt$Chr
    pnt0$Ref_pos  = pnt$Ref_pos
    pnt0$Phys_pos = pnt$Phys_pos
    pnt0$Delta = pnt$Delta
    pnt0$Copy_number = pnt$Copy_number
    pnt0$MalfunctionedByPointMut  =  pnt$MalfunctionedByPointMut
    
    pnt0$mut_order  =  pnt$mut_order
    
    pnt_clones <<- c( pnt_clones, pnt0 )
    
    pnt1  =  copy_pnt_no_mutation( pnt0 )
    pnt_clones <<- c( pnt_clones, pnt1 )
    
    return( pnt0 )
}

# the function to choose probability of CNA mutation for several genes
get_u_cna <- function( genes, dupOrdel ){
    # input: genes, u_del or u_dup for oncogenes and suppressors
    u = NULL
    for (gene in genes ) {
        os = onco$onsp[ which( onco$name == gene) ]
        u1 = ifelse( dupOrdel == 'dup', 
                     ifelse( os == 'o', uo_dup, us_dup ),   # for duplication
                     ifelse( os == 'o', uo_del, us_del ) )  # for deletion
        u = c( u, u1 )
    }
    
    return( max( u ) )
}

### Generation CNA mutation info

### Generation point mutation info
get_cna_mutation <- function( onco1, dupOrdel, gm_1_2 ){
    ### The function to get position of point mutation related to gm (gene_map) info for each chr
    prntl  =  sample( c(1,2), size = 1, replace = TRUE, prob = c( sum(onco1$cds_1), sum(onco1$cds_2) ))
    lambda =  ifelse( dupOrdel == 'dup', lambda_dup, lambda_del )
    l_cna  =  round( rexp(10, 1/lambda) ) +1 # rpois( n = 1, lambda ) + 1   # length of CNA
    
    gm     =  gm_1_2[[ prntl ]]
    
    w1  =  sample( 1:length( gm$Start ), size = 1, prob = gm$Len )  # choose the row in gene_map
    Chr = gm$Chr[ w1 ]
    w   =  which( gm$Chr == Chr )    ### Area of checking for CNA
    pos_st =  sample( gm$Start[w1]:gm$End[w1], 1, replace=FALSE)   ### position of start CNA
    pos_end  = pos_st  +  l_cna
    
    w3  =  which( gm$Start < pos_end & gm$Chr == Chr )
    w_cna  = w3[ which(w3 >= w1)]    ### rows of CNA in gene_map
    start_end  =  c(0,0)
    start_end[1]  =  pos_st
    start_end[2]  =  ifelse( pos_end < gm$End[ max(w_cna) ], pos_end, gm$End[ max(w_cna) ] )
    genes   =   unique( gm$Gene[w_cna] )   ### Genes of CNA at same Chr 

    return( list( prntl, Chr, genes, start_end, w_cna ) )
}

# function to get order of mutation for all types
mixed_mut_order   <-   function( clone1 ) {
    # return the order, type and number of mutation (del, dup or point in order of appearance)
    order_tp_num  <-  data.frame( order = NULL, type = NULL, ID = NULL )
    i  =  0
    if ( length( clone1$PointMut_ID ) > 0  & clone1$PointMut_ID  != 0 ){
        for (i in 1:length( clone1$PointMut_ID )) {
            order_tp_num[i,'type']   =  'pnt'
            order_tp_num[i,'ID']     =  as.numeric( clone1$PointMut_ID[i] )
            order_tp_num[i,'order']  =  pnt_clones[[ order_tp_num[i,'ID'] ]]$mut_order   ### pn[order_tp_num[i,'ID'], 'mut_order']
        }
    }    
    
    if ( length( clone1$CNA_ID ) > 0 & clone1$CNA_ID  != 0 ){
        for (j in 1:length( clone1$CNA_ID ) ) {
            order_tp_num[j+i,'ID']      =   clone1$CNA_ID[ j ]
            cn1   =   cna_clones[[ order_tp_num[ j+i, 'ID' ] ]]
            order_tp_num[j+i,'order']   =   cn1$mut_order ###  cn[order_tp_num[ j+i, 'ID' ], 'mut_order']
            order_tp_num[j+i,'type']    =   cn1$dupOrdel  ###  as.character( cn[order_tp_num[ j+i, 'ID' ], 'dupOrdel'] )
        }
    }
    
    if ( length(order_tp_num)  !=  0 ){
        order_tp_num  <-  order_tp_num[ order( order_tp_num$order), ]
        row.names(order_tp_num) <- 1:length( order_tp_num[,'type'])
    } else order_tp_num = NULL
    
    return( order_tp_num )
}

generate_cna  <-  function( prntl, genes, start_end, onco1, dupOrdel ) {
    ### The function to generate object of point mutation cna
    cna0 = CNA_Mutations$new()
    cna0$CNA_ID =  ifelse( length(cna_clones) == 0, 1, 
                                cna_clones[[ length(cna_clones) ]]$CNA_ID + 1 ) 
    cna0$Parental_1or2  =  prntl   
    cna0$dupOrdel = dupOrdel 
    cna0$Chr = gene_map$Chr[ which( gene_map$Gene == genes[1] )[1] ]
    cna0$Ref_start  = start_end[1]
    cna0$Ref_end    = start_end[2]
    cna0$Gene_names = genes
    u = get_u_cna( genes, dupOrdel )
    cna0$MalfunctionedByCNA  =  ifelse( ( runif(1) < u ), TRUE, FALSE )
    
    mut_order  <<-  mut_order  +  1
    cna0$mut_order  =  mut_order
    
    cna_clones <<- c( cna_clones, cna0 )
    
    return( cna0 )
}

### The function to make gene map for the clone1 (add the mutations)
modify_gene_map  <-  function( clone1 , onco1 ){
  # index    =   which( sapply( onco_clones , FUN = function(x)  x$id == clone1$id) )  # index of clone1 in clones list
  # onco1  =  onco_clones[[ index ]]
  gm1    =  gm2  =  gene_map        # locations for the first and second chromosomes 
  gm1$pnts  =  gm2$pnts  = ''
  mixed_order  =  mixed_mut_order( clone1 = clone1 )
  if ( is.null(mixed_order) ) return( list( gm1, gm2 ) )
  
  for (l in 1:nrow( mixed_order ) ) {
      cn1  =  NULL
      if ( mixed_order[l, 'type'] == 'del' ) {
          
          cn1        =  cna_clones[[  mixed_order$ID[l] ]]   # cn[ mixed_order$ID[l], ]
          Ref_start  =  cn1$Ref_start
          Ref_end    =  cn1$Ref_end
          Chr        =  cn1$Chr
          
          # change gene map only one of two chromosomes: 1 or 2 
          ifelse( cn1$Parental_1or2 == 1, gm <- gm1, gm <- gm2 )  
          
          gm  =  add_deletion( gm = gm, Ref_start = Ref_start, 
                               Ref_end = Ref_end, Chr = Chr )
          
          ### come back to gene map for certain chromosome 1 or 2
          ifelse( cn1$Parental_1or2 == 1, gm1 <- gm, gm2 <- gm )
          rm( gm )
      }
      
      if ( mixed_order[l, 'type'] == 'dup' ) {
          
          cn1        =  cna_clones[[  mixed_order$ID[l] ]]   ###  cn[ mixed_order$ID[l], ]
          Ref_start  =  cn1$Ref_start
          Ref_end    =  cn1$Ref_end
          Chr        =  cn1$Chr
          
          # change gene map only one of two chromosomes: 1 or 2 
          ifelse( cn1$Parental_1or2 == 1, gm <- gm1, gm <- gm2 ) 
          
          gm  =  add_duplication( gm = gm, Ref_start = Ref_start, 
                                  Ref_end = Ref_end, Chr = Chr )
          
          ### come back to gene map for the certain chromosome 1 or 2
          ifelse( cn1$Parental_1or2 == 1, gm1 <- gm, gm2 <- gm )
          rm( gm )          
      } 
      
      if ( mixed_order[l, 'type'] == 'pnt' ) {
          
          pn1        =  pnt_clones[[ mixed_order$ID[l] ]]   ###  pn[ mixed_order$ID[l], ]
          pos_pnt1   =  pn1$Ref_pos
          Chr        =  pn1$Chr
          
          # change gene map only one of two chromosomes: 1 or 2 
          ifelse( pn1$Parental_1or2 == 1, gm <- gm1, gm <- gm2 ) 
          
          gm  =  add_pnt_mutation( gm = gm, pos_pnt = pos_pnt1, Chr = Chr )
          
          ### come back to gene map for the certain chromosome 1 or 2
          ifelse( pn1$Parental_1or2 == 1, gm1 <- gm, gm2 <- gm )
          rm( gm ) 
          
      }  
    
  }
  
  return( list( gm1, gm2 ) )
}

### The function to add point mutation to gene map (chromosomal location data frame)
add_pnt_mutation   <-  function( gm = gm, pos_pnt = pos_pnt1, Chr = Chr ){
    
    w = which( gm$Chr == Chr & gm$Start <= pos_pnt  & gm$End >= pos_pnt )
    
    if ( length(w) != 1 ) {
        print( w )
        print( pos_pnt )
        print( gm )
        stop( 'Point mutation does not meet into location ranges' ) 
    } else {
        gm[w, 'pnts']  =  ifelse( gm[w, 'pnts'] ==  '', 
                                  gm[w, 'pnts']  <-  as.character(pos_pnt),  
                                  gm[w, 'pnts']  <-  paste(gm[w, 'pnts'], pos_pnt, sep = ',') )
    }
    
    return( gm )
}

### The function to add deletion to gene map (chromosomal location data frame)
add_deletion  <-  function( gm, Ref_start, Ref_end, Chr ){
    ### Change gene_map with CNA deletion:
    if (Ref_end < Ref_start)  stop( 'End should be larger or equal then Start for CNA' )
    if (Ref_start < min(gm$Start[gm$Chr == Chr]) ) stop( 'The CNA should be inside of the cromosomal locations' )
    w1  =  max( which( gm$Start <= Ref_start  &  gm$Chr == Chr ) )
    w2  =  max( which( gm$Start <= Ref_end    &  gm$Chr == Chr ) )
    
    if ( w1 == w2 ) {
        gm  =  rbind( gm[1:w1,], gm[ w1:nrow(gm), ] )  # duplicate the row w1=w2
        w2  =  w1 + 1
    }
    ### Correction of rows w1 and w2:
    if ( gm[ w1, 'Start']   <   Ref_start ){  # to delete part of row or leave as before
        gm[ w1, 'End']    =   ifelse( Ref_start <=  gm[ w1, 'End'], Ref_start - 1, gm[ w1, 'End']  )
        gm[ w1, 'Len']    =   gm[ w1, 'End']  -  gm[ w1, 'Start']    +  1
        gm[ w1, 'pnts']   =   check_pnts( gm[ w1,  ] )   
    } else {
        w1  =  w1  -  1   # to delete whole row
    }
    
    if ( gm[ w2, 'End']   >   Ref_end ){   # to delete part of row or leave as before
        gm[ w2, 'Start']    =   Ref_end + 1 
        gm[ w2, 'Len']      =   gm[ w2, 'End']  -  gm[ w2, 'Start']  +  1 
        gm[ w2, 'pnts']     =   check_pnts( gm[ w2,  ] )
    } else {
        w2  =  w2  +  1    # to delete whole row 
    }
    
    ### delete part of gm related to deletion
    if ( w2 - w1 > 1 )      gm   =    gm[ -c( (w1+1):(w2-1) ), ]  
    
    # To subtract the delta from positions
    w    =  which( gm$Start > Ref_end & gm$Chr == Chr )
    dlt  =  Ref_end  -  Ref_start  +  1
    if ( length(w) > 0 ){
        gm[w, 'Start']    =  gm[w, 'Start']  -  dlt
        gm[w, 'End']      =  gm[w, 'End']    -  dlt
        gm[w, 'pnts']     =  sapply( w, FUN = function(x)  pnts_add_dlt( gm[ x, ], dlt  =  -dlt)  ) 
    }
    
    return( gm )
}

### The function to add duplication to gene map (chromosomal location data frame)
add_duplication  <-  function( gm, Ref_start, Ref_end, Chr ){

    if (Ref_end < Ref_start)  stop( 'End should be larger or equal then Start for CNA' )
    if (Ref_start < min(gm$Start[gm$Chr == Chr]) ) stop( 'The CNA should be inside of the cromosomal locations' )
    dlt        =  Ref_end  -  Ref_start  +  1  # delta for all next chromosomal positions
    
    ### Change gene_map with CNA duplication:
    w1  =  max( which( gm$Start <= Ref_start  &  gm$Chr == Chr ) )
    w2  =  max( which( gm$Start <= Ref_end    &  gm$Chr == Chr ) )
    
    # ADD rows for duplication 
    if ( (w2+1) <= nrow(gm) ){
        gm  =  rbind( gm[1:w2, ], gm[w1:w2, ] , gm[(w2+1):nrow(gm), ])
    } else {
        gm  =  rbind( gm[1:w2, ], gm[w1:w2, ] ) 
    }
    
    w3  =  w2 + 1

    w5  =  max( which( gm$Chr  ==  Chr ) )
    # change the final position BEFORE duplication
    gm[w2, 'End']    =  ifelse( gm[w2, 'End'] > Ref_end, Ref_end, gm[w2, 'End'] )
    gm[w2, 'Len']    =  gm[w2, 'End'] - gm[w2, 'Start'] + 1
    gm[w2, 'pnts']   =  check_pnts( gm[ w2,  ] )
    
    # Add delta to all positions after duplication
    gm[w3:w5, 'Start']  =  gm[w3:w5, 'Start'] + dlt
    gm[w3:w5, 'End']    =  gm[w3:w5, 'End']   + dlt
    gm[w3:w5, 'pnts']   =  sapply(w3:w5, FUN = function(x) pnts_add_dlt( gm[ x, ], dlt )  ) 
    
    # Correct Start of w3 row if it's necessary 
    if ( gm[w3, 'End']  >= (Ref_start + dlt) ){
        gm[w3, 'Start']  =  Ref_start + dlt
        gm[w3, 'Len']    =  gm[w3, 'End'] - gm[w3, 'Start'] + 1
        gm[w3, 'pnts']   =  check_pnts( gm[ w3,  ] )
    } else { # delete the row if duplication started from outside of gene positions
        gm  =  gm[-w3,]
    }
    
    return( gm )
}

## to check add_del ... add_dupl, please, use daff library and:
## diff_data( gm_ref, gm, show_unchanged_columns = TRUE, always_show_order = TRUE )

### The function to subtract delta from position of point mutations
pnts_add_dlt  <-  function( gm_w1 , dlt, lst = NULL ){
    ### Return the pnts - dlt for one row of gene map
    if ( is.null(gm_w1) )  stop( 'The input is null' )
    pnts = gm_w1$pnts 
    if ( length(pnts)  ==  0 )  stop( 'The input is empty' )
    if ( pnts == '' ) return( '' )
    pnts  =  as.numeric( unlist( strsplit( pnts, split = ',') ) )
    if ( !is.numeric( pnts ) ) stop( 'Incorrect format of points mutation: should be numeric')
    pnts  =  pnts  +  dlt
    pnts  =  paste( as.character( pnts ) , collapse = ',') 
    
    return( pnts )
}

### The function to check what pnts do fall into the range?
check_pnts  <-  function( gm_w1 ){
    ### Return the pnts which fall into the range
    if ( is.null(gm_w1) )  stop( 'The input is null' )
    pnts = gm_w1$pnts 
    if ( length(pnts)  ==  0 )  stop( 'The input is empty' )
    if ( pnts == '' ) return( '' )
    pnts  =  as.numeric( unlist( strsplit( pnts, split = ',') ) )
    if ( !is.numeric( pnts ) ) stop( 'Incorrect format of points mutation: should be numeric')
    
    check_log  =  sapply( pnts, FUN = function( x ) ( x <= gm_w1$End & x >= gm_w1$Start ) )
    pnts  =  paste( as.character( pnts[ check_log ] ) , collapse = ',') 
    pnts  =  ifelse( length(pnts) == 0, '', pnts )
    
    return( pnts )
}
## for tests:  
## Prepare the input and then make test:
## data.frame( test = sapply( X = 1:12 , FUN = function( x ) check_pnts( gm_w1 = gm1[ x, ]) ) )
## data.frame( test = sapply(1:12 , FUN = function( x ) pnts_add_dlt( gm1[ x, ], dlt  =  -dlt ) ) )


# The function to get length of CDS and of genes from gm (gene_map) and related probabilities:
get_cds_rna  <-  function( gm ){
    
    name0  =  unique( gm$Gene )
    rna0  =  cds0  =  NULL
    for ( i in 1:length(name0) ) {
      w     =  which( gm$Gene == name0[ i ] )
      cds0  =  c( cds0, sum( gm$End[w]  -  gm$Start[w] + 1 ) )
      rna0  =  c( rna0, max( gm$End[w] ) - min( gm$Start[w]) + 1  )
    }
    
    sum_prob   =  sum( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) )
    prob       =  c(   m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) ) / sum_prob
    p0         =   (1-m0)**sum(cds0) * (1-m_dup)**sum(rna0) * (1-m_del)**sum(rna0)
    
    return( list( names = name0, CDS = cds0, RNA = rna0, PROB = prob, SUM = sum_prob, P0 = p0 ) )

}

### IV)  MODEL: Functions --------------------------------------------------------


# aggregate
sum_cell <- function(env, clones) {
  if (length(clones) > 0) {
    avg = apply(matrix(unlist(lapply(clones, sum_mutation)),ncol=length(clones)),1,sum)  #  /length(clones)
    env$c = avg[1] / (env$N + env$M)
    env$d = avg[2] / (env$N + env$M)
    env$i = avg[3] / (env$N + env$M)
    env$a = avg[4] / (env$N + env$M)
    env$k = avg[5] / (env$N + env$M)
    env$E = avg[6] / (env$N + env$M)
    env$Nmax = avg[7] / (env$N + env$M)
    env$im = avg[8] / (env$N + env$M)
    env$Ha = avg[9] / (env$N + env$M)
    env$Him = avg[10] / (env$N + env$M)
    env$Hi = avg[11] / (env$N + env$M)
    env$Hb = avg[12] / (env$N + env$M)
    env$Hd = avg[13] / (env$N + env$M)
    env$type = avg[14] / (env$N + env$M)
    env$mutden = avg[15] / (env$N + env$M)
  } else {
    env$M = 0
    env$N = 0
    env$c = 0
    env$d = 0
    env$i = 0
    env$a = 0
    env$k = 0
    env$E = 0
    env$Nmax = 0
    env$im = 0
    env$Ha = 0
    env$Him = 0
    env$Hi = 0
    env$Hb = 0
    env$Hd = 0
    env$type = 0
  }
  # env$posdriver = rep("", length(onco$name))
  # env$pospasngr = rep("", length(onco$name))
}

sum_mutation <- function(clone1) {
  return(c(clone1$c*clone1$N_cells,      clone1$d*clone1$N_cells,    clone1$i*clone1$N_cells, 
           clone1$a*clone1$N_cells,      clone1$k*clone1$N_cells,    clone1$E*clone1$N_cells,
           clone1$Nmax*clone1$N_cells,   clone1$im*clone1$N_cells,   clone1$Ha*clone1$N_cells, 
           clone1$Him*clone1$N_cells,    clone1$Hi*clone1$N_cells,   clone1$Hb*clone1$N_cells, 
           clone1$Hd*clone1$N_cells,     ifelse(clone1$invasion,1,0), 
           clone1$mutden*clone1$N_cells)) # , 
  # clone1$N_cells * ifelse(clone1$invasion,0,1),   clone1$N_cells * ifelse(clone1$invasion,1,0) ))
  
  #           clone1$gene*clone1$N_cells) )
}

# To calculate N and M numbers - normal and metastasis cells
sum_N_M <- function(env, clones) {
  if (length(clones) > 0) {
    avg = apply(matrix(unlist(lapply(clones, number_N_M)),ncol=length(clones)),1,sum)  #  /length(clones)
    env$N = avg[1]
    env$M = avg[2]
    return(env$N + env$M)
  }
}

number_N_M <- function(clone1) {
  return( c(clone1$N_cells * ifelse(clone1$invasion,0,1),   clone1$N_cells * ifelse(clone1$invasion,1,0) ))
}


# trial function
# The resul is a number of new cells, if N_New < 0 it means that the number is decreased.
trial <- function( clone1, onco1 ) { 
    
    # trial for Environmental death of cell 
    N_die  =  calc_binom( 1, clone1$N_cells, clone1$k )   # The number of cells to die due to the Environmental death of cells in clone
    
    # Apoptosis trial
    N_die  =  calc_binom( 1, clone1$N_cells, clone1$a ) 
    
    # invasion / metastasis trial 
    if (clone1$im > 0) {
        if (!clone1$invasion) {
            N_die = calc_binom( 1, clone1$N_cells, ( 1 - clone1$im ) )
            if ( model_name == 'proportional_metastatic') {
                clone1$invasion = ifelse( clone1$im > 0, TRUE, FALSE )  # condition here is related to the model
            } else {
                clone1$invasion = ifelse( clone1$im != 1, TRUE, FALSE )  
            }
        }
    }
    
    # The new number of cells in the clone: 
    clone1$N_cells = ifelse( (clone1$N_cells - N_die) > 0, (clone1$N_cells - N_die) , 0)
    
    N_new = clone1$N_cells   # the initial number to split / before trial / - all cells "want" to split
    
    # Fragmentation restriction trial
    if (clone1$c > 50) {
        N_new = calc_binom(1, N_new, (1 - clone1$i))
    }
    
    # Divide trial
    N_new = calc_binom( 1, N_new, clone1$d )   # The ! FINAL ! number of cells to split, this number is the output value of function "trial"
    
    if ( clone1$N_cells > 0 ) clone1$c = clone1$c + N_new / clone1$N_cells  # How to calculate the counter of division ? - As an average value of the counters for all cells !  
    
    clone1$N_cells = clone1$N_cells + N_new       # The number of cells are increased due to the splitting 
    
    N_new_clones = 0             # The number of new clones arising from clone1
    
    # p = clone1$m * sum(onco$cds) 
    # if (p < 1) N_new_clones = calc_binom(1, 2 * N_new, p ) else N_new_clones = 2 * N_new
    
    N_new_clones = calc_binom( 1, 2 * N_new, 1 - (onco1$p0_1 + onco1$p0_2 ) / 2  )   # 1 and 2 chr
    
    clone1$N_cells = ifelse( ( clone1$N_cells - N_new_clones ) > 0, ( clone1$N_cells - N_new_clones ) , 0 )
    
    return( N_new_clones )
}

# mutagenesis trial
trial_mutagenesis <- function( clone1, num_mut, onco1 ) {
  
                            # num_mut is a number of mutations in this NEW clone1
                            # length of onco - the number of genes
                            # onco1 is onco related to new clone1
    
    ### if (sum(mut1) == 0) stop("The mutation is zero, that is incorrect", call. = TRUE)    # We have known that mutation must occur 
    type_events  =  sample( x = c('point_mut', 'del', 'dup' ), size = num_mut,
                                replace = TRUE, prob = (onco1$prob_1 + onco1$prob_2) / 2 )
    gm  =  modify_gene_map( clone1 , onco1 )
    # print( gm[[1]] )
    # print( gm[[2]] )
    for ( t in type_events ) {
        ### For each type of mutation to generate event of mutation
        if (t == 'point_mut') {
            pm  =  get_point_mutation( onco1, gm )
            prntl = unlist( pm[[1]] )
            gene  = unlist( pm[[2]] )
            pos   = unlist( pm[[3]] )
            Chr   = unlist( pm[[4]] )
            pnt0 = generate_pnt( prntl, gene, pos, onco1, Chr )
            if ( (clone1$PointMut_ID == 0)[1] ) {
                        id   =  pnt_clones[[ length(pnt_clones) ]]$PointMut_ID
                      } else  id   =  c( clone1$PointMut_ID, pnt_clones[[ length(pnt_clones) ]]$PointMut_ID )
            clone1$PointMut_ID  =  id   # pnt_clone is generated in the generate_pnt function
            
            if ( pnt0$MalfunctionedByPointMut ){ 
                clone1$gene[ which( onco$name == gene ) ] = 1
            } else {
                clone1$pasgene[ which( onco$name == gene ) ] = 1
            }
            gm[[ prntl ]]  =  add_pnt_mutation( gm = gm[[ prntl ]], pos_pnt = pos, Chr = Chr )
        }
        
        if (t == 'dup' | t == 'del') {
            cna_mut = get_cna_mutation( onco1, dupOrdel = t , gm = gm)
            prntl =  unlist( cna_mut[[1]] )
            Chr   =  unlist( cna_mut[[2]] )
            genes =  unlist( cna_mut[[3]] )
            start_end   = unlist( cna_mut[[4]] )
            cna0 = generate_cna( prntl, genes, start_end, onco1, t )
            if ( (clone1$CNA_ID == 0)[1] ) {
                id   =  cna_clones[[ length(cna_clones) ]]$CNA_ID
            } else  id   =  c( clone1$CNA_ID, cna_clones[[ length(cna_clones) ]]$CNA_ID )
            clone1$CNA_ID  =  id 
            
            if ( cna0$MalfunctionedByCNA ){ 
                clone1$gene[ sapply(genes, FUN = function(x) which(onco$name == x) ) ] = 1
            } else {
                clone1$pasgene[ sapply(genes, FUN = function(x) which(onco$name == x) ) ] = 1
            }
            
            ### Change the gene_map related chromosome
            ifelse( t == 'dup', 
                    gm[[ prntl ]]  <-  add_duplication( gm = gm[[ prntl ]], Ref_start = start_end[1], Ref_end = start_end[2], Chr = Chr ),
                    gm[[ prntl ]]  <-     add_deletion( gm = gm[[ prntl ]], Ref_start = start_end[1], Ref_end = start_end[2], Chr = Chr )  ) 
            
            ### Check what point mutations match into the CNA 
            sp   = FALSE
            sp_A = FALSE
            if ( clone1$PointMut_ID != 0 ){
                sp = sapply( clone1$PointMut_ID , FUN = function( x )  {
                              chk_pnt_mut( pnt1  =  pnt_clones[[ x ]], Ref_start = start_end[1], 
                                           Ref_end = start_end[2], Chr = Chr, prntl  =  prntl )
                          })
                ### to check original allele A do/don't match into CNA:
                prntl_inv  =  ifelse( prntl  ==  1, 2, 1 )
                sp_A = sapply( clone1$PointMut_ID , FUN = function( x )  {
                  chk_pnt_mut( pnt1  =  pnt_clones[[ x ]], Ref_start = start_end[1], 
                               Ref_end = start_end[2], Chr = Chr, prntl  =  prntl_inv )
                })
                
                
            }
            
            if ( any( sp ) | any( sp_A ) ){
                ### Before changing point mutation for NEW clone
                ### we have to copy it and avoid any changes in OTHER (parental) clones
                
                ### Copy pnt mutations matched into CNA
                
                ### circle and function to copy pnt for a NEW clone
                for( q in clone1$PointMut_ID[ sp | sp_A ] ){
                    # pnt_clone is generated in the generate_to_copy_pnt function
                    pnt0 = generate_to_copy_pnt( pnt = pnt_clones[[ q ]] )
                    x  =  which( clone1$PointMut_ID == q )
                    clone1$PointMut_ID[ x ]  =  pnt0$PointMut_ID
                }
                ### end of the new part of the code
                
                sapply( clone1$PointMut_ID[ sp ], 
                        FUN = function( x ) change_pnt_by_cna( pnt1  =  pnt_clones[[x]], start_end, t )  )
                sapply( clone1$PointMut_ID[ sp_A ], 
                        FUN = function( x ) change_allele_A_by_cna( pnt1  =  pnt_clones[[x+1]], start_end, t )  )
                
            }
            
        }
    }
    
    onco_update( onco1, gm )

}

init_pnt_clones   <- function( clones, onco_clones ) {
    
    if ( is.null( clones ) )  return( NULL )
    
    for( i in 1:length( clones ) ){
        clone1  =  clones[[ i ]]
        onco1   =  onco_clones[[ i ]]
        
        genes  =  onco1$name[ clone1$gene == 1 ]
        gm  =  modify_gene_map( clone1 , onco1 )
        
        if ( length( genes ) == 0 ) return( stop('Length of mutated genes should be non-zero ') )
        
        for( gene in genes ){
            
            pm  =   get_point_mutation_for_gene( onco1, gm_1_2 = gm, gene )   #  get_point_mutation( onco1, gm )
            prntl = unlist( pm[[1]] )
            # gene  = unlist( pm[[2]] )
            pos   = unlist( pm[[3]] )
            Chr   = unlist( pm[[4]] )
            pnt0 = generate_pnt( prntl, gene, pos, onco1, Chr, mutation = TRUE )
            
            ### Add pnt mutation ID to a clone:
            if ( (clone1$PointMut_ID == 0)[1] ) {
              id   =  pnt_clones[[ length(pnt_clones) ]]$PointMut_ID
            } else  id   =  c( clone1$PointMut_ID, pnt_clones[[ length(pnt_clones) ]]$PointMut_ID )
            clone1$PointMut_ID  =  id
        }
    }
    
}

### The function to change the point mutation due to CNA:
change_pnt_by_cna  <-  function( pnt1, start_end, t ) {
    
    ntrs  =  intersect( which( pnt1$Phys_pos >= start_end[1] ), 
                        which( pnt1$Phys_pos <= start_end[2] )  )
    cf    =  length( ntrs )  ### coefficient for CN  
    
    if ( pnt1$Copy_number != 0 )  {
        pn_cn  =  pnt1$Copy_number + cf * ifelse(t == 'dup', 1, -1 )
        pnt1$Copy_number  =  ifelse( pn_cn >= 0, pn_cn, 0 )
    }
    
    pos_pnt       =   pnt1$Phys_pos[ ntrs ]
    
    dlt  =  ifelse( t == 'dup', start_end[2]  -  start_end[1], start_end[1]  -  start_end[2] )
    pnt1$Phys_pos   =  c( pnt1$Phys_pos, pos_pnt + dlt )
    pnt1$Delta      =    pnt1$Phys_pos  -  pnt1$Ref_pos

}

### The function to change copy number of the allele A 
###     of the point mutation at the allele B due to CNA:
change_allele_A_by_cna  <-  function( pnt1, start_end, t ) {
    if ( pnt1$Copy_number != 0 )  {
        pn_cn  =  pnt1$Copy_number + ifelse(t == 'dup', 1, -1 )
        pnt1$Copy_number  =  ifelse( pn_cn >= 0, pn_cn, 0 )
    }
}

### The function to update onco1 after mutation ( used in trial_mutagenesis )
onco_update  <-  function( onco1, gm ){
    
    lst1  =  get_cds_rna( gm[[1]] )
    rd1   =  as.integer( sapply( onco1$name, FUN = function(x) which( x  ==  lst1[[1]]) ) )
    
    lst2  =  get_cds_rna( gm[[2]] )   
    rd2   =  as.integer( sapply( onco1$name, FUN = function(x) which( x  ==  lst2[[1]]) ) )
    
    # change the onco1 related to new gene_map:
    onco1$cds_1  =  lst1$CDS[ rd1 ]
    onco1$cds_2  =  lst2$CDS[ rd2 ]
    
    onco1$rna_1  =  lst1$RNA[ rd1 ]
    onco1$rna_2  =  lst2$RNA[ rd2 ]
    
    onco1$prob_1 =  lst1$PROB
    onco1$prob_2 =  lst2$PROB
    
    onco1$sum_prob_1  =  lst1$SUM
    onco1$sum_prob_2  =  lst2$SUM
    
    onco1$p0_1   =  lst1$P0
    onco1$p0_2   =  lst2$P0
    
    return( onco1 )
}

### The function to check point mutations match or don't match into duplication or deletion
chk_pnt_mut  <-  function( pnt1 , Ref_start, Ref_end, Chr, prntl ){
    
    for( X in pnt1$Phys_pos ){
        if ( pnt1$Chr == Chr &  X <= Ref_end & X >= Ref_start & pnt1$Parental_1or2 == prntl ) {
            return( TRUE )}
    }
    return( FALSE )
}

# to make one copy for clone1 in clone_init function
clone_copy <- function(clone1) {
    env$last_id = env$last_id + 1
    
    return(clone$new(id=env$last_id, parent=clone1$id, c=clone1$c, d=clone1$d, 
                        i=clone1$i, mutden=clone1$mutden, a=clone1$a, 
                        k=clone1$k, E=clone1$E, Nmax=clone1$Nmax, 
                        gene=clone1$gene, pasgene=clone1$pasgene, 
                        PointMut_ID = clone1$PointMut_ID, CNA_ID = clone1$CNA_ID,
                        # posdriver=clone1$posdriver, pospasngr=clone1$pospasngr, 
                        invasion=clone1$invasion, s=clone1$s, birthday=env$T)) #, 
                        # len = clone1$len ))
}

# to make one copy for onco1 in init_onco_clones function
onco_copy <- function( onco1 ){
    
    onco2 = oncogene$new(id = onco1$id, name = onco1$name, onsp = onco1$onsp, len = onco1$len,
                         cds_1 = onco1$cds_1, cds_2 = onco1$cds_2, 
                         rna_1 = onco1$rna_1, rna_2 = onco1$rna_2,
                         p0_1 = onco1$p0_1, p0_2 = onco1$p0_2, 
                         prob_1 = onco1$prob_1, prob_2 = onco1$prob_2,
                         sum_prob_1 = onco1$sum_prob_1, sum_prob_2 = onco1$sum_prob_2  )
    return( onco2 )
}

# write log file
write_log <- function(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, 
                      E0, F0, m0, uo, us, s=10, k, 
                      m_dup, m_del, lambda_dup, lambda_del, # CNA parameters
                      uo_dup, us_dup, uo_del, us_del,       # CNA parameters
                      censore_n, censore_t, d0, Compaction_factor, model_name, time_stop ) {
    data <- c("genefile", "clonefile", "geneoutfile", "cloneoutfile", "logoutfile",
              "E", "F", "m0", "uo", "us", "s", "k", 
              "m_dup", "m_del", "lambda_dup", "lambda_del",
              "uo_dup", "us_dup", "uo_del", "us_del",
              "censore_n", "censore_t", "d0", 'Compaction_factor', 'model_name', 'time_stop')
    data <- rbind( data, c(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, 
                             E0, F0, m0, uo, us, s, k, 
                             m_dup, m_del, lambda_dup, lambda_del, # CNA parameters
                             uo_dup, us_dup, uo_del, us_del,       # CNA parameters
                             censore_n, censore_t, d0, Compaction_factor, model_name, time_stop ) )
    write(data, logoutfile, ncolumn=2, sep="\t")
}


write_geneout <- function(outfile, hall) {
    data <- c(onco$name[hall$Ha], onco$name[hall$Hi], onco$name[hall$Hd], onco$name[hall$Hb], onco$name[hall$Him])
    data <- rbind(data, c(rep("apoptosis", length(onco$name[hall$Ha])),
                          rep("immortalization", length(onco$name[hall$Hi])),
                          rep("growth|anti-growth", length(onco$name[hall$Hd])),
                          rep("angiogenesis", length(onco$name[hall$Hb])),
                          rep("invasion", length(onco$name[hall$Him]))))
    data <- rbind(data, c(hall$Ha_w, hall$Hi_w, hall$Hd_w, hall$Hb_w, hall$Him_w))
    data <- rbind(data, c(onco$onsp[hall$Ha], onco$onsp[hall$Hi], onco$onsp[hall$Hd], onco$onsp[hall$Hb], 
                          onco$onsp[hall$Him]))
    write(data, outfile, ncolumn=4, sep="\t")
}

write_header <- function(outfile, env, onco) {
    header <- c('Time', 'N_cells', 'AvgOrIndx', 'ID', 'ParentID', 'Birth_time', 'c', 'd', 'i', 'im', 'a',
                'k', 'E', 'N', 'Nmax', 'M', 'Ha', 'Him', 'Hi', 'Hd', 'Hb', 'type', 'mut_den',
                # paste("PosDriver:", onco$name, sep=""), paste("PosPasngr:", onco$name, sep="") )      #   , 'Clone number', 'Passengers Clone number', 'Mix Clone number')
                'driver_genes', 'passenger_genes', 
                'PointMut_ID', 'CNA_ID' , 'onco_ID', 
                paste('CDS_', onco$name, sep=''), paste('Len_', onco$name, sep=''), 
                'p0', 'prob_point_mut', 'prob_del', 'prob_dup' )
    write(header, outfile, append=FALSE, ncolumn=length(header), sep="\t")
}

write_cloneout <- function( outfile, env, clones, isFirst, onco_clones ) {
    data <- c(env$T, '-', 'avg', '-', '-', '-', env$c, env$d, env$i, env$im, env$a, env$k, env$E, env$N,
              env$Nmax, env$M, env$Ha, env$Him, env$Hi, env$Hd, env$Hb, env$type, env$mutden,
              rep('-',17) )
    
    write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
       
    if (length(clones) > 0 & isFirst) {
        for (i in 1:length(clones)) {
            clone1 = clones[[i]]
            onco1  = onco_clones[[i]]
            data <- c(env$T, clone1$N_cells, i, clone1$id, clone1$parent, clone1$birthday, clone1$c, clone1$d, 
                      clone1$i, clone1$im, clone1$a, clone1$k, clone1$E, env$N, clone1$Nmax, env$M,
                      clone1$Ha, clone1$Him, clone1$Hi, clone1$Hd, clone1$Hb, ifelse(clone1$invasion,1,0), 
                      clone1$mutden, 
                      paste(clone1$gene, collapse =  ' '), paste(clone1$pasgene, collapse =  ' '),
                      paste(clone1$PointMut_ID, collapse = ', '), paste(clone1$CNA_ID, collapse = ', '), 
                      onco1$id, onco1$cds_1, onco1$rna_1, onco1$p0_1, onco1$prob_1 )     #  clone1$posdriver, clone1$pospasngr)           #   , vc, pasvc, mixvc)
            write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
        }
    }
}

write_weights <- function(outfile, hall) {
    #data <- c("Hallmarks", "Designation", onco$name)
    data <- data.frame( "Gene" = onco$name)   
    data$Gene <-   as.character(data$Gene)
    
    w <- rep(0.0, onco$len)
    for (j in 1:onco$len) { if ( length( which(j==hall$Ha) ) > 0  ) w[j] = hall$Ha_w[which(j==hall$Ha)]  }
    data <- cbind(data, "Apoptosis ($H_a$)" = w)
    
    w <- rep(0.0, onco$len)
    for (j in 1:onco$len) { if ( length( which(j==hall$Hb) ) > 0  ) w[j] = hall$Hb_w[which(j==hall$Hb)]  }
    data <- cbind(data, "Angiogenesis ($H_b$)" = w)
    
    w <- rep(0.0, onco$len)
    for (j in 1:onco$len) { if ( length( which(j==hall$Hd) ) > 0  ) w[j] = hall$Hd_w[which(j==hall$Hd)]  }
    data <- cbind(data, "Growth / Anti-growth ($H_d$)" = w)
    
    w <- rep(0.0, onco$len)
    for (j in 1:onco$len) { if ( length( which(j==hall$Hi) ) > 0  ) w[j] = hall$Hi_w[which(j==hall$Hi)]  }
    data <- cbind(data, "Immortalization ($H_i$)" = w)
    
    w <- rep(0.0, onco$len)
    for (j in 1:onco$len) { if ( length( which(j==hall$Him) ) > 0  ) w[j] = hall$Him_w[which(j==hall$Him)]  }
    data <- cbind(data, "Invasion / Metastasis ($H_{im}$)" = w)
    
    write.table(data, outfile, sep="\t", row.names = FALSE)
}

# write the point mutation info for all clones for all time steps:
write_pnt_clones <- function(pnt_clones, file_out = 'Output/point_mutations.txt'){
    # if ( is.null(pnt_clones) ) return( print( 'Empty data.' ) )
    pn  <-  NULL
    if ( !is.null( pnt_clones ) ){
        for (i in 1:length(pnt_clones)) {
            pnt1 <-  unlist( pnt_clones[[i]] )
            pn1  <-  safe_pnt_mut( pnt1 )     ### pnt1$safe()
            pn   <-  rbind( pn, pn1)
        }
    }
    
    write.table( pn, file = file_out, append=FALSE, sep="\t", row.names = FALSE)
}

# initial clone setting
init_clones <- function(clonefile, clone1) {
    mpos <- regexpr("\\.", clonefile)[1]
    if (mpos != -1) {
        name <- substr(clonefile, 1, mpos - 1)
    } else {
        name <- clonefile
    }
    clones = NULL
    n <- as.numeric(name)
    if (!is.na(n) && is.numeric(n)) {
        factor = n / sum(clone1$m*onco$cds_1)
        f2 = 1.0
        while (TRUE) {
            if (sum(floor(clone1$m*onco$cds_1*factor*f2 + 0.5)) >= n) {
                break
            }
            f2 = f2 + 0.1
        }
        nums = floor(clone1$m*onco$cds_1*factor*f2 + 0.5)
        clones = NULL
        for (i in 1:n) {
            clones = c(clones, clone_copy(clone1))
        }
        pos = 0
        for (i in 1:length(nums)) {
            if (nums[i] > 0) {
                for (j in 1:nums[i]) {
                    if (pos + j <= n) {
                        clones[[pos + j]]$gene[i] = 1
                    }
                }
                pos = pos + nums[i]
            }
        }
    } else {
        data = read.table(clonefile, sep="\t")
        n <- nrow(data)
        
        for (i in 1:n) {
            clone2 = clone_copy(clone1)
            p <- match(onco$name, str_trim(strsplit(as.character(data[i,2]),",")[[1]]))
            clone2$gene[seq(1,length(onco$name))[!is.na(p)]] = 1
            clone2$N_cells = as.numeric(data[i,3])
            clones = c(clones, clone2)
        }
    }
    for (i in 1:n) {
        clones[[i]]$id = i
        clones[[i]]$parent = 0
        clones[[i]]$birthday = 0
        # clones[[i]]$posdriver = ifelse(clones[[i]]$gene == 1,
        #                               paste(ceiling(runif(onco$len)*onco$cds),"0",sep = ":"),
        #                               clones[[i]]$posdriver)
        clones[[i]]$calcMutden()
        clones[[i]]$calcApoptosis()
    }
    env$last_id = n
    return(as.list(clones))
}

### initial onco settings for all clones (onco_clones)
init_onco_clones <- function( onco1, clones ) {
    # initialization of onco_clones in order to assign an onco to each clone
    onco_clones = NULL
    for ( i in 1:length( clones ) ) {
        clone1 = clones[[i]]
        onco_clone2 = onco_copy( onco1 )
        onco_clone2$id = clone1$id
        onco_clones = c( onco_clones, onco_clone2 )
    }
    
    return( as.list( onco_clones ) )
}

### Function to calculate binominal distribution including BIG NUMBERS like 10^12 and more using approximation with normal distribution 
calc_binom <- function(tr,n,p){
    if (n*p < 10^8){
        ou <- rbinom(tr,n,p) 
    } else {
        m <- n * p
        s <- sqrt(  n*p* (1-p)  )
        ou <- rnorm(tr,mean = m, sd = s)
    }
    
    return(  round( ou )  )
}

### V)   MODEL: main function 'model' --------------------------------------------

model <- function(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, 
                  E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0) {
    write_log(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, 
              E0, F0, m0, uo, us, s=10, k=k0, 
              m_dup, m_del, lambda_dup, lambda_del, # CNA parameters
              uo_dup, us_dup, uo_del, us_del,       # CNA parameters
              censore_n, censore_t, d0, Compaction_factor, model_name, time_stop)   # write input parameters
    onco = oncogene$new()        # make the vector onco about the hallmarks
    onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
    hall = hallmark$new()        # make a vector hall with hallmarks parameters
    hall$read( genefile, onco$name, normalization = TRUE )     # read from the genefile - 'gene_hallmarks.txt'
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
    write_geneout(geneoutfile, hall)                  # write the geneout.txt file with initial hallmarks 
    write_weights("Output/Weights.txt", hall)                 # write the weights of genes for hallmarks 
    write_header( cloneoutfile, env, onco )                   # 
    
    cells_number <- sum_N_M(env, clones)                 # to calculate cells numbers - N,M 
    init_pnt_clones( clones, onco_clones )              # initialization of pnt_clones for point mutations
    
    lapply(clones,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    hall$updateEnviron(env, clones)                     # make averaging for cells 
    isFirst = TRUE
    write_cloneout( cloneoutfile, env, clones, isFirst, onco_clones )     #  write initial clones

    print( paste0("The probability of an absence of the mutations is p0 = ", as.character(onco$p0_1) )) 
    time_start  =  Sys.time()
    time_current  =  Sys.time()
    while(length(clones) > 0 && censore_n > cells_number && 
          env$T < censore_t  && 
          ( as.numeric( time_current - time_start ) < time_stop ) ){

        k_old = length(clones)          # the number of clones from last step
        
        clones_new <- NULL 
        onco_clones_new <- NULL
        
        N_clones_new = unlist( mapply( trial, clones, onco_clones ) ) 
        
        survived_clones = NULL
        
        for (i in 1:k_old) {
            if (N_clones_new[i] > 0) { 
                for (j in 1:N_clones_new[i])  {
                    clones_new = c(clones_new,clone_copy(clones[[i]]) ) 
                    onco_clones_new = c(onco_clones_new, onco_copy(onco_clones[[i]]))
                    onco_clones_new[[length(onco_clones_new)]]$id = clones_new[[length(clones_new)]]$id
                }
            }
            
            # To delete the clones with N_cells == 0 because they are died 
            if (clones[[i]]$N_cells == 0 )  survived_clones = c(survived_clones, FALSE)  else survived_clones = c(survived_clones, TRUE )   
        }
        
        
        # The number of mutations for each NEW clone
        N_new <- length(clones_new)
        
        if ( N_new > 0) {
            num_mut <- numeric(0)
            sm <- unlist( lapply(onco_clones_new, function(x) (x$sum_prob_1+x$sum_prob_2)/2 ) ) # sum of prob for new clones
            num_mut <- rztpois(N_new, sm ) # Numbers of mutations for each new clone
            # To apply the mutagenesis only to new clones with a number of mutations: 
            for ( nn in 1:N_new )  {
                trial_mutagenesis( clones_new[[nn]], num_mut[nn], onco_clones_new[[nn]]  )
            }
        }
        
        # the new generation = the survived clones + new_clones 
        clones = c(clones[survived_clones],clones_new)
        onco_clones = c(onco_clones[survived_clones],onco_clones_new)
        
        cells_number <- sum_N_M(env, clones)                 # to calculate cells numbers - N,M for next step
        lapply(clones,update_Hallmarks) 
        hall$updateEnviron(env, clones)                      # to average probabilities and hallmarks
        
        env$T = env$T + 1                                    # to next step
        
        write_cloneout( cloneoutfile, env, clones, isFirst, onco_clones )
        #print(c(env$T,env$N,env$M,env$last_id, length(clones), "N_clones_new = ", N_clones_new))
        
        time_current  =  Sys.time()

    }
    
    # write_pnt_clones( pnt_clones, file = 'Output/point_mutations.txt' )
    # write_pnt_clones( cna_clones, file = 'Output/CNA_mutations.txt' )
    
    return( list( clones , onco_clones ) )
}

