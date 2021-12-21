

# Libraries ---------------------------------------------------------------

library(stringr)
library(ape)
# library(ggplot2)
# library(ggtree)
options(stringsAsFactors=FALSE)

par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5), tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))


# Analyze data block ------------------------------------------------------

# function to read file
read_file  <-  function( file_name = '', stringsAsFactors = FALSE ){
    if ( file.size( file_name )  < 10 ) return( NULL )
    return( read.table( file = file_name, stringsAsFactors  =  stringsAsFactors ,
                        sep="\t", header = TRUE ))
}

### The function to get data about last simulation
get_flow_data <- function(cloneoutfile, genefile ) {

    # Get data about onco and hallmarks
    onco   <<- oncogene$new()        # make the vector onco about the hallmarks
    onco$read( genefile )          # read the input info to the onco from genefile - 'gene_cds2.txt'
    hall   <<- hallmark$new()        # make a vector hall with hallmarks parameters
    hall$read( genefile, onco$name )     # read from the genefile - 'gene_cds2.txt'
  
    ## Get data from output file    
    data_out  <-  read_file( file_name = cloneoutfile )   # read.table( cloneoutfile, sep="\t", header = TRUE )
    data_out[is.na(data_out)]  <-  ""
    
    names(data_out)
    
    # average data
    data_avg  <<-  data_out[ which( data_out$AvgOrIndx  ==  "avg"), ]
    
    # data without averaging - flow data
    data_flow  <<-  data_out[ which( !data_out$AvgOrIndx == "avg" ), ]
    
    # the data of the last time step 
    time_max   <<-  max( data_flow$Time )
    data_last  <<-  data_flow[ which( data_flow$Time  ==  time_max ), ]
    rm( data_out )
    
    cna_mut  <<-  read_file( file_name = 'Output/CNA_mutations.txt' )
      # read.table( file = 'Output/CNA_mutations.txt', sep = '\t', header = TRUE )
    pnt_mut  <<-  read_file( file_name = 'Output/point_mutations.txt' )
      #read.table( file = 'Output/point_mutations.txt', sep = '\t', header = TRUE )
    
}  

### the function to plot main data
plot_average_simulation_data  <-  function(){
    # plot number of cells:
  
    flnm =  readline(prompt=" Press Enter to skip saving to a file or enter name of file to save data ")
    if ( flnm != '' ){
        safe_pdf = TRUE
    } else safe_pdf = FALSE
    
    g_range_y  =  range( 0, data_avg$N + 1 , data_avg$M + 1 )
    g_range_x  =  range( min( data_avg$Time ), max( data_avg$Time ) )
    plot_2D_lines( x = data_avg$Time, data_avg, 
                   nl = c('N', 'M'), xr = g_range_x, yr = g_range_y,
                   legend_names = c('Primary tumor', 'Metastatic'),
                   names = c( 'Time step', 'Number of cells' ) ,
                   safe_pdf  =  safe_pdf, 
                   filename = paste0( flnm, '_NM.pdf' )  )
    rl =  readline(prompt="This is a plot for Numbers of primary tumor and Metastasis cells - Press Enter  ")
    
 
    g_range_y  =  range( 0, 1 )
    plot_2D_lines( x = data_avg$Time, data_avg, 
                   nl = c('a', 'd', 'i', 'im', 'k'), xr = g_range_x, yr = g_range_y,
                   legend_names = '',
                   names = c( 'Time step', 'Probabilities' ) ,
                   safe_pdf  =  safe_pdf, 
                   filename = paste0( flnm, '_probabilities.pdf' )  )
    rl =  readline(prompt="This is a plot for average probabilities - Press Enter  ")
    
    
    g_range_y  =  range( 0, max( data_avg[c('Ha', 'Hd', 'Hi', 'Him', 'Hb')] ) )
    plot_2D_lines( x = data_avg$Time, data_avg, 
                   nl = c('Ha', 'Hd', 'Hi', 'Him', 'Hb'), xr = g_range_x, yr = g_range_y,
                   legend_names = '',
                   names = c( 'Time step', 'Hallmarks' ) ,
                   safe_pdf  =  safe_pdf, 
                   filename = paste0( flnm, '_hallmarks.pdf' )  )
    
    rl =  readline(prompt="This is a plot for average values of Hallmarks - Press Enter")
    
    
    g_range_y  =  range( 0, max( data_avg[c('mut_den')] ) + 0.05 )
    plot_2D_lines( x = data_avg$Time, data_avg, 
                   nl = c('mut_den'), xr = g_range_x, yr = g_range_y,
                   legend_names = '',
                   names = c( 'Time step', 'Mutation rate' ) ,
                   safe_pdf  =  safe_pdf, 
                   filename = paste0( flnm, '_mutation_rate.pdf' )  )
    
    rl =  readline(prompt="This is a plot for average value of mutation rate - Press Enter")
    
    
    
    
    
}


### the function to get order of genes' dysfunction:

get_order_of_genes_dysfunction  <-  function(){
    
    # as.numeric( unlist( str_split( data_last[ , 'PointMut_ID' ], pattern = ',' ) ) )
    
    # sum( as.numeric( unlist( str_split( data_last[ 1, 'driver_genes' ], pattern = ' ' ) ) ) ) != 0 
    
    ch = sapply(X = 1:nrow( data_last ), FUN = function( x )  sum( as.numeric( unlist( str_split( data_last[ x , 'driver_genes' ], pattern = ' ' ) ) ) ) != 0 )
    genes_dysfunction  =  data_last[ ch, c('N_cells', 'ID', 'ParentID', 'Birth_time',
                                           'type', 'mut_den', 'driver_genes', 'passenger_genes',
                                           'PointMut_ID', 'CNA_ID' )  ]
    if ( nrow( genes_dysfunction ) > 0  ) {
        for (i in 1:nrow( genes_dysfunction ) ){
            
            drivers  =  as.numeric( unlist( str_split( genes_dysfunction$driver_genes[ i ] , pattern = ' ') )  )
            # nm       =  onco$name[ as.logical(drivers) ]
            
            ### Point mutations:
            id_pnt_mut    =  as.numeric( unlist( str_split( genes_dysfunction$PointMut_ID[ i ] , pattern = ',') ) )
            if ( id_pnt_mut[1] != 0 ){
                r = as.vector( sapply( X = id_pnt_mut, FUN = function( x ) which( pnt_mut$PointMut_ID ==   x  ) ) )
            
                l = pnt_mut[ r , 'MalfunctionedByPointMut' ]
                l[ is.na( l ) ] = FALSE 
                nm_pnt  =  pnt_mut[ r, 'Gene_name'][ l ]
                order_pnt  =  pnt_mut[ r, 'mut_order'][ l ]
            } else{
                nm_pnt  =  NULL
                order_pnt  =  NULL
            }
            
            ### CNA mutations:
            id_cna   =  as.numeric( unlist( str_split( genes_dysfunction$CNA_ID[ i ] , pattern = ',') ) )
            if ( id_cna[1] != 0 ){
                r_cna    =  as.vector( sapply( X = id_cna, FUN = function( x ) which( cna_mut$CNA_ID ==   x  ) ) )
                l_cna    =  cna_mut[ r_cna , 'MalfunctionedByCNA' ]
                l_cna[ is.na( l_cna ) ] = FALSE 
                nm_cna   =  cna_mut[ r_cna, 'Gene_names' ][ l_cna ]
                order_cna  =  cna_mut[ r_cna, 'mut_order'][ l_cna ]
            } else {
                nm_cna    =  NULL
                order_cna =  NULL
            }
            
            # Combine together:
            nm  =  c( nm_pnt, nm_cna )
            rdr = c( order_pnt, order_cna )
            genes_dysfunction[ i, 'order' ]  =  paste( nm[ order( rdr ) ], collapse = ' -> ')
        }
    } else { return( NULL ) }
    
    write.table( genes_dysfunction, file = './Output/order_genes_dysfunction.txt', 
               sep = '\t', row.names = FALSE, col.names = TRUE, append = FALSE )
    cat('Order of genes dysfunction saved to the file in Output directory')
    
    return( genes_dysfunction )
}

get_VAF  <-  function(){
    
    pnt_mut_A  <<-  pnt_mut[ which(  is.na( pnt_mut$MalfunctionedByPointMut ) ) , ]
    pnt_mut    <<-  pnt_mut[ which( !is.na( pnt_mut$MalfunctionedByPointMut ) ) , ]
    
    ids  =  str_split( data_last$PointMut_ID, pattern = ',' )
    ids =  sapply( X = 1:nrow( data_last), FUN = function(x) as.numeric( ids[[ x ]] ) )
    
    nqu  =  unique( unlist( ids ) )
    if ( nqu[1] == 0 ) nqu = nqu[ -1 ]
    
    # start VAF 
    VAF  =  NULL
    
    if ( length( nqu ) == 0 ) return( NULL )
    
    for( j in 1:length( nqu ) ){
        wc  =  sapply( X = 1:length( ids ), FUN = function( x ) is.element( nqu[j] , ids[[ x ]] ) )
        
        VAF_1  =  pnt_mut[ which( pnt_mut$PointMut_ID == nqu[ j ] ) , ]
        
        ### number of cells N - primary cells, M - metastasis cells
        if ( any( data_last[ which( wc ), 'type' ] == 0 ) ){
            VAF_1$N  =  sum( as.numeric( data_last[ which( wc & data_last$type ==  0 ),  'N_cells'  ] ) )
        } else VAF_1$N  =  0
        
        if ( any( data_last[ which( wc), 'type' ] == 1 ) ){
            VAF_1$M  =  sum( as.numeric( data_last[ which( wc & data_last$type ==  1 ),  'N_cells'  ] ) )
        } else VAF_1$M  =  0
        
        # add copy number of original allele A:
        VAF_A  =  pnt_mut_A[ which( pnt_mut_A$PointMut_ID == nqu[ j ] ) , ]
        
        VAF_1$Copy_number_A  =  VAF_A$Copy_number
        
        VAF  =  rbind( VAF, VAF_1)
    }
    
    if ( any( data_last[ , 'type' ] == 0 ) ) {
        VAF$N_total = sum( as.numeric( data_last[ which( data_last$type == 0 ), 'N_cells' ] ) )  
    } else  VAF$N_total  =  0 
    
    if ( any( data_last[ , 'type' ] == 1 ) ) {
        VAF$M_total = sum( as.numeric( data_last[ which( data_last$type == 1 ), 'N_cells' ] ) )  
    } else  VAF$M_total  =  0
    
    ### Save VAF to the file:
    write.table( VAF,file = "Output/VAF_data.txt", append = FALSE, row.names = FALSE, sep="\t" )
    print("VAF data for allele B and A is saved to the file `Output/VAF_data.txt` ")
    return( VAF )
}

