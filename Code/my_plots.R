

# Libraries ---------------------------------------------------------------

library( randomcoloR )
# library(fishplot)


library(stringr)
library(ape)
library(ggplot2)
# library(ggtree)



# Define default parameters of the plots ----------------------------------

par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5), tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))

cl = randomColor(count = 20, 
                 hue = c(" ", "random", "red", "orange", "yellow", 
                         "green", "blue", "purple", "pink", "monochrome")[1], 
                 luminosity = c(" ", "random", "light", "bright", "dark")[5] )

basic_cls  =  c( 'blue3', 'darkmagenta', 'red', 'green4', 
                 'darkorange', 'steelblue1' )


# Plot 2D --------------------------------------------------------------------

### to plot 2D figure of points y = y(x)
plot_2D   <-  function( x, y, names = c( 'X', 'Y' ), pch = 18, col = 'blue', cex = 1.2, 
                        xr = c(-10,10), yr = c(-10,10), 
                        safe_pdf = FALSE, filename = './plot.pdf' ){
    rp = 1 
    if ( safe_pdf )    {
        pdf( filename, width = 8, height = 8 )
        rp = 2
    }
    for( i in 1:rp ){
        par( mgp = c(2.2, 0.5, 0), font.axis = 2, font.lab = 2 )
        plot( x, y, pch = pch, xlab = names[1], xlim = xr, ylim = yr,
              ylab = names[2], axes = FALSE, cex.lab = 1.5, col = col,
              cex = cex )
        
        axis( 1, font = 2, tck = 0.03, cex.axis = 1.4 )
        axis( 2, font = 2, tck = 0.03, cex.axis = 1.4)    
        axis( 3, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE )
        axis( 4, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE ) 
        
        if ( safe_pdf && i == 1 )      dev.off( )
    }
}

### to plot 2D figure of lines  yi = DF[ , nl[i] ), i - index
plot_2D_lines   <-  function( x, DF, nl = 1:2, names = c( 'X', 'Y' ), 
                              legend_names = '',
                               col = basic_cls, cex = 1.2, lwd = 2.0, 
                              lt = c(1:6), xr = c(-10,10), yr = c(-10,10), 
                        safe_pdf = FALSE, filename = './plot.pdf', 
                        type = 'l' , logscale = '' , draw_key  =  TRUE ){
    rp = 1 
    if ( safe_pdf )    {
        pdf( filename, width = 8, height = 8 )
        rp = 2
    }
    for( i in 1:rp ){
        par( mgp = c(2.2, 0.5, 0), font.axis = 2, font.lab = 2 )
        ### Plot the first line:
        y = DF[, nl[1] ]
        plot( x, y, xlab = names[1], xlim = xr, ylim = yr, 
              ylab = names[2], axes = FALSE, cex.lab = 1.5, col = col[1],
              lwd = lwd, lty = lt[1], 
              type = type, log = logscale )
        
        axis( 1, font = 2, tck = 0.03, cex.axis = 1.4 )
        axis( 2, font = 2, tck = 0.03, cex.axis = 1.4)    
        axis( 3, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE )
        axis( 4, font = 2, tck = 0.03, cex.axis = 1.4, labels = FALSE ) 
        
        ### More than 1 line: 
        if ( length( nl ) > 1 ){
            for( j in 2:length(nl) ){
                st = nl[j]
                lines( x, DF[, st ], lwd = lwd, col = col[j], lty = lt[j] )
            }
        }
        if ( draw_key ){
            key = names( DF[ nl ] )
            if( legend_names[1] != '') key = legend_names 
            legend( x = 'bottom', legend = key,  
                    horiz = TRUE, xpd = TRUE,  inset = c(0, 1.03), box.col = "white", 
                    lty = lt[ 1:length(nl) ], col = col[ 1:length(nl) ]  )
        }
        if ( safe_pdf & i == 1 )      dev.off( )
    }
}


plot_order_dysfunction  <-  function( rdr_dysf , pos = c(0,100), 
                                      logscale = 'y', cex = 1 ){
    tbl_rdr_dysf  =  aggregate( N_cells ~ order, data = rdr_dysf, FUN = sum )
    tbl_rdr_dysf  =  tbl_rdr_dysf[ order( tbl_rdr_dysf$N_cells, decreasing = T ), ]
    
    if ( logscale == '' ) cfcnt = 1.05 else cfcnt = 10.5
    
    plot_2D_lines( x = 1:nrow(tbl_rdr_dysf), DF = tbl_rdr_dysf, nl = 2 ,
                   names = c( 'Index', 'Number of cells' ),
                   yr = c(1, max( tbl_rdr_dysf$N_cells ) * cfcnt ),
                   xr = c(0.1, nrow( tbl_rdr_dysf )+2 ), 
                   type = 's', logscale = logscale ) 
    txt = NULL
    for( i in 1:nrow( tbl_rdr_dysf) ) {
        txt  = paste( txt, paste( i, tbl_rdr_dysf$order[ i ] ) , '\n', collapse = '   ') 
    }
    text( x = pos[1], y = pos[2], 
          labels = txt , pos = 4, cex = cex )
}




# Plot clones evolution ---------------------------------------------------

# Make a large number of colors
gen_colors  <-  function(nm = 12){
    # nm is a number of colors    
    w <- (nm^(1/3)) %/% 1 +1
    
    st <-  w^3 %/% nm
    
    sq <- seq(0,1-1/w,1/w)
    
    cr <- 1:nm
    
    l <- 0
    R <- 1:(w^3)
    G <- R
    B <- R
    
    for (i in 1:w) {
        for (j in 1:w) {
            for (k in 1:w) {
                l <- l+1
                R[l] <- sq[i]
                G[l] <- sq[j]
                B[l] <- sq[k]
            } 
        }  
    }
    
    # seq(1,w^3,st) # is consequence of each color to make a high diversity of colors
    jColor <- data.frame( number = 1:length( seq( 1,w^3, st ) ),
                          color  = rgb( R[seq( 1, w^3, st ) ], G[seq( 1, w^3, st)], 
                                        B[seq( 1, w^3, st ) ] ), stringsAsFactors = FALSE )
    
    return(jColor)
    
}

plot_clone_evolution  <-  function( threshold = 0.05, lwd = 2.0 ){

    clones_flow  =  data_flow[ ,c('Time', 'ID', 'ParentID', 'Birth_time', 'N_cells' ) ]
    Nmax  =  max( clones_flow$N_cells )
    time_max  =  max( clones_flow$Time )
    Nthreshold  =  round( Nmax * threshold )
    
    # delete clones with number of cells less than Nthreshold 
    w = sapply( X = ( max(clones_flow$ID[ which( clones_flow$Time == 0 ) ]) + 1 ):( max(clones_flow$ID ) ), 
                FUN = function( x ) { 
                    if ( max( clones_flow$N_cells[ which( clones_flow$ID == x ) ] ) > Nthreshold  ) 
                        return( x ) 
                    else 
                        return(NULL) 
                    }
                )
    w  =  unlist( w )
    w  =  c( clones_flow$ID[ which( clones_flow$Time == 0 ) ], w )
    
    DF  =  list( )
    for ( i in 1:length( w ) ){
        ss  =  which( clones_flow$ID == w[ i ] ) 
        DF[[ i ]]  =  data.frame( x = clones_flow$Time[ ss ], y = clones_flow$N_cells[ ss ] )
    }
    
    plot_2D_lines( x = DF[[ 1 ]]$x, DF = DF[[ 1 ]], nl = 2, names = c( 'Time step', 'Number of cells'),
                    xr = c( 1, time_max+7 ), yr = c( 1, Nmax ), draw_key = FALSE )
    # clrs  =  gen_colors( nm = length( w ) )
    clrs  =  randomColor(count = length( w ), 
                     hue = c(" ", "random", "red", "orange", "yellow", 
                             "green", "blue", "purple", "pink", "monochrome")[1], 
                     luminosity = c(" ", "random", "light", "bright", "dark")[5] )
    if ( length( w ) > 1 ){
        for( i in 2:length( w ) ){
                    lines( x = DF[[ i ]]$x, y = DF[[ i ]]$y, lwd = lwd, col = clrs[ i ] ) # clrs[ i, 'color']  )
        }
    }
    
}





