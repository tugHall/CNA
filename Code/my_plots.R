

# Libraries ---------------------------------------------------------------

library( randomcoloR )
# library(fishplot)


library(stringr)
# library(ape)
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

#' Function to plot 2D figure of points y = y(x)
#'
#' @param x Input data for axes X
#' @param y Input data for axes Y
#' @param names Vector of two characters with names for X and Y axes
#' @param pch Parameter pch for plot function
#' @param col Colors of points
#' @param cex Parameter cex for plot function
#' @param xr Range for X
#' @param yr Range for Y
#' @param safe_pdf Indicator to save plot to a file or not
#' @param filename Name of file to save plot if safe_pdf == TRUE
#'
#' @return NULL, making 2D plot using points
#' @export
#'
#' @examples plot_2D( x=-5:5, y=-3:7 )
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


#' Function to plot 2D figure of lines  yi = DF[ , nl[i] ), i - index
#'
#' @param x Input data for axes X
#' @param DF data.frame with data to plot
#' @param nl indexes of columns in DF to plot
#' @param names Vector of two characters with names for X and Y axes
#' @param legend_names Name of legend
#' @param col Vector of colors for lines
#' @param cex Parameter cex for plot function
#' @param lwd Vector of width of lines
#' @param lt Vector of types of lines
#' @param xr Range for X
#' @param yr Range for Y
#' @param safe_pdf Indicator to save plot to a file or not
#' @param filename Name of file to save plot if safe_pdf == TRUE
#' @param type Parameter type in plot function
#' @param logscale Parameter logscale in plot function
#' @param draw_key Indicator to draw key or not
#'
#' @return NULL, making 2D plot using lines
#' @export
#'
#' @examples plot_2D_lines( x = DF[, 3 ], DF, nl = 1:2 )
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


#' Function to plot order of genes dysfunction as a step function with number of cells related to each order
#'
#' @param rdr_dysf Order of genes dysfunction as a data.frame
#' @param pos Coordinates of list of order of genes dysfunction
#' @param logscale Parameter logscale for plot function
#' @param cex Parameter cex for plot function
#'
#' @return NULL, making plot with step function of order of genes' dysfunction
#' @export
#'
#' @examples plot_order_dysfunction( rdr_dysf )
plot_order_dysfunction  <-  function( rdr_dysf , pos = c(0,100),
                                      logscale = 'y', cex = 1 ){
    tbl_rdr_dysf  =  aggregate( N_cells ~ order, data = rdr_dysf, FUN = sum )
    tbl_rdr_dysf  =  tbl_rdr_dysf[ order( tbl_rdr_dysf$N_cells, decreasing = T ), ]

    if ( logscale == '' ) cfcnt = 1.05 else cfcnt = 10.5

    plot_2D_lines( x = 1:nrow(tbl_rdr_dysf), DF = tbl_rdr_dysf, nl = 2 ,
                   names = c( 'Index', 'Number of cells' ),
                   yr = c(1, max( tbl_rdr_dysf$N_cells ) * cfcnt ),
                   xr = c(0.1, round( nrow( tbl_rdr_dysf )+5, digits = -1) ),
                   type = 's', logscale = logscale,  draw_key  =  FALSE )
    txt = NULL
    for( i in 1:nrow( tbl_rdr_dysf) ) {
        txt  = paste( txt, paste( i, tbl_rdr_dysf$order[ i ] ) , '\n', collapse = '   ')
    }
    text( x = pos[1], y = pos[2],
          labels = txt , pos = 4, cex = cex )
}


# Plot clones evolution ---------------------------------------------------

#' Function to make a large number of colors
#'
#' @param nm Number of colors
#'
#' @return Vector of colors with length more than nm
#' @export
#'
#' @examples gen_colors( nm = 120 )
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


#' Function to plot clone evolution
#'
#' @param threshold Vector two numbers from 0 to 1 to show clones with relative final numbers of cells in the range of threshold
#' @param lwd Line width in the plot function
#' @param hue Parameter hue in the function randomColor from library randomcoloR
#' @param luminosity Parameter luminosity in the function randomColor from library randomcoloR
#' @param yr Range for Y axes
#' @param add_initial Indicator to add or do not add initial clones to plot
#' @param log_scale Indicator to use log_scale or not for Y axes
#'
#' @return NULL, making plot with clones evolution
#' @export
#'
#' @examples     plot_clone_evolution( threshold = c(0.01, 1 ), add_initial = TRUE, log_scale = FALSE )
#' @examples     plot_clone_evolution( threshold = c(0, 0.01 ), add_initial = FALSE, log_scale = TRUE )
plot_clone_evolution  <-  function( threshold = c(0.05,1.0), lwd = 2.0,
                                    hue = c(" ", "random", "red", "orange", "yellow",
                                            "green", "blue", "purple", "pink", "monochrome")[1],
                                    luminosity = c(" ", "random", "light", "bright", "dark")[5] ,
                                    yr = NA , add_initial = TRUE, log_scale = FALSE ){

    clones_flow  =  data_flow[ ,c('Time', 'ID', 'ParentID', 'Birth_time', 'N_cells' ) ]
    Nmax  =  max( clones_flow$N_cells )
    time_max  =  max( clones_flow$Time )
    Nthreshold  =  round( Nmax * threshold[1] )
    N_max  =  round( Nmax * threshold[2] )

    # delete clones with number of cells less than Nthreshold
    w = sapply( X = ( max(clones_flow$ID[ which( clones_flow$Time == 0 ) ]) + 1 ):( max(clones_flow$ID ) ),
                FUN = function( x ) {
                    if ( max( clones_flow$N_cells[ which( clones_flow$ID == x ) ] ) > Nthreshold &
                         max( clones_flow$N_cells[ which( clones_flow$ID == x ) ] ) < N_max )
                        return( x )
                    else
                        return(NULL)
                    }
                )
    w  =  unlist( w )
    if ( add_initial )  w  =  c( clones_flow$ID[ which( clones_flow$Time == 0 ) ], w )

    # clrs  =  gen_colors( nm = length( w ) )
    clrs  =  randomColor(count = length( w ),
                         hue = hue,
                         luminosity = luminosity )


    DF  =  list( )
    for ( i in 1:length( w ) ){
        ss  =  which( clones_flow$ID == w[ i ] )
        DF[[ i ]]  =  data.frame( x = clones_flow$Time[ ss ], y = clones_flow$N_cells[ ss ] )
    }
    if ( is.na(yr) ) yr = c( 1, N_max )
    plot_2D_lines( x = DF[[ 1 ]]$x, DF = DF[[ 1 ]], nl = 2, names = c( 'Time step', 'Number of cells'),
                    xr = c( 1, time_max+5 ), yr = yr, draw_key = FALSE,
                   logscale = ifelse( log_scale, 'y', '' ) )

    if ( length( w ) > 1 ){
        for( i in 2:length( w ) ){
                    lines( x = DF[[ i ]]$x, y = DF[[ i ]]$y, lwd = lwd, col = clrs[ i ] ) # clrs[ i, 'color']  )
        }
    }
    title( main = paste0('Number of shown clones is ', length( w ), ' from ',  max( clones_flow$ID ), ' clones' ) )
}





