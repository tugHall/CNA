

# Libraries ---------------------------------------------------------------

library( randomcoloR )
library(fishplot)


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
                        safe_pdf = FALSE, filename = './plot.pdf' ){
    rp = 1 
    if ( safe_pdf )    {
        pdf( filename, width = 8, height = 8 )
        rp = 2
    }
    for( i in 1:rp ){
        par( mgp = c(2.2, 0.5, 0), font.axis = 2, font.lab = 2 )
        ### Plot the first line:
        y = DF[, nl[1] ]
        plot( x, y, type = 'l', xlab = names[1], xlim = xr, ylim = yr,
              ylab = names[2], axes = FALSE, cex.lab = 1.5, col = col[1],
              lwd = lwd, lty = lt[1] )
        
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
        key = names( DF[ nl ] )
        if( legend_names[1] != '') key = legend_names 
        legend( x = 'bottom', legend = key,  
                horiz = TRUE, xpd = TRUE,  inset = c(0, 1.03), box.col = "white", 
                lty = lt[ 1:length(nl) ], col = col[ 1:length(nl) ]  )
        if ( safe_pdf & i == 1 )      dev.off( )
    }
}



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




### Plot with circles: 
if ( FALSE ){
    cl = 'purple'    
    text(-5 + dltx, 5 + dlty, "A", cex = cx, col = 'black')
    text( 5 + dltx, 5 + dlty, "B", cex = cx, col = 'black')
    text(-5 + dltx,-5 + dlty, "C", cex = cx, col = 'black')
    text( 5 + dltx,-5 + dlty, "D", cex = cx, col = 'black')
    
    
    plot_2D( DF, sAB, sAC, sBD, sCD, names = c( 'X', 'Y' ),
             xr = c(-10,10), yr = c(-10,10))
    draw.arc(x=dltx,y=dlty,radius = 7, angle1 = 50*pi/180, 
             angle2 = 130*pi/180, col = cl, lwd = 1.6 )
    draw.arc(x=dltx,y=dlty,radius = 7.0, angle1 = 146*pi/180, 
             angle2 = 216*pi/180, col = cl, lwd = 1.6 )
    draw.arc(x=dltx,y=dlty,radius = 7.0, angle1 = 235*pi/180, 
             angle2 = 312*pi/180, col = cl, lwd = 1.6 )
    draw.arc(x=dltx,y=dlty,radius = 7.0, angle1 = 322*pi/180, 
             angle2 = 390*pi/180, col = cl , lwd = 1.6)
    
    clt = 'purple'
    text( -1, 9, round(sAB,3), cex = cx/1.5, col = clt)
    text( -8.7, 1, round(sAC,3), cex = cx/1.5, col = clt)
    text( 0, -7.6, round(sCD,3), cex = cx/1.5, col = clt)
    text( 7.2, 2, '0.00', cex = cx/1.5, col = clt)

    for( ngl in c(50, 130, 146, 216, 235, 312, 322, 390 ) ){
        draw.radial.line(start = 6.8, end = 7.2, center = c( dltx, dlty ),
                         angle = ngl*pi/180, col = cl, lwd = 1.6)
    }

}




# Plot clones -------------------------------------------------------------

if ( FALSE ){
    
    #provide a list of timepoints to plot
    #You may need to add interpolated points to end up with the desired
    #visualization. Example here was actually sampled at days 0 and 150
    timepoints=c(0,30,75,150, 200, 230, 270, 300 )     
    
    ### An example of the fish plot:
    #provide a matrix with the fraction of each population
    #present at each timepoint
    frac.table = matrix(
        c(53.79, 25, 13, 00,
          02, 20, 02, 00,
          02, 10, 02, 00,
          40, 05, 45, 07,
          25, 00, 45, 25.9, 
          08, 00, 02, 40, 
          03, 00, 01, 60, 
          07, 00, 03, 90),
        ncol=length(timepoints))
    
    #provide a vector listing each clone's parent
    #(0 indicates no parent)
    parents = c(0,0,0,0)
    
    #create a fish object
    fish = createFishObject(frac.table,parents,timepoints=timepoints)
    
    #calculate the layout of the drawing
    fish = layoutClones(fish)
    
    #draw the plot, using the splining method (recommended)
    #and providing both timepoints to label and a plot title
    fishPlot(fish,shape = "spline",title = "Fishplot of the clone evolution",
             title.btm = "Title bottom",
             cex.title = 1.5, vlines = timepoints, 
             cex.vlab = 1.2, vlab = timepoints)
    
    ##-------------------------------------------
    ##panel A
    timepoints=c(0,30,200,423)
    parents = c(0,1,1,3)
    frac.table = matrix(
        c(100, 38, 24, 00,
          002, 00, 01, 00,
          002, 00, 01, 01,
          100, 00, 98, 30),
        ncol=length(timepoints))
    
    fish = createFishObject(frac.table,parents,timepoints=timepoints)
    fish = layoutClones(fish, separate.independent.clones=F)
    fishPlot(fish, shape="spline", vlines=c(0,423), vlab=c(0,423), title="Sample 150288", cex.title=0.9, cex.vlab=0.8)
    ##panel label
    text(-100,130,"A",xpd=T,font=2)
    
    
    ##panel C
    
    timepoints=c(0,34,69,187,334,505,530)
    parents = c(0,1,1,1,3,4,0)
    frac.table = matrix(
        c(90, 28,     2, 60, 0,     2, 8,
          1,   0,   0.1, 00, 0,     0, 1,
          3,   0,   2.5, 00, 0,     0, 1,
          1,   0,   0.9, 00, 0,     0, 10,
          3,   0,   0.9, 00, 0.1,   0, 20,
          80,  0,    76, 00, 60,    0, 15,
          0.1, 0, 0.005, 00, 0.001, 0, 0),
        ncol=7)
    
    fish = createFishObject(frac.table,parents,timepoints=timepoints)
    fish = layoutClones(fish, separate.independent.clones=T)
    fish = setCol(fish,c("#888888", "#EF0000", "#8FFF40", "#FF6000", "#50FFAF", "#FFCF00", "#0070FF"))
    vlines=c(0,34,69,187,334,505,530,650,750)
    
    fishPlot(fish, shape="spline", vlines=vlines, vlab=vlines, title.btm="AML31",cex.vlab=0.9)
    
    ##panel label
    text(-125,130,"C",xpd=T,font=2)
    

}
