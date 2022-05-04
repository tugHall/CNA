


library(readtext)

get_functions_names  <-  function( file = './Code/Functions_read_maps.R' ){    
    txt <- readtext( file )[,2]
    
    # cat( txt )
    
    st <- str_locate_all( txt, '\\{')[[1]][ , 1 ]
    
    end <- str_locate_all( txt, '\\}')[[1]][ , 1 ]
    
    sr  <-  sort( c( st, end) )
    
    m_st  <-  match(st,  sr)
    m_end <-  match(end, sr)
    
    
    for( i in 1:length( st ) ){
        if ( i == 1 ){
            st_true <- st[1]
            next
        }
        if ( m_st[i-1] == m_st[i] - 1 ) next
        st_true = c( st_true, st[i] )
    }
    
    end_true  <- NULL
    for( i in 1:length( end ) ){
        if ( i == length( end )){
            end_true = c( end_true, end[i] )
            break
        }
        if ( m_end[i] == m_end[i+1] - 1 ) next
        end_true = c( end_true, end[i] )
    }
    
    for( j in 1:length( st_true ) ){
        
        if ( j == 1 ) {
            txt_fncts  <- substr( txt, 1, st_true[ 1 ] - 1 )
            next
        }
        txt_fncts  <- paste0( txt_fncts, substr( txt, end_true[ j-1 ] + 1 , st_true[ j ] - 1 ) )
    }
    txt_fncts  <- paste0( txt_fncts, substr( txt, end_true[ length( st_true ) ] + 1 , length( txt ) ) )
    
    
    # cat( txt_fncts )
    return( txt_fncts )
}    

txt <- get_functions_names( file = './Code/Functions_read_maps.R' )


txt <- get_functions_names( file = './Code/tugHall_3.0_functions copy.R' )


