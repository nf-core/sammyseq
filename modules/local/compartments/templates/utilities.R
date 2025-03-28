# Note: parameter names with "_" are from calder

## Note: possible riscale [0,1] as for hic
read.marks<-function(dir.file,name.file,bin.size,select,bin.select){

	cat("\n",name.file,"(",select,")",": read.. ")


	TR.chr<-get(load(paste0(dir.file,name.file)))
	TR.chr.sel<-data.table::data.table(do.call("cbind",(lapply(TR.chr[[select]],function(x) x$score))))


	## Ids compatible with calder (starts are 0-based)
	# all TR object have the same length
	Id<-( (start(TR.chr[[1]][[1]])-1  ) / as.numeric(bin.size) ) +1

	## Common bins between HiC and SAMMY/Chip-seq
	cat("select.. ")
	ind.take<-match(bin.select,Id)
	TR.chr.sel<-TR.chr.sel[ind.take,]
	TR.chr.sel$Id<-bin.select

	cat("\n")
	return(TR.chr.sel)
}

compute.fastcor.calder<-function(A,cor.cor=TRUE,trans.atanh=TRUE,const=1+1E-7){

	cat("\t1.cor.. ")
	cA<-CALDER::fast_cor(A)

	if(cor.cor){
		cat("2.cor.. ")
		ccA<- CALDER::fast_cor(cA)
	} else { ccA<-cA }

	if(trans.atanh){
		cat("inv.hyper.tangent.. ")
		accA<- atanh( ccA / const)
	} else { accA<- ccA }

	cat("\n")
	return(accA)

}

# A: should be the ccaA_oe_compressed_log
get.blocks.calder<-function(A,bin.size,chr,window.sizes = 3){

	p.th <- ifelse(as.numeric(bin.size) < 40000, 0.05, 1)

	# change the indices to be sure of the consistences of boundary predictions
	info.index<-data.table(Id=rownames(A),index=paste0(1:nrow(A)))
	rownames(A)<-colnames(A)<-1:nrow(A)

	#chr_name = paste0("chr", chr)
	TD.out<- CALDER::generate_compartments_bed(input_mat=A,chr=gsub("chr","",chr),p_thresh=p.th, bin_size=as.numeric(bin.size), window.sizes = window.sizes, out_file_name=NULL, stat_window_size = NULL )



	blocks<-lapply(1:nrow(TD.out$domain),function(i,D,info) {
		x<-D[i,]
		cl.id<-info$Id[x$from.id:x$to.id]

		return(cl.id)
	},D=TD.out$domain,info=info.index)
	names(blocks)<-1:length(blocks)

	return(blocks)

}

# A: should be the contact/distance matrix
cor.trend.blocks<-function(A,blocks,lag=4,trans.atanh=TRUE,scale=TRUE,const=1+1E-7,metric="mean"){

	# Summarize by blocks (it takes the rownames of the blocks)
	cat("\tSummarize by blocks.. ")
	B <- CALDER::HighResolution2Low_k_rectangle(mat=A,row_split=blocks,col_split=blocks,sum_or_mean =metric)
	rownames(B)<-colnames(B)<-names(blocks)

	# Compute enrichment trends at different lags
	cat("Trend.. ")
	n.block<-length(blocks)
	T.lags <- lapply( 1:lag, function(v,mat,n) {
		n<-nrow(B)
		1 * (mat[, -(1:v)] > mat[, - n - 1 + (v:1)])
	},mat=B,n=n.block)
	T <- do.call(cbind, T.lags)


	cat("Corr.. ")
	cT<- CALDER::fast_cor(t(T))

	if(trans.atanh){
		cat("Trans atanh.. ")
		acT<-atanh(cT/const)
	} else { acT<-cT }

	if(scale){
		cat("Scale..\n")
		acT.scaled<-scale(acT)
	} else {  acT.scaled<-acT }
	cat("\n")

	return(acT.scaled)

}

#distinguish_uniques è necessario

distinguish_uniques <- function( mat, sum = 1e-15 ){

    umat <- unique( mat )

    mat_length <- dim( mat )[ 1 ]
    umat_length <- dim( umat )[ 1 ]

    if( mat_length == umat_length ){

        return( mat )

    } else{

        not_unique_indices <- which( ! rownames( mat ) %in% rownames( umat ) )

        mat[ not_unique_indices, ] <- mat[ not_unique_indices, ] + sum

        mat <- distinguish_uniques( mat )

    }

    return( mat )

}

distinguish_uniques__in_vect  <- function( vect ){

    mat <- as.matrix( vect )
    rownames( mat ) <- seq( length( vect ) )

    umat <- distinguish_uniques( mat )

    uvect <- as.vector( umat )

    return( uvect )

}

adjust_hs <- function( l_r_h ){

    hs = sapply( l_r_h, function( v ) v$h )
    all_names = sapply( l_r_h, function( v ) paste0( collapse = "_",
        sort( c( v$l, v$r ) ) ) )
    r_names = sapply( l_r_h, function( v ) paste0( collapse = "_",
        sort( c( v$r ) ) ) )
    sizes = sapply( l_r_h, function( v ) length( v$l ) + length( v$r ) )

    hs = hs + sizes * 1e-07

    hs <- distinguish_uniques__in_vect( hs )

    l_b = 2
    r_b = which( r_names[ 1 ] == all_names )
    l_h = hs[ l_b ]
    r_h = hs[ r_b ]
    max_h = max( l_h, r_h )
    hs_new = mean( sort( hs, decreasing = TRUE )[ 2:3 ] )
    hs[ l_b ] = ifelse( l_h > r_h, max_h, hs_new )
    hs[ r_b ] = ifelse( r_h > l_h, max_h, hs_new )
    if( any( duplicated( hs ) ) )
        stop( "ERROR: DUPLICATED HEIGHTS exist in bisecting_kmeans" )

    return( hs )

}

bisecting_kmeans <- function( data ){

    dist_mat = as.matrix( stats::dist( data ) )
    indices = 1:nrow( data )
    l_r_h <<- list()

    get_h <- function( l_indices, r_indices ){

        combined_indices = c( l_indices, r_indices )
        idx <- as.matrix( expand.grid( combined_indices, combined_indices ) )
        max( dist_mat[ idx ] )

    }

    get_sub_tree <- function( indices ){

        n_nodes = length( indices )

        if( n_nodes == 1 ){

            h = NULL
            return()

        }

        if( n_nodes == 2 ){

            cluster = c( 1, 2 )

        } else{

            cluster = my_kmeans( x = data[ indices, ], centers =  2 )$cluster

        }

        l_indices = indices[ cluster == 1 ]
        r_indices = indices[ cluster == 2 ]
        h = get_h( l_indices, r_indices )
        l_r_h <<- c(
            l_r_h,
            list( list( l = l_indices, r = r_indices, h = h  ) )
        )

        l_branch = get_sub_tree( l_indices )
        r_branch = get_sub_tree( r_indices )

    }

    get_sub_tree( indices )
    hs = adjust_hs( l_r_h )
    r_hs = rank( hs )

    for( i in 1:length( l_r_h ) ){

        name = r_hs[ i ]
        names( name ) = paste0(
            collapse = "_",
            sort( c( l_r_h[[ i ]]$l, l_r_h[[i]]$r ) ) )
        l_r_h[[ i ]]$name = name

    }

    pos_names = sapply( l_r_h, function( v ) v$name )
    neg_names = -( 1:length( indices ) )
    names( neg_names ) = 1:length( indices )
    all_names = c( pos_names, neg_names )

    for(i in 1:length( l_r_h ) ){

        l_r_h[[ i ]]$l_name = unname(
            all_names[paste0( l_r_h[[ i ]]$l,
            collapse = "_" ) ] )
        l_r_h[[ i ]]$r_name = unname(
            all_names[ paste0( l_r_h[[ i ]]$r, collapse = "_" )])

    }

    merge_height = data.frame( l = sapply( l_r_h, function( v ) v$l_name ),
                              r = sapply( l_r_h, function( v ) v$r_name ), h = hs )

    merge_height = merge_height[ order( merge_height$h ), ]
    rownames( merge_height ) = NULL
    data_tmp = cbind( c( 0, 0, 1, 1 ), c( 0, 1, 1, 0 ) )
    hc = hclust( stats::dist( data_tmp ), "com" )
    hc$merge = as.matrix( unname( merge_height[ , 1:2 ] ) )
    hc$height = merge_height$h
    hc$labels = 1:length( indices )
    den <- as.dendrogram( hc )
    hc_r <- as.hclust( reorder( den, 1:length( indices ) ) )
    hc_r$method = "complete"
    hc_r$dist.method = "euclidean"
    l_r_h <<- list()
    rm( l_r_h )

    return( hc_r )

}

get.subcompartment.calder <- function( T, blocks, chr, genes_gr, bins_gr, n.comp = 10, const.comp = 5 ){


    ## Take first 10 principal component
    PC.comp <- CALDER::get_PCs( T, which = 1:n.comp )
    PC.comp[ , 2:n.comp ] <- PC.comp[ , 2:n.comp ] / const.comp

    ## First PCA should have the positive values in gene dense regions
    PC.comp[ , 1 ] <- set_sign_from_genedens( pc1 = PC.comp[ , 1 ], genes_gr = genes_gr, bins_gr = bins_gr, blocks = blocks )

    ## Distinguish identical lines, if there are, adding a not significant number
    PC.comp <- distinguish_uniques( PC.comp )

    ## Complete k(=2)-iterative clustering, with eucledean distance return a hclust dendogram object)
    H.k2 <- bisecting_kmeans( PC.comp )

    ## Reorder blocks using the first projected linear component
    ## Non-linear projection using the first two components
    new.pc1 <- CALDER::project_to_major_axis(PC.comp)
    ord.block<-CALDER::get_best_reorder(hc_hybrid_x_pro=H.k2, x_pro=new.pc1$x_pro)

    H.k2.ord <- dendextend::rotate(x=H.k2, order=ord.block)

    ## vector of
    AB.sub<-CALDER::get_cluser_levels(H.k2.ord, k_clusters=Inf, balanced_4_clusters=FALSE)$cluster_labels

    AB.sub.dt<-data.table::data.table(
                               chr=chr,
                               bin=unlist(blocks),
                               block=rep(names(blocks),sapply(blocks, length)),
                               sub=rep(AB.sub[names(blocks)],sapply(blocks, length)) )

    AB.sub.dt$sub.2<-substr(AB.sub.dt$sub,start=1,stop=1)
    AB.sub.dt$sub.4<-substr(AB.sub.dt$sub,start=1,stop=3)
    AB.sub.dt$sub.8<-substr(AB.sub.dt$sub,start=1,stop=5)

    info.pca<-data.table::data.table(
                              block=names(blocks),
                              sub=AB.sub[names(blocks)],
                              sub.2=substr(AB.sub[names(blocks)],start=1,stop=1),
                              sub.4=substr(AB.sub[names(blocks)],start=1,stop=3),
                              sub.8=substr(AB.sub[names(blocks)],start=1,stop=5),
                              pc1=PC.comp[,1],pc2=PC.comp[,2],
                              new.pc1=new.pc1$x_pro
                          )

    return(list(Bin=AB.sub.dt,Block=info.pca,Dendro=H.k2.ord))

}
