##############################################################
# ELISA_FUNCTIONS.R
#####################################################

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

##############################################################
# SAMMY_SUBCOMPARTMENTS.R
####################################################

## FUNCTIONS
### It takes tracks in input and calculate the euclidean distance matrix
my_read.SAMMY.calder <- function( tracks, track_names, bins_gr,  keeping_bins = "all", metric = "euclidean", cores = 4 ){

    track_matrix_info <- make_tracks_matrix(
        tracks = tracks,
        track_names = track_names,
        bins_gr = bins_gr,
        keeping_bins = keeping_bins,
        cores = cores
    )

    bws_dtable <- track_matrix_info[[ "bws_dtable" ]]
    keeping_bins <- track_matrix_info[[ "keeping_bins" ]]

    ## Annotate the bins with 0 coverage in all fractions
    ## They will be removed in all other samples
    bws_df <- as.data.frame( bws_dtable )
    rownames( bws_df ) <- as.character( keeping_bins )

    removing_bins1 <- rownames( bws_df[ ( rowSums( bws_dtable ) == 0 ), ] )
    print( "Bins with no coverage annotated" )

    ## Calculate eucledean distance between pairs of points (i.e., bins)
    ## Each point is define in the n-dimensional space, where n is 3,4, or 6 based on the number of fractions or Chip-seq experiments
    dist_mat <- as.matrix( dist( bws_dtable, method = metric ) )
    rownames( dist_mat ) <- colnames( dist_mat ) <- keeping_bins

    print( "Distance matrix made" )

    return( list( dist_mat = dist_mat, removing_bins1 = removing_bins1 ) )

}

### Function to extract gr object from subcompartment obj for selected subcompartment call
subgr_extractor <- function( subcompartment_bin, bins_gr, keeping_bins1, sub_colors, sublevel = "sub.8" ){

    subcomps_vect <- as.data.frame( subcompartment_bin )[ , sublevel ]
    subcolor_vect <- sub_colors[ subcomps_vect ]
    subcomp_df <- as.data.frame( cbind( subcomps_vect, subcolor_vect ) )

    ## bins_df <- as.data.frame( bins_gr[ keeping_bins1 ] )
    bins_df <- as.data.frame( bins_gr )

    subcomp_gr <- makeGRangesFromDataFrame( as.data.frame( cbind( bins_df, subcomp_df ) ), keep.extra.columns = TRUE )

    return( subcomp_gr )

}

### Function to make annotrack obj from subcompartment obj for selected subcompartment call
subanno_maker <- function( subcomp_gr, patient ){

    subcomp_anno <- AnnotationTrack(
        subcomp_gr,
        name = patient,
        stacking = "dense",
        showFeatureId = FALSE,
        id = subcomp_gr$subcomps_vect,
        fill= subcomp_gr$subcolor_vect,
        col = "transparent"
    )

    return( subcomp_anno )

}

removing_sammynocov_bins <- function( keeping_bins1, sammy_dist_objs, patients ){

    all_removing_bins1 <- c()
        for( patient in patients ){

            all_removing_bins1 <- c(
                all_removing_bins1,
                sammy_dist_objs[[ patient ]]
            )

        }
        all_removing_bins1 <- as.numeric( unique( all_removing_bins1 ) )

        keeping_bins1 <- keeping_bins1[ !( keeping_bins1 %in% all_removing_bins1 ) ]
        bins_gr <- bins_gr[ keeping_bins1 ]

        return( bins_gr )

}

### Wrapper to call subcompartments and return objects containing all the informations
call_subcompartments_sammy <- function( patients, tracks_db, bins_gr, subs_file, binsize, chr, genes_gr, keeping_bins1 = "all", sublevel = "sub.8", sub_colors = c( "B.2.2" = "#4575b4", "B.2.1" = "#74add1", "B.1.2" = "#abd9e9", "B.1.1" = "#e0f3f8", "A.1.1" = "#fee090", "A.1.2" = "#fdae61", "A.2.1" = "#f46d43", "A.2.2" = "#d73027" ), cores = 4, n.comp = 10, const.comp = 5 ){


    ## If a list of bins to analyzed has not been passed, use all genes in bins_gr
    if( keeping_bins1[ 1 ] == "all" ){

        keeping_bins1 <- seq( 1, length( bins_gr ) )
        print( "No bin removed" )

    }

    ## Proceed with compartment calculation
    print( "Calling subcompartments" )

    sammy_dist_objs <- mclapply( patients, mc.cores = cores, function( patient ){

        print( paste0( "Analysing: ", patient ) )

        ### Load files from a database containing for each row patient_name, fraction and file path
        sammy_files <- tracks_db[
            which( tracks_db$Patient_name == patient ),
            "File" ]
        names( sammy_files ) <- tracks_db[
            which( tracks_db$Patient_name == patient ),
            "Fraction" ]

        print( "Got file info" )

        ### Make the distance matrix
        sammy_distobj_file <- paste0( patient, "_distance-matrix___", chr, '_', binsize, ".Rdata" )

        if( !file.exists( sammy_distobj_file ) ){

            sammy_dist_obj <- my_read.SAMMY.calder(
                tracks = sammy_files,
                track_names = names( sammy_files ),
                bins_gr = bins_gr,
                keeping_bins = keeping_bins1,
                cores = cores
            )

            print( "Saving matrix..." )
            save( sammy_dist_obj, file = sammy_distobj_file )
            print( "Matrix saved" )

        } else{

            print( "Distance matrix already exists" )
            load( sammy_distobj_file )

        }

        return( sammy_dist_obj[[ "removing_bins1" ]] )

    })

    ## Make a list of bins with no coverage in at list one sample
    bins_gr <- removing_sammynocov_bins(
        keeping_bins1,
        sammy_dist_objs,
        patients
    )

    ## Calculate sub compartments
    sub_objs <- mclapply( patients, mc.cores = cores, function( patient ){

        print( paste( "Analysing patient", patient ) )

        ### Load the previously created distance matrix
        sammy_distobj_file <- paste0( patient, "_distance-matrix___", chr, '_', binsize, ".Rdata" )
        load( sammy_distobj_file )

        ### Remove from matrix bins with no coverage in all fractions in at least one sample
        sammy_dist_fullmat <- sammy_dist_obj[[ "dist_mat" ]]
        sammy_dist_mat <- sammy_dist_fullmat[ keeping_bins1, keeping_bins1 ]
        print( "Removed from the analysis bin with no coverage in all fraction in at least one patient" )

        rm( sammy_dist_obj )
        rm( sammy_dist_fullmat )

        ### Make the correlation matrix
        sammy_corrmat_file <- paste0( patient, "_corr-matrix___", chr, '_', binsize, ".Rdata" )

        if( !file.exists( sammy_corrmat_file ) ){

            sammy_corr_mat <- compute.fastcor.calder(
                A = sammy_dist_mat,
                cor.cor = TRUE,
                trans.atanh = TRUE
            )

            print( "Correlation matrix made" )

            print( "Saving matrix..." )
            save( sammy_corr_mat, file = sammy_corrmat_file )
            print( "Correlation matrix saved" )

        } else{

            print( "Correlation matrix already exists" )
            load( sammy_corrmat_file )

        }

        #### Get the Calder blocks
        sammy_blocks <- get.blocks.calder(
            A = sammy_corr_mat,
            bin.size = binsize,
            chr = chr
        )

        print( "Blocks calculated" )

        sammy_blockstrend_file <- paste0( patient, "_blocks-trend___", chr, '_', binsize, ".Rdata" )

        if( !file.exists( sammy_blockstrend_file ) ){

            sammy_blocks_trend <- cor.trend.blocks(
                A = sammy_dist_mat,
                blocks = sammy_blocks,
                lag = 4,
                trans.atanh = TRUE,
                scale = TRUE,
                metric = "mean"
            )

            print( "Saving matrix..." )
            save( sammy_blocks_trend, file = sammy_blockstrend_file )
            print( "Correlation matrix saved" )


        } else{

            print( "Blocks trend file exists" )
            load( sammy_blockstrend_file )

        }

        print( "Trend per block calculated" )

        ### Call subcompartments
        subcompartment_obj <- get.subcompartment.calder(
            T = sammy_blocks_trend,
            blocks = sammy_blocks,
            chr = chr,
            genes_gr = genes_gr,
            bins_gr = bins_gr,
            n.comp = n.comp,
            const.comp = const.comp
        )


        print( "Calculated subcompartments" )

        ## Transfrom subcompartment in a GRanges object to plot it with Givz
        subcomp_gr <- subgr_extractor(
            subcompartment_bin = subcompartment_obj$Bin,
            bins_gr = bins_gr,
            keeping_bins1 = keeping_bins1,
            sub_colors = sub_colors,
            sublevel = sublevel
        )

        ## Subcompartment object
        subcomp_anno <- subanno_maker( subcomp_gr, patient )

            return(

                list(
                    sammy_blocks = sammy_blocks,
                    sammy_blocks_trend = sammy_blocks_trend,
                    subcompartment_obj = subcompartment_obj,
                    annotrack = subcomp_anno,
                    gr = subcomp_gr
                )

            )

        })
        names( sub_objs ) <- patients

    save( sub_objs, file = subs_file )

    return( sub_objs )

}
##############################################################
# UTILITIES.R
##############################################################

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


### Import a set of tracks and arrange in a matrix (dtable, columns are the tracks and rows are the genomic bins)
make_tracks_matrix <- function( tracks, track_names, bins_gr, keeping_bins = "all", cores = 4 ){

    genome <- as.character( genome( bins_gr ) )

    ## Import tracks
    bws <- import_and_rebin__bw(
        files = tracks,
        bin_list = bins_gr,
        genome = genome,
        names = track_names,
        cores = cores
    )
    print( "Tracks imported" )

    ## Keep bins not having NA row in Hi-C
    after_bins_selected <- bins_selector( bws, track_names, keeping_bins, cores )

    bws <- after_bins_selected[[ "bws" ]]
    keeping_bins <- after_bins_selected[[ "keeping_bins" ]]

    ## Make the matrix (Elisa's code adapted)
    ### Extract scores and arrange them in a data.table obj
    score_list <- mclapply( track_names, mc.cores = cores, function( name ){

        bw <- bws[[ name ]]
        return( score( bw ) )

    })

    names( score_list ) <- track_names

    scores_df <- do.call( "cbind", score_list )
    bws_dtable <- data.table( scores_df )

    ### Transform the keeping_bin from 1-based to 0-based (Calder works on 0-based bin list)
    keeping_bins0 <- keeping_bins - 1

    print( "Files ready for calculating the euclidean distance" )

    return( list( bws_dtable = bws_dtable, keeping_bins = keeping_bins0 ) )

}


import_and_rebin__bw <- function( files, bin_list, genome, names, cores = 10 ){

    bws <- mclapply( files, mc.cores = cores, function( file ){

        ## Import the bigwig as RleLists
        bwR <- import.bw( file, as = "RleList" )

        ## Sort the "bin_list" seqlevels names to make them coincide with the bigwig imported
        bin_list__names <- seqlevels( bin_list )
        bwR <- bwR[ bin_list__names ]

        ## Rebin the imported bigwig according to the previous calculated bins
        bw <- binnedAverage( bins = bin_list, numvar = bwR, varname = "score" )

        ## Add extra information to the rebinned bigwig
        genome( bw ) <- genome

        return( bw )

    })
    names( bws ) <- names

    return( bws )

}

### Select the bin to discard from a set of tracks imported in a GRange object list
bins_selector <- function( bws, track_names, keeping_bins, cores = 4 ){

    if( keeping_bins[ 1 ] == "all" ){

        keeping_bins <- seq( 1, length( bins_gr ) )
        print( "No bin removed" )

    } else if( is.numeric( keeping_bins ) ){

        tmp_bws <- mclapply( track_names, mc.cores = cores, function( name ){

            bw_df <- as.data.frame( bws[[ name ]] )
            smallbw_df <- bw_df[ keeping_bins, ]

            return( makeGRangesFromDataFrame( smallbw_df, keep.extra.columns = TRUE ) )

        })

        names( tmp_bws ) <- track_names

        bws <- tmp_bws

        print( "Keeping bin list updated" )

    } else{

        print( "Wrong keeping bins list, is not numeric" )

    }

    return( list( bws = bws, keeping_bins = keeping_bins ) )

}


### Function to orient according to standard the first principal component defining the compartments
set_sign_from_genedens <- function( pc1, genes_gr, bins_gr, blocks ){

    ## Get bin coordinates
    ### Make a database to associate bins used for the analysis to compartment block
    blocks_dblist <- lapply( names( blocks ), function( nblock ){

        block_db <- cbind(
            as.numeric( nblock ),
            as.numeric( unlist( blocks[ nblock ] ) )
        )

        return( block_db )

    })
    blocks_db <- as.data.frame( do.call( rbind, blocks_dblist ) )
    names( blocks_db ) <- c( "block", "nbin" )

    ### Add to the database the information for each bin of pca value and if it is A or B according to the sign automatically calculated
    blocks_db$pc1 <- as.numeric( pc1[ blocks_db$block ] )
    blocks_db$fakecomp <- ifelse( blocks_db$pc1 > 0, 'A', 'B' )

    ### Merge the bin pca info to the bin genomic coordinates info
    bins_df <- as.data.frame( bins_gr )

    blocks_df <- cbind( bins_df, blocks_db )
    blocks_gr <- makeGRangesFromDataFrame( blocks_df, keep.extra.columns = TRUE )

    ## Calculate gene density for positive and negative bins
    ### Calculate the genes per bin
    blocks_df$ngenes <- countOverlaps( blocks_gr, genes_gr, type = "any", ignore.strand	= TRUE )

    ### Calculate the number of bins corrisponding to A and to B
    AB_nbins <- table( blocks_df$fakecomp )
    A_nbins <- as.numeric( AB_nbins[ 1 ] )
    B_nbins <- as.numeric( AB_nbins[ 2 ] )

    ### Count the genes in compartment A and B
    A_genes <- sum( blocks_df[ which( blocks_df$fakecomp == 'A' ), "ngenes" ] )
    B_genes <- sum( blocks_df[ which( blocks_df$fakecomp == 'B' ), "ngenes" ] )

    ### Calculate the gene density for A and B
    A_gendens <- A_genes / A_nbins
    B_gendens <- B_genes / B_nbins

    ## Decide to flip the sign or not
    if( A_gendens > B_gendens ){

        pc1_correct_sign <- pc1

    } else if( A_gendens < B_gendens ){

        pc1_correct_sign <- pc1 * -1

    } else{

        print( "Error! A and B compartments have exactly the same gene density!" )
        return( "Error" )

    }

    return( pc1_correct_sign )

}

######################
# BED AND BEDGRAPHS
######################


# Function to convert hex colors to RGB
rgb_str <- function(hex) {
  rgb_col <- paste(as.vector(col2rgb(hex)), collapse = ",")
  return(rgb_col)
}

# Function to generate TSV and BED files
generate_files <- function(sub_objs, chr) {
  for (ctrl in names(sub_objs)) {
    df_tp <- as.data.frame(sub_objs[[ctrl]][["gr"]])
    df_tp_chronly <- df_tp[df_tp$seqnames == chr,]

    # Generate the compartments TSV file
    write.table(df_tp_chronly,
                paste0(ctrl, "_", chr, "_compartments_merged.tsv"),
                sep = "\t",
                row.names = FALSE)

    # Generate the eigenvectors TSV file
    prvbin <- as.data.frame(sub_objs[[ctrl]][["subcompartment_obj"]][["Bin"]])
    prvblock <- as.data.frame(sub_objs[[ctrl]][["subcompartment_obj"]][["Block"]])

    df_eigenvect <- data.frame()
    for (i in seq_len(nrow(prvbin))) {
      partdf <- prvblock[prvblock$block == as.integer(prvbin$block[i]), c("block", "pc1")]
      df_eigenvect <- rbind(df_eigenvect, partdf)
    }

    df_tp_chronly_eigenvect <- cbind(df_tp_chronly, df_eigenvect)
    write.table(df_tp_chronly_eigenvect,
                paste0(ctrl, "_", chr, "_compartments_eigenvector.tsv"),
                sep = "\t",
                row.names = FALSE)

    # Generate the BED file
    df_tp_chronly$strand <- gsub("\\*", "\\.", df_tp_chronly$strand)
    df_tp_chronly$zero <- 0
    bed_data <- df_tp_chronly[, c('seqnames', 'start', 'end', 'subcomps_vect', 'zero', 'strand', 'start', 'end', 'subcolor_vect')]

    # Modify colors for A and B compartments
    bed_data$subcolor_vect <- ifelse(substr(bed_data$subcomps_vect, 1, 1) == "A", "69,117,180", "207,207,207")

    bed_file <- paste0(ctrl, "_", chr, "_compartments.bed")
    header_bedfile <- paste0('track name="', ctrl, '" description="', ctrl, ' (Emission ordered)" visibility=1 itemRgb="On"')

    writeLines(header_bedfile, bed_file)
    write.table(bed_data,
                bed_file,
                append = TRUE,
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)

    # Generate the bedGraph file for eigenvectors
    bedgraph_data <- df_tp_chronly_eigenvect[, c('seqnames', 'start', 'end', 'pc1')]
    bedgraph_file <- paste0(ctrl, "_", chr, "_comp_eigenvector.bedgraph")
    header_bedgraph <- paste0('track type=bedGraph name="', ctrl, '_eigenvector" description="', ctrl, ' eigenvector" visibility=full color=200,100,0 altColor=0,100,200 priority=20')

    writeLines(header_bedgraph, bedgraph_file)
    write.table(bedgraph_data,
                bedgraph_file,
                append = TRUE,
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
  }
}
