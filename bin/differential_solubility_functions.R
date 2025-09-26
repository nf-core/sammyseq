#!/usr/bin/env Rscript
# differential_solubility_functions.R

#####################################################################
## IMPORT AND REBIN BIGWIG FUNCTION
#####################################################################

import_and_rebin__bw <- function(files, bin_list, names, cores = 1, genome = NULL) {
    bws <- parallel::mclapply(files, mc.cores = cores, function(file) {
    
    bwR <- rtracklayer::import(file, format = "BigWig", as = "RleList")

    bin_names <- GenomeInfoDb::seqlevels(bin_list)
    bw_names  <- names(bwR)
    chr_order <- sapply(paste0("^", bw_names, "$"), function(chr) grep(chr, bin_names))
    bins_for_bw <- bin_list
    GenomeInfoDb::seqlevels(bins_for_bw) <- bin_names[as.vector(unlist(chr_order))]

    bw <- GenomicRanges::binnedAverage(
    bins    = bins_for_bw,
    numvar  = bwR[ GenomeInfoDb::seqlevels(bins_for_bw) ],
    varname = "score"
    )

    # CHANGED -> set genome tag only if provided (we will provide the genome parameter later...)
    if (!is.null(genome)) GenomeInfoDb::genome(bw) <- genome

    bw
    })
    names(bws) <- names
    bws
}

#####################################################################
## CHECK SIGN FUNCTION
#####################################################################


check_sign<- function(x,meann){
            if ( sign(x) == sign(meann) ){
                return("constant_solubility")
        } else if (sign(x) > sign(meann)) {
        
            return("shift_increase")
        }else if (sign(x) < sign(meann)) {
        
            return("shift_decrease")
        }
}

#####################################################################
## CONFIDENCE INTERVAL FUNCTION
#####################################################################
##standard error

confidence_interval <- function(vector, nm="prove") {
    # Standard deviation of sample
    vec_serr <- sd(vector)/sqrt(length(vector))
    vec_serr2x<- vec_serr*2
    # Sample size
    n <- length(vector)
    # Mean of sample
    vec_mean <- mean(vector)
    nm_confint_low<-paste0(nm,"_serrx2_lower")
    nm_confint_up<-paste0(nm,"_serrx2_upper")
    name_mean <- paste0(nm,"_mean")
    name_serr <- paste0(nm,"_serrX2")
    result <- c(nm_confint_low = vec_mean - vec_serr2x,
                nm_confint_up = vec_mean + vec_serr2x, 
                name_mean = vec_mean, 
                name_serr = vec_serr2x
                )
    names(result) <- c(nm_confint_low,nm_confint_up,name_mean,name_serr)
    return(result)
}

#####################################################################
## RANGE CHECK FUNCTION
#####################################################################

is_in_serrx2_range_and_shift <-function(vector,
                                        xgroup_name,
                                        y_name
                                        ){
    #meann,lower,upper,mean_tp_name,whereis_tp_name){
    xgroup_lower_bound_confint_name<-paste0(xgroup_name,"_serrx2_lower")
    xgroup_upper_bound_confint_name<-paste0(xgroup_name,"_serrx2_upper")
    xgroup_mean_name<-paste0(xgroup_name,"_mean")

    xgroup_lower_bound_confint<-vector[xgroup_lower_bound_confint_name][1]
    xgroup_upper_bound_confint<-vector[xgroup_upper_bound_confint_name][1]
    
    yname_val<-vector[y_name][1]
    xgroup_mean<-vector[xgroup_mean_name][1]

    confint_check <-''
    mean_startsign<- ''
    mean_sign <-''
    shift_solubility<- ''
    
    col_is_in_confint_name<-paste0(y_name,"_ov_check")
    col_whereis_tp_name<- paste0(y_name,"_ov_specs")
    shift_solubility_name<- paste0(y_name, "_sol_shift")
    mean_startsign_name<- paste0(xgroup_name, "_mean_startsign")
    ##########check if both or one mean is in the confint of the other mean
    if ( between(yname_val, xgroup_lower_bound_confint, xgroup_upper_bound_confint) 
        ) {
        confint_check <- sign(yname_val)
        mean_sign <- "nodiff"
        shift_solubility <- "constant_solubility"
        mean_startsign<- sign(xgroup_mean)
        ##########check if confint are not overlapping, the ygroup timepoint confint is lower than xgroup
        }else if (yname_val < xgroup_lower_bound_confint ) {
        
        confint_check <- sign(yname_val)

        #mean_sign <- paste0(y_name ,"_lower_than_",xgroup_name)
        mean_sign <- "lower"
        shift_solubility <- check_sign(yname_val,xgroup_mean)
        mean_startsign<- sign(xgroup_mean)
        
        ##########check if confint are not overlapping, the xgroup timepoint confint is lower than ygroup
        }else if (xgroup_upper_bound_confint < yname_val) {
        
        confint_check <-sign(yname_val)
        mean_sign <- "higher"
        shift_solubility <- check_sign(yname_val,xgroup_mean)
        mean_startsign<- sign(xgroup_mean)
        }

        else{
    
        confint_check <- sign(yname_val)
        mean_sign <- "no_idea"
        shift_solubility <- check_sign(yname_val,xgroup_mean)
        mean_startsign<- sign(xgroup_mean)
    }
    
    result <- c(confint_check,mean_sign,shift_solubility,mean_startsign)
    
    names(result) <- c(col_is_in_confint_name,col_whereis_tp_name,shift_solubility_name,mean_startsign_name)
    return(result)
}   

#####################################################################
## MAIN BINS SELECTOR FUNCTION
#####################################################################
#(ex list_groups)
#list_groups<- lapply(1:ncol(pr),Bins_selector(x,allmixeddf_grobj=allmixeddf_S2svsS3_grobj,fr1="S2S", fr2="S3"))

#Bins_selector <- function(combination, allmixeddf_grobj, fraction1 = "S2S", fraction2 = "S3") { ##aggiunto ths alla funzione
Bins_selector <- function(combination, allmixeddf_grobj, fraction1 = "S2S", fraction2 = "S3", ths = 0.1) {
    
    cat("Running Bins_selector for combination:", combination, "\n")
    
    # Define groups to compare and select their samples
    x <- get(pr[, combination][1])
    sample_typex <- pr[, combination][1]
    sample_typey <- pr[, combination][2]
    y <- get(pr[, combination][2])
    
    cat("Groups:", paste(x, collapse = ", "), "vs", paste(y, collapse = ", "), "\n")
    
    # Define constraints for name uniqueness
    xgroup <- gsub("_.*", "", x[1], perl = TRUE)
    ygroup <- gsub("_.*", "", y[1], perl = TRUE)
    
    cat("Group names:", xgroup, "vs", ygroup, "\n")

    ##new_selection<-NULL
    new_selection <- allmixeddf_grobj ## Assigned
    mcols(new_selection) <- mcols(new_selection)[c(x, y)]

    # Calculations   
    # Confidence intervals and means
    confint_mean_sd_first <- apply(as.matrix(mcols(new_selection)[c(x)]), 1, confidence_interval, nm = xgroup) 
    confint_mean_sd_second <- apply(as.matrix(mcols(new_selection)[c(y)]), 1, confidence_interval, nm = ygroup)     
    delta <- confint_mean_sd_first[1, ] - confint_mean_sd_second[1, ]
    res <- cbind(t(confint_mean_sd_first), t(confint_mean_sd_second), delta) 
    df_toadd1 <- do.call("cbind", as.data.frame(res))
    mcols(new_selection) <- cbind(mcols(new_selection), df_toadd1)

    # Range analysis forward comparison
    range_analysis <- mclapply(1:length(y), mc.cores = 1, FUN = function(n) {        
        y_name <- y[n]
        z <- apply(as.matrix(mcols(new_selection)[c(paste0(xgroup, "_serrx2_lower"), paste0(xgroup, "_serrx2_upper"), paste0(xgroup, "_mean"), y_name)]), 1, 
            is_in_serrx2_range_and_shift, # invece di is_in_sdx2_range_and_shift,
            xgroup_name = xgroup,
            y_name = y_name 
        )
        return(as.data.frame(t(z)))
    })

    df_toadd <- do.call("cbind", range_analysis)
    mcols(new_selection) <- cbind(mcols(new_selection), df_toadd)
    
    # Range analysis reverse comparison   
    range_analysis_rev <- mclapply(1:length(x), mc.cores = 1, FUN = function(n) {
        x_name <- x[n]
        z <- apply(as.matrix(mcols(new_selection)[c(paste0(ygroup, "_serrx2_lower"), paste0(ygroup, "_serrx2_upper"), paste0(ygroup, "_mean"), x_name)]), 1, 
            is_in_serrx2_range_and_shift, # invece di is_in_sdx2_range_and_shift,
            xgroup_name = ygroup,
            y_name = x_name
        )
        return(as.data.frame(t(z)))
    })

    df_toadd2 <- do.call("cbind", range_analysis_rev)
    mcols(new_selection) <- cbind(mcols(new_selection), df_toadd2)
    
    # Select bins
    #   ths <- 0.1
    prvdf <- as.data.frame(new_selection)
    
    # Out of range check
    prvdftest <- prvdf[abs(prvdf[paste0(xgroup, "_mean")]) >= ths, ]
    pprvlow <- prvdftest[paste0(y, "_ov_specs")] == "lower"   # <-- cambiato da _shift a _ov_specs
    pprvhigh <- prvdftest[paste0(y, "_ov_specs")] == "higher" # <-- cambiato da _shift a _ov_specs
    
    # Sum up by group 1
    prvdftest[, paste0(ygroup, "_sign_SUM")] <- rowSums(sapply(prvdftest[, paste0(y, "_ov_check")], as.numeric))
    prvdftest$ovlow <- apply(pprvlow, 1, sum) * -1
    prvdftest$ovvhigh <- apply(pprvhigh, 1, sum)     

    # Commutative group testing
    pprvlow_X <- prvdftest[paste0(x, "_ov_specs")] == "lower"   #  <-- cambiato da _shift a _ov_specs
    pprvhigh_X <- prvdftest[paste0(x, "_ov_specs")] == "higher" #  <-- cambiato da _shift a _ov_specs

    # Sum up by group 2
    ## PRIMA ERA prvdftest[, paste0(xgroup, "_sign_SUM")] <- rowSums(sapply(prvdftest[, paste0(x, "_sign")], as.numeric))
    prvdftest[, paste0(xgroup, "_sign_SUM")] <- rowSums(sapply(prvdftest[, paste0(x, "_ov_check")], as.numeric))
    prvdftest$ovlow_X <- apply(pprvlow_X, 1, sum) * -1
    prvdftest$ovvhigh_X <- apply(pprvhigh_X, 1, sum)
    prvdftest_gr <- makeGRangesFromDataFrame(prvdftest, keep.extra.columns = TRUE)
    
    # INFORMATIVE over THS BINS SELECTION
    # Separate bins according to group X start sign
    a <- paste0(xgroup, "_mean")
    column <- which(names(prvdftest_gr@elementMetadata@listData) == a)
    startmeanpos <- prvdftest_gr[prvdftest_gr@elementMetadata[[column]] >= 0]  
    startmeanneg <- prvdftest_gr[prvdftest_gr@elementMetadata[[column]] < 0]  
    
    # Select coherent bins with value out of the IC range of the comparison group for all the values 
    ovvhighconservedpos <- startmeanpos[startmeanpos$ovvhigh == length(y) & startmeanpos$ovlow_X == -length(x)] 
    ovlowconservedpos <- startmeanpos[startmeanpos$ovlow == -length(y) & startmeanpos$ovvhigh_X == length(x)]
    ovvhighconservedneg <- startmeanneg[startmeanneg$ovvhigh == length(y) & startmeanneg$ovlow_X == -length(x)] 
    ovlowconservedneg <- startmeanneg[startmeanneg$ovlow == -length(y) & startmeanneg$ovvhigh_X == length(x)]
    
    # Make a list of bins to save and analyze 
    #(GENERALIZED for different comparisons in nextflow)

    list_ofbins_to_save_and_analyse <- setNames(
        list(
        ovvhighconservedpos,
        ovlowconservedpos,
        ovvhighconservedneg,
        ovlowconservedneg
        ),
        c(
        paste0(fraction1, "_up"),
        paste0(fraction1, "_down"),
        paste0(fraction2, "_up"),
        paste0(fraction2, "_down")
        )
    )

    
    # Numeric coding for the groups
    mcols(ovvhighconservedpos)[[paste0(ygroup, "_vs_", xgroup)]] <- rep(2, length(ovvhighconservedpos))
    mcols(ovlowconservedpos)[[paste0(ygroup, "_vs_", xgroup)]] <- rep(1, length(ovlowconservedpos))
    mcols(ovvhighconservedneg)[[paste0(ygroup, "_vs_", xgroup)]] <- rep(-1, length(ovvhighconservedneg))    
    mcols(ovlowconservedneg)[[paste0(ygroup, "_vs_", xgroup)]] <- rep(-2, length(ovlowconservedneg))
    
    # Prepare results list
    x <- list()
    x[[paste0(ygroup, "_vs_", xgroup, "_all_shifting_bins")]] <- c(
        ovvhighconservedpos,                                                   
        ovlowconservedpos,
        ovvhighconservedneg,
        ovlowconservedneg
    )    
    x[[names(list_ofbins_to_save_and_analyse[1])]] <- list_ofbins_to_save_and_analyse[[1]]
    x[[names(list_ofbins_to_save_and_analyse[2])]] <- list_ofbins_to_save_and_analyse[[2]]
    x[[names(list_ofbins_to_save_and_analyse[3])]] <- list_ofbins_to_save_and_analyse[[3]]
    x[[names(list_ofbins_to_save_and_analyse[4])]] <- list_ofbins_to_save_and_analyse[[4]]
    x[[paste0(ygroup, "_allgr_", xgroup)]] <- prvdftest_gr  

    names(x) <- c(
        paste0(ygroup, "_vs_", xgroup, "_all_shifting_bins"),
        paste0(ygroup, "_", names(list_ofbins_to_save_and_analyse[1]), "_", xgroup),
        paste0(ygroup, "_", names(list_ofbins_to_save_and_analyse[2]), "_", xgroup),
        paste0(ygroup, "_", names(list_ofbins_to_save_and_analyse[3]), "_", xgroup),
        paste0(ygroup, "_", names(list_ofbins_to_save_and_analyse[4]), "_", xgroup),
        paste0(ygroup, "_allgr_", xgroup)
    )
    
    cat("Results:", paste0(names(x[paste0(ygroup, "_vs_", xgroup, "_all_shifting_bins")]), "_", length(x[[paste0(ygroup, "_vs_", xgroup, "_all_shifting_bins")]]), "_bins"), "\n")
    return(x)
}