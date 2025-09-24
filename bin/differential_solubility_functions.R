#####analisi

##definisce trend del campione rispetto a media di gruppo di confronto
check_sign<- function(x,meann){
            if ( sign(x) == sign(meann) ){
                return("constant_solubility")
        } else if (sign(x) > sign(meann)) {
        
            return("shift_increase")
        }else if (sign(x) < sign(meann)) {
        
            return("shift_decrease")
        }
}

###calcola stdev,sd2x,intervallo di confidenza (interval non è usata) su sd2x, media,delta per ciascun gruppo
confidence_interval <- function(vector, nm="prove") {
    vec_sd <- sd(vector)
    vec_sd2x<- vec_sd*2
    vec_mean <- mean(vector)
    nm_confint_low<-paste0(nm,"_sdx2_lower")
    nm_confint_up<-paste0(nm,"_sdx2_upper")
    name_mean <- paste0(nm,"_mean")
    name_sd <- paste0(nm,"_sdX2")
    name_mean_sign <- paste0(nm,"_mean_sign")
        result <- c(
                vec_mean = vec_mean,
                vec_mean_sign = sign(vec_mean),
                vec_sd2x = vec_sd2x,
                vec_confint_low = vec_mean - vec_sd2x,
                vec_confint_up = vec_mean + vec_sd2x
                )
    names(result) <- c(name_mean,name_mean_sign,name_sd,nm_confint_low,nm_confint_up)
    return(result)
}



###La funzione valuta se le medie candono nell'int di conf dell'altro gruppo,definiscono se ci sia uno shift di solubilità e il verso e registra
il segno della media e del valore testato

is_in_sdx2_range_and_shift_eva<-function(vector,
                                xgroup_name,
                                y_name
                                ){
    xgroup_lower_bound_confint_name<-paste0(xgroup_name,"_sdx2_lower")
    xgroup_upper_bound_confint_name<-paste0(xgroup_name,"_sdx2_upper")
    xgroup_mean_name<-paste0(xgroup_name,"_mean")
    y_sign_name<-paste0(y_name,"_sign")
    y_shift_name<- paste0(y_name,"_shift")

    xgroup_lower_bound_confint<-vector[xgroup_lower_bound_confint_name][1]
    xgroup_upper_bound_confint<-vector[xgroup_upper_bound_confint_name][1]
    yname_val<-vector[y_name][1]
    y_sign <-''
    y_shift <-''
    ##########check if both or one mean is in the confint of the other mean
    if ( between(yname_val, xgroup_lower_bound_confint, xgroup_upper_bound_confint) 
        ) {
        y_sign <- sign(yname_val)
        y_shift <- "nodiff"
    ##########check if confint are not overlapping, the ygroup timepoint confint is lower than xgroup
        }else if (yname_val < xgroup_lower_bound_confint ) {  
        y_sign <- sign(yname_val)
        y_shift <- "lower"
    ##########check if confint are not overlapping, the xgroup timepoint confint is lower than ygroup
        }else if (xgroup_upper_bound_confint < yname_val) { 
        y_sign <-sign(yname_val)
        y_shift <- "higher"
        }
        else{
        y_sign <- sign(yname_val)
        y_shift <- "no_idea"
    }
    
    result <- c(y_sign,y_shift)   
    names(result) <- c(y_sign_name,y_shift_name)
    return(result)
}



###Ths selection clean
##add normalized ratio table    

#list_groups<- lapply(1:ncol(pr),Bins_selector(x,allmixeddf_grobj=allmixeddf_S2svsS3_grobj,fr1="S2S", fr2="S3"))


Bins_selector=function(combination,allmixeddf_grobj=allmixeddf_S2svsS3_grobj,fraction1=S2S, fraction2=S3) {
##define groups to compare and select their samples
##test combination<-1
    x<-get(pr[,combination][1])
    sample_typex<- pr[,combination][1]
    sample_typey<- pr[,combination][2]
    y<-get(pr[,combination][2])
    print(c(x,y))
##definire vincoli per unicità del nome e se repX va obbligatoriamente finale sep da . o _
    xgroup<-gsub("_.*", "", x[1],perl = TRUE)
    ygroup<-gsub("_.*", "", y[1],perl = TRUE)
    print(xgroup)
    print(ygroup)

    new_selection<-NULL
    new_selection <-allmixeddf_grobj  ###da assegnare
    mcols(new_selection) <- mcols(new_selection)[c(x,y)]

##Calculations   

#   1)calcolo int di confidenza, media e stdev per ogni gruppo
    confint_mean_sd_first<- apply(as.matrix(mcols(new_selection)[c(x)]),1, confidence_interval,nm=xgroup) 
    confint_mean_sd_second<- apply(as.matrix(mcols(new_selection)[c(y)]),1, confidence_interval,nm=ygroup)     
    delta<-confint_mean_sd_first [1,]-confint_mean_sd_second[1,]
    res<-cbind(t(confint_mean_sd_first),t(confint_mean_sd_second),delta) 
    df_toadd1<- do.call("cbind", as.data.frame(res))
    mcols(new_selection)<- cbind(mcols(new_selection),df_toadd1)

###########Range analysis fw 
    range_analysis<- mclapply(1: length(y) , mc.cores=1 , FUN=function(n) {        
        y_name<- y[n]
        z<-apply(as.matrix(mcols(new_selection)[c(paste0(xgroup,"_sdx2_lower"),paste0(xgroup,"_sdx2_upper"), paste0(xgroup,"_mean"), y_name)]), 1, 
            is_in_sdx2_range_and_shift,
            xgroup_name= xgroup,
            y_name= y_name
            )
            return(as.data.frame(t(z)))
            })

    df_toadd<- do.call("cbind", range_analysis)
    mcols(new_selection)<- cbind(mcols(new_selection),df_toadd)
    
    
###########Range analysis reverse comparison   

    range_analysis_rev<- mclapply(1: length(x) , mc.cores=1 , FUN=function(n) {
        
        x_name<- x[n]
        z<-apply(as.matrix(mcols(new_selection)[c(paste0(ygroup,"_sdx2_lower"),paste0(ygroup,"_sdx2_upper"), paste0(ygroup,"_mean"), x_name)]), 1, 
            is_in_sdx2_range_and_shift,
            xgroup_name= ygroup,
            y_name= x_name
            )
            return(as.data.frame(t(z)))
            })

    df_toadd2<- do.call("cbind", range_analysis_rev)
    mcols(new_selection)<- cbind(mcols(new_selection),df_toadd2)
    

####Select bins
#seleziona i bin le cui medie variano oltre ad un det ths e li divide in base a alterazione del segno (lower or higher),alterazione solubilità     
    ths<-0.3
    prvdf<- as.data.frame(new_selection)
    #### out of range check
    prvdftest<-prvdf[abs(prvdf[paste0("delta","")])>=ths,]
    pprvlow<- prvdftest[paste0(y,"_shift")]=="lower"
    pprvhigh<- prvdftest[paste0(y,"_shift")]=="higher"
    ####sum up by group 1
    prvdftest[,paste0(ygroup,"_sign_SUM")] <- rowSums(sapply(prvdftest[,paste0(y,"_sign")],as.numeric))
    prvdftest$ovlow<-apply(pprvlow,1,sum)*-1
    prvdftest$ovvhigh<-apply(pprvhigh,1,sum)     
    
    ##################################################################################################
    ### commutative group testing
    #### out of range check
    pprvlow_X<- prvdftest[paste0(x,"_shift")]=="lower"
    pprvhigh_X<- prvdftest[paste0(x,"_shift")]=="higher"

    ####sum up by group 2
    prvdftest[,paste0(xgroup,"_sign_SUM")] <- rowSums(sapply(prvdftest[,paste0(x,"_sign")],as.numeric))
    ###count characterisic
    prvdftest$ovlow_X<-apply(pprvlow_X,1,sum)*-1
    prvdftest$ovvhigh_X<-apply(pprvhigh_X,1,sum)
    prvdftest_gr<-makeGRangesFromDataFrame(prvdftest,keep.extra.columns = T)
        
    

    ###INFORMATIVE over THS BINS SELECTION
    #separate bins according to group X start sign
    a<- paste0(xgroup,"_mean")
    column<- which(names(prvdftest_gr@elementMetadata@listData)==a)
    startmeanpos<-prvdftest_gr[prvdftest_gr@elementMetadata[[column]] >= 0]  
    startmeanneg<-prvdftest_gr[prvdftest_gr@elementMetadata[[column]] < 0 ]  
    ###   select coherent bins with value out of the IC range of the comparison group for all the values 
    ovvhighconservedpos <- startmeanpos[startmeanpos$ovvhigh == length(y) & startmeanpos$ovlow_X == - length(x)] 
    ovlowconservedpos <-  startmeanpos[startmeanpos$ovlow == -length(y)  & startmeanpos$ovvhigh_X == length(x)]
    ovvhighconservedneg <- startmeanneg[startmeanneg$ovvhigh == length(y)& startmeanneg$ovlow_X == - length(x)] 
    ovlowconservedneg <-  startmeanneg[startmeanneg$ovlow == -length(y)& startmeanneg$ovvhigh_X == length(x)] 
    ####make a list of which I need to do the same and make a for loop
    ####names need generalizazion
    list_ofbins_to_save_and_analyse<- list("S2S_up"=ovvhighconservedpos,
                                            "S2S_down"=ovlowconservedpos,
                                            "S3_up"=ovvhighconservedneg,
                                            "S3_down"=ovlowconservedneg )

    
    print(names(list_ofbins_to_save_and_analyse[1]))
    ###numeric coding for the groups

    ###+2
    mcols(ovvhighconservedpos )[[paste0(ygroup,"_vs_",xgroup)]]<- rep(2,length(ovvhighconservedpos))
    ####+1
    mcols(ovlowconservedpos)[[paste0(ygroup,"_vs_",xgroup)]] <- rep(1,length(ovlowconservedpos))
    #### -1
    mcols(ovvhighconservedneg)[[paste0(ygroup,"_vs_",xgroup)]] <- rep(-1,length(ovvhighconservedneg))    
    ###-2
    mcols(ovlowconservedneg)[[paste0(ygroup,"_vs_",xgroup)]] <- rep(-2,length(ovlowconservedneg))

    
    

    x<- list()
    x[[paste0(ygroup,"_vs_",xgroup,"_all_shifting_bins")]] <-c(ovvhighconservedpos,                                                   
                                                                ovlowconservedpos,
                                                                ovvhighconservedneg,
                                                                ovlowconservedneg )    
    x[[names(list_ofbins_to_save_and_analyse[1])]] <-list_ofbins_to_save_and_analyse[[1]]
    x[[names(list_ofbins_to_save_and_analyse[2])]] <-list_ofbins_to_save_and_analyse[[2]]
    x[[names(list_ofbins_to_save_and_analyse[3])]] <-list_ofbins_to_save_and_analyse[[3]]
    x[[names(list_ofbins_to_save_and_analyse[4])]] <-list_ofbins_to_save_and_analyse[[4]]
          
    x[[paste0(ygroup,"_allgr_",xgroup)]] <-prvdftest_gr  

     names(x)<-c(paste0(ygroup,"_vs_",xgroup,"_all_shifting_bins"),
                 paste0(ygroup,"_",names(list_ofbins_to_save_and_analyse[1]),"_",xgroup),
                 paste0(ygroup,"_",names(list_ofbins_to_save_and_analyse[2]),"_",xgroup),
                 paste0(ygroup,"_",names(list_ofbins_to_save_and_analyse[3]),"_",xgroup),
                 paste0(ygroup,"_",names(list_ofbins_to_save_and_analyse[4]),"_",xgroup),
                 paste0(ygroup,"_allgr_",xgroup)
                )
    
    print(paste0(names(x[paste0(ygroup,"_vs_",xgroup,"_all_shifting_bins")]),"_",length(x[[paste0(ygroup,"_vs_",xgroup,"_all_shifting_bins")]]),"_bins"))
    return (x)
})