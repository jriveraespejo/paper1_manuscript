require(docstring)


entropy_func = function( N_i ){
  #' Shannon's (1948) entropy function
  #'
  #' @param N_i number of occurrences
  #'
  #' @return entropy value
  #' @export
  #'
  #' @examples
  
  pi = N_i / sum(N_i)
  H = -1 * ( sum( ( pi * log(pi, base=2) ) ) / log( sum(N_i), base=2) )
  return( H )
}






file_id = function(chains_path, model_int){
  #' file_id
  #' 
  #' it detects all stanfit object of interest in a folder
  #' it only applies for a study simulation with multiple conditions
  #'
  #' @param chains_path location of csv files corresponding to stanfit objects
  #' @param model_int character VECTOR with the names of the models of interest
  #'
  #' @return name list
  #' @export
  #'
  #' @examples
  
  # # test
  # chains_path = model_out
  # model_int = model_nam
  
  # packages
  require(stringr)
  
  # list all files
  chains_list = list.files( chains_path )
  
  # identify files of interest
  idx1 = str_detect(chains_list, '.csv')
  idx2 = str_detect(chains_list, paste0(model_int, '-[:digit:]', '.csv') )
  chains_list = chains_list[idx1 & idx2]
  
  return(chains_list)
}






plot_compare = function(waic_obj, psis_obj, ns=1, dM=T){
  #' Plot of model comparison
  #'
  #' @param waic_obj waic comparison object
  #' @param psis_obj psis comparison objects
  #' @param ns confidence interval threshold, default ns=1
  #' @param dM plot difference in criteria, default dm=T
  #'
  #' @return plot
  #' @export
  #'
  #' @examples
  
  # # test
  # waic_obj=RQ3_WAIC
  # psis_obj=RQ3_PSIS
  # ns=1
  # dM=T
  
  
  # extra calculations
  waic_obj$WAIC_lower = with( waic_obj, WAIC - ns*SE )
  waic_obj$WAIC_upper = with( waic_obj, WAIC + ns*SE )
  waic_obj$dWAIC_lower = with( waic_obj, dWAIC - ns*dSE )
  waic_obj$dWAIC_upper = with( waic_obj, dWAIC + ns*dSE )
  waic_obj$model = row.names(waic_obj)
  idx = which( names(waic_obj) %in% c('SE','dSE') )
  names(waic_obj)[idx] = c('SE_WAIC','dSE_WAIC')
  
  psis_obj$PSIS_lower = with( psis_obj, PSIS - ns*SE )
  psis_obj$PSIS_upper = with( psis_obj, PSIS + ns*SE )
  psis_obj$dPSIS_lower = with( psis_obj, dPSIS - ns*dSE )
  psis_obj$dPSIS_upper = with( psis_obj, dPSIS + ns*dSE )
  psis_obj$model = row.names(psis_obj)
  idx = which( names(psis_obj) %in% c('SE','dSE') )
  names(psis_obj)[idx] = c('SE_PSIS','dSE_PSIS')
  
  dA = merge(waic_obj, psis_obj[,-1], by='model')
  
  
  # selecting appropriate data (sorted)
  if(!dM){
    var_plot = c('DIC','WAIC','WAIC_lower','WAIC_upper',
                 'PSIS','PSIS_lower','PSIS_upper')
    dT = dA[,var_plot]
    dT = dT[ order(dT$WAIC, decreasing=T), ]
    
    x_lab = 'deviance'
    v_grid = min(dT$WAIC)
    leg = c( 'DIC','WAIC','PSIS' )
    leg_col = c( 'black', rgb(0,0,0,0.5), rgb(0,0,1,0.5) )
    point_pch = c(1,19,19)
    point_jif = c(0,-0.2,0.2)
    
  } else{
    var_plot = c('dWAIC','dWAIC_lower','dWAIC_upper',
                 'dPSIS','dPSIS_lower','dPSIS_upper')
    dT = dA[,var_plot]
    dT = dT[ order(dT$dWAIC, decreasing=T), ]
    
    x_lab = 'difference in deviance'
    v_grid = min(dT$dWAIC)
    leg = c('dWAIC','dPSIS')
    leg_col = c( rgb(0,0,0,0.5), rgb(0,0,1,0.5) )
    point_pch = c(19,19)
    point_jif = c(-0.2,0.2)
    
  }
   
  
  # plot range
  x_lim = range( dT, na.rm=T )
  if( max(x_lim) < 0 ){
    x_lim[ x_lim==max(x_lim) ] = 0
  } else if( min(x_lim) > 0 ){
    x_lim[ x_lim==min(x_lim) ] = 0
  }
  
  
  # empty plot
  par0 = par()
  
  par( mar=c(5.1,5.1,4.1,2.1) )
  plot( NULL, xlim=x_lim, ylim=c(0, nrow(dT)+1),
        xlab=x_lab, ylab='', yaxt='n' )
  axis(side=2, at=1:nrow(dT), labels=dA$model, las=2)
  abline( v=v_grid, h=1:nrow(dT), lty=2, col=rgb(0,0,0,0.3) )
  legend('top', horiz=T, bty='n', legend=leg, fill=leg_col)
  
  for(i in 1:length(leg)){
    points( dT[,leg[i]], (1:nrow(dT))+point_jif[i], 
            pch=point_pch[i], col=leg_col[i] )
  }
  
  start = ifelse( length(leg)>2, 2, 1 )
  for(i in start:length(leg)){
    var_int = paste0(leg[i], c('_lower','_upper'))
    for(j in 1:nrow(dT)){
      lines( x=dT[j, var_int], y=rep(j+point_jif[i],2), 
             lw=2, col=leg_col[i] )
    }
  }
  
  par( mar=par0$mar )
  
}



# setMethod("plot_compare" , "compareIC" , 
#           function( x, y, xlim, SE=TRUE, dSE=TRUE, weights=FALSE,...){
#             dev_in <- x[[1]] - x[[5]]*2 # criterion - penalty2
#             dev_out <- x[[1]]
#             if ( !is.null(x[['SE']]) ) devSE <- x[['SE']]
#             dev_out_lower <- dev_out - devSE
#             dev_out_upper <- dev_out + devSE
#             if ( weights==TRUE ) {
#               dev_in <- ICweights(dev_in)
#               dev_out <- ICweights(dev_out)
#               dev_out_lower <- ICweights(dev_out_lower)
#               dev_out_upper <- ICweights(dev_out_upper)
#             }
#             n <- length(dev_in)
#             if ( missing(xlim) ) {
#               xlim <- c(min(dev_in),max(dev_out))
#               if ( SE==TRUE & !is.null(x[['SE']]) ) {
#                 xlim[1] <- min(dev_in,dev_out_lower)
#                 xlim[2] <- max(dev_out_upper)
#               }
#             }
#             main <- colnames(x)[1]
#             set_nice_margins()
#             dotchart( dev_in[n:1] , labels=rownames(x)[n:1] , xlab="deviance" , pch=16 , xlim=xlim , ... )
#             points( dev_out[n:1] , 1:n )
#             mtext(main)
#             # standard errors
#             if ( !is.null(x[['SE']]) & SE==TRUE ) {
#               for ( i in 1:n ) {
#                 lines( c(dev_out_lower[i],dev_out_upper[i]) , rep(n+1-i,2) , lwd=0.75 )
#               }
#             }
#             if ( !all(is.na(x@dSE)) & dSE==TRUE ) {
#               # plot differences and stderr of differences
#               dcol <- col.alpha("black",0.5)
#               abline( v=dev_out[1] , lwd=0.5 , col=dcol )
#               diff_dev_lower <- dev_out - x$dSE
#               diff_dev_upper <- dev_out + x$dSE
#               if ( weights==TRUE ) {
#                 diff_dev_lower <- ICweights(diff_dev_lower)
#                 diff_dev_upper <- ICweights(diff_dev_upper)
#               }
#               for ( i in 2:n ) {
#                 points( dev_out[i] , n+2-i-0.5 , cex=0.5 , pch=2 , col=dcol )
#                 lines( c(diff_dev_lower[i],diff_dev_upper[i]) , rep(n+2-i-0.5,2) , lwd=0.5 , col=dcol )
#               }
#             }
#           }
# )


pred_obs = function(d, stanfit_obj){
  #' Observations prediction function
  #'
  #' @param d original data
  #' @param stanfit_obj stanfit object 
  #'
  #' @return model's prediction samples
  #' @export
  #'
  #' @examples
  
  # # test
  # d = data_H
  # stanfit_obj=model01
  
  # packages
  require(rethinking)
  require(stringr)
  
  # posterior distribution
  post = extract.samples(stanfit_obj)
  model = str_sub(stanfit_obj@model_name, start=1, end=7 )
  # str(post)
  
  # model groups
  model_normal_g1 = paste0( 'model0', 1:3 )
  model_normal_g2 = paste0( 'model0', 4:6 )
  model_beta_g1 = paste0( 'model0', 7:9 )
  model_beta_g2 = paste0( 'model', 10:12 )
  
  
  # data prediction
  H_var = matrix(NA, nrow=dim(post$mu)[1], ncol=dim(post$mu)[2])
  # n=1
  if( model %in% model_normal_g1 ){
    for(n in 1:nrow(d)){
      samples = nrow(H_var)
      mu = with(post, mu[, d$cid[n]] - b_i[, d$bid[n]] )
      sd = with(post, s_w)
      H_var[,n] = rnorm( n=samples, mean=mu, sd=sd )
    }
  } else if( model %in% model_normal_g2 ){
    for(n in 1:nrow(d)){
      samples = nrow(H_var)
      mu = with(post, mu[, d$cid[n]] - post$b_i[, d$bid[n]] )
      sd = with(post, s_w[, d$cid[n]] )
      H_var[,n] = rnorm( n=samples, mean=mu, sd=sd )
    }
  } else if( model %in% model_beta_g1 ){
    for(n in 1:nrow(d)){
      samples = nrow(H_var)
      mu = with(post, inv_logit( -SI[, d$cid[n]] ) )
      Mw = with(post, Mw )
      H_var[,n] = rbeta2( n=samples, prob=mu, theta=Mw )
    }
  } else if( model %in% model_beta_g2 ){
    for(n in 1:nrow(d)){
      samples = nrow(H_var)
      mu = with(post, inv_logit( -SI[, d$cid[n]] ) )
      Mw = with(post, Mw[, d$cid[n]] )
      H_var[,n] = rbeta2( n=samples, prob=mu, theta=Mw )
    }
  }
  # str(H_var)
  # Consideration:
  #   DO NOT include blocks in the calculation of the speakers, because these 
  #   blocks will never repeat themselves. When you consider blocks the predictions
  #   of both models are quite similar, which is not possible considering the
  #   abysmal difference in their log-likelihoods.
  
  
  # return object
  return( data.frame(H_var) )
  
  
}




pred_speaker = function(d, stanfit_obj, p=0.95, raw=F){
  #' speaker prediction function
  #'
  #' @param d original data
  #' @param stanfit_obj stanfit object 
  #' @param p significance level, default p=0.95
  #' @param raw raw sample data?, default raw=F
  #' 
  #' @return if raw=F returns precis_obj, otherwise list(precis_obj, raw_data)
  #' @export
  #'
  #' @examples
  
  # # test
  # d=data_H
  # stanfit_obj=model04
  # p=0.95
  # raw=T
  
  # packages
  require(rethinking)
  
  # speaker prediction
  pred = pred_obs(d=d, stanfit_obj=stanfit_obj)
  
  raw_sample = list()
  pred_precis = c()
  for(i in 1:max(d$cid)){
    
    idx = which(d$cid==i)
    pred_join = unlist( pred[,idx] )
    attr(pred_join, 'names') = NULL
    raw_sample[[i]] = pred_join
    
    pred_mom = precis(pred_join, prob=p)
    pred_mom = pred_mom[,-5]
    pred_mom = cbind( pred_mom, t( HPDI(pred_join, prob=p) ) )
    names(pred_mom)[3:6] = c('CI_lower','CI_upper','HPDI_lower','HPDI_upper')
    rownames(pred_mom) = paste0('speaker',i)
    
    pred_precis = rbind(pred_precis, pred_mom)
    
  }
  
  # return(object)
  if(raw){
    return( list(precis=pred_precis, samples=raw_sample) )
  } else{
    return( pred_precis )  
  }
  
}





plot_speaker = function(d, stanfit_obj1, stanfit_obj2, 
                        p=0.95, decreasing=T, 
                        leg=c('model1','model2')){
  #' speakers prediction plot function
  #'
  #' @param d original data
  #' @param stanfit_obj1 stanfit object, model 1
  #' @param stanfit_obj2 stanfit object, model 2
  #' @param p significance level, default p=0.95
  #' @param decreasing plot in decreasing order?, default decreasing=T
  #' @param leg character vector for legend, default leg=c('model1','model2')
  #' 
  #' @return model's prediction plots
  #' @export
  #'
  #' @examples
  
  # # test
  # d=data_H
  # stanfit_obj1=model06
  # stanfit_obj2=model12
  # p=0.95
  # decreasing=F
  # leg=c('normal','beta-proportion')
  
  
  # packages
  require(rethinking)
  
  # to sort by intelligibility
  pred1 = pred_speaker(d=d, stanfit_obj=stanfit_obj1)
  pred2 = pred_speaker(d=d, stanfit_obj=stanfit_obj2)
  
  speakers = order(pred2$mean, decreasing=decreasing)
  
  # sort true data
  data_mom = c()
  data_agg = c()
  for(i in 1:length(speakers)){
    idx = which(d$cid==speakers[i])
    dm = d[idx, c('cid','Hwsib')]
    dm$cid_mod = i
    
    data_agg = rbind(data_agg, sapply(dm, mean) )  
    data_mom = rbind(data_mom, dm)
  }
  
  
  # data plot
  y_lim = range( data.frame(pred1[,5:6], pred2[,5:6]))
  y_lim[2] = y_lim[2] + 0.2
  
  plot( data_mom$cid_mod, data_mom$Hwsib, 
        xaxt='n', ylim=y_lim,
        xlab='speakers', ylab='entropy',
        pch=19, col=rgb(0,0,0,0.15) )
  # points( data_agg[,3], data_agg[,2], col='red', pch=19)
  axis(side=1, at=1:length(speakers), labels=speakers, las=2)
  abline( h=c(0,1), lty=2, col=rgb(0,0,0, 0.3))
  
  
  # predictions
  col_sel = rep(NA,2)
  for(i in 1:2){
    
    if(i==1){
      pred = pred1
      col_sel[i] = col.alpha(rethink_palette[2], 1) 
      jitt = -0.2
      
    } else{
      pred = pred2
      col_sel[i] = col.alpha(rethink_palette[1], 1)
      jitt = +0.2
      
    }
    
    points( 1:length(speakers)+jitt, pred$mean[speakers], pch=19, 
            col=col_sel[i] )
    for(j in 1:length(speakers)){
      lines( x=rep(j+jitt,2), y=pred[speakers[j], c('HPDI_lower','HPDI_upper')],
             col=col_sel[i])
    }
    
  }
  
  legend('top', horiz=T, bty='n', legend=c('entropy',leg),
         fill=c('black', rethink_palette[2], rethink_palette[1]) )
  
}


pred_SI = function(d, stanfit_obj, p=0.95){
  #' Speech intelligibility per speaker
  #'
  #' @param d original data
  #' @param stanfit_obj stanfit object
  #' @param p significance level, default p=0.95
  #' 
  #' @return speech intelligibility plot
  #' @export
  #'
  #' @examples
  
  # # test
  # d=data_H
  # stanfit_obj=model07
  # p=0.95
  
  # package requirements
  require(stringr)
  require(rethinking)
  
  # posterior
  post = extract.samples(stanfit_obj)
  model = str_sub(stanfit_obj@model_name, start=1, end=7 )
  # str(post)
  
  # model groups
  model_normal = paste0( 'model0', 1:6 )
  if(model %in% model_normal){
    stop( 'Only works with models 7 to 12' )
  }
  
  # calculations
  post = post$SI
  SI = precis( data.frame(post), prob=p) #, pars=est_par
  SI = SI[,-5]
  hpdi_res = coda::HPDinterval(coda::as.mcmc(post), prob=p)
  SI = cbind( SI, hpdi_res) 
  names(SI)[3:6] = c('CI_lower','CI_upper','HPDI_lower','HPDI_upper')
  
  # join true data
  d_uniq = unique( d[ c('cid','HS','A','Am') ] ) 
  SI = cbind(d_uniq, SI)
  rownames(SI) = paste0( ifelse(1:nrow(SI) < 10, 'speaker0', 'speaker'), 1:nrow(SI))
  
  # return object
  return( round(SI, 3) )  
  
}



contrast_SI = function(d, stanfit_obj, speakers=1:32, p=0.95, raw=F){
  #' Speech intelligibility per speaker
  #'
  #' @param d original data
  #' @param stanfit_obj stanfit object
  #' @param speakers selection of speakers, default speakers=1:32
  #' @param p significance level, default p=0.95
  #' @param raw return raw data, default raw=T
  #' 
  #' @return speech intelligibility contrast data
  #' @export
  #'
  #' @examples
  
  # # test
  # d=data_H
  # stanfit_obj=model07
  # speakers=20
  # p=0.95
  # raw=T
  
  # package requirements
  require(stringr)
  require(rethinking)
  
  # posterior
  post = extract.samples(stanfit_obj)
  model = str_sub(stanfit_obj@model_name, start=1, end=7 )
  # str(post)
  
  # model groups
  model_normal = paste0( 'model0', 1:6 )
  if(model %in% model_normal){
    stop( 'Only works with models 7 to 12' )
  }
  
  # calculations
  post = post$SI
  # str(post)
  
  SI_raw = c()
  SI_final = c()
  SI_names = c()
  # i=1; j=2
  for( i in 1:dim(post)[2]){
    for( j in 1:dim(post)[2]){
      if(i<j){
        
        # name
        speaker_i = paste0( ifelse(i<10, 'speaker0', 'speaker'), i)
        speaker_j = paste0( ifelse(j<10, 'speaker0', 'speaker'), j)
        SI_names = c(SI_names, paste0( speaker_j, '-', speaker_i ))
        
        # raw data
        contr_res = post[,j] - post[,i]
        SI_raw = c( SI_raw, list(contr_res) )
        
        SI_contr = precis( data.frame(contr_res), prob=p) #, pars=est_par
        SI_contr = SI_contr[,-5]
        hpdi_res = coda::HPDinterval(coda::as.mcmc(contr_res), prob=p)
        SI_contr = cbind( SI_contr, hpdi_res) 
        names(SI_contr)[3:6] = c('CI_lower','CI_upper','HPDI_lower','HPDI_upper')
        
        # concatenate data
        SI_final = rbind(SI_final, SI_contr)
        
      } 
    }
  }
  names(SI_raw) = SI_names
  rownames(SI_final) = SI_names
  # str(SI_raw)
  
  # select speakers
  speakers = paste0( ifelse(speakers<10, 'speaker0', 'speaker'), speakers)
  if( length(speakers)==1 ){
    idx = str_detect( rownames(SI_final), speakers) 
  } else{
    idx = t( sapply( rownames(SI_final), str_detect, speakers) )
    idx = as.logical(rowSums(idx))
  }
  SI_raw = SI_raw[idx]
  SI_final = SI_final[idx, ]
  
  # return object
  if(raw){
    return( list(SI_contr=round(SI_final, 3), SI_raw=SI_raw) )
  } else{
  return( round(SI_final, 3) )  
  }
  
}



plot_SI = function(d, stanfit_obj, p=0.95, decreasing=F){
  #' Speech intelligibility plot per speaker
  #'
  #' @param d original data
  #' @param stanfit_obj stanfit object
  #' @param p significance level, default p=0.95
  #' @param decreasing plot in decreasing order?, default decreasing=T
  #' 
  #' @return speech intelligibility plot
  #' @export
  #'
  #' @examples
  
  # # test
  # d=data_H
  # stanfit_obj=model07
  # p=0.95
  # decreasing=T
  
  # generating SI
  SI = pred_SI(d=d, stanfit_obj=stanfit_obj, p=p)
  SI = SI[order(SI$mean, decreasing=decreasing),]
  
  
  # plot
  y_lim = with(SI, c( min(HPDI_lower)-0.5, max(HPDI_upper)+0.5 ) )
  
  plot( NULL, # <3> 
        xlim=c(1,32), ylim=y_lim, 
        xaxt='n', xlab='speakers', 
        ylab='Potential intelligibility' )
  abline( h=c(-3,-2,-1,0,1,2,3), lty=2, col=rgb(0,0,0,0.3) )
  
  axis( side=1, at=1:nrow(SI),
        labels=SI$cid, las=2, cex.axis=0.8 )
  
  for( i in 1:nrow(SI) ){ 
    points( i, SI$mean[i], # <4>
            pch=19, col=rgb(0,0,0,0.5) )
    
    lines( col=rgb(0,0,0,0.5), # <5>
           x=rep( i, each=2 ),
           y=c( SI$HPDI_lower[i],
                SI$HPDI_upper[i] ) )  
    
  }
  
}




binning = function(y, min_y, max_y, n_bins, dens=F){
  #' Binning function
  #'
  #' @param y data to bin
  #' @param min_y minimum value observed in y
  #' @param max_y maximum value observed in y
  #' @param n_bins number of bins
  #' @param dens in density?, default=F
  #'
  #' @return
  #' @export
  #'
  #' @examples
  
  # # test
  # y=dm
  # min_y=0
  # max_y=1
  # n_bins=20
  # dens=T
  
  # binning
  by_prog = round( (max_y-min_y)/(n_bins-1), 2)
  cuts = seq(from=min_y, to=max_y, by=by_prog)
  x = cut(y, cuts)
  x = table(x)
  if(dens){
    x = x/max(x) 
  }
  names(x) = round( cuts[ 1:(length(cuts)-1) ], 2)  
  
  # return object
  return(x)  
  
}





pred_speaker_pairs = function(speakers, d, stanfit_obj, 
                              p=0.95, nbins=30, col_string){
  #' Prediction of speaker pairs
  #'
  #' @param speakers pairs of speakers in one vector 
  #' @param d original data
  #' @param stanfit_obj stanfit object 
  #' @param p significance level, default p=0.95
  #' @param nbins binning of data for plot, default nbins=30
  #' @param col_string string of colors to use
  #'
  #' @return (3,3) plot of pairs of speakers
  #' @export
  #'
  #' @examples
  
  # # test
  # speakers = c(20,31, 24,30, 11,6)
  # d=data_H
  # stanfit_obj=model10
  # p=0.95
  # nbins=20
  
  # predict obs.
  pred = pred_speaker(d=d, stanfit_obj=stanfit_obj, p=p, raw=T)
  model = str_sub(stanfit_obj@model_name, start=1, end=7 )
  
  par(mfrow=c(2,3))
  # j=1
  for( j in 1:length(speakers) ){
    
    # data binning and plot  
    idx = d$cid==speakers[j]
    hs = unique( d$HS[idx] )
    dm = d$Hwsib[idx]
    dat = binning( y=dm, min_y=0, max_y=1, n_bins=nbins, dens=T )
    
    plot(dat, ylab="Frequency-Density", ylim=c(-0.15,max(dat)), xaxt='n',
         xlim=c(-0.05,1.05), xlab='entropy', col=rgb(0,0,0,0.6) )
    abline( h=0, col='gray' )
    abline( v=c(0,1), lty=2, col=rgb(0,0,0,0.3) )
    axis( side=1, at=as.numeric(names(dat)),
          labels=names(dat), las=2, cex.axis=0.8 )
    mtext( text=paste0('speaker ', speakers[j]), side=3, adj=0, cex=1.1)
    
    
    # prediction points and intervals
    points(pred$precis$mean[ speakers[j] ], -0.1, pch=19, col=col_string[hs] )
    lines(x=pred$precis[ speakers[j], c('HPDI_lower','HPDI_upper')], 
          y=rep(-0.1,2), col=col_string[hs])
    
    
    # samples binning and plot
    mom = binning( y=pred$samples[[ speakers[j] ]], 
                   min_y=0, max_y=1, n_bins=nbins, dens=T )
    rownames(mom) = as.character( as.numeric(rownames(mom)) + 0.012 )
    lines( mom, ylab="Frequency", ylim=c(-0.1,max(dat)), col=col_string[hs])
    
    if(j==1){
      legend('topright', legend=c('entropy', model), 
             fill=c(rgb(0,0,0,0.6), col_string ), bty='n' )
    }
    
  }
  
  par(mfrow=c(1,1))
  
}




hpdi = function(stanfit_obj, p=0.95) {
  #' Highest Percentile Density Interval (HPDI)
  #'
  #' @param stanfit_obj stanfit object
  #' @param p significance level, default p=0.95
  #'
  #' @return HPDI interval values
  #' @export
  #'
  #' @examples
  
  # # test
  # stanfit_obj=model06
  # p=0.95
  
  # packages
  require(rstan)
  require(runjags)
  
  # converting 
  samples = As.mcmc.list(stanfit_obj) # to coda mcmc
  samples = combine.mcmc(samples) # combine mcmc
  # str(samples)
  
  # calculating HPDI
  hpdi_res = coda::HPDinterval(samples, prob=p)
  # str(hpdi_res)
  
  return(hpdi_res)
  
}




par_recovery = function(stanfit_obj, est_par, p=0.95){
  #' Parameter recovery function
  #'
  #' @param stanfit_obj stanfit object
  #' @param est_par character vector with the names of parameters of interest
  #' @param p significance level, default p=0.95
  #'
  #' @return data.frame of sttatistics for parameters
  #' @export
  #'
  #' @examples
  
  # # test
  # stanfit_obj=model10
  # est_par='SI'
  # true_par=NULL
  # p=0.95
  
  # packages
  require(rethinking)
  
  
  # get the point estimates
  params = precis(stanfit_obj, depth=5, prob=p) #, pars=est_par
  
  # HPDI
  hpdi_res = hpdi( stanfit_obj, p=p)
  rem = which( str_detect( row.names(hpdi_res), 'log_lik') |
                 str_detect( row.names(hpdi_res), 'lp__') )
  hpdi_res = hpdi_res[-rem,]
  params = cbind( params, hpdi_res)
  
  # sort again
  params = params[,c(1:4,7:8,5:6)]
  names(params)[3:6] = c('CI_lower','CI_upper','HPDI_lower','HPDI_upper')
  
  
  # selecting parameters
  idx_par = c()
  for(j in 1:length(est_par)){
    idx_par = c(idx_par, which( str_detect( row.names(params), paste0('^',est_par[j],'[:punct:]*') )))
  }
  params = params[idx_par,]
  params = round( params, 3)
  
  
  # return object
  return(params)
  
}



contrast_intel = function(stanfit_obj, p=0.95, rope=c(-0.05,0.05)){
  #' Contrasts of HS parameters
  #'
  #' @param stanfit_obj stanfit object
  #' @param p significance level, default p=0.95
  #' @param rope region of practical equivalence, default rope=c(-0.05,0.05) 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  
  # # test
  # stanfit_obj=model12
  # p=0.95
  # rope=c(-0.05,0.05)
  
  # contrast
  post = extract.samples( stanfit_obj )
  # str(post)
  
  idx1 = as.logical( sum( names(post) == 'aHS' ) )
  idx2 = as.logical( sum( names(post) == 'bAmHS' ) )
  
  if(idx1 & idx2){
    post = with(post, data.frame(aHS[,2]-aHS[,1], bAmHS[,2]-bAmHS[,1]) )  
    colnames(post) = c('aHS[2]-aHS[1]', 'bAmHS[2]-bAmHS[1]')
  } else if( idx1 & !idx2){
    post = with(post, data.frame(aHS[,2]-aHS[,1]) )
    colnames(post) = 'aHS[2]-aHS[1]'
  } else{
    stop('there are no parameter to contrast')
  }
  
  contrs = precis(post, depth=2, prob=0.95) #, pars=est_par
  contrs = contrs[,-5]
  contrs = cbind( contrs, coda::HPDinterval(coda::as.mcmc(post), prob=p) )
  names(contrs)[3:6] = c('CI_lower','CI_upper','HPDI_lower','HPDI_upper')
  
  contrs$P_ROPE_lower = colSums(post < rope[1])/nrow(post)
  contrs$P_ROPE_upper = colSums(post > rope[2])/nrow(post)
  
  # return object
  return( round(contrs, 3) )

}




pred_intel = function(d, stanfit_obj, p=0.95, ns=500, seed=NULL){
  #' Model prediction for intelligibility
  #'
  #' @param d original data
  #' @param stanfit_obj stanfit object
  #' @param p significance level, default p=0.95
  #' @param ns number of samples to plot, default ns=500
  #' @param seed seed for replicatiom, defaul seed=NULL
  #'
  #' @return plot of intelligibility prediction
  #' @export
  #'
  #' @examples
  
  # # test
  # d=data_H
  # stanfit_obj=model10
  # p=0.95
  # ns=500
  # seed=NULL
  
  # packages
  require(RColorBrewer)
  require(dplyr)
  
  col_string = rethink_palette[c(3,5)] 
  
  # model name
  model = str_sub(stanfit_obj@model_name, start=1, end=7 )
  if(model %in% paste0('model0',1:6)){
    stop('plot is not available for these models')
  }
  
  # SI estimates
  SI = pred_SI(d=d, stanfit_obj=stanfit_obj, p=p)
  idx = with( SI, order( HS, Am, mean, decreasing=F ) )
  SI = SI[idx, ]
  
  # parameter estimates
  post = extract.samples(stanfit_obj)
  # str(post)
  
  # identification of parameters
  idx1 = as.logical( sum( names(post)=='aHS') )
  idx2 = as.logical( sum( names(post)=='bAmHS' ) )
  
 
  if(!is.null(seed)){
    set.seed(seed)
  }
  S = sample( 1:dim(post$mu)[1], ns, replace=F) 
  
  # hearing status
  HS = unique(SI$HS)
  
  # plot
  par( mfrow=c(1,2) )
  # hs=1
  for( hs in HS){
    
    # selecting observed data
    idx = SI$HS==hs
    SI_mom = SI[idx, ]
    
    # creating missing data
    idx = !( 0:max(SI$Am) %in% unique(SI_mom$Am) )
    A_miss = (0:max(SI$Am))[idx]
    A = A_miss + min(SI$A)
    
    dm = data.frame( matrix(NA, nrow=length(A_miss), ncol=ncol(SI_mom)) )
    names(dm) = names(SI_mom)
    dm[,c('HS','A','Am')] = data.frame(hs,A,A_miss)
    SI_mom = rbind( SI_mom, dm )
    SI_mom = SI_mom[order(SI_mom$Am),]
    
    # calculating regression line
    # idx = which( row.names(params_obj) %in% c('m_i','m_u') )
    if(idx1 & idx2){
      int = post$aHS[, hs]  
      slop = post$bAmHS[, hs]
    } else if(idx1 & !idx2){
      int = post$aHS[, hs]
      slop = post$bAm
    } else{
      int = post$a
      slop = rep(0.00, length(post$a))
    }
    
    # plot
    plot( SI_mom[,c('Am','mean')], 
          pch=19, col=col_string[hs], ylim=c(-2.5,4), xaxt='n',
          xlab='chronological age (months)', ylab='potential intelligibility' )
    axis( side=1, at=SI_mom$Am, 
          labels=SI_mom$A, las=2, cex.axis=0.8 )
    
    
    for( i in 1:nrow(SI_mom) ){ 
      lines( col=col_string[hs], # <5>
             x=rep( SI_mom$Am[i], each=2 ),
             y=c( SI_mom$HPDI_lower[i],
                  SI_mom$HPDI_upper[i] ) )  
    }
    
    
    for(s in S){
      abline( a=int[s], b=slop[s], lty=1, col=col.alpha(col_string[hs],0.03) )  
    }
    abline( a=mean(int), b=mean(slop), lty=2, col=col_string[hs], lwd=2.5 )
    
    
    with(SI_mom, text( x=Am-1.2, y=mean, cid, cex=0.8 ) )
    
    
    if( idx1 & idx2 ){
      t_lab = paste0('aHS[', hs, '] = ', round(mean(int),2), '\n ', 
                     'bAmHS[',hs,'] = ', round(mean(slop),2)) 
      text( x=27, y=-0.5, t_lab )
    } else if(idx1 & !idx2){
      t_lab = paste0('aHS[', hs, '] = ', round(mean(int),2), '\n ', 
                     'bAm = ', round(mean(slop),2))
      text( x=27, y=-0.5, t_lab )
    } else{
      t_lab = paste0('a = ', round( mean(int), 2), '\n ', 
                     'bAm = ', round( mean(slop), 2))
      text( x=27, y=3.5, t_lab )
    }
    
    
    if( hs==1 ){
      legend('topleft', legend=c('NH','HI/CI'), fill=col_string, bty='n' )
    }
    
  }
  
  par( mfrow=c(1,1) )
  
  
}




plot_outlier = function(d, stanfit_obj){
  #' Prediction of outliers
  #'
  #' @param d original data
  #' @param stanfit_obj stanfit object
  #'
  #' @return plot of outliers
  #' @export
  #'
  #' @examples
  
  # # test
  # d=data_H
  # stanfit_obj = model07
  
  # model name
  model = str_sub(stanfit_obj@model_name, start=1, end=7 )
  
  # params
  WAIC_E = WAIC(stanfit_obj, pointwise=T)
  PSIS_E = PSIS(stanfit_obj, pointwise=T)
  ps = with( d, paste0(cid, ',', uid ) )
  
  plot( PSIS_E$k, WAIC_E$penalty, 
        col=col.alpha(rangi2, 0.2), pch=19, lwd=2, xlim=c(0,1.5),
        xlab="PSIS Pareto k", ylab="WAIC penalty" )
  abline(v=0.5, lty=2)
  abline(v=0.7, lty=2, lwd=2)
  mtext( text=model, side=3, adj=0, cex=1.1)
  
  idx = which(PSIS_E$k >= 0.5)
  if( !length(idx)==0 ){
    text( x=PSIS_E$k[idx]-0.025, y=WAIC_E$penalty[idx]-0.025, ps[idx] )  
  }
  
}




trace_plot = function(stan_object, pars) {
  #' Trace plot
  #'
  #' @param stan_object stanfit object
  #' @param pars parameters of interest
  #'
  #' @return trace plots
  #' @export
  #'
  #' @examples
  
  # # test
  # stan_object = res
  # pars = 'mu_a'
  
  # packages 
  require(RColorBrewer)
  require(rethinking)
  
  # posterior
  post = rstan::extract(stan_object, pars=pars, permuted=FALSE)
  # str(post)
  
  # parameters
  n_chains = dim(post)[2]
  chain.cols = rep_len(rethink_palette, n_chains)
  wstart = 1
  wend = dim(post)[1]
  ylim = range(post[wstart:wend, , ])
  ytick = (ylim[2] - ylim[1])/6
  yaxis = round( seq(ylim[1], ylim[2], by=ytick), 2)
  neff = summary(stan_object)$summary[, "n_eff"]
  neff_use <- neff[names(neff) == pars]
  
  # plot
  plot(NULL, type="l", xlim=c(wstart, wend), ylim=ylim,
       xlab="", ylab="", axes=F)
  box(bty="l")
  axis(side=1, at=seq(0, wend, by=100))
  axis(side=2, at=yaxis, las=1 )
  #mtext(paste("n_eff =", round(neff_use, 0)), 3, adj = 1, cex = 1.1)
  mtext(pars, 3, adj = 0, cex=1.1)
  for(c in 1:n_chains){
    lines(1:wend, post[, c, ], col=chain.cols[c], lwd = 0.5)
  }
  
}




trank_plot = function(stan_object, pars, wide=50){
  #' Trace rank plot (trank plot)
  #'
  #' @param stan_object stanfit object
  #' @param pars parameters of interest
  #' @param wide number of iterations considered, default wide=50
  #'
  #' @return trank plots
  #' @export
  #'
  #' @examples
  
  # # test
  # stan_object = res
  # pars = 'mu_a'
  # wide=50
  
  # for colors
  require(RColorBrewer)
  require(rethinking)
  
  # posterior
  post = rstan::extract(stan_object, pars=pars, permuted=FALSE)
  # str(post)
  
  # parameters
  n_chains = dim(post)[2]
  chain.cols = rep_len(rethink_palette, n_chains)
  wstart = 1
  wend = dim(post)[1]
  neff = summary(stan_object)$summary[, "n_eff"]
  neff_use <- neff[names(neff) == pars]
  
  # rank calculation
  ranks = list()
  xrange = rep(1:wide, each=2)
  yrange = vector('list', n_chains)
  for(c in 1:n_chains){
    ranks[[c]] = rank( post[1:(wide+1), c, ] )
    y_ran = c()
    for(i in 2:(wide+1)){
      y_ran = c(y_ran, c( ranks[[c]][i-1], ranks[[c]][i] ) )
    }
    yrange[[c]] = y_ran
  }
  
  
  # plot
  plot(NULL, type='l', xlim=c(0, wide+1), ylim=c(0, wide+1),
       xlab="", ylab="", xaxt ="n", yaxt ="n", axes=F)
  box(bty="l")
  # mtext(paste("n_eff =", round(neff_use, 0)), 3, adj = 1, cex = 0.9)
  # mtext(pars, 3, adj = 0, cex = 1)
  for(c in 1:n_chains){
    lines(xrange, yrange[[c]], col=chain.cols[c], lwd=1.5)  
  }
  
}



acf_plot = function(stan_object, pars){
  #' Autocorrelation function plot (acf plot)
  #'
  #' @param stan_object stanfit object
  #' @param pars parameters of interest
  #'
  #' @return acf plots
  #' @export
  #'
  #' @examples
  
  # # test
  # stan_object = res
  # pars = 'mu_a'
  
  # posterior
  post = rstan::extract(stan_object, pars=pars, permuted=FALSE)
  # str(post)
  
  # plot
  acf( post[, 1, ], main='', xlab='', ylab='', mar = c(0, 0, 0, 0) )
  # mtext(paste("n_eff =", round(neff_use, 0)), 3, adj = 1, cex = 0.9)
  # mtext(pars, 3, adj = 0, cex = 1)
  
}




tri_plot = function(stan_object, pars){
  #' Trace, trank and acf plots
  #'
  #' @param stan_object stanfit object
  #' @param pars parameters of interest
  #'
  #' @return trace, trank and acf plots
  #' @export
  #'
  #' @examples
  
  # # test
  # stan_object = stan_model
  # pars = paste0('m_b[',1:5,']')
  
  # figure parameters
  opar = par()
  
  # ensure there is only 5 paramaters
  if(length(pars)>5){
    pars = pars[1:5]
  }
  
  # plot
  par(mfrow=c(length(pars), 3), mar=c(3,3.5,1.5,1)+0.1)
  
  for(i in 1:length(pars)){
    trace_plot(stan_object, pars=pars[i]) 
    trank_plot(stan_object, pars=pars[i])
    acf_plot(stan_object, pars=pars[i])
  }
  
  par(mfrow=c(1,1), mar=opar$mar)
  
}





dens_plot = function(stanfit_obj, pars, p=0.95){
  #' Density plots with HPID
  #'
  #' @param stanfit_obj stanfit object
  #' @param pars parameters of interest
  #' @param p significance level, default p=0.95
  #'
  #' @return parameters' density plot with HPDI
  #' @export
  #'
  #' @examples
  
  # # test
  # stanfit_obj = model12
  # pars = 'Mw'
  # p=0.95
  
  require(rethinking)
  
  # getting names 
  post = extract.samples(stanfit_obj)
  
  idx = names(post) %in% pars 
  post = post[idx]
  # str(post)
  
  post_list = list()
  nam = c()
  # i=2
  for( i in 1:length(post) ){
    check = is.na( dim(post[[i]])[2] )
    if( check ){
      nam = c( nam, names(post)[i] )
      post_list = c( post_list, post[i] )
    } else{
      dim2 = dim(post[[i]])[2]
      # j=1
      for( j in 1:dim2 ){
        nam = c( nam, paste0( names(post)[i], '[', j, ']' ) )
        post_list = c( post_list, list( post[[i]][,j] ) )
      }
    }
  }
  names(post_list) = nam
  # str(post_list)

    
  # plot
  nplot = min( length(post_list), 6)
  par(mfrow=c(2,3))
  for(i in 1:nplot){
    dens( x=post_list[[i]], show.HPDI=p, 
          xlab='', col=rgb(0,0,0,0.5) )
    abline( v=mean(post_list[[i]]), lty=2, col='black' )
    mtext( text=nam[i], side=3, adj=0, cex=1.1)
    if(i ==1){
      legend( 'topright', legend=c('mean', paste0(p*100,'% HPDI') ),
              col=c('black','gray'), fill=c('black','gray'), 
              lty=c(2,1), bty='n' )
    }
  }
  par(mfrow=c(1,1))
  
}

