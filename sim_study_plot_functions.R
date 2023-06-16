



make_rocs = function(P, gt, gt.po, name = "TAP"){
  if(length(P) == 0){
    return(NULL)
  }
  rocs = list()
  less = 0
  for(rc in 1:length(P)){
    # rc = 1
    if(is.null(P[[rc]])){
      less = less + 1
    } else{
      gt.mat = gt[[rc]]
      pred.mat = P[[rc]]
      for(t in 1:nrow(gt.mat)){
        pred.mat[t,t] = NA
        gt.mat[t,t] = NA
      }
      gt.po.mat = gt.po[[rc]]
      gt.po.vec = as.vector(gt.po.mat)
      gt.vec = as.vector(gt.mat)
      pred.vec = as.vector(pred.mat)
      pred.vec[which(is.na(gt.vec))] = NA
      gt.vec[which(!is.na(gt.po.vec))] = NA
      gt.vec = gt.vec[which(!is.na(gt.vec))]
      pred.vec[which(!is.na(gt.po.vec))] = NA
      pred.vec = pred.vec[which(!is.na(pred.vec))]
      
      theroc = pROC::roc(formula = gt.vec ~ pred.vec,
                         quiet = TRUE)
      rocs[[rc - less]] = theroc
      names(rocs)[rc - less] = paste(name , rc - less)
    }
  }
  return(rocs)
}


# For the labels: It is important that the rocs are named via: "LABEL SOMETHING" with a blank space as seperator.

plot_roc = function(rocs, labels, add.points = NULL, alp = 1,
                    title = "ROC curves"){
  
  nbr.labels = length(labels)
  indices.included = NULL
  indices = list()
  nbr.roc.pts = rep(0, times = nbr.labels)
  labels.list = stringr::str_split(names(rocs), " ")
  labels.list = lapply(labels.list, function(l){l[1]})
  
  for(t in 1:nbr.labels){
    #t = 1
    indices[[t]] = which(labels.list == labels[t])
    subset.rocs = rocs[indices[[t]]]
    length(subset.rocs$`TAP 1`$sensitivities)
    nbr.roc.pts[t] = sum(unlist(lapply(subset.rocs, function(l){length(l$sensitivities)})))
    indices.included = c(indices.included, indices[[t]])
  }
  rocs = rocs[indices.included]
  my.colors = c("#006600","#33CC33","#66FF00","#CCFF00","#FFFF00","#FFCC33","#FF9900",
                "#FF3300","#CC3300","#990000","#CC3366","#FF6699","#FF33CC","#993366",
                "#996699","#660066","#6600CC","#3300FF","#000099","#0033CC","#3366FF",
                "#6699FF","#66CCFF","#33CCFF","#0099CC","#339999","#00CCCC","#00FFFF",
                "#33FFCC","#00CC99","#00FF99","#CCCCCC","#999999","#333333","#99CC99",
                "#669966")
  if(is.null(add.points)){
    nbr.pt.labels = 0
  } else{
    nbr.pt.labels = length(add.points)
  }
  nbr.tot.labels = nbr.labels + nbr.pt.labels
  col.indices = floor(((0:(nbr.tot.labels-1))*(length(my.colors)-1)/(nbr.tot.labels-1)) + 1)
  clrs = NULL
  actual.alpha = alp
  for(nu in 1:nbr.labels){
    clrs = c( clrs, rep(alpha(my.colors[col.indices[nu]],actual.alpha), times = nbr.roc.pts[nu]))
  }
  # View(names(rocs))
  if(is.null(add.points)){
    
    ggroc(rocs, aes = "alpha", color = clrs, legacy.axes = TRUE) +
      theme_bw() +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "black", alpha = .5) +
      labs(title = title,
           x = "False Positive Rate",
           y = "True Positive Rate") +
      scale_alpha_manual(values = rep(actual.alpha, times = 200),
                         limits = c(labels),
                         breaks = c(labels)) +
      guides(alpha = guide_legend(title = "Methods",
                                  override.aes = list(color = my.colors[col.indices])))
    
  } else{
    pt.frame = data.frame("Methods" = rep(names(add.points)[1], times = nrow(add.points[[1]])),
                          "FPR" = add.points[[1]][,1],
                          "TPR" = add.points[[1]][,2])
    if(nbr.pt.labels >= 2){
      for(k in 2:nbr.pt.labels){
        pt.frame.add = data.frame("Methods" = rep(names(add.points)[k], times = nrow(add.points[[k]])),
                                  "FPR" = add.points[[k]][,1],
                                  "TPR" = add.points[[k]][,2])
        pt.frame = rbind(pt.frame, pt.frame.add)
      }
    }
    
    ggroc(rocs, aes = "alpha", color = clrs, legacy.axes = TRUE) +
      theme_bw() +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "black", alpha = .5) +
      labs(title = "ROC curves",
           x = "False Positive Rate",
           y = "True Positive Rate") +
      scale_alpha_manual(values = rep(actual.alpha, times = 200),
                         limits = c(labels),
                         breaks = c(labels)) +
      geom_point( data = pt.frame, aes(x = FPR, y = TPR, color = Methods),
                  shape = 9,
                  alpha = actual.alpha , size = 2.5) +
      scale_color_manual( values = my.colors[col.indices[(nbr.labels + 1): (nbr.labels+nbr.pt.labels)]]) +
      guides(alpha = guide_legend(title = "ROC Methods",
                                  override.aes = list(color = my.colors[col.indices[1:nbr.labels]])))
  }
}


get_mean_roc = function(rocs, pt.dense = 1001){
  if(is.null(rocs)){
    return(NULL)
  }
  corrd.mat = matrix(0, nrow = pt.dense, ncol = length(rocs) )
  pts = seq(from = 0, to = 1, length.out = pt.dense)
  for(t in 1:length(rocs)){
    corrd.mat[,t] = as.vector(unlist(pROC::coords(roc = rocs[[t]], x = pts, input = "specificity", ret = "sensitivity")))
  }
  mean.curve = rocs[[1]]
  mean.curve$specificities = c(0,pts,1)
  mean.curve$thresholds = c(0, (1/(2*pt.dense)), pts[2:(pt.dense-1)], 1 - (1/(2*pt.dense)) , 1)
  mean.curve$sensitivities = c(1,rowMeans(corrd.mat),0)
  return(mean.curve)
}


make_auc_comp = function(auc_table, title = "Comparison of AUCs", x_axis = "p", x_limits = c(0,1)){
  
  ggplot(data = auc_table, aes(y = AUC, x = as.factor(auc_table[,3]), color = Method, fill= Method)) +
    geom_boxplot() +
    theme_bw() +
    labs(title = title,
         x = x_axis,
         y = "AUC") +
    scale_y_continuous(limits = x_limits) +
    geom_hline(yintercept = 0.5)
  
}


make_time_comp = function(time_table, title = "Comparison of AUCs", x_axis = "p", x_limits = c(0,1)){
  
  
  
  ggplot(data = time_table, aes(y = Time, x = p, color = Method, fill= Method)) +
    # geom_boxplot() +
    stat_summary(fun = mean, geom = "line") +
    theme_bw() +
    labs(title = title,
         x = x_axis,
         y = "Time in s") +
    scale_y_continuous(trans = "log1p")
  # scale_y_continuous(limits = x_limits) +
  # geom_hline(yintercept = 0.5)
  
}



# pred_graphs = P.pc
# gt_graphs = G.As
# observed_graphs = G.Apos
# diag.diff = NULL


get_rocpoints = function(pred_graphs, gt_graphs, observed_graphs, diag.diff = NULL){
  if(length(pred_graphs) == 0){
    return(NULL)
  }
  tim = length(pred_graphs)
  for(cc in 1:tim){
    # cc = 2
    pred_graph = pred_graphs[[cc]]
    gt_graph = gt_graphs[[cc]]
    chk = observed_graphs[[cc]][,1]
    chk[1] = observed_graphs[[cc]][1,2]
    int_indices = which(!is.na(chk))
    
    gt_graph[int_indices,] = NA
    pred_graph[int_indices,] = NA
    gt.mat = gt_graph
    pred.mat = pred_graph
    if(is.null(diag.diff)){
      for(t in 1:nrow(gt.mat)){
        pred.mat[t,t] = NA
        gt.mat[t,t] = NA
      }
    } else{
      for(t in 1:nrow(gt.mat)){
        pred.mat[t,diag.diff[t]] = NA
        gt.mat[t,diag.diff[t]] = NA
      }
    }
    gt.vec = as.vector(gt.mat)
    pred.vec = as.vector(pred.mat)
    pred.vec[which(is.na(gt.vec))] = NA
    gt.vec = gt.vec[which(!is.na(gt.vec))]
    pred.vec = pred.vec[which(!is.na(pred.vec))]
    fpr = length(which((gt.vec == 0) & (pred.vec != gt.vec)))/length(which(gt.vec == 0))
    tpr = length(which((gt.vec == 1) & (pred.vec == gt.vec)))/length(which(gt.vec == 1))
    if(cc == 1){
      output = matrix(c(fpr,tpr), ncol = 2, nrow = 1)
    } else{
      output = rbind(output , c(fpr,tpr))
    }
  }
  return(output)
}


make_list_avoiding_NULL = function(initial.list){
  list.len = length(initial.list)
  to.remove = 0
  counter = 1
  for(t in 1:list.len){
    if(is.null(initial.list[[t]])){
      to.remove[counter] = t
      counter = counter + 1
    }
  }
  if(length(to.remove) > 1){
    initial.list = initial.list[-to.remove]
  } else {if(to.remove != 0){
    initial.list = initial.list[-to.remove]
  }
  }
  output = list(list = initial.list,
                to.remove = to.remove)
  return(output)
}


