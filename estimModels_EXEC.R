source("FLSA.R"); source("utils.R"); library(Rdimtools)


cv_evalTopics_FLSA = function(K=5,cvs=NULL,Xw=NULL,dtm_matrix=NULL,words_dic=NULL,noc=NULL,J=NULL,Utype=NULL){
  Cvout = matrix(nrow=9,ncol=K)
  for(k in 1:K){
    cat("\n FLSA: cv no.: ",k)
    train_iid = cvs$subsets[cvs$which==k]; test_iid = cvs$subsets[cvs$which!=k]
    
    train_model = FLSA_core(X=Xw[train_iid,],N = length(train_iid),M = NCOL(Xw),K = noc,type = Utype,J_svd = J)
    fcm_test <- Matrix::crossprod(Matrix::Matrix(dtm_matrix[test_iid,]),Matrix::Matrix(dtm_matrix[test_iid,]))
    pdT_test <- FLSA_core(X=Xw[test_iid,],N = length(test_iid),M = NCOL(Xw),K = noc,type = "svd")$PD_T
    PW_T <- train_model$PW_T
    PW_T_sort = mapply(function(j)names(PW_T[order(PW_T[,j],decreasing = TRUE),j][1:30]),1:noc) #most frequent terms for each topic
    Cvout[1:3,k] = coherence(dtm_data=dtm_matrix[test_iid,],fcm_data=fcm_test,topic_words_dist=t(train_model$PW_T),top_terms_topics=PW_T_sort,
                             average=TRUE,smth=1,probcoh = FALSE,
                             metrics_text2vec = c("mean_logratio","mean_pmi","mean_npmi"))[1:3]
    
    rpca_out <- rpca::rpca(M = pdT_test,term.delta = 1e-05)
    cl <- rl <- stb1 <- stb2 <- stb3 <- NA
    if(rpca_out$convergence$converged){
      
      cl <- 1-norm(rpca_out$S,type="F")/norm(rpca_out$L+rpca_out$S,type="F")
      rl <- 1-length(rpca_out$L.svd$d)/NCOL(rpca_out$L)
      
      pp_raw <- apply(rpca_out$L,1,sum)/apply(rpca_out$L+rpca_out$S,1,sum)*100
      out_smooth <- smooth.spline(pp_raw,df = median(pp_raw))
      pp_smooth <- out_smooth$y
      
      out_deriv <- predict(out_smooth, 1:length(pp_raw), deriv = 1)
      stb1 <- mean(abs(out_deriv$y - mean(out_deriv$y)))
      stb2 <- sum(diff(sign(out_deriv$y))!=0)
      stb3 <- mean(diff(sign(out_deriv$y)))
      
      Cvout[4:8,k] <- c(cl,rl,stb1,stb2,stb3)
    }
    Cvout[9,k] <- mean_inf.rm(mapply(function(j)dist_from_corpus(PW_T=PW_T,dtm_data=dtm_matrix[test_iid,]),1:noc))
    
  }
  return(Cvout)
}


args <- commandArgs(trailingOnly = TRUE); for (arg in args) {eval(parse(text = arg))} #to retrieve input args
#Input: noc, K

load("casestudy_data.RData")
dtm_matrix <- dtm_matrix_filtered

cvs <- cvTools::cvFolds(NROW(dtm_matrix),K,type="random")[c(4,5)]


out_J2_std <- list()
out_J2_std[[1]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = dtm_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 2,Utype = "svd")},error=function(e){return(NA)})
out_J2_std[[2]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = dtm_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 2,Utype = "nnmf")},error=function(e){return(NA)})
out_J2_std[[3]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = dtm_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 2,Utype = "le")},error=function(e){return(NA)})

out_J3_std <- list()
out_J3_std[[1]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = dtm_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 3,Utype = "svd")},error=function(e){return(NA)})
out_J3_std[[2]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = dtm_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 3,Utype = "nnmf")},error=function(e){return(NA)})
out_J3_std[[3]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = dtm_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 3,Utype = "le")},error=function(e){return(NA)})

out_J4_std <- list()
out_J4_std[[1]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = dtm_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 4,Utype = "svd")},error=function(e){return(NA)})
out_J4_std[[2]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = dtm_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 4,Utype = "nnmf")},error=function(e){return(NA)})
out_J4_std[[3]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = dtm_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 4,Utype = "le")},error=function(e){return(NA)})


out_J2_pcp <- list()
out_J2_pcp[[1]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = L_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 2,Utype = "svd")},error=function(e){return(NA)})
out_J2_pcp[[2]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = L_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 2,Utype = "nnmf")},error=function(e){return(NA)})
out_J2_pcp[[3]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = L_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 2,Utype = "le")},error=function(e){return(NA)})

out_J3_pcp <- list()
out_J3_pcp[[1]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = L_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 3,Utype = "svd")},error=function(e){return(NA)})
out_J3_pcp[[2]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = L_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 3,Utype = "nnmf")},error=function(e){return(NA)})
out_J3_pcp[[3]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = L_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 3,Utype = "le")},error=function(e){return(NA)})

out_J4_pcp <- list()
out_J4_pcp[[1]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = L_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 4,Utype = "svd")},error=function(e){return(NA)})
out_J4_pcp[[2]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = L_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 4,Utype = "nnmf")},error=function(e){return(NA)})
out_J4_pcp[[3]] = tryCatch({cv_evalTopics_FLSA(K = K,cvs = cvs,Xw = L_matrix,dtm_matrix = dtm_matrix,words_dic = words_dic,noc = noc,J = 4,Utype = "le")},error=function(e){return(NA)})

out_std <- list('J2'=out_J2_std,'J3'=out_J3_std,'J4'=out_J4_std)
out_pcp <- list('J2'=out_J2_pcp,'J3'=out_J3_pcp,'J4'=out_J4_pcp)

saveRDS(out_std,file = paste0("results/dataout_std_",noc,".rds"))
saveRDS(out_pcp,file = paste0("results/dataout_pcp_",noc,".rds"))


