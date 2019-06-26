


#' Title
#'
#' @param s location coordinate matrix (generally a N by 2 matrix)
#' @param phi_mat basis function matrix
#'
#' @return B matrix in thesis p18
B_gen <- function(s, phi_mat){
        
        s_dim = ncol(s)
        if(s_dim==1){
                g <- function(r){(r^3)/12}
                G <- g(fields::rdist(s))
        }
        if(s_dim==2){
                g <- function(r){(r^2)*log(r)/(16*pi)}
                G <- g(fields::rdist(s)+diag(0.0001,nrow=nrow(s)))
        }
        E <- cbind(rep(1,nrow(s)),s)
        GE = cbind(rbind(G,t(E)),rbind(E,matrix(0,nrow=ncol(E),ncol=ncol(E))))
        phi_aug = rbind(phi_mat, matrix(0, nrow=(nrow(GE)-nrow(phi_mat)), ncol=ncol(phi_mat)))
        B = solve(GE, phi_aug)
}


#generate phi(S0), in thesis p18
new_loc_phi <- function(s_new, s_old, B){
        s_dim = ncol(s_new)
        if(s_dim==1){
                g <- function(r){(r^3)/12}
        }
        if(s_dim==2){
                g <- function(r){(r^2)*log(r)/(16*pi)}
        }
        dist_mat = fields::rdist(s_new, s_old)
        zero_index = which(dist_mat==0)
        G_new <- g(fields::rdist(s_new, s_old))
        G_new[zero_index]=0
        E_new <- cbind(rep(1,nrow(s_new)), s_new)
        new_phi <- cbind(G_new, E_new) %*% B
}

# generate Ca, aka Sigma_a, in thesis p10 
V_mat_calc <- function(A, lambda_mat){
        lambda_mat - A %*% lambda_mat %*% t(A)
}

# generate lambda matrix by A & Ca
lambda_mat_calc <- function(A, Ca){
        k = nrow(A)
        vec = solve(diag(k^2)-kronecker(A, A)) %*% as.vector(Ca)
        matrix(vec, nrow=k)
}

# generate data in thesis p23~24
data_gen_1d_2 <- function(s_old_mat, s_new_mat, k, n,
                          lambda_mat, sigma_sq, A, V, 
                          normalize=FALSE){ 
        
        #--------------------------
        #Step 1: Initial setting
        #--------------------------
        
        s <- s_old_mat
        p = nrow(s)
        d = ncol(s)
        
        #--------------------------
        #Step 2: generate phi & phi_mat_new
        #--------------------------
        
        if(d==1){
                phi_1 <- exp(-s^2)/norm(as.matrix(exp(-s^2)),type="F")
                phi_2 <- s*exp(-s^2)/norm(as.matrix(s*exp(-s^2)), "F")
                phi_mat <- cbind(phi_1, phi_2)
        }
        if(d==2){
                phi_1 <- exp(-(s[,1]^2+s[,2]^2))/norm(as.matrix(exp(-(s[,1]^2+s[,2]^2))), "F")
                phi_2 <- s[,1]*s[,2]*exp(-(s[,1]^2+s[,2]^2))/norm(as.matrix(s[,1]*s[,2]*exp(-(s[,1]^2+s[,2]^2))), "F")
                phi_mat <- cbind(phi_1, phi_2)
        }
        phi_mat_new <- new_loc_phi(s_new_mat, s, B=B_gen(s, phi_mat))
        
        
        #--------------------------
        #Step 3: generate xi matrix
        #--------------------------
        
        xsi_1 <- t(mvtnorm::rmvnorm(n=1, mean=c(0,0), sigma=lambda_mat))
        eta <- mvtnorm::rmvnorm(n=n, mean=c(0,0), sigma=V)
        xsi <- matrix(0,nrow=n,ncol=k)
        for(i in 1:n){
                if(i==1){xsi[i,] <- t(xsi_1)}
                if(i!=1){
                        xsi[i,] <- A %*% xsi_1 + eta[i,]
                        xsi_1 <- matrix(xsi[i,],nrow=k)
                }
                
        }
        if(normalize==T){xsi = scale(xsi)}
        
        #--------------------------
        #Step 4: generate Y matrix
        #--------------------------
        
        epsilon <- matrix(rnorm(sd=sqrt(sigma_sq),n*p), nrow=n, ncol=p)
        y_mat <- xsi %*% t(phi_mat) + epsilon
        
        
        #--------------------------
        #Step 5: output
        #--------------------------
        
        list(y_mat=y_mat, xsi_mat=xsi, epsilon_mat=epsilon,
             location=s, phi_mat=phi_mat, 
             location_new=s_new_mat, phi_mat_new=phi_mat_new)
}



# ---------------------------------------------

# our method (parameter estimation via EM)
EM_func3 = function(s, y_mat, y_mat_new, sigma2_eps=1, itermax = 30, 
                    tol = 0.001, tau, k){
        
        # generate phi_int
        yty <- t(y_mat) %*% y_mat
        T=nrow(y_mat)
        N=ncol(y_mat)
        init <- RSpectra::eigs_sym(yty, k, which="LM")$vectors # fast eigen-decomposition
        
        #omega matrix
        s_dim = ncol(s)
        if(s_dim==1){
                g <- function(r){(r^3)/12}
                G <- g(fields::rdist(s))
        }
        if(s_dim==2){
                g <- function(r){(r^2)*log(r)/(16*pi)}
                G <- g(fields::rdist(s)+diag(0.0001,nrow=nrow(s)))
        }
        E <- cbind(rep(1,nrow(s)),s)
        G_inv <- solve(G)
        G_inv_E <- G_inv %*% E
        omega = G_inv - G_inv_E %*% solve(t(E)%*%G_inv_E, t(G_inv_E))
        
        xi_mat = y_mat %*% init
        Cov0 = crossprod(xi_mat)/nrow(xi_mat)
        C0inv = solve(Cov0)
        Cov1 = crossprod(xi_mat[-1,], xi_mat[-nrow(xi_mat),])/(nrow(xi_mat)-1)
        Ceta = Cov0 - Cov1 %*% C0inv %*% t(Cov1)
        Ceta = as.matrix(Matrix::forceSymmetric(Ceta))
        Mest = Cov1%*%C0inv
        H = init
        muinit = matrix(0,k,1)
        
        # for proj.grad.descent setting
        tol.m = .Machine$double.eps
        tol1=1e-8
        tol2=1e-6
        
        EM_nonfixedPhi_rcpp2(y_mat=y_mat, y_mat_new=y_mat_new, 
                             Phi=init, mu0 = rep(0,k),
                             Cov0 = Cov0,  
                             A=Mest, Ca=Ceta, sigma2_eps=1,
                             itermax=30, tol=tol,
                             tau=tau, omega=omega, tol_m=tol.m, tol1=tol1, tol2=tol2)
}

temp_spat_cv_final3 <- function(s, y_mat, sigma2_eps=1, itermax = 30, 
                                tol = 0.001, tau, k){
        
        p=ncol(y_mat)
        s_dim = ncol(s)
        fold_index <- rep(1:5, each=ceiling(nrow(y_mat)/5), length=nrow(y_mat))
        
        k_result <- numeric(length=length(k))
        
        if(length(tau)>1){
                print("Start tau tuning...")
                tau_list = list()
                k_like <- c()
                for(k_index in 1:length(k)){
                        k_select = k[k_index]
                        tau_cv_mat <- matrix(0, nrow=5, ncol=length(tau))
                        colnames(tau_cv_mat) <- as.character(round(tau,2))
                        for(kk in 1:5){
                                train_cv <- y_mat[which(fold_index!=kk),]
                                test_cv <- y_mat[which(fold_index==kk),]
                                
                                n_test <- nrow(test_cv)
                                n_cv = nrow(train_cv)
                                p_cv = ncol(train_cv)
                                
                                tau_cv_mat[kk,] = foreach::foreach (i = 1:length(tau), .combine='c') %dopar% {
                                        
                                        tau.gd <- tau[i]
                                        result = EM_func3(s=s, y_mat=train_cv, 
                                                          y_mat_new=test_cv, 
                                                          sigma2_eps=1, 
                                                          itermax = itermax, 
                                                          tol = tol, 
                                                          tau=tau[i], 
                                                          k=k_select)
                                        result$test_like
                                }
                                print(paste("Tune tau Fold- ", kk, " of 5 is done!!"))
                        }
                        tau_cv_result = colMeans(tau_cv_mat)
                        stau <- tau[which.min(tau_cv_result)]
                        k_like[k_index] <- min(tau_cv_result)
                        
                        process_list = list(k=k_select, stau=stau, tau=tau, 
                                            tau_cv_score=k_like[k_index],
                                            tau_cv_mat=tau_cv_mat, 
                                            tau_cv_result=tau_cv_result)
                        k_name<-paste('k_', k_select, sep = '')
                        tau_list[[k_name]] <-  process_list
                        print(paste("For k = ",k_select,", tau_select = ", stau, sep=""))
                }
                k_best_em <- k[which.min(k_like)]
                print(paste("k_select = ", k_best_em, sep=""))
                k_best_name<-paste('k_', k_best_em, sep = '')
                k_best_tau = tau_list[[k_best_name]]
                return(list(k_best=k_best_em, 
                            tau_best=k_best_tau$stau,
                            tau_list=tau_list, 
                            k_process = k_like,
                            best_list = k_best_tau))
        }else{
                
                if(length(k)==1){
                        print("Start fitting model...")
                        result = EM_func3(s=s, y_mat=y_mat, 
                                          y_mat_new=y_mat, 
                                          sigma2_eps=1, 
                                          itermax = itermax, 
                                          tol = tol, 
                                          tau=tau, 
                                          k=k)
                        return(list(k_best=k, 
                                    tau_best=tau,
                                    tau_list=list(), 
                                    k_process = result$test_like,
                                    best_list = result))
                        
                        
                }else{
                        print("Start k tuning...")
                        like_cv_mat <- matrix(0, nrow=5, ncol=length(k))
                        colnames(like_cv_mat) <- as.character(k)
                        for(mm in 1:5){
                                train_cv <- y_mat[which(fold_index!=mm),]
                                test_cv <- y_mat[which(fold_index==mm),]
                                like_cv_mat[mm,] = foreach (kkk = 1:length(k), .combine='c') %dopar% {
                                        result = EM_func3(s=s, y_mat=train_cv, 
                                                          y_mat_new=test_cv, 
                                                          sigma2_eps=1, 
                                                          itermax = itermax, 
                                                          tol = tol, 
                                                          tau=tau, 
                                                          k=k[kkk])
                                        result$test_like
                                }
                                print(paste("Tune k Fold- ", kk, " of 5 is done!!"))
                        }
                        like_cv_result = as.numeric(colMeans(like_cv_mat))
                        k_best_em <- k[which.min(like_cv_result)]
                        k_best_name<-paste('k_', k_best_em, sep = '')
                        print(paste("k_select = ", k_best_em, sep=""))
                        return(list(k_best=k_best_em, 
                                    tau_best=tau,
                                    tau_list=list(), 
                                    k_process = like_cv_result,
                                    best_list = list()))
                }
        }
}


pred_func = function(new_phi, y_mat_new, EM_result){
        
        n_train = ncol(EM_result$smooth_result$xi_smooth)
        xi_result = Kfilter_rcpp(y_mat = y_mat_new,
                                 Phi = EM_result$Phi,
                                 mu0 = EM_result$smooth_result$xi_smooth[,n_train],
                                 Cov0 = EM_result$smooth_result$Cov_filter[,,n_train],
                                 A = EM_result$A,Ca = EM_result$Ca,
                                 sigma2_eps = EM_result$sigma2_eps)
        y_spat <- t(new_phi %*% xi_result$xi_filter)
        y_temp <- t(new_phi %*% xi_result$xi_pred)
        
        Cov_filter_list = purrr::array_branch(xi_result$Cov_filter,margin = 3)
        Cov_pred_list = purrr::array_branch(xi_result$Cov_pred,margin = 3)
        
        y_spat_se =  t(purrr::map_dfc(Cov_filter_list, function(x) 
                sqrt(diag(new_phi %*% x %*% t(new_phi))+EM_result$sigma2_eps)))
        y_spat_up = as.matrix(y_spat +1.96*y_spat_se)
        y_spat_lw = as.matrix(y_spat -1.96*y_spat_se)
        
        y_temp_se =  t(purrr::map_dfc(Cov_pred_list, function(x) 
                sqrt(diag(new_phi %*% x %*% t(new_phi))+EM_result$sigma2_eps)))
        y_temp_up = as.matrix(y_temp +1.96*y_temp_se)
        y_temp_lw = as.matrix(y_temp -1.96*y_temp_se)
        
        lambda = lambda_mat_calc(A=EM_result$A, Ca=EM_result$Ca)
        lag0_cov_est = new_phi%*%lambda%*%t(new_phi)
        lag1_cov_est = new_phi%*%EM_result$A%*%lambda%*%t(new_phi)
        
        
        list(new_phi=new_phi, xi_result=xi_result,
             y_spat=y_spat, y_spat_se=y_spat_se, 
             y_spat_up=y_spat_up, y_spat_lw=y_spat_lw,
             y_temp=y_temp, y_temp_se=y_temp_se, 
             y_temp_up=y_temp_up, y_temp_lw=y_temp_lw,
             lag0_cov_est=lag0_cov_est,lag1_cov_est=lag1_cov_est)

}

performance_func = function(y_new_true,y_spat, y_temp,
                            lag0_cov_true,lag0_cov_est,
                            lag1_cov_true,lag1_cov_est){
        
        N = ncol(y_new_true)
        T = nrow(y_new_true)
        spat_acc = sum((y_new_true-y_spat)^2)/(T*N)
        temp_acc = sum((y_new_true-y_temp)^2)/(T*N)
        lag0_cov_acc = sum((lag0_cov_true-lag0_cov_est)^2)/(N*N)
        lag1_cov_acc = sum((lag1_cov_true-lag1_cov_est)^2)/(N*N)
        
        list(spat_acc=spat_acc,
             temp_acc=temp_acc,
             lag0_cov_acc=lag0_cov_acc,
             lag1_cov_acc=lag1_cov_acc)
        
}

temp_spat_sim_em3 = function(s, y_mat, y_mat_new, sigma2_eps=1, itermax = 30, 
                             tol = 0.001, tau, k,
                             s_new, y_new_true){
        EM_result = EM_func3(s, y_mat, y_mat_new, sigma2_eps=1, itermax = 30, 
                             tol = 0.001, tau, k)
        new_phi <- new_loc_phi(s_new=s_new, 
                               s_old = s, 
                               B = B_gen(s, EM_result$Phi))
        pred_result = pred_func(new_phi, y_mat_new, EM_result)
        performance = performance_func(y_new_true = signal_mat_new,
                                       y_spat = pred_result$y_spat,
                                       y_temp = pred_result$y_temp,
                                       lag0_cov_true = lag0_cov_true,
                                       lag0_cov_est = pred_result$lag0_cov_est,
                                       lag1_cov_true = lag1_cov_true,
                                       lag1_cov_est = pred_result$lag1_cov_est)
        list(EM_result = EM_result,
             pred_result=pred_result,
             performance=performance)
}



# --------------------------------------------------------
# model 1
model1_func = function(s, s_new, y_mat,y_mat_new,
                       xsi_mat, phi_mat, 
                       new_xsi_mat, phi_mat_new, lambda_mat, A){
        
        exp_cov_func = function(sigma2, alpha, h){
                sigma2*exp(-1*h/(2*alpha))
        }
        
        # mle
        cost_func = function(parameters){
                phi = parameters[1]
                sigma_eta = parameters[2]
                alpha = parameters[3]
                sigma2 = parameters[4]
                
                w_mat = y_mat[-1,]-phi*y_mat[-nrow(y_mat),]
                w_vec = as.vector(t(w_mat))
                
                block_cov_mat = exp_cov_func(sigma2=sigma_eta,
                                             alpha=alpha,
                                             fields::rdist(s))
                block_cov_mat = block_cov_mat + 
                        diag(sigma2*(1+phi^2),nrow(block_cov_mat))
                inv_block_cov_mat = solve(block_cov_mat)
                inv_data_cov_mat = Matrix::bdiag(replicate(n = (nrow(y_mat)-1),inv_block_cov_mat, simplify=F))
                neg_two_like = (nrow(y_mat)-1)*log(det(block_cov_mat)) + 
                        as.numeric(t(matrix(w_vec)) %*% inv_data_cov_mat %*% matrix(w_vec))
                neg_two_like
        }
        trytry = optim(fn=cost_func, par=c(0.2,0.1,0.1,0.1),
                       method="L-BFGS-B",
                       lower=c(-0.95,0.01,0.01,0.01),
                       upper=c(0.95,4,10,4))
        phi = trytry$par[1]
        sigma_eta = trytry$par[2]
        alpha = trytry$par[3]
        sigma2_est = trytry$par[4]
        
        # Ceta <- exp_cov_func(sigma2=trytry$par[2],
        #                      alpha=trytry$par[3],
        #                      fields::rdist(rbind(s,s_new)))
        # A_est = phi*diag(nrow=(nrow(s)+nrow(s_new)))
        # Phi_est = cbind(diag(nrow=nrow(s)),
        #                 matrix(0, nrow=nrow(s),ncol=nrow(s_new)))
        # lambda_est = solve((diag(nrow(s)+nrow(s_new))-
        #                             diag(phi^2, nrow=(nrow(s)+nrow(s_new)))),
        #                    Ceta)
        s_all = rbind(s, s_new)
        which_duplicate = which(duplicated(s_all)==T)
        s_all_clean = matrix(s_all[-which_duplicate,])
        
        Ceta <- exp_cov_func(sigma2=trytry$par[2],
                             alpha=trytry$par[3],
                             fields::rdist(s_all_clean))
        Ceta = (Ceta+t(Ceta))/2
        
        A_est = phi*diag(nrow=(nrow(s_all_clean)))
        Phi_est = cbind(diag(nrow=nrow(s)),
                        matrix(0, nrow=nrow(s),ncol=nrow(s_new)-2))
        lambda_est = solve((diag(nrow(s_all_clean))-
                                    diag(phi^2, nrow=nrow(s_all_clean))),
                           Ceta)
        
        
        # prediction ------------------------------------------
        y_all = rbind(y_mat, y_mat_new)
        xi_result = Kfilter_rcpp(y_mat = y_all,
                                 Phi = Phi_est,
                                 mu0 = rep(0,(nrow(s_all_clean))),
                                 Cov0 = lambda_est,
                                 A = A_est,Ca = Ceta,
                                 sigma2_eps = sigma2_est)
        y_spat <- t(xi_result$xi_filter[c(1,51:198,50),(nrow(y_mat)+1):nrow(y_all)])
        y_temp <- t(xi_result$xi_pred[c(1,51:198,50),(nrow(y_mat)+1):nrow(y_all)])
        
        
        
        # result
        y_new_true = new_xsi_mat %*% t(phi_mat_new)
        lag0_cov_true = phi_mat_new %*% lambda_mat %*% t(phi_mat_new)
        lag1_cov_true = phi_mat_new %*% A %*% lambda_mat %*% t(phi_mat_new)
        
        C_eta_test = exp_cov_func(sigma2=sigma_eta,
                                  alpha=alpha,
                                  fields::rdist(s_new))
        lag0_cov_est = solve(diag(nrow(s_new))-
                                        diag(phi^2, nrow=nrow(s_new)),
                                C_eta_test)
        lag1_cov_est = diag(phi,nrow=nrow(lag0_cov_est))%*%lag0_cov_est
        
        performance = performance_func(y_new_true = y_new_true,
                                       y_spat = y_spat,
                                       y_temp = y_temp,
                                       lag0_cov_true = lag0_cov_true,
                                       lag0_cov_est = lag0_cov_est,
                                       lag1_cov_true = lag1_cov_true,
                                       lag1_cov_est = lag1_cov_est)
        
        
        # check
        signal = xsi_mat %*% t(phi_mat)
        obtain_phi = function(vec){
                tryCatch({
                        arima(vec, c(1,0,0))$coef["ar1"]
                },
                
                warning = function(msg) {
                        arima(vec, c(1,0,0))$coef["ar1"]
                },
                
                error = function(msg) {
                        return(-2)
                }
                )
        }
        obtain_arma_phi = function(vec){
                tryCatch({
                        arima(vec, c(1,0,1))$coef["ar1"]
                },
                
                warning = function(msg) {
                        arima(vec, c(1,0,1))$coef["ar1"]
                },
                
                error = function(msg) {
                        return(-2)
                }
                )
                
        }
        obtain_arma_theta = function(vec){
                tryCatch({
                        arima(vec, c(1,0,1))$coef["ma1"]
                },
                
                warning = function(msg) {
                        arima(vec, c(1,0,1))$coef["ma1"]
                },
                
                error = function(msg) {
                        return(-2)
                }
                )
        }
        signal_ar1_result = apply(signal, 2, obtain_phi)
        signal_ar1_mean = mean(signal_ar1_result[which(signal_ar1_result>-2)])
        data_ar1_result = apply(y_mat, 2, obtain_phi)
        data_ar1_mean = mean(data_ar1_result[which(data_ar1_result>-2)])
        data_arma_phi_result = apply(y_mat, 2, obtain_arma_phi)
        data_arma_phi_mean = mean(data_arma_phi_result[which(data_arma_phi_result>-2)])
        data_arma_theta_result = apply(y_mat, 2, obtain_arma_theta)
        data_arma_theta_mean = mean(data_arma_theta_result[which(data_arma_theta_result>-2)])
        
        
        
        list(y_spat=y_spat, y_temp=y_temp, 
             lag0_cov_est=lag0_cov_est,
             lag1_cov_est=lag1_cov_est,
             xi_result=xi_result,
             performance=performance,
             phi=phi, sigma_eta=sigma_eta, alpha=alpha, sigma2_est=sigma2_est,
             Phi_est=Phi_est, A_est=A_est, Ceta=Ceta, 
             lambda_est=lambda_est, 
             signal_ar1_result=signal_ar1_result,
             signal_ar1_mean = signal_ar1_mean,
             data_ar1_result=data_ar1_result,
             data_ar1_mean=data_ar1_mean,
             data_arma_phi_result=data_arma_phi_result,
             data_arma_phi_mean=data_arma_phi_mean,
             data_arma_theta_result=data_arma_theta_result,
             data_arma_theta_mean=data_arma_theta_mean)
}

# --------------------------------------------------------
# model 2
auto_gau_basis_1d = function(s_mat, num_center,scale=1.5){
        
        s_max =max(s_mat)
        s_min = min(s_mat)
        
        
        if(num_center == 1){
                s_center = matrix(seq(s_min,s_max,length=(num_center+2)))
                s_center = matrix(s_center[-c(1,nrow(s_center)),])
                delta = abs(s_max - s_min)
                sigma2 = (delta/scale)^2
        }else{
                s_center = matrix(seq(s_min,s_max,length=(num_center)))
                delta = abs(s_center[1,1]-s_center[2,1])
                sigma2 = (delta/scale)^2
        }
        out = matrix(0,nrow=nrow(s_mat),ncol=num_center)
        for(i in 1:nrow(s_center)){
                s_center_mat = matrix(rep(s_center[i,],nrow(s_mat)),
                                      nrow=nrow(s_mat),byrow = T)
                out[,i]=exp(-rowSums((s_mat-s_center_mat)^2)/(2*sigma2))
        }
        #print(sigma2)
        out
}

EM_func_model2 = function(s, y_mat, y_mat_new, sigma2_eps=1, itermax = 30, 
                          tol = 0.001, k, phi_mat){
        
        yty <- t(y_mat) %*% y_mat
        T=nrow(y_mat)
        N=ncol(y_mat)
        init = phi_mat
        
        # initial 
        xi_mat = y_mat %*% init
        Cov0 = crossprod(xi_mat)/nrow(xi_mat)
        C0inv = solve(Cov0)
        Cov1 = crossprod(xi_mat[-1,], xi_mat[-nrow(xi_mat),])/(nrow(xi_mat)-1)
        Ceta = Cov0 - Cov1 %*% C0inv %*% t(Cov1)
        Ceta = as.matrix(Matrix::forceSymmetric(Ceta))
        Mest = Cov1%*%C0inv
        H = init
        muinit = matrix(0,k,1)
        # proj.grad.descent 
        tol.m = .Machine$double.eps
        tol1=1e-8
        tol2=1e-6
        
        EM_fixedPhi_rcpp(y_mat = y_mat, y_mat_new = y_mat_new,
                         Phi = init, mu0 = muinit, 
                         Cov0 = Cov0,A = Mest,Ca = Ceta,
                         sigma2_eps = 1,itermax = itermax,tol = tol)
        
        
}

temp_spat_sim_em_model2 <- function(s, s_new, y_mat, y_mat_new, 
                                    bas_mat, bas_mat_new,
                                    new_xsi_mat, phi_mat_new, lambda_mat, A,
                                    k, sigma2_eps=1, itermax = 30, tol = 0.001){
        
        
        EM_result = EM_func_model2(s, y_mat, y_mat_new, 
                                     sigma2_eps=sigma2_eps, 
                                     itermax = itermax, 
                                     tol = tol, k=k, phi_mat=bas_mat)
        
        pred_result = pred_func(new_phi=bas_mat_new, y_mat_new, EM_result)
        
        performance = performance_func(y_new_true = signal_mat_new,
                                       y_spat = pred_result$y_spat,
                                       y_temp = pred_result$y_temp,
                                       lag0_cov_true = lag0_cov_true,
                                       lag0_cov_est = pred_result$lag0_cov_est,
                                       lag1_cov_true = lag1_cov_true,
                                       lag1_cov_est = pred_result$lag1_cov_est)
        list(EM_result = EM_result,
             pred_result=pred_result,
             performance=performance)
        
}
