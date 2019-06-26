


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





