


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
