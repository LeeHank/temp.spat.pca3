#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double f_new(arma::mat phi, arma::mat Ceps_inv, arma::mat yxi,arma::mat S11, double tau, arma::mat omega) {
        double result = -2*trace(Ceps_inv*yxi*trans(phi))+trace(Ceps_inv*phi*S11*trans(phi))+tau*trace(trans(phi)*omega*phi);
        return result;
}

// [[Rcpp::export]]
double f_new2(arma::mat y_mat, arma::mat phi, arma::mat Ceps_inv, arma::mat yxi,arma::mat S11, double tau, arma::mat omega) {
        int T=y_mat.n_rows;
        int N=y_mat.n_cols;
        double tr = arma::trace(trans(y_mat)*y_mat-yxi*trans(phi)-phi*trans(yxi)+phi*S11*trans(phi));
        double result = T*N*log(tr/(T*N))+tau*trace(trans(phi)*omega*phi);
        return result;
}

// [[Rcpp::export]]
arma::mat f_grad(arma::mat phi,
                 arma::mat Ceps_inv,
                 arma::mat yxi,
                 arma::mat S11,
                 double tau,
                 arma::mat omega){
        return -2*Ceps_inv*yxi+Ceps_inv*phi*(S11+trans(S11))+2*tau*omega*phi;
}

// [[Rcpp::export]]
arma::mat f_grad2(arma::mat y_mat,
                  arma::mat phi,
                  arma::mat Ceps_inv,
                  arma::mat yxi,
                  arma::mat S11,
                  double tau,
                  arma::mat omega){
        int T=y_mat.n_rows;
        int N=y_mat.n_cols;
        double tr = arma::trace(trans(y_mat)*y_mat-yxi*trans(phi)-phi*trans(yxi)+phi*S11*trans(phi));
        return (T*N/tr)*(-2*yxi+2*phi*S11)+2*tau*omega*phi;
}

// [[Rcpp::export]]
arma::mat line_search(arma::mat init, arma::mat z2, double learning_rate){
        arma::mat z3 = init + learning_rate * z2;
        arma::mat U;
        arma::vec d;
        arma::mat V;
        svd_econ(U, d, V, z3);
        arma::mat m_next = U*trans(V);
        return m_next;
}


// [[Rcpp::export]]
arma::mat find_phi(arma::mat init, arma::mat Ceps_inv, arma::mat yxi, arma::mat S11, double tau, arma::mat omega,
                   double tol, double tol1, double tol2, bool verbose){
        
        arma::mat phi_est = init;
        int N = init.n_rows;
        arma::mat I_N = arma::eye<arma::mat>(N, N);
        double current;
        for(int i=0; i<1000; ++i){
                arma::mat z1 = -1*f_grad(phi_est, Ceps_inv, yxi, S11, tau, omega);
                arma::mat z2 = 0.5*phi_est*(trans(phi_est)*z1 - trans(z1)*phi_est)+
                        (I_N - (phi_est*trans(phi_est)))*z1;
                double x = 2;
                current = f_new(phi_est, Ceps_inv, yxi, S11, tau, omega);
                arma::mat m_res;
                for(int ll=0; ll<50;++ll){
                        m_res = line_search(phi_est, z2, x);
                        double res = f_new(m_res, Ceps_inv, yxi, S11, tau, omega);
                        if(res>current){
                                x = 1/(pow(2,ll));
                        }else break;
                }
                arma::mat m_next = line_search(phi_est, z2, x);
                
                if(verbose){
                        printf("** i = %d \n", i);
                        printf("** current loss = %f \n", current);
                }
                
                arma::vec state = vectorise(arma::abs(phi_est-m_next)/(arma::abs(phi_est)+tol2));
                if(arma::all( state < tol1)){break;}
                
                phi_est = m_next;
                
        }
        return phi_est;
}

// [[Rcpp::export]]
arma::mat find_phi2(arma::mat y_mat, arma::mat init, arma::mat Ceps_inv, arma::mat yxi, arma::mat S11, double tau, arma::mat omega,
                    double tol, double tol1, double tol2, bool verbose){
        
        arma::mat phi_est = init;
        int N = init.n_rows;
        arma::mat I_N = arma::eye<arma::mat>(N, N);
        //double current, update_loss;
        for(int i=0; i<200; ++i){
                arma::mat z1 = -1*f_grad2(y_mat, phi_est, Ceps_inv, yxi, S11, tau, omega);
                arma::mat z2 = 0.5*phi_est*(trans(phi_est)*z1 - trans(z1)*phi_est)+
                        (I_N - (phi_est*trans(phi_est)))*z1;
                double x = 2;
                double current = f_new2(y_mat, phi_est, Ceps_inv, yxi, S11, tau, omega);
                arma::mat m_res;
                for(int ll=0; ll<50;++ll){
                        m_res = line_search(phi_est, z2, x);
                        double res = f_new2(y_mat, m_res, Ceps_inv, yxi, S11, tau, omega);
                        if(res>current){
                                x = 1/(pow(2,ll));
                        }else break;
                }
                arma::mat m_next = line_search(phi_est, z2, x);
                double update_loss = f_new2(y_mat, m_next, Ceps_inv, yxi, S11, tau, omega);
                if(verbose){
                        printf("** i = %d \n", i);
                        printf("** current loss = %f \n", current);
                        printf("** update loss = %f \n", update_loss);
                }
                
                //arma::vec state = vectorise(arma::abs(phi_est-m_next)/(arma::abs(phi_est)+tol2));
                //if(arma::all( state < tol1)){break;}
                if(fabs(update_loss - current) < 0.001 ){
                        //printf("diff = %f \n", fabs(update_loss - current));
                        break;
                };
                
                
                phi_est = m_next;
                
        }
        return phi_est;
}


using namespace arma;
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
Rcpp::List Kfilter_rcpp(arma::mat y_mat, arma::mat Phi, 
                        arma::vec mu0,arma::mat Cov0, 
                        arma::mat A, arma::mat Ca, double sigma2_eps){
        
        //declare
        int dim_K = A.n_cols;
        int N = y_mat.n_cols;
        int T = y_mat.n_rows;
        
        arma::cube Cov_pred, Cov_filter, sig_t;
        arma::mat xi_pred, xi_filter, Ceps, innov,  siginv, K; 
        arma::mat I_N = arma::eye(N,N);
        arma::mat I_K = arma::eye(dim_K,dim_K);
        double like=0;
        
        //initialize
        xi_pred.zeros(dim_K,T);
        xi_filter.zeros(dim_K,T);
        innov.zeros(N,T);
        Cov_pred.zeros(dim_K,dim_K,T);
        Cov_filter.zeros(dim_K,dim_K,T);
        sig_t.zeros(N,N,T);
        siginv.zeros(N,N);
        K.zeros(dim_K,N);
        Ceps = sigma2_eps*I_N;
        
        //i=0
        xi_pred.col(0) = A * mu0;
        Cov_pred.slice(0) = A * Cov0 * trans(A) + Ca;
        sig_t.slice(0) = Phi * Cov_pred.slice(0) * trans(Phi) + Ceps;
        sig_t.slice(0) = (sig_t.slice(0)+trans(sig_t.slice(0)))/2;
        siginv = inv(sig_t.slice(0));
        K = Cov_pred.slice(0) * trans(Phi) * siginv;
        innov.col(0) = trans(y_mat.row(0)) - Phi * xi_pred.col(0);
        
        xi_filter.col(0) = xi_pred.col(0) + K * innov.col(0);
        Cov_filter.slice(0) =  Cov_pred.slice(0) - K * Phi * Cov_pred.slice(0);
        //sigmat = as.matrix(sig[, , 1], nrow = qdim, ncol = qdim)
        like = log(det(sig_t.slice(0))) + trace(trans(innov.col(0)) * siginv * innov.col(0));
        
        //filtering
        //printf("...Filtering...\n");
        for (int i=1; i<T; i++) {
                xi_pred.col(i) = A * xi_filter.col(i-1);
                Cov_pred.slice(i) = A * Cov_filter.slice(i-1) * trans(A) + Ca;
                sig_t.slice(i) = Phi * Cov_pred.slice(i) * trans(Phi) + Ceps;
                siginv = inv(sig_t.slice(i));
                K = Cov_pred.slice(i) * trans(Phi) * siginv;
                innov.col(i) = trans(y_mat.row(i)) - Phi * xi_pred.col(i);
                xi_filter.col(i) = xi_pred.col(i) + K * innov.col(i);
                Cov_filter.slice(i) =  Cov_pred.slice(i) - K * Phi * Cov_pred.slice(i);
                
                //sigmat = as.matrix(sig[, , i], nrow = qdim, ncol = qdim)
                like = like + log(det(sig_t.slice(i))) + 
                        trace(trans(innov.col(i)) * siginv * innov.col(i));
        }
        return Rcpp::List::create(Named("xi_pred")=xi_pred,
                                  Named("Cov_pred")=Cov_pred,
                                  Named("xi_filter")=xi_filter,
                                  Named("Cov_filter")=Cov_filter,
                                  Named("like")=like,
                                  Named("innov")=innov,
                                  Named("sig")=sig_t,
                                  Named("K_T")=K);
}

using namespace arma;
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List Ksmooth_rcpp(arma::mat y_mat, arma::mat Phi, 
                        arma::vec mu0,arma::mat Cov0, 
                        arma::mat A, arma::mat Ca, double sigma2_eps){
        
        //declare
        Rcpp::List kf;
        int dim_K = A.n_cols;
        int N = y_mat.n_cols;
        int T = y_mat.n_rows;
        
        arma::mat xi_smooth, J0, Cov_0_smooth;
        arma::cube Cov_smooth, J;
        arma::vec xi_0_smooth;
        
        //initialized
        xi_smooth.zeros(dim_K,T);
        Cov_smooth.zeros(dim_K,dim_K,T);
        J.zeros(dim_K,dim_K,T);
        
        kf = Kfilter_rcpp(y_mat, Phi, mu0, Cov0,A,Ca,sigma2_eps);
        arma::mat xi_pred=kf["xi_pred"];
        arma::cube Cov_pred=kf["Cov_pred"];
        arma::mat xi_filter=kf["xi_filter"];
        arma::cube Cov_filter=kf["Cov_filter"];
        arma::mat K_T = kf["K_T"];
        
        xi_smooth.col(T-1) = xi_filter.col(T-1);
        Cov_smooth.slice(T-1) = Cov_filter.slice(T-1);
        //printf("...Smoothing...\n");
        for (int i=(T-2); i>=0;i--) {
                J.slice(i) = Cov_filter.slice(i) * trans(A) * inv(Cov_pred.slice(i+1));
                xi_smooth.col(i) = (xi_filter.col(i) + J.slice(i) * (xi_smooth.col(i+1) - xi_pred.col(i+1)));
                Cov_smooth.slice(i) = Cov_filter.slice(i) + 
                        J.slice(i) *(Cov_smooth.slice(i+1) - Cov_pred.slice(i+1))*trans(J.slice(i));
                
        }
        J0 = Cov0 * trans(A) * inv(Cov_pred.slice(0));
        xi_0_smooth = mu0 + J0 * (xi_smooth.col(0)-xi_pred.col(0));
        Cov_0_smooth = Cov0 + J0*(Cov_smooth.slice(0) - Cov_pred.slice(0))*trans(J0);
        
        return Rcpp::List::create(Named("xi_smooth")=xi_smooth,
                                  Named("Cov_smooth")=Cov_smooth,
                                  Named("xi_0_smooth")=xi_0_smooth,
                                  Named("Cov_0_smooth")=Cov_0_smooth,
                                  Named("J0")=J0,
                                  Named("J")=J,
                                  Named("xi_pred")=xi_pred,
                                  Named("Cov_pred")=Cov_pred,
                                  Named("xi_filter")=xi_filter,
                                  Named("Cov_filter")=Cov_filter,
                                  Named("like")=kf["like"],
                                                  Named("innov")=kf["innov"],
                                                                   Named("sig")=kf["sig"],
                                                                                  Named("K_T")=K_T);
}

using namespace arma;
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List EM_fixedPhi_rcpp(arma::mat y_mat, arma::mat y_mat_new, arma::mat Phi, 
                            arma::vec mu0,arma::mat Cov0, 
                            arma::mat A, arma::mat Ca, double sigma2_eps,
                            int itermax, double tol){
        
        
        
        //double tau, double tol_m, double tol1, double tol2
        
        
        //declare
        Rcpp::List ks, y_new_filter;
        int dim_K = A.n_cols;
        int N = y_mat.n_cols;
        int T = y_mat.n_rows;
        
        double test_like,current_diff;
        arma::mat Ceps_nondiag, S00, S11, S10, R;
        arma::vec u;
        arma::mat I_N = arma::eye(N,N);
        arma::mat I_K = arma::eye(dim_K,dim_K);
        arma::cube Cov_cross;
        //arma::vec twonegloglik=zeros<vec>(itermax);
        bool endloop = false;
        
        Rcpp::NumericVector twonegloglik(0);
        
        
        //initialized 
        Cov_cross.zeros(dim_K,dim_K,T);
        
        
        for (int iter=0; iter<(itermax); iter++){
                // Filtering & Smoothing
                ks = Ksmooth_rcpp(y_mat, Phi, mu0, Cov0,A,Ca,sigma2_eps);
                arma::mat xi_pred=ks["xi_pred"];
                arma::cube Cov_pred=ks["Cov_pred"];
                arma::mat xi_filter=ks["xi_filter"];
                arma::cube Cov_filter=ks["Cov_filter"];
                arma::mat xi_smooth=ks["xi_smooth"];
                arma::cube Cov_smooth=ks["Cov_smooth"];
                arma::mat K_T = ks["K_T"];
                arma::cube J = ks["J"];
                arma::mat J0 = ks["J0"];
                arma::mat Cov_0_smooth = ks["Cov_0_smooth"];
                arma::vec xi_0_smooth = ks["xi_0_smooth"];
                double current_like = ks["like"];
                
                //check convergence
                twonegloglik.push_back(current_like);
                //twonegloglik(iter)=current_like;
                printf("iter = %d twonegloglike = %f1 \n", iter, current_like);
                if (iter == (itermax-1)) {
                        printf("Maximum iterations reached\n");
                        endloop = true;
                }
                if (iter>0){
                        current_diff = (twonegloglik(iter-1)-twonegloglik(iter))/twonegloglik(iter-1);
                        if(current_diff<0){current_diff = -current_diff;};
                        if(current_diff < tol){
                                printf("Tolerance reached\n");
                                endloop = true;
                        }
                }
                if (endloop) {
                        
                        y_new_filter = Kfilter_rcpp(y_mat_new, 
                                                    Phi, 
                                                    xi_0_smooth,
                                                    Cov_0_smooth,
                                                    A, 
                                                    Ca, 
                                                    sigma2_eps);
                        test_like = y_new_filter["like"];
                        return Rcpp::List::create(Named("Phi")=Phi,
                                                  Named("A")=A,
                                                  Named("Ca")=Ca,
                                                  Named("sigma2_eps")=sigma2_eps,
                                                  Named("Ceps_nondiag")=Ceps_nondiag,
                                                  Named("smooth_result")=ks,
                                                  Named("twonegloglik")=twonegloglik,
                                                  Named("test_like")=test_like,
                                                  Named("tol")=tol);        
                        
                }
                
                //printf("...Computing cross-covariances...\n");
                Cov_cross.slice(T-1) <- (I_K - K_T*Phi) *A*Cov_filter.slice(T-2);
                for (int i=(T-1); i>1; i--) {
                        Cov_cross.slice(i-1) = Cov_filter.slice(i-1) * trans(J.slice(i-2)) + 
                                J.slice(i-1) * (Cov_cross.slice(i) - A * Cov_filter.slice(i-1)) * 
                                trans(J.slice(i-2));
                }
                Cov_cross.slice(0) = Cov_filter.slice(0) * trans(J0) + 
                        J.slice(0) * (Cov_cross.slice(1) - A * Cov_filter.slice(0))*
                        trans(J0);
                S11 = sum(Cov_smooth.slices(0,(T-1)), 2);
                S11 = S11 + xi_smooth*trans(xi_smooth);
                S00 = sum(Cov_smooth.slices(0,(T-2)), 2);
                S00 = S00 + xi_smooth.cols(0,(T-2))*trans(xi_smooth.cols(0,(T-2)));
                S00 = S00 + Cov_0_smooth+xi_0_smooth*trans(xi_0_smooth);
                S10 = sum(Cov_cross.slices(0,(T-1)),2);
                S10 = S10 + xi_smooth.col(0)* trans(xi_0_smooth)+
                        xi_smooth.cols(1,(T-1))*trans(xi_smooth.cols(0,(T-2)));
                
                //S11 = xi_smooth.col(0) * trans(xi_smooth.col(0)) + Cov_smooth.slice(0);
                //S10 = xi_smooth.col(0) * trans(xi_0_smooth) + Cov_cross.slice(0);
                //S00 = xi_0_smooth * trans(xi_0_smooth) + Cov_0_smooth;
                //u = trans(y_mat.row(0)) - Phi * xi_smooth.col(0);
                //R = u * trans(u) + Phi * Cov_smooth.slice(0)*trans(Phi);
                //for(int i=1; i<T; i++){
                //        S11 = S11 + xi_smooth.col(i)*trans(xi_smooth.col(i))+Cov_smooth.slice(i);
                //        S10 = S10 + xi_smooth.col(i)*trans(xi_smooth.col(i-1))+Cov_cross.slice(i);
                //        S00 = S00 + xi_smooth.col(i-1)*trans(xi_smooth.col(i-1))+Cov_smooth.slice(i-1);
                //        u = trans(y_mat.row(i)) - Phi * xi_smooth.col(i);
                //        R = R + u * trans(u) + Phi * Cov_smooth.slice(i)*trans(Phi);
                //}
                
                //printf("...M-step...\n");
                A = S10 * inv(S00);
                Ca = (S11 - A*trans(S10))/T;
                Ca = (Ca+trans(Ca))/2;
                R = trans(y_mat)*y_mat-trans(y_mat)*trans(xi_smooth)*trans(Phi)-Phi*xi_smooth*y_mat+Phi*S11*trans(Phi);
                Ceps_nondiag = R/T;
                sigma2_eps = trace(R)/(T*N);
                mu0 = xi_0_smooth;
                Cov0 = Cov_0_smooth;
        }
}


using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List EM_nonfixedPhi_rcpp(arma::mat y_mat, arma::mat y_mat_new, arma::mat Phi, 
                               arma::vec mu0,arma::mat Cov0, 
                               arma::mat A, arma::mat Ca, double sigma2_eps,
                               int itermax, double tol,
                               double tau, arma::mat omega, double tol_m, double tol1, double tol2){
        
        //declare
        Rcpp::List ks, y_new_filter;
        int dim_K = A.n_cols;
        int N = y_mat.n_cols;
        int T = y_mat.n_rows;
        
        double test_like,current_diff;
        arma::mat Ceps_nondiag, S00, S11, S10, R, Qeps, yxi;
        arma::vec u;
        arma::mat I_N = arma::eye(N,N);
        arma::mat I_K = arma::eye(dim_K,dim_K);
        arma::cube Cov_cross;
        //arma::vec twonegloglik=zeros<vec>(itermax);
        bool endloop = false;
        
        Rcpp::NumericVector twonegloglik(0);
        
        
        //initialized 
        Cov_cross.zeros(dim_K,dim_K,T);
        
        
        for (int iter=0; iter<(itermax); iter++){
                // Filtering & Smoothing
                ks = Ksmooth_rcpp(y_mat, Phi, mu0, Cov0,A,Ca,sigma2_eps);
                arma::mat xi_pred=ks["xi_pred"];
                arma::cube Cov_pred=ks["Cov_pred"];
                arma::mat xi_filter=ks["xi_filter"];
                arma::cube Cov_filter=ks["Cov_filter"];
                arma::mat xi_smooth=ks["xi_smooth"];
                arma::cube Cov_smooth=ks["Cov_smooth"];
                arma::mat K_T = ks["K_T"];
                arma::cube J = ks["J"];
                arma::mat J0 = ks["J0"];
                arma::mat Cov_0_smooth = ks["Cov_0_smooth"];
                arma::vec xi_0_smooth = ks["xi_0_smooth"];
                double current_like = (double)ks["like"];
                current_like = current_like + tau*trace(trans(Phi)*omega*Phi) ;
                
                //check convergence
                twonegloglik.push_back(current_like);
                //twonegloglik(iter)=current_like;
                printf("iter = %d twonegloglike = %f1 \n", iter, current_like);
                if (iter == (itermax-1)) {
                        printf("Maximum iterations reached\n");
                        endloop = true;
                }
                if (iter>0){
                        current_diff = (twonegloglik(iter-1)-twonegloglik(iter))/abs(twonegloglik(iter-1));
                        if(current_diff<0){
                                current_diff = -current_diff;
                                endloop = true;
                        };
                        if(current_diff < tol){
                                printf("Tolerance reached\n");
                                endloop = true;
                        }
                }
                if (endloop) {
                        
                        y_new_filter = Kfilter_rcpp(y_mat_new, 
                                                    Phi, 
                                                    xi_0_smooth,
                                                    Cov_0_smooth,
                                                    A, 
                                                    Ca, 
                                                    sigma2_eps);
                        test_like = y_new_filter["like"];
                        return Rcpp::List::create(Named("Phi")=Phi,
                                                  Named("A")=A,
                                                  Named("Ca")=Ca,
                                                  Named("sigma2_eps")=sigma2_eps,
                                                  Named("Ceps_nondiag")=Ceps_nondiag,
                                                  Named("smooth_result")=ks,
                                                  Named("twonegloglik")=twonegloglik,
                                                  Named("test_like")=test_like,
                                                  Named("tol")=tol);        
                        
                }
                
                //printf("...Computing cross-covariances...\n");
                Cov_cross.slice(T-1) <- (I_K - K_T*Phi) *A*Cov_filter.slice(T-2);
                for (int i=(T-1); i>1; i--) {
                        Cov_cross.slice(i-1) = Cov_filter.slice(i-1) * trans(J.slice(i-2)) + 
                                J.slice(i-1) * (Cov_cross.slice(i) - A * Cov_filter.slice(i-1)) * 
                                trans(J.slice(i-2));
                }
                Cov_cross.slice(0) = Cov_filter.slice(0) * trans(J0) + 
                        J.slice(0) * (Cov_cross.slice(1) - A * Cov_filter.slice(0))*
                        trans(J0);
                S11 = sum(Cov_smooth.slices(0,(T-1)), 2);
                S11 = S11 + xi_smooth*trans(xi_smooth);
                S00 = sum(Cov_smooth.slices(0,(T-2)), 2);
                S00 = S00 + xi_smooth.cols(0,(T-2))*trans(xi_smooth.cols(0,(T-2)));
                S00 = S00 + Cov_0_smooth+xi_0_smooth*trans(xi_0_smooth);
                S10 = sum(Cov_cross.slices(0,(T-1)),2);
                S10 = S10 + xi_smooth.col(0)* trans(xi_0_smooth)+
                        xi_smooth.cols(1,(T-1))*trans(xi_smooth.cols(0,(T-2)));
                
                //S11 = xi_smooth.col(0) * trans(xi_smooth.col(0)) + Cov_smooth.slice(0);
                //S10 = xi_smooth.col(0) * trans(xi_0_smooth) + Cov_cross.slice(0);
                //S00 = xi_0_smooth * trans(xi_0_smooth) + Cov_0_smooth;
                //u = trans(y_mat.row(0)) - Phi * xi_smooth.col(0);
                //R = u * trans(u) + Phi * Cov_smooth.slice(0)*trans(Phi);
                //for(int i=1; i<T; i++){
                //      S11 = S11 + xi_smooth.col(i)*trans(xi_smooth.col(i))+Cov_smooth.slice(i);
                //        S10 = S10 + xi_smooth.col(i)*trans(xi_smooth.col(i-1))+Cov_cross.slice(i);
                //        S00 = S00 + xi_smooth.col(i-1)*trans(xi_smooth.col(i-1))+Cov_smooth.slice(i-1);
                //        u = trans(y_mat.row(i)) - Phi * xi_smooth.col(i);
                //        R = R + u * trans(u) + Phi * Cov_smooth.slice(i)*trans(Phi);
                //}
                
                //printf("...M-step...\n");
                A = S10 * inv(S00);
                Ca = (S11 - A*trans(S10))/T;
                Ca = (Ca+trans(Ca))/2;
                R = trans(y_mat)*y_mat-trans(y_mat)*trans(xi_smooth)*trans(Phi)-Phi*xi_smooth*y_mat+Phi*S11*trans(Phi);
                Ceps_nondiag = R/T;
                sigma2_eps = trace(R)/(T*N);
                Qeps = (1/sigma2_eps)*I_N;
                mu0 = xi_0_smooth;
                Cov0 = Cov_0_smooth;
                yxi = trans(y_mat)*trans(xi_smooth);
                Phi = find_phi(Phi, Qeps, yxi, S11, tau, omega, tol_m, tol1, tol2, false);
        }
}

using namespace arma;
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List EM_nonfixedPhi_rcpp2(arma::mat y_mat, arma::mat y_mat_new, arma::mat Phi, 
                                arma::vec mu0,arma::mat Cov0, 
                                arma::mat A, arma::mat Ca, double sigma2_eps,
                                int itermax, double tol,
                                double tau, arma::mat omega, double tol_m, double tol1, double tol2){
        
        //declare
        Rcpp::List ks, y_new_filter;
        int dim_K = A.n_cols;
        int N = y_mat.n_cols;
        int T = y_mat.n_rows;
        
        double test_like,current_diff;
        arma::mat Ceps_nondiag, S00, S11, S10, R, Qeps, yxi;
        arma::vec u;
        arma::mat I_N = arma::eye(N,N);
        arma::mat I_K = arma::eye(dim_K,dim_K);
        arma::cube Cov_cross;
        //arma::vec twonegloglik=zeros<vec>(itermax);
        bool endloop = false;
        
        Rcpp::NumericVector twonegloglik(0);
        
        
        //initialized 
        Cov_cross.zeros(dim_K,dim_K,T);
        
        
        for (int iter=0; iter<(itermax); iter++){
                // Filtering & Smoothing
                ks = Ksmooth_rcpp(y_mat, Phi, mu0, Cov0,A,Ca,sigma2_eps);
                arma::mat xi_pred=ks["xi_pred"];
                arma::cube Cov_pred=ks["Cov_pred"];
                arma::mat xi_filter=ks["xi_filter"];
                arma::cube Cov_filter=ks["Cov_filter"];
                arma::mat xi_smooth=ks["xi_smooth"];
                arma::cube Cov_smooth=ks["Cov_smooth"];
                arma::mat K_T = ks["K_T"];
                arma::cube J = ks["J"];
                arma::mat J0 = ks["J0"];
                arma::mat Cov_0_smooth = ks["Cov_0_smooth"];
                arma::vec xi_0_smooth = ks["xi_0_smooth"];
                double current_like = (double)ks["like"];
                current_like = current_like + tau*trace(trans(Phi)*omega*Phi) ;
                
                //check convergence
                twonegloglik.push_back(current_like);
                //twonegloglik(iter)=current_like;
                printf("iter = %d twonegloglike = %f1 \n", iter, current_like);
                if (iter == (itermax-1)) {
                        printf("Maximum iterations reached\n");
                        endloop = true;
                }
                if (iter>0){
                        current_diff = (twonegloglik(iter-1)-twonegloglik(iter))/abs(twonegloglik(iter-1));
                        if(current_diff<0){
                                current_diff = -current_diff;
                                endloop = true;
                        };
                        if(current_diff < tol){
                                printf("Tolerance reached\n");
                                endloop = true;
                        }
                }
                if (endloop) {
                        
                        y_new_filter = Kfilter_rcpp(y_mat_new, 
                                                    Phi, 
                                                    xi_0_smooth,
                                                    Cov_0_smooth,
                                                    A, 
                                                    Ca, 
                                                    sigma2_eps);
                        test_like = y_new_filter["like"];
                        return Rcpp::List::create(Named("Phi")=Phi,
                                                  Named("A")=A,
                                                  Named("Ca")=Ca,
                                                  Named("sigma2_eps")=sigma2_eps,
                                                  Named("Ceps_nondiag")=Ceps_nondiag,
                                                  Named("smooth_result")=ks,
                                                  Named("twonegloglik")=twonegloglik,
                                                  Named("test_like")=test_like,
                                                  Named("tol")=tol);        
                        
                }
                
                //printf("...Computing cross-covariances...\n");
                Cov_cross.slice(T-1) <- (I_K - K_T*Phi) *A*Cov_filter.slice(T-2);
                for (int i=(T-1); i>1; i--) {
                        Cov_cross.slice(i-1) = Cov_filter.slice(i-1) * trans(J.slice(i-2)) + 
                                J.slice(i-1) * (Cov_cross.slice(i) - A * Cov_filter.slice(i-1)) * 
                                trans(J.slice(i-2));
                }
                Cov_cross.slice(0) = Cov_filter.slice(0) * trans(J0) + 
                        J.slice(0) * (Cov_cross.slice(1) - A * Cov_filter.slice(0))*
                        trans(J0);
                S11 = sum(Cov_smooth.slices(0,(T-1)), 2);
                S11 = S11 + xi_smooth*trans(xi_smooth);
                S00 = sum(Cov_smooth.slices(0,(T-2)), 2);
                S00 = S00 + xi_smooth.cols(0,(T-2))*trans(xi_smooth.cols(0,(T-2)));
                S00 = S00 + Cov_0_smooth+xi_0_smooth*trans(xi_0_smooth);
                S10 = sum(Cov_cross.slices(0,(T-1)),2);
                S10 = S10 + xi_smooth.col(0)* trans(xi_0_smooth)+
                        xi_smooth.cols(1,(T-1))*trans(xi_smooth.cols(0,(T-2)));
                
                //S11 = xi_smooth.col(0) * trans(xi_smooth.col(0)) + Cov_smooth.slice(0);
                //S10 = xi_smooth.col(0) * trans(xi_0_smooth) + Cov_cross.slice(0);
                //S00 = xi_0_smooth * trans(xi_0_smooth) + Cov_0_smooth;
                //u = trans(y_mat.row(0)) - Phi * xi_smooth.col(0);
                //R = u * trans(u) + Phi * Cov_smooth.slice(0)*trans(Phi);
                //for(int i=1; i<T; i++){
                //      S11 = S11 + xi_smooth.col(i)*trans(xi_smooth.col(i))+Cov_smooth.slice(i);
                //        S10 = S10 + xi_smooth.col(i)*trans(xi_smooth.col(i-1))+Cov_cross.slice(i);
                //        S00 = S00 + xi_smooth.col(i-1)*trans(xi_smooth.col(i-1))+Cov_smooth.slice(i-1);
                //        u = trans(y_mat.row(i)) - Phi * xi_smooth.col(i);
                //        R = R + u * trans(u) + Phi * Cov_smooth.slice(i)*trans(Phi);
                //}
                
                //printf("...M-step...\n");
                A = S10 * inv(S00);
                Ca = (S11 - A*trans(S10))/T;
                Ca = (Ca+trans(Ca))/2;
                //R = trans(y_mat)*y_mat-trans(y_mat)*trans(xi_smooth)*trans(Phi)-Phi*xi_smooth*y_mat+Phi*S11*trans(Phi);
                //Ceps_nondiag = R/T;
                //sigma2_eps = trace(R)/(T*N);
                Qeps = (1/sigma2_eps)*I_N;
                mu0 = xi_0_smooth;
                Cov0 = Cov_0_smooth;
                yxi = trans(y_mat)*trans(xi_smooth);
                //Phi = find_phi2(Phi, Qeps, yxi, S11, tau, omega, tol_m, tol1, tol2, true);
                
                // 先解Phi再解sigma^2
                Phi = find_phi2(y_mat, Phi, Qeps, yxi, S11, tau, omega, tol_m, tol1, tol2, false);
                R = trans(y_mat)*y_mat-trans(y_mat)*trans(xi_smooth)*trans(Phi)-Phi*xi_smooth*y_mat+Phi*S11*trans(Phi);
                Ceps_nondiag = R/T;
                sigma2_eps = trace(R)/(T*N);
                Qeps = (1/sigma2_eps)*I_N;
                
        }
}
