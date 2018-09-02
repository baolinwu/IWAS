#' Transcriptome informed adaptive variant-set association test using GWAS summary data
#'
#' Take input of summary Z-statistics, variant LD matrix, transcript derived variant weight.
#' 
#' @param  Z summary Z-statistics for a set of variants
#' @param  R the estimated variant LD matrix
#' @param  W variant weight derived from outside reference panel of transcript and SNP data
#' @param  rho sequence of weights assigned to the MetaXcan (summary statistics based PrediXcan) test 
#' @return
#' \describe{
#'   \item{p.value}{ p-values for: AT (adaptive test); QT (quadratic test); MT (MetaXcan based test) }
#'   \item{pvals}{ vector of p-values for all weighted tests }
#'   \item{rho.est}{ estimated optimal \eqn{\rho} value }
#' }
#' @export
#' @references
#' Gamazon, E.R., Wheeler, H.E., Shah, K.P. et al. (2015) A gene-based association method for mapping traits using reference transcriptome data. Nat. Genet. 47, 1091–1098.
#'
#' Gusev, A., Ko, A., Shi, H., et al. (2016) Integrative approaches for large-scale transcriptome-wide association studies. Nat. Genet. 48, 245–252. 
#' 
#' Wu,B., Guan,W., Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. Annals of human genetics, 80(2), 123-135.
#'
#' Wu,B., Guo,B. and Liu,N. (2018) A powerful and efficient statistical method for transcriptome-wide association test using GWAS summary data. \emph{Bioinformatics}, under revision.
#' @examples
#' library(IWAS)
#' ## pairwise cor = 0.2
#' m = 20
#' R = cor(matrix(rnorm(100*m),100,m)*sqrt(0.8)+rnorm(100)*sqrt(0.2))
#' Z = rnorm(m)*sqrt(0.8)+rnorm(1)*sqrt(0.2); W = rnorm(m)
#' TSATS(Z,R,W)
#' Z[1:2] = Z[1:2] + 1*W[1:2]
#' TSATS(Z,R,W)
TSATS <- function(Z,R,W, rho=c((0:5/10)^2,0.5,1)){
  M = length(Z)
  ##
  Zw = Z*W; Rw = t(R*W)*W
  eR = eigen(Rw,sym=TRUE);  lamR = abs(eR$val); eta = colSums(eR$vec)*sqrt(lamR)
  R1 = sum(eta^2); R2 = sum(eta^2*lamR)
  c2 = outer(eta,eta)
  Lamq = eigen(diag(lamR) - R2/R1^2*c2, symmetric=TRUE,only.values=TRUE)$val
  ## MT: PXU  (S-PrediXcan)
  Qb = sum(Zw)^2
  pvalb = pchisq(Qb/R1, 1,lower=FALSE)
  ## QT: PXV 
  Qv = sum(Zw^2)
  pvalv = KATpval(Qv,lamR)
  ## AT: PXA
  L = length(rho)
  L1 = L-1; rho1 = rho[-L]
  Qw = (1-rho)*Qv + rho*Qb
  pval = rep(1, L)
  pval[1] = pvalv; pval[L] = pvalb
  Lamk = vector('list', L)
  Lamk[[L]] = R1;  Lamk[[1]] = lamR
  for(k in 2:L1){
    mk = rho[k]*c2;  diag(mk) = diag(mk) + (1-rho[k])*lamR
    aak = zapsmall( abs(eigen(mk, sym=TRUE, only.val=TRUE)$val) )
    Lamk[[k]] = aak[aak>0]
    pval[k] = KATpval(Qw[k],Lamk[[k]])
  }
  minP = min(pval)
  qval = rep(0,L1)
  for(k in 1:L1) qval[k] = Liu0.qval(minP, Lamk[[k]])
  Rs = rowSums(Rw); R3 = sum(Rs*colSums(Rw*Rs))
  lam = eigen(Rw-outer(Rs,Rs)/R1,sym=TRUE,only.val=TRUE)$val
  tauk = (1-rho1)*R2/R1 + rho1*R1;  vp2 = 4*(R3/R1-R2^2/R1^2)
  MuQ = sum(lam);  VarQ = sum(lam^2)*2
  sd1 = sqrt(VarQ)/sqrt(VarQ+vp2)
  q1 = qchisq(minP,1,lower=FALSE)
  katint = function(xpar){
    eta1 = sapply(xpar, function(eta0) min((qval-tauk*eta0)/(1-rho1)))
    x = (eta1-MuQ)*sd1 + MuQ
    KATpval(x,lam)*dchisq(xpar,1)
  }
  prec = 1e-4
  p.value = try({ minP + integrate(katint, 0,q1,  subdivisions=1e3,abs.tol=minP*prec)$val }, silent=TRUE)
  while(class(p.value)=='try-error'){
    p.value = try({ minP + integrate(katint, 0,q1, abs.tol=minP*prec)$val }, silent=TRUE)
    prec = prec*2
  }
  p.value = c(A=p.value, V=pvalv, U=pvalb)
  names(p.value) = c('AT', 'QT', 'MT')
  return(list(p.value=p.value, pval=pval, rho.est=rho[which.min(pval)]) )
}

