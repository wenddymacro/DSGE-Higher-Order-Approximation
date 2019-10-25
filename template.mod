% Higher order approximation:
% (1) reducing size of shocks
% (2)pruning
% (3)NLMA¡ªnonlinear maving aerage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the template is based on the model in Sims' "Using Dynare"
% @Wenddy Xu,CIMERS, xuweny87@163.com

var y I k a c w R r;
varexo e;
parameters alpha beta delta rho sigma sigmae;

alpha=1/3;
beta=0.99;
delta=0.025;
rho=0.95;
sigma=1;
sigmae=1;

model;

exp(c)^(-sigma)=beta*exp(c(1))^(-sigma)*(alpha*exp(a(1))*exp(k)^(alpha-1)+(1-delta));

exp(y) = exp(a)*exp(k(-1))^(alpha);
exp(k) = exp(a)*exp(k(-1))^(alpha)-exp(c) + (1-delta)*exp(k(-1));
a = rho*a(-1) + e;
exp(y) = exp(c) + exp(I);
exp(c)^(-sigma) = beta*exp(c(+1))^(-sigma)*(1+r);
exp(R) = alpha*exp(a)*exp(k(-1))^(alpha-1);
exp(w) = (1-alpha)*exp(a)*exp(k(-1))^(alpha);

end;

initval;
k = log(30);
y = log(3);
c = log(2.5);
I = log(0.5);
a = 0;
r = (1/beta)-1;
R = log((1/beta)-(1-delta));
w = log(1);
end;

shocks;
var e = sigmae^2;
end;

steady;

stoch_simul(order=2,irf=40,pruning);

%nlma_theoretical_moments = nlma_th_moments(M_,oo_,options_,var_list_);

%options_.irf = 40;

%nlma_irf = nlma_irf(M_,options_,var_list_);



