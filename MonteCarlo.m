%% Monte Carlo Simulation to price a European Put option


% Variables
r = 0.03;
del = 0.05;
a = 0.3;

T = 1;
%Stock Price at time 0
S0 = 100;
%Strike price
K1 = 120;
K2 = 80;
%Number of MC simulations
N_mc = 100000;
% Number of time steps
M = 2000;
tau = linspace(0,T,M+1);
dtau = T/M;
% dtau = T/M;

% For Antithetic Sampling method generating two samples 
%Vector to store stock at time points on path for two different samples
Sim_sam_1 = zeros(M+1,1);
Sim_sam_1(1) = S0;

Sim_sam_2 = zeros(M+1,1);
Sim_sam_2(1) = S0;

%Vector to store our final stock prices at Maturity
Y1 = zeros(N_mc,1);
Y2 = zeros(N_mc,1);

%Vector to store final discounted payoff
simulated_put_payoffs_K1 = zeros(N_mc,1);
simulated_put_payoffs_K2 = zeros(N_mc,1);

%Vector to store local volatility for both samples

sigma1 = zeros(M+1,1);
sigma2 = zeros(M+1,1);

for i=1:N_mc

    %entire path
    for m=1:M
        Z = randn; %standard normal distribution mean 0, sigma 1
        
%       1st sample 
      sigma1(m) = 0.25*exp(-tau(m)) * (100/Sim_sam_1(m))^a;
      
      Sim_sam_1(m+1) = Sim_sam_1(m) * exp((r-del-0.5*sigma1(m)^2)*dtau + sigma1(m)*sqrt(dtau)*Z);
        
%    2nd sample 
      sigma2(m) = 0.25*exp(-tau(m)) * (100/Sim_sam_2(m))^a;

      Sim_sam_2(m+1) = Sim_sam_2(m) * exp((r-del-0.5*sigma2(m)^2)*dtau + sigma2(m)*sqrt(dtau)*(-Z));
        
    end
    Z = randn; %standard normal distribution mean 0, sigma 1
    WT = sqrt(T) * Z;
    %end point only
    %Y1(i) = S0*exp((r - del - 0.5*sigma1(M)^2)*T + sigma1(M)*WT);
    %Y2(i) = S0*exp((r - del - 0.5*sigma2(M)^2)*T + sigma2(M)*(-WT));
    
    simulated_put_payoffsK1_1 = exp(-r*T)*max(K1-Sim_sam_1(M),0);
    simulated_put_payoffsK1_2 = exp(-r*T)*max(K1-Sim_sam_2(M),0);
    
    simulated_put_payoffsK2_1 = exp(-r*T)*max(K2-Sim_sam_1(M),0);
    simulated_put_payoffsK2_2 = exp(-r*T)*max(K2-Sim_sam_2(M),0);

    
    simulated_put_payoffs_K1(i) = (simulated_put_payoffsK1_1+simulated_put_payoffsK1_2)/2;
    
    simulated_put_payoffs_K2(i) = (simulated_put_payoffsK2_1+simulated_put_payoffsK2_2)/2;

    
end


fprintf("\n Monte Carlo Simulation European Put Price \n")

fprintf("\n Price at K = 120 \n")
put_mc_K1 = mean(simulated_put_payoffs_K1);
disp(put_mc_K1)

fprintf("Price at K = 80 \n")
put_mc_K2 = mean(simulated_put_payoffs_K2);
disp(put_mc_K2)

