function sim_result = shock_response(P0, R0, A0, c_init, fp, fr, F0, amin, kmax, asym)
% Function to simulate shock response and cascade starting with a change (decline) in production dP of a particular product in the target country
% imput: P0, R0 and A0 are the initial production vector, reserve vector and trade matrix
%        c_init is the ID of the initial shocked country
%        fp and fr are the shock intensity and reserve release fraction
%        F0 is the substitution fraction matrix, where element is substitution fraction of one country-product (row) for another domestic product (column)
%        amin is the propagation threshold
%        kmax is the max number of iterations
%        if asym = TRUE, countries with depleted available reserves cannot increase exports
% ouput: Changes in core variables, including consumption C, reserves R, net supply S, trade A, and substitution O
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
nc = length(P0);
% Check input validity
if ~isvector(P0) || ~isvector(R0) || ~ismatrix(A0) || ~ismatrix(F0)
     error('P0 must be a vector, R0 must be a vector, F0 must be a matrix, A0 must be a matrix.');
end
if length(R0) ~= length(P0) ||any(size(A0) ~= length(P0)) || any(size(F0) ~= length(P0))
    error('Dimensions of P0, R0, F0 and A0 do not match.');
end
if ~all(sum(F0, 1) >= 0 & sum(F0, 1) <= 1)
    error('The sum of the substitution coefficients Fs for each crop must be between 0 and 1.');
end
if ~ismember(c_init, 1:nc)
    error('c_init must be one of the initial countries.');
end
if ~all(fp >= 0 & fp <= 1)
    error('fp must be between 0 and 1.');
end
if ~isequal(size(fp, 2), 1)
    error('The number of columns in the fp must be 1.');
end
if ~all(fr >= 0 & fr <= 1)
    error('fr must be between 0 and 1.');
end
if ~all(amin >= 0 & amin <= 1)
    error('amin must be between 0 and 1.');
end

% Initialize variables
S = zeros(nc, kmax);
S(:, 1) = P0 + sum(A0, 1)' - sum(A0, 2);
C = zeros(nc, kmax);
C(:, 1) = S(:, 1);
dC = zeros(nc, kmax);
R = zeros(nc, kmax);
R(:, 1) = fr .* R0;
dR = zeros(nc, kmax);
A = zeros(nc, nc, kmax);    
A(:, :, 1) = A0;
dA = zeros(nc, nc, kmax);
O = zeros(nc, nc, kmax);       
dO = zeros(nc, nc, kmax);   
Fs = zeros(nc, nc, kmax);    
Fs(:, :, 1) = F0;
blocked = false(nc, 1);

% Initial shock (dS = drop in supply)
dP = zeros(nc, 1);
dP(c_init) = -fp .* P0(c_init);
dS = dP;
shocked = zeros(nc, kmax); 
shocked(:, 1) = dS;

% Main loop
k=1; 
equilibrium = 1; %1=equilibrium��0=no-equilibrium

while ~ all(dS == 0) 
    % If the maximum number of iterations is reached, terminate the run
    if k > kmax
        disp('No equilibrium reached after maximum number of iterations.');
        equilibrium = 0;    
        break;
    end
    
    % Shock first absorbed by reserves
    dR(:, k) = max(dS, -R(:, k));
    R(:, k+1) = R(:, k) + dR(:, k);
    if asym
        blocked( R(:, k+1) == 0) = true;
    end
    
    res_shock = dS - dR(:, k);
    
    % Shocks propagated through cross-product substitution
    Ak = A(:, :, k);
    vol_forSubs = R(:, k+1) + sum(Ak, 2) + sum(Ak(~blocked, :), 1)';   %substitution capacity limit
    
    [dO(:, :, k), dOvol, Fs(:, :, k+1)] = substitution_strategy(res_shock, vol_forSubs, Fs(:, :, k), S(:,k), amin);
    
    res_shock = res_shock + dOvol;

    % Shocks propagated through trade
    [dA(:, :, k), dTvol] = trade_strategy(A(:, :, k), res_shock, S(:,k), amin, blocked);
    
    % Shock absorbed by consumption
    dC(:, k) = dS - dR(:, k) + dOvol - dTvol;  
    
    % Update values for next iteration
    % Decrease in net supply corresponds to shock absorbed internally (R, C)
    A_temp = A(:, :, k) + dA(:, :, k);
    A_temp(A_temp < 0) = 0;                     
    A(:, :, k+1) = A_temp;

    C(:, k+1) = C(:, k) + dC(:, k);
    O(:, :, k+1) = O(:, :, k) + dO(:, :, k);   
    S(:, k+1) = S(:, k) + dR(:, k) + dC(:, k); 
    
    % calculate next dS
    dS = sum(dA(:, :, k), 1)' - sum(dA(:, :, k), 2) - sum(dO(:, :, k), 2) + dTvol;
    k = k+1;
    shocked(:,k) = dS;
end
sim_result = struct('dP', dP, 'dR', dR(:,1:k-1), 'dC', dC(:,1:k-1), 'dA', dA(:, :, 1:k-1), 'dO', dO(:, :, 1:k-1),...
                'R', R(:,1:k), 'C', C(:,1:k), 'A', A(:,:,1:k), 'O', O(:,:,1:k), 'S',S(:,1:k), 'Fs', Fs(:,:,1:k),...
                'shocked', shocked(:,1:k), 'blocked', blocked, 'equilibrium', equilibrium);
end


function [dO, dOvol, F] = substitution_strategy(res_shock, vol_forSubs, F, Sk, amin)
% Function to calculate the relative change in each substitution link
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
res_shock(abs(res_shock) < amin .* Sk) = 0;   
dO = - F .* res_shock';  
D = sum(dO, 2);        
prop_subs = min(D, vol_forSubs)./D;   
prop_subs(isnan(prop_subs)) = 0; 
dO = dO.*prop_subs;
dOvol = sum(dO, 1)';   
F(dOvol > 0, :) = 0;  
end


function [dA, dTvol] = trade_strategy(Ak, res_shock, Sk, amin, blocked)
% Function to calculate the relative change in each trade link
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
nc = length(Sk);
Tvol = sum(Ak, 2) + sum(Ak(~blocked, :), 1)';  
dTvol = zeros(nc, 1);    
dTvol(abs(res_shock) >= amin .* Sk) = res_shock(abs(res_shock) >= amin .* Sk);    

dTvolE = (sum(Ak, 2) ./ Tvol) .* dTvol;
dTvolE(isnan(dTvolE)) = 0;   
dTvolI = (sum(Ak(~blocked, :), 1)' ./ Tvol) .* dTvol;
dTvolI(isnan(dTvolI)) = 0; 

dA = zeros(nc, nc);
dTvol = zeros(nc,1);

if ~all(abs(dTvolE) < 1) || ~all(abs(dTvolI) < 1 )
    [dA, dE_E, dE_I] =  dEalloc(dTvolE, dTvolI, blocked, Ak);  
    dTvol = sum(dE_E, 2) - sum(dE_I, 1)';   
end
end

function [dA, dE_E, dE_I] = dEalloc(dTvolE, dTvolI, blocked, Ak)
    n = length(blocked);
    GE = Ak;
    GI = Ak .* repmat(~blocked, 1, n);
    GE(isnan(GE)) = 0; 
    GI(isnan(GI)) = 0; 
    
    prop_E = GE ./ sum(GE,2);
    prop_E(isnan(prop_E)) = 0; 
    prop_E(dTvolE == 0,:) = 0;
    
    prop_I = GI ./ sum(GI,1);
    prop_I(isnan(prop_I)) = 0; 
    prop_I(:, dTvolI == 0) = 0;
    
    dE_E = prop_E .* dTvolE;
    dE_I = prop_I .* dTvolI';
    
    dE_E = max(dE_E, -Ak);
    dE_I = min(-dE_I, Ak);
    
    dA = dE_E + dE_I;
end