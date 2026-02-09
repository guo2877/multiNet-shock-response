function diag_result = sim_diagnostics(sim_result, N, asym)
% function takes the output of shock_response and tests whether it respects mass balance equation
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
% If no tolerance error is specified, the default value will be used
if nargin < 3
    asym = 1e-5;
end

% Extract output results
dP = sim_result.dP;
dS0 = sum(dP);
dR = sim_result.dR;
dC = sim_result.dC;
A = sim_result.A;
O = sim_result.O;
S0 = sim_result.S(:,1);

% Test 1: Net supply change dP + dI - dE + dOI - dOE for each country-product must equal dR + dC
dE_All = A(:,:,size(A,3)) - A(:,:,1);
dO_All = O(:,:,size(O,3)) - O(:,:,1);
dS_byc = dP + sum(dE_All, 1)' - sum(dE_All, 2) + sum(dO_All, 1)' - sum(dO_All, 2);
error = (dS_byc - sum(dR,2) - sum(dC,2)) ./ (dS_byc + S0);
error(isinf(error)) = 0;         
discrep = abs(error) > asym;
discrep_countries = find(discrep, 1);
if ~isempty(discrep_countries)
    %disp('Net supply change does not match dC + dR for each country-product');
    %diag_result = discrep_countries;
    diag_result = 1;
    return;
end

% Test 2: Net supply change dP + dI - dE for each country must equal dR + dC
for  l = 1:length(S0)/N
    if l == 1
        dE_All_c = dE_All(1:N,1:N);
        dP_c = dP(1:N,:);
        dR_c = dR(1:N,:);
        dC_c = dC(1:N,:);
        S0_c = S0(1:N);
    else
        dE_All_c = dE_All_c + dE_All(1+N*(l-1):N*l,1+N*(l-1):N*l);
        dP_c = dP_c + dP(1+N*(l-1):N*l,:);
        dR_c = dR_c + dR(1+N*(l-1):N*l,:);
        dC_c = dC_c + dC(1+N*(l-1):N*l,:);
        S0_c = S0_c + S0(1+N*(l-1):N*l);
    end
end
dS_byc_c = dP_c + sum(dE_All_c, 1)' - sum(dE_All_c, 2);
error = (dS_byc_c - sum(dR_c,2) - sum(dC_c,2)) ./ (dS_byc_c + S0_c);
error(isinf(error)) = 0;          
discrep = abs(error) > asym;
discrep_countries = find(discrep, 1);
if ~isempty(discrep_countries)
    %disp('Net supply change does not match dC + dR for each country');
    %diag_result = discrep_countries;
    diag_result = 3;
    return;
end

% Test 3: Total dR + dC must match initial shock
discrep = (sum(dR(:)) + sum(dC(:))) / dS0 - 1;
if abs(discrep) > asym
    %diag_result = sprintf('Total dC + dR does not match initial shock. Relative difference: %f', discrep);
    diag_result = 1;
    return;
end

%diag_result = 'All tests passed.';
diag_result = 4;
end