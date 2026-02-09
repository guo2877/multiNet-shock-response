%The script validates mass balance equations and calculates shock impacts (consumption deficits)

%% Set parameters
Products_Names ={'Rice','Wheat'};  %analyze products
year = 2017; %analyze data year
year_id = year - 1993 + 1;  

shocked_product = ['Rice'];  %shocked product name

% Input parameter value interval
fps = 1.0;   
frs = 0.5;           
fss = [0.2];          
amin = 0.00001;     
kmax = 100;          
asym = true;                   

for i = 1:size(Products_Names, 2)
    if i == 1
        net_name = [Products_Names{i}];   %multilayer network name
    else
        net_name = [net_name,'_',Products_Names{i}];
    end
end
data_name = [net_name,'_Data'];   %multilayer network data name

%% Load data
load(['./input/',data_name,'.mat'],data_name);
eval(['CNS = {',data_name ,'.CName};']) %country name
eval(['FS = {',data_name ,'.scMatrix};']) %initial substitution fraction matrix

CNAME = CNS{year_id};  F0 = FS{year_id};

M = size(Products_Names, 2);  %number of product types 
N = length(CNAME)/M;   %number of countries

%% Shock protocol
for i = 1:size(Products_Names, 2)
   if strcmp(shocked_product, [Products_Names{i}])
       country_id = [1:171] + 171*(i-1);  %exhaustive shock protocol: shocks are applied separately to each country in the shocked layer
       %country_id = 162 + 171*(i-1);       %single shock protocol: shock to one country£¨cid, 162 = ukr) in the shocked layer
   end 
end

%% Validate and calculate on two-layer network 
for p = 1:length(fps) 
    for r = 1:length(frs)
        for s = 1:length(fss)
            Parameters_Name = ['fp',num2str(fps(p)*100),'_fr',num2str(frs(r)*100),'_fs',num2str(fss(s)*100)];
            sim_results_Name=['sim_results_',Parameters_Name];
            diag_results_Name=['diag_results_',Parameters_Name];
            impact_results_Name=['impact_results_',Parameters_Name];
            load(['./output/shock_results/',net_name,'/',shocked_product,'/',sim_results_Name,'.mat'],sim_results_Name);
            eval(['sim_results = ', sim_results_Name,';']);
            
            diag_results = diag_analysis(sim_results, country_id, CNAME, N);  
            eval([diag_results_Name,'= diag_results;']);
            clear diag_results
            
            impact_results = impact_analysis(sim_results, country_id, CNAME, N, M); 
            eval([impact_results_Name,'= impact_results;']);
            clear impact_results
            
            save(['./output/analysis_results/',net_name,'/',shocked_product,'/',diag_results_Name,'.mat'],diag_results_Name);
            save(['./output/analysis_results/',net_name,'/',shocked_product,'/',impact_results_Name,'.mat'],impact_results_Name);
            disp(['Calculation complete: two-layer_analysis_results_',Parameters_Name]);
            eval(['clear ',diag_results_Name,';']); 
            eval(['clear ',impact_results_Name,';']); 
            clear sim_results
        end
    end
end

%% Validate and calculate on three-layer network with chain structure 
% for p = 1:length(fps) 
%     for r = 1:length(frs)
%         for s = 1:length(fss)
%             fs = fss(s);       
%             fs1 = fs;        
%             fss2 = 0.1:0.1:(1 - fs1);
%             for s2 = 1:length(fss2)
%                 fs2 = fss2(s2);
%                 Parameters_Name = ['fp',num2str(fps(p)*100),'_fr',num2str(frs(r)*100),'_fs',num2str(fs1*100),'_fs',num2str(fs2*100)];
%                 sim_results_Name=['sim_results_',Parameters_Name];
%                 diag_results_Name=['diag_results_',Parameters_Name];
%                 impact_results_Name=['impact_results_',Parameters_Name];
%                 load(['./output/shock_results/',net_name,'/',shocked_product,'/',sim_results_Name,'.mat'],sim_results_Name);
%                 eval(['sim_results = ', sim_results_Name,';']);
%             
%                 diag_results = diag_analysis(sim_results, country_id, CNAME, N);  
%                 eval([diag_results_Name,'= diag_results;']);
%                 clear diag_results
%             
%                 impact_results = impact_analysis(sim_results, country_id, CNAME, N, M); 
%                 eval([impact_results_Name,'= impact_results;']);
%                 clear impact_results
%             
%                 save(['./output/analysis_results/',net_name,'/',shocked_product,'/',diag_results_Name,'.mat'],diag_results_Name);
%                 save(['./output/analysis_results/',net_name,'/',shocked_product,'/',impact_results_Name,'.mat'],impact_results_Name);
%                 disp(['Calculation complete: three-layer-chain_analysis_results_',Parameters_Name]);
%                 eval(['clear ',diag_results_Name,';']); 
%                 eval(['clear ',impact_results_Name,';']); 
%                 clear sim_results
%             end
%         end
%     end
% end

%% Validate and calculate on three-layer network with triadic structure
% for p = 1:length(fps) 
%     for r = 1:length(frs)
%         for s = 1:length(fss)
%             fs = fss(s);          
%             fss1 = 0.1:0.1:(fs-0.1);   
%             for s1 = 1:length(fss1)
%                 fs1 = fss1(s1);  
%                 fs2 = fs - fs1;  
%                 Parameters_Name = ['fp',num2str(fps(p)*100),'_fr',num2str(frs(r)*100),'_fs',num2str(fs1*100),'_fs',num2str(fs2*100)];
%                 sim_results_Name=['sim_results_',Parameters_Name];
%                 diag_results_Name=['diag_results_',Parameters_Name];
%                 impact_results_Name=['impact_results_',Parameters_Name];
%                 load(['./output/shock_results/',net_name,'/',shocked_product,'/',sim_results_Name,'.mat'],sim_results_Name);
%                 eval(['sim_results = ', sim_results_Name,';']);
%             
%                 diag_results = diag_analysis(sim_results, country_id, CNAME, N);  
%                 eval([diag_results_Name,'= diag_results;']);
%                 clear diag_results
%             
%                 impact_results = impact_analysis(sim_results, country_id, CNAME, N, M); 
%                 eval([impact_results_Name,'= impact_results;']);
%                 clear impact_results
%             
%                 save(['./output/analysis_results/',net_name,'/',shocked_product,'/',diag_results_Name,'.mat'],diag_results_Name);
%                 save(['./output/analysis_results/',net_name,'/',shocked_product,'/',impact_results_Name,'.mat'],impact_results_Name);
%                 disp(['Calculation complete: three-layer-triangular_analysis_results_',Parameters_Name]);
%                 eval(['clear ',diag_results_Name,';']); 
%                 eval(['clear ',impact_results_Name,';']); 
%                 clear sim_results
%             end
%         end
%     end
% end

%% Validate and calculate on three-layer network with triangular structure
% for p = 1:length(fps) 
%     for r = 1:length(frs)
%         for s = 1:length(fss)
%             fs = fss(s);          
%             fss1 = 0.1:0.1:(fs-0.1);   
%             for s1 = 1:length(fss1)
%                 fs1 = fss1(s1);  
%                 fs2 = fs - fs1;  
%                 fss3 = 0.1:0.1:min(1-fs1,1-fs2);  
%                 for s3 = 1:length(fss3)
%                      fs3 = fss3(s3);
%                      Parameters_Name = ['fp',num2str(fps(p)*100),'_fr',num2str(frs(r)*100),'_fs',num2str(fs1*100),'_fs',num2str(fs2*100),'_fs',num2str(fs3*100)];
%                      sim_results_Name=['sim_results_',Parameters_Name];
%                      diag_results_Name=['diag_results_',Parameters_Name];
%                      impact_results_Name=['impact_results_',Parameters_Name];
%                      load(['./output/shock_results/',net_name,'/',shocked_product,'/',sim_results_Name,'.mat'],sim_results_Name);
%                      eval(['sim_results = ', sim_results_Name,';']);
%             
%                      diag_results = diag_analysis(sim_results, country_id, CNAME, N);  
%                      eval([diag_results_Name,'= diag_results;']);
%                      clear diag_results
%             
%                      impact_results = impact_analysis(sim_results, country_id, CNAME, N, M); 
%                      eval([impact_results_Name,'= impact_results;']);
%                      clear impact_results
%             
%                      save(['./output/analysis_results/',net_name,'/',shocked_product,'/',diag_results_Name,'.mat'],diag_results_Name);
%                      save(['./output/analysis_results/',net_name,'/',shocked_product,'/',impact_results_Name,'.mat'],impact_results_Name);
%                      disp(['Calculation complete: three-layer-triadic_analysis_results_',Parameters_Name]);
%                      eval(['clear ',diag_results_Name,';']); 
%                      eval(['clear ',impact_results_Name,';']); 
%                      clear sim_results
%                 end
%             end
%         end
%     end
% end


%% Verification and computation functions
function diag_results = diag_analysis(sim_results, country_id, CNAME, N)
%Function to validate mass balance equations
equ_num = 0;    %Number of simulates reaching equilibrium within kmax iterations
equ_res = zeros(1,length(country_id));   %Equilibrium status record for each simulatek: 1 = equilibrium, 0 = no-equilibrium.
diag_res = zeros(1,length(country_id));  %Mass balance compliance record for each simulatek: 4 = full compliance.
for c = 1:length(country_id)
    c_init = country_id(c);     
    sim_res = sim_results.(CNAME{c_init});
                        
    equ_num = equ_num + sim_res.equilibrium;
    equ_res(c) = sim_res.equilibrium;
    diag_res(c) = sim_diagnostics(sim_res,N);
                        
    clear sim_res
end
diag_num = sum(diag_res == 4);   %Number of simulations compliant with mass balance equations
if isequal(equ_num, diag_num)   %Equality of both indicates correct results: reaching equilibrium within kmax and complying with mass balance equations, with a value of 1.
    iscor_res = 1;
else
    iscor_res = 0;
end
diag_results = struct('equ_num', equ_num, 'equ_result', equ_res, 'diag_num', diag_num, 'diag_result', diag_res, 'iscor_res', iscor_res);
end


function impact_results = impact_analysis(sim_results, country_id, CNAME, N, M)  
%Function to calculate shock impact.
delta_C = zeros(length(country_id), length(CNAME));   %consumption changes for individual country-product, with rows representing different initial shock countries and columns representing country-products across layers
nonTrade_country = cell(1,M+1);

%impact of single shocks
for c = 1:length(country_id)
    c_init = country_id(c);     
    sim_res = sim_results.(CNAME{c_init});
    
    S = sim_res.S;
    dC = sim_res.dC;
    dC(abs(dC) < 1) = 0;
    delta_C(c, :) = sum(dC, 2)';  
    C0 = S(:,1)';  %initial consumption of country-products.
    
    A_SUM = zeros(length(1:N), length(1:N));
    if c == 1
        A = sim_res.A;
        A0 =  A(:,:,1);
        for m = 1:(M+1)
           if m < M+1
               index = [1:171] + 171*(m-1);
               A_L = A0(index,index);
               flowSum = sum(A_L, 2) + sum(A_L, 1)';
               nonTrade_country{m} = find(flowSum == 0)'; 
               A_SUM  = A_SUM + A_L;  
           else
               flowSum = sum(A_SUM, 2) + sum(A_SUM, 1)';
               nonTrade_country{m} = find(flowSum == 0)'; 
           end
        end     
    end
    
    delta_C_SUM = zeros(1, length(1:N));
    C0_SUM = zeros(1, length(1:N));
    D_CP = zeros(length(1:N), length(1:M));   %deficit of country-products, with rows representing different countries and columns representing different product layers
    D_L = zeros(1, length(1:M));   %avenge deficit of product layers
    U_L = zeros(1, length(1:M));   %deficit unevenness of product layers
    for m = 1:(M+1)
        if m < M+1
            index = [1:171] + 171*(m-1);
            D_CP_L = (delta_C(c, index)./C0(1, index))';        %relative consumption changes
            D_CP_L(isnan(D_CP_L)) = 0; 
            D_CP(:,m) = D_CP_L;
            
            D_CP_L(nonTrade_country{m},1) = NaN;
            D_L(1,m) = nanmean(D_CP_L);
            U_L(1,m) = nanstd(D_CP_L);
            
            delta_C_SUM = delta_C_SUM + delta_C(c, index);
            C0_SUM = C0_SUM +  C0(1, index);
        else
            D_C= (delta_C_SUM ./ C0_SUM)';
            D_C(isnan(D_C)) = 0; 
            
            D_C_N = D_C;
            D_C_N(nonTrade_country{m},1) = NaN;
            D_N = nanmean(D_C_N);  %avenge deficit of overall network
            U_N = nanstd(D_C_N);   %deficit unevenness of overall network
        end
    end
    
    shock_impact = struct('deficit_forCountryProduct', D_CP, 'deficit_forCountry', D_C, ...
       'deficit_forProduct', D_L, 'unevenness_forProduct', U_L, ...
       'deficit_forNet', D_N, 'unevenness_forNet', U_N);
   
    singleShock_impact.(CNAME{c_init}) = shock_impact;
end

%impact of multiple shocks
D_CP_all = zeros(length(1:N), length(1:M));
D_L_all = zeros(1, length(1:M));
U_L_all = zeros(1, length(1:M));
delta_C_all = zeros(length(country_id), length(1:N));
for m = 1:(M+1)
    if m < M+1
        index = [1:171] + 171*(m-1);
        D_CP_L_all =  (sum(delta_C(:, index), 1) ./ C0(1, index))';
        D_CP_L_all(isnan(D_CP_L_all)) = 0; 
        D_CP_all(:,m) = D_CP_L_all;
        
        D_CP_L_all(nonTrade_country{m},1) = NaN;
        D_L_all(1,m) = nanmean(D_CP_L_all);
        U_L_all(1,m) = nanstd(D_CP_L_all);
        
        delta_C_all = delta_C_all + delta_C(:, index);
    else
        D_C_all= (sum(delta_C_all, 1) ./ C0_SUM)';
        D_C_all(isnan(D_C_all)) = 0; 
            
        D_C_all_N = D_C_all;
        D_C_all_N(nonTrade_country{m},1) = NaN;
        D_N_all = nanmean(D_C_all_N);
        U_N_all = nanstd(D_C_all_N);
    end
end

multiShock_impact = struct('deficit_forCountryProduct',D_CP_all, 'deficit_forCountry', D_C_all, ...
       'deficit_forProduct', D_L_all, 'unevenness_forProduct', U_L_all, ...
       'deficit_forNet', D_N_all, 'unevenness_forNet', U_N_all);

impact_results = struct('singleShock_impact', singleShock_impact,...
    'multiShock_impact', multiShock_impact);
end
