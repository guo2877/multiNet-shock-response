% script simulates shock responses and their cascades on multilayer networks 
% under different parameter configurations

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
eval(['CIS = {',data_name ,'.CId};'])  %country id
eval(['CNS = {',data_name ,'.CName};']) %country name
eval(['PS = {',data_name ,'.Prod};'])  %initial production vector 
eval(['RS = {',data_name ,'.endStock};']) %initial reserve vector
eval(['AS = {',data_name ,'.tradeMatrix};']) %initial trade matrix
eval(['FS = {',data_name ,'.scMatrix};']) %initial substitution fraction matrix

CID = CIS{year_id}; CNAME = CNS{year_id};
P0 = PS{year_id}; R0 = RS{year_id}; A0 = AS{year_id}; F0 = FS{year_id}; 

M = size(Products_Names, 2);   %number of product types 
N = length(CNAME)/M;  %number of countries

%% Shock protocol
for i = 1:size(Products_Names, 2)
   if strcmp(shocked_product, [Products_Names{i}])
       country_id = [1:N] + N*(i-1);  %exhaustive shock protocol: shocks are applied separately to each country in the shocked layer
       %country_id = 162 + N*(i-1);   %single shock protocol: shock to one country£¨cid, 162 = ukr) in the shocked layer
   end 
end

%% Simulate shock on two-layer network 
for p = 1:length(fps) 
    for r = 1:length(frs)
        for s = 1:length(fss)
            Parameters_Name = ['fp',num2str(fps(p)*100),'_fr',num2str(frs(r)*100),'_fs',num2str(fss(s)*100)];
            sim_results_Name=['sim_results_',Parameters_Name];
            eval([sim_results_Name,' = struct();']);
            for c = 1:length(country_id)
                c_init = country_id(c);    
                fp = fps(p);               
                fr = frs(r);               
                fs = fss(s);                
                F = F0 .* fs;  % set substitution fraction matrix
                sim_result = shock_response(P0, R0, A0, c_init, fp, fr, F, amin, kmax, asym); 
                eval([sim_results_Name,'.(CNAME{c_init})= sim_result;']);
                clear sim_result
            end
            save(['./output/shock_results/',net_name,'/',shocked_product,'/',sim_results_Name,'.mat'],sim_results_Name, '-v7.3');
            disp(['Processing complete: two-layer_sim_results_',Parameters_Name]);
            eval(['clear ',sim_results_Name,';']);
        end
    end
end

%% Simulate shock on three-layer network with chain structure
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
%                 eval([sim_results_Name,' = struct();']);
%                 for c = 1:length(country_id)
%                      c_init = country_id(c);    
%                      fp = fps(p);                
%                      fr = frs(r); 
%                      % set substitution fraction matrix
%                      F = F0;
%                      F(1:171,172:342) = F0(1:171,172:342) .* fs1;
%                      F(172:342,1:171) = F0(172:342,1:171) .* fs1;
%                      F(172:342,343:513) = F0(172:342,343:513) .* fs2;
%                      F(343:513,172:342) = F0(343:513,172:342) .* fs2;
%                      sim_result = shock_response(P0, R0, A0, c_init, fp, fr, F, amin, kmax, asym);
%                      eval([sim_results_Name,'.(CNAME{c_init})= sim_result;']);
%                      clear sim_result
%                 end
%                 save(['./output/shock_results/',net_name,'/',shocked_product,'/',sim_results_Name,'.mat'],sim_results_Name, '-v7.3');
%                 disp(['Processing complete: three-layer-chain_sim_results_',Parameters_Name]);
%                 eval(['clear ',sim_results_Name,';']); 
%             end
%         end
%     end
% end


%% Simulate shock on three-layer network with triadic structure
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
%                 eval([sim_results_Name,' = struct();']);
%                 for c = 1:length(country_id)
%                      c_init = country_id(c);    
%                      fp = fps(p);                
%                      fr = frs(r); 
%                      % set substitution fraction matrix
%                      F = F0;
%                      F(1:171,172:342) = F0(1:171,172:342) .* fs1;
%                      F(172:342,1:171) = F0(172:342,1:171) .* fs1;
%                      F(172:342,343:513) = F0(172:342,343:513) .* fs2;
%                      F(343:513,172:342) = F0(343:513,172:342) .* fs2;
%                      sim_result = shock_response(P0, R0, A0, c_init, fp, fr, F, amin, kmax, asym);
%                      eval([sim_results_Name,'.(CNAME{c_init})= sim_result;']);
%                      clear sim_result
%                 end
%                 save(['./output/shock_results/',net_name,'/',shocked_product,'/',sim_results_Name,'.mat'],sim_results_Name, '-v7.3');
%                 disp(['Processing complete: three-layer-triadic_sim_results_',Parameters_Name]);
%                 eval(['clear ',sim_results_Name,';']); 
%             end
%         end
%     end
% end


%% Simulate shock on three-layer network with triangular structure
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
%                      eval([sim_results_Name,' = struct();']);
%                      for c = 1:length(country_id)
%                         c_init = country_id(c);    
%                         fp = fps(p);                
%                         fr = frs(r); 
%                         % set substitution fraction matrix
%                         F = F0;
%                         F(1:171,172:342) = F0(1:171,172:342) .* fs1;
%                         F(172:342,1:171) = F0(172:342,1:171) .* fs1;
%                         F(1:171,343:513) = F0(1:171,343:513) .* fs2;
%                         F(343:513,1:171) = F0(343:513,1:171) .* fs2;
%                         F(172:342,343:513) = F0(172:342,343:513) .* fs3;
%                         F(343:513,172:342) = F0(343:513,172:342) .* fs3;
%                         sim_result = shock_response(P0, R0, A0, c_init, fp, fr, F, amin, kmax, asym);
%                         eval([sim_results_Name,'.(CNAME{c_init})= sim_result;']);
%                         clear sim_result
%                      end
%                      save(['./output/shock_results/',net_name,'/',shocked_product,'/',sim_results_Name,'.mat'],sim_results_Name, '-v7.3');
%                      disp(['Processing complete: three-layer-triangular_sim_results_',Parameters_Name]);
%                      eval(['clear ',sim_results_Name,';']); 
%                 end
%             end
%         end
%     end
% end
