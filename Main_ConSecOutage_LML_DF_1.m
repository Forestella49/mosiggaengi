%===================================================================================
% PHY-Security for Time-Selective Fading Channels : Simulation
% Connection Outage & Secrecy Outage Prob. for Alamouti 2x1 system: LML, ZF Combination 
% Copyright: Prof. Seong Ho Chae / 2022.03.16 
%===================================================================================

clear all; close all; clc;

%% simulation setting
% OnSave = true;
% Es_dBm_min = 20;
% Es_dBm_max = 50;
% Es_dBm_stepsize = 5;

Es_dBm = 40;   %

rho_B = 0.95; % Time Correlation of Bob's Channel
rho_E = 0.95; % Time Correlation of Eve's Channel

Noise_Var_B = 1;  % Bob's Noise Variance
Noise_Var_E = 1;  % Eve's Noise Variance

n = 100; % stepsize

R_T_min = 0;
R_T_max = 8;
R_T_stepsize = (R_T_max - R_T_min)/(n-1);

R_E_min = 0;
R_E_max = 8;
R_E_stepsize = (R_E_max - R_E_min)/(n-1);

R_T = [R_T_min : R_T_stepsize : R_T_max]; % Target Secrecy Rate [bps/Hz]
R_E = [R_E_min : R_E_stepsize : R_E_max]; 

% R_S = R_T - R_E;

R_T_th = 2.^(R_T) -1;
R_E_th = 2.^(R_E) -1;

iter = 10^4;  %iteration

nData = length(R_T);

Num_ConOutage_LML = zeros(1,nData);
Num_ConOutage_ZF  = zeros(1,nData);
Num_ConOutage_DF  = zeros(1,nData);
Num_SecOutage_LML = zeros(1,nData);
Num_SecOutage_ZF  = zeros(1,nData);
Num_SecOutage_DF  = zeros(1,nData);

parfor k = 1:nData 
    

    R_T_sample = R_T(k);
    R_E_sample = R_E(k);
    
    disp([]);
    
    SystemParam = {iter, R_T_sample, R_E_sample, rho_B, rho_E, Es_dBm, Noise_Var_B, Noise_Var_E};
    SystemOut = ConSecOutageProb_LML_ZF_DF( SystemParam );
    Num_ConOutage_LML(1,k) = [SystemOut{1}];
    Num_ConOutage_ZF(1,k)  = [SystemOut{2}];
    Num_ConOutage_DF(1,k)  = [SystemOut{3}];
    Num_SecOutage_LML(1,k)  = [SystemOut{4}];
    Num_SecOutage_ZF(1,k)   = [SystemOut{5}];
    Num_SecOutage_DF(1,k)   = [SystemOut{6}];
%     R_S(1,k)                = [SystemOut{7}];
    
    disp(['E_s(dBm) = ',num2str(R_T_sample),...
        ', Num_ConOutage_LML = ',num2str(SystemOut{1}), ', Num_ConOutage_ZF = ',num2str(SystemOut{2}), ', Num_ConOutage_DF = ',num2str(SystemOut{3}),...
        ', Num_SecOutage_LML = ',num2str(SystemOut{4}), ', Num_SecOutage_ZF = ',num2str(SystemOut{5}), ', Num_SecOutage_DF = ',num2str(SystemOut{6}),...
%         ', R_S = ',num2str(SystemOut{5}) ...
        ]);
    disp('   '); 

end
%% Rs * Num_ConOutage * Num_SecOutage

Num_Bob_Pr_LML = [];
Num_Bob_Pr_ZF = []; 
Num_Bob_Pr_DF = []; 
Num_Eve_Pr_LML = [];
Num_Eve_Pr_ZF = [];
Num_Eve_Pr_DF = [];

for  Bob = 1:1:n
   for Eve = 1:1:n
       if Bob > Eve
           R_S(Eve,Bob) = R_T(Bob) - R_E(Eve);           
       else
           R_S(Eve,Bob) = 0;
       end       
   end
   Num_Bob_Pr_LML = [Num_Bob_Pr_LML; Num_ConOutage_LML];
   Num_Bob_Pr_DF = [Num_Bob_Pr_DF; Num_ConOutage_DF];

   Num_Eve_Pr_LML = [Num_Eve_Pr_LML transpose(Num_SecOutage_LML)];
   Num_Eve_Pr_DF = [Num_Eve_Pr_DF transpose(Num_SecOutage_DF)];

   Num_Bob_Pr_ZF = [Num_Bob_Pr_ZF; Num_ConOutage_ZF];
   Num_Eve_Pr_ZF = [Num_Eve_Pr_ZF transpose(Num_SecOutage_ZF)];
end

Num_Bob_LML_Eve_LML = R_S.*(1-Num_Bob_Pr_LML).*(1-Num_Eve_Pr_LML);
Num_Bob_LML_Eve_ZF = R_S.*(1-Num_Bob_Pr_LML).*(1-Num_Eve_Pr_ZF);
Num_Bob_LML_Eve_DF = R_S.*(1-Num_Bob_Pr_LML).*(1-Num_Eve_Pr_DF);


Num_Bob_ZF_Eve_LML = R_S.*(1-Num_Bob_Pr_ZF).*(1-Num_Eve_Pr_LML);
Num_Bob_ZF_Eve_ZF = R_S.*(1-Num_Bob_Pr_ZF).*(1-Num_Eve_Pr_ZF);
Num_Bob_ZF_Eve_DF = R_S.*(1-Num_Bob_Pr_ZF).*(1-Num_Eve_Pr_DF);

Num_Bob_DF_Eve_LML = R_S.*(1-Num_Bob_Pr_DF).*(1-Num_Eve_Pr_LML);
Num_Bob_DF_Eve_ZF = R_S.*(1-Num_Bob_Pr_DF).*(1-Num_Eve_Pr_ZF);
Num_Bob_DF_Eve_DF = R_S.*(1-Num_Bob_Pr_DF).*(1-Num_Eve_Pr_DF);

[R_T_LML_LML,R_E_LML_LML] = find(Num_Bob_LML_Eve_LML==max(max(Num_Bob_LML_Eve_LML)));
[R_T_LML_DF,R_E_LML_DF] = find(Num_Bob_LML_Eve_DF==max(max(Num_Bob_LML_Eve_DF)));
[R_T_LML_ZF,R_E_LML_ZF] = find(Num_Bob_LML_Eve_ZF==max(max(Num_Bob_LML_Eve_ZF)));

[R_T_ZF_LML,R_E_ZF_LML] = find(Num_Bob_ZF_Eve_LML==max(max(Num_Bob_ZF_Eve_LML)));
[R_T_ZF_ZF,R_E_ZF_ZF] = find(Num_Bob_ZF_Eve_ZF==max(max(Num_Bob_ZF_Eve_ZF)));
[R_T_ZF_DF,R_E_ZF_DF] = find(Num_Bob_ZF_Eve_DF==max(max(Num_Bob_ZF_Eve_DF)));

[R_T_DF_LML,R_E_DF_LML] = find(Num_Bob_DF_Eve_LML==max(max(Num_Bob_DF_Eve_LML)));
[R_T_DF_ZF,R_E_DF_ZF] = find(Num_Bob_DF_Eve_ZF==max(max(Num_Bob_DF_Eve_ZF)));
[R_T_DF_DF,R_E_DF_DF] = find(Num_Bob_DF_Eve_DF==max(max(Num_Bob_DF_Eve_DF)));



% Num_Eve_DF_Bob_DF = R_S .* (1-Num_ConOutage_LML) .* (1-Num_SecOutage_LML);
%% 
close all;
MS = 10; LW = 1.5;

f = figure('Name','Num Eve_LML * Bob_LML','NumberTitle','off');
% f.Position(3:4) = [420 840];
f.Position(3:4) = [840 840];
tiledlayout(3,3);

ax1 = nexttile;
contour(ax1, R_T,R_E, Num_Bob_LML_Eve_LML); hold on;
p1 = plot(R_E(R_E_LML_LML),R_T(R_T_LML_LML),'r*');
legend(p1,{'MAX'},'Location','northwest')
colorbar
grid on; box on;
xlabel(ax1,'{\it R}_T[bps/Hz]');
ylabel(ax1,'{\it R}_E[bps/Hz]');

ax2 = nexttile;
contour(ax2, R_T,R_E, Num_Bob_LML_Eve_DF); hold on;
p1 = plot(R_E(R_E_LML_DF),R_T(R_T_LML_DF),'r*');
legend(p1,{'MAX'},'Location','northwest')
colorbar
grid on; box on;
xlabel(ax2,'{\it R}_T[bps/Hz]');
ylabel(ax2,'{\it R}_E[bps/Hz]');

ax3 = nexttile;
contour(ax3, R_T,R_E, Num_Bob_LML_Eve_ZF); hold on;
p1 = plot(R_E(R_E_LML_ZF),R_T(R_T_LML_ZF),'r*');
legend(p1,{'MAX'},'Location','northwest')
colorbar
grid on; box on;
xlabel(ax3,'{\it R}_T[bps/Hz]');
ylabel(ax3,'{\it R}_E[bps/Hz]');

ax4 = nexttile;
contour(ax4, R_T,R_E, Num_Bob_ZF_Eve_LML); hold on;
p1 = plot(R_E(R_E_ZF_LML),R_T(R_T_ZF_LML),'r*');
legend(p1,{'MAX'},'Location','northwest')
colorbar
grid on; box on;
xlabel(ax4,'{\it R}_T[bps/Hz]');
ylabel(ax4,'{\it R}_E[bps/Hz]');

ax5 = nexttile;
contour(ax5, R_T,R_E, Num_Bob_ZF_Eve_ZF); hold on;
p1 = plot(R_E(R_E_ZF_ZF),R_T(R_T_ZF_ZF),'r*');
legend(p1,{'MAX'},'Location','northwest')
colorbar
grid on; box on;
xlabel(ax5,'{\it R}_T[bps/Hz]');
ylabel(ax5,'{\it R}_E[bps/Hz]');

ax6 = nexttile;
contour(ax6, R_T,R_E, Num_Bob_ZF_Eve_DF); hold on;
p1 = plot(R_E(R_E_ZF_DF),R_T(R_T_ZF_DF),'r*');
legend(p1,{'MAX'},'Location','northwest')
colorbar
grid on; box on;
xlabel(ax6,'{\it R}_T[bps/Hz]');
ylabel(ax6,'{\it R}_E[bps/Hz]');

ax7 = nexttile;
contour(ax7, R_T,R_E, Num_Bob_DF_Eve_LML); hold on;
p1 = plot(R_E(R_E_DF_LML),R_T(R_T_DF_LML),'r*');
legend(p1,{'MAX'},'Location','northwest')
colorbar
grid on; box on;
xlabel(ax7,'{\it R}_T[bps/Hz]');
ylabel(ax7,'{\it R}_E[bps/Hz]');

ax8 = nexttile;
contour(ax8, R_T,R_E, Num_Bob_DF_Eve_ZF); hold on;
p1 = plot(R_E(R_E_DF_ZF),R_T(R_T_DF_ZF),'r*');
legend(p1,{'MAX'},'Location','northwest')
colorbar
grid on; box on;
xlabel(ax8,'{\it R}_T[bps/Hz]');
ylabel(ax8,'{\it R}_E[bps/Hz]');

ax9 = nexttile;
contour(ax9, R_T,R_E, Num_Bob_DF_Eve_DF); hold on;
p1 = plot(R_E(R_E_DF_DF),R_T(R_T_DF_DF),'r*');
legend(p1,{'MAX'},'Location','northwest')
colorbar
grid on; box on;
xlabel(ax9,'{\it R}_T[bps/Hz]');
ylabel(ax9,'{\it R}_E[bps/Hz]');

title(ax1,'LML-LML')
title(ax2,'LML-ZF')
title(ax3,'LML-DF')

title(ax4,'ZF-LML')
title(ax5,'ZF-ZF')
title(ax6,'ZF-DF')

title(ax7,'DF-LML')
title(ax8,'DF-ZF')
title(ax9,'DF-DF')


