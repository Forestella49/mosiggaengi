%==========================================================================
% PHY-Security of Alamouti Code for 2x1 system 
% LML-LML, LML-DF, DF-LML, DF-DF
% Last updated : 2022. 10. 04.
%==========================================================================


function SystemOut = ConSecOutageProb_LML_ZF_DF( SystemParam )


    dBToLinear = @(x)( 10.^(x./10) );    

    % SystemParam : iter, R_T, R_E, rho_B, rho_E, Es_dBm_sample, Noise_Var_B, Noise_Var_E
    
    iter = SystemParam{1};
    R_T = SystemParam{2};
    R_E = SystemParam{3};  
    rho_B = SystemParam{4};
    rho_E = SystemParam{5};
    Es_dBm = SystemParam{6};
    Noise_Var_B = SystemParam{7}; % Bob's Noise Variance
    Noise_Var_E = SystemParam{8}; % Eve's Noise Variance
    
    R_S = R_T - R_E; 
    Es = 10^(Es_dBm/10)/1000; % (dBm to watts)

    sigma_E = sqrt(Noise_Var_E); % Eve's Noise Standard Deviation
    sigma_B = sqrt(Noise_Var_B); % Bob's Noise Standard Deviation

    avg_SNR_B = (2*Es) / (sigma_B)^2;
    avg_SNR_E = (2*Es) / (sigma_E)^2;

    %=======================================================================================
    % avg_SNR_B / avg_SNR_E = (2Es / sigma_B^2)/(2Es / sigma_E^2) = sigma_E^2 / sigama_B^2
    %=======================================================================================     
    %==========================================================================
    % Correlated Channel Generation
    % Paper Title: "A Generalized Linear Quasi-ML Decoder of OSTBCs for WIreless
    %             Communications Over Time-Selective Fading Channels"
    %             -> h_i[n] = rho * h_i[n-1] + v_i[n], for i = 1,2
    %             where |rho|^2 + sigma_v^2 = 1,  v_i ~ CN(0,sigma_v^2)
    %==========================================================================
    Noise_Var_E_v = 1 - rho_E ^2; 
    Noise_Var_B_v = 1 - rho_B ^2;

    count_LML_LML = zeros(1,iter);
    count_LML_DF = zeros(1,iter);
    count_DF_LML = zeros(1,iter);
    count_DF_DF = zeros(1,iter);
        
    for k = 1:iter
    
        %if (rem(k,10000)==0)
        %    disp(['iteration = ', num2str(k)]);
        %end
      
   %=====================================
   % Time Selective Channel Generation
   %=====================================
    H_B1 = (randn(2,1)+sqrt(-1)*randn(2,1)) / sqrt(2);
    H_B2 = rho_B.*H_B1 + (randn(2,1)+sqrt(-1)*randn(2,1)).* sqrt(Noise_Var_B_v)./sqrt(2);
    H_B3 = [H_B1,H_B2]; % Each column : 1st time slot, 2nd time slot
    H_B3_JML = [H_B1,H_B1]; % Each column : 1st time slot, 2nd time slot
    H_B = [H_B3(1,1), H_B3(2,1); conj(H_B3(2,2)), - conj(H_B3(1,2))];
    G_B = H_B'*H_B; 
    H_B_JML = [H_B3_JML(1,1), H_B3_JML(2,1); conj(H_B3_JML(2,2)), -conj(H_B3_JML(1,2))]; 
    G_B_JML= H_B_JML'*H_B_JML;
      
    H_E1 = (randn(2,1)+sqrt(-1)*randn(2,1)) / sqrt(2); % 1st time slot
    H_E2 = rho_E.*H_E1 +  (randn(2,1)+sqrt(-1)*randn(2,1)).* sqrt(Noise_Var_E_v)./sqrt(2); % 2nd time slot
    H_E3 = [H_E1,H_E2]; % Each column : 1st time slot, 2nd time slot
    H_E3_JML = [H_E1,H_E1]; % Each column : 1st time slot, 2nd time slot
    H_E = [H_E3(1,1), H_E3(2,1); conj(H_E3(2,2)), - conj(H_E3(1,2))];
    G_E = H_E'*H_E; 
    H_E_JML = [H_E3_JML(1,1), H_E3_JML(2,1); conj(H_E3_JML(2,2)), -conj(H_E3_JML(1,2))]; 
    G_E_JML= H_E_JML'*H_E_JML;
        
    %=========================================
    % Effective received SNR for each antenna
    %=========================================     
    % LML at Bob
    B_Eff_SNR_Ant1_LML = avg_SNR_B / ((1-rho_B^2)*avg_SNR_B + 2) * G_B(1,1);
    B_Eff_SNR_Ant2_LML = avg_SNR_B / ((1-rho_B^2)*avg_SNR_B + 2) * G_B(2,2);
   
    % LML at Eve
    E_Eff_SNR_Ant1_LML = avg_SNR_E / ((1-rho_E^2)*avg_SNR_E + 2) * G_E(1,1);
    E_Eff_SNR_Ant2_LML = avg_SNR_E / ((1-rho_E^2)*avg_SNR_E + 2) * G_E(2,2);

    % ZF at Bob
    B_zeta = sqrt(G_B(1,1)*G_B(2,2) - abs(G_B(1,2))^2);   
    B_Eff_SNR_Ant1_ZF = B_zeta^2 / (2*G_B(2,2)) * avg_SNR_B;
    B_Eff_SNR_Ant2_ZF = B_zeta^2 / (2*G_B(1,1)) * avg_SNR_B;
    
    % ZF at Eve
    E_zeta = sqrt(G_E(1,1)*G_E(2,2) - abs(G_E(1,2))^2);
    E_Eff_SNR_Ant1_ZF = E_zeta^2 / (2*G_E(2,2)) * avg_SNR_E;
    E_Eff_SNR_Ant2_ZF = E_zeta^2 / (2*G_E(1,1)) * avg_SNR_E;
    
    % DF at Bob
    B_zeta = sqrt(G_B(1,1)*G_B(2,2) - abs(G_B(1,2))^2);   
    B_Eff_SNR_Ant1_DF = B_zeta^2 / (2*G_B(2,2)) * avg_SNR_B; % ZF
    B_Eff_SNR_Ant2_DF = avg_SNR_B / 2 * G_B_JML(2,2); % JML
    
    % DF at Eve
    E_zeta = sqrt(G_E(1,1)*G_E(2,2) - abs(G_E(1,2))^2);
    E_Eff_SNR_Ant1_DF = E_zeta^2 / (2*G_E(2,2)) * avg_SNR_E; % ZF
    E_Eff_SNR_Ant2_DF = avg_SNR_E / 2 * G_E_JML(2,2); % JML
      
    C_B_LML = 1/2*( log2(1+B_Eff_SNR_Ant1_LML) + log2(1+B_Eff_SNR_Ant2_LML) ); % Capacity of LML at Bob
    C_B_ZF  = 1/2*( log2(1+B_Eff_SNR_Ant1_ZF)  + log2(1+B_Eff_SNR_Ant2_ZF)  ); % Capacity of ZF at Bob
    C_B_DF  = 1/2*( log2(1+B_Eff_SNR_Ant1_DF)  + log2(1+B_Eff_SNR_Ant2_DF)  ); % Capacity of DF at Bob
    
    C_E_LML = 1/2*( log2(1+E_Eff_SNR_Ant1_LML) + log2(1+E_Eff_SNR_Ant2_LML) ); % Capacity of LML at Eve
    C_E_ZF  = 1/2*( log2(1+E_Eff_SNR_Ant1_ZF)  + log2(1+E_Eff_SNR_Ant2_ZF)  ); % Capacity of ZF at Eve
    C_E_DF  = 1/2*( log2(1+E_Eff_SNR_Ant1_DF)  + log2(1+E_Eff_SNR_Ant2_DF)  ); % Capacity of DF at Eve
    
    count_ConOutage_LML(k) = (C_B_LML <= R_T);
    count_ConOutage_ZF(k)  = (C_B_ZF <= R_T);
    count_ConOutage_DF(k)  = (C_B_DF <= R_T);

    count_SecOutage_LML(k) = (C_E_LML > R_E);
    count_SecOutage_ZF(k)  = (C_E_ZF > R_E);
    count_SecOutage_DF(k)  = (C_E_DF > R_E);
    
    
    end
    
    ConOutage_LML = sum(count_ConOutage_LML)/iter;
    ConOutage_ZF  = sum(count_ConOutage_ZF) /iter;
    ConOutage_DF  = sum(count_ConOutage_DF) /iter;

    SecOutage_LML  = sum(count_SecOutage_LML) /iter;
    SecOutage_ZF   = sum(count_SecOutage_ZF) /iter;
    SecOutage_DF   = sum(count_SecOutage_DF) /iter;
    
    SystemOut = {ConOutage_LML, ConOutage_ZF, ConOutage_DF, SecOutage_LML, SecOutage_ZF, SecOutage_DF, R_S};
    
end
     