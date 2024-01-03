%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code for Zhou, et al.

% Based on the code structure from Kwong et al., PNAS (2015), 112(41):12627.
% and Chan et al., Nature nanotechnology, 2020, 15(9): 792-800.

clear;
close all;
clc;

M_ratio = 0.1;                              % macrophage membrane ratio

%initial ode parameters
options = odeset('AbsTol',1e-6,'RelTol',1e-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Models the first 12 hours after SEDOM administration. This could be
% changed by replacing 720 with any number of minites.
%
% The SEDOM dose was calculated by function SEconc(), which is included at
% the end of this code. Unless otherwise stated, all concentration units
% are um.
%
[t,y] = ode15s(@(t,y) NumSimulation(t,y,M_ratio),0:0.1:720, [SEconc(M_ratio), 0, 0, 0, 0, 0, 0], options); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(t,y(:,7),'r');
ylabel('urineR(mM)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = NumSimulation(t,y,M_ratio)

%%% Clearing rates %%%
%liver clearing rates
l_NP_liClear = 2.82e-2;                     % m^-1 Reflect clearance rate of Nanoparticle from liver (measured).
%kidney clearing rates
l_R_kiClear = 0.166;                        % m^-1 Urinary filtration of Reporter (measured).

%%% Urine parameters %%%
a_urine = 1.6/24/60;                        % urine generation rate 1.6 mL in 24 hours.
V_urine = 0.2;                              % initial urine volume 200 uL.

%%% Transport rates %%%
k_NP_tumor = 8e-03;                         % m^-1 Transport of Nanoparticle into tumor (manually fitted).
k_R_tumor = 0.1;                            % m^-1 Transport of Reporter out of tumor (manually fitted).
k_NP_liver = 2e-03;                         % m^-1 Transport of Nanoparticle into liver (manually fitted).
k_R_liver = 0.1;                            % m^-1 Transport of Reporter out of liver (manually fitted).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Enzyme cleavage kinetics and partition coefficients %%%
% Tumor-SE3
kcat_tumor = 2.267 *60e-3;                  % Tumor esterase cleavage kcat (measured).
Km_tumor = 0.5201;                          % Tumor esterase cleavage Km (measured).
E_tumor = 1;                                % Tumor esterase content, measured with Kcat that Kat*content = 2.267 in 10 mg/mL tissue lysis buffer.

% Liver-SE3
kcat_liver = 1.288 *60e-3;                  % Liver esterase cleavage kcat (measured).
Km_liver = 9.696;                           % Liver esterase cleavage Km (measured).
E_liver = 50;                               % Liver esterase content, measured with Kcat that Kat*content = 1.288 in 10 mg/mL tissue lysis buffer.


% Tumor-SE1
% kcat_tumor = 0.8087 *60e-3;               % Tumor esterase cleavage kcat (measured).
% Km_tumor = 3.272;                         % Tumor esterase cleavage Km (measured).
% E_tumor = 1;                              % Tumor esterase content, measured with Kcat that Kat*content = 2.267 in 10 mg/mL tissue lysis buffer.
 
% Liver-SE1
% kcat_liver = 0.4873 *60e-3;               % Liver esterase cleavage kcat (measured).
% Km_liver = 6.494;                         % Liver esterase cleavage Km (measured).
% E_liver = 50;                             % Liver esterase content, measured with Kcat that Kat*content = 1.288 in 10 mg/mL tissue lysis buffer.


% Tumor-SE2
% kcat_tumor = 1.457 *60e-3;                % Tumor esterase cleavage kcat (measured).
% Km_tumor = 1.600;                         % Tumor esterase cleavage Km (measured).
% E_tumor = 1;                              % Tumor esterase content, measured with Kcat that Kat*content = 2.267 in 10 mg/mL tissue lysis buffer.

% Liver-SE2
% kcat_liver = 0.6087 *60e-3;               % Liver esterase cleavage kcat (measured).
% Km_liver = 6.773;                         % Liver esterase cleavage Km (measured).
% E_liver = 50;                             % Liver esterase content, measured with Kcat that Kat*content = 1.288 in 10 mg/mL tissue lysis buffer.


% Tumor-SE4
% kcat_tumor = 2.706 *60e-3;                % Tumor esterase cleavage kcat (measured).
% Km_tumor = 2.887;                         % Tumor esterase cleavage Km (measured).
% E_tumor = 1;                              % Tumor esterase content, measured with Kcat that Kat*content = 2.267 in 10 mg/mL tissue lysis buffer.

% Liver-SE4
% kcat_liver = 1.019 *60e-3;                % Liver esterase cleavage kcat (measured).
% Km_liver = 5.384;                         % Liver esterase cleavage Km (measured).
% E_liver = 50;                             % Liver esterase content, measured with Kcat that Kat*content = 1.288 in 10 mg/mL tissue lysis buffer.
               

B_tumor = 0.597;                            % Tumor targeting of nanoparticles (calculated with measurements).
%B_tumor = memratio(M_ratio);               % When calculating different cell membrane contents, SEDOM has different targeting coefficients for tumors, and the coefficient could be calculated by this function.
B_liver = 1-B_tumor;                        % Liver targeting of nanoparticles (calculated with measurements).

%Plasma-N.S.
kcat_plasma = kcat_liver/28;                % Plasma esterase cleavage kcat (manually fitted).
Km_plasma = Km_liver*12;                    % Plasma esterase cleavage Km (manually fitted).
E_plasma = 10;                              % Plasma esterase content (manually fitted).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% State Conditions %%%
C_NP_plasma = y(1);                         % Concentration of nanoparticle in plasma
C_NP_tumor = y(2);                          % Concentration of nanoparticle in tumor
C_NP_liver = y(3);                          % Concnetration of nanoparticle in liver
C_R_tumor = y(4);                           % Concentration of reporter in tumor
C_R_liver = y(5);                           % Concentration of reporter in liver
C_R_plasma = y(6);                          % Concentration of reporter in plasma
C_R_urine = y(7);                           % Concentration of reporter in urine

%%% Equation S1 %%%
dC_NP_plasmadt = - l_NP_liClear*C_NP_plasma - k_NP_tumor*B_tumor*(C_NP_plasma-C_NP_tumor) - k_NP_liver*B_liver*(C_NP_plasma-C_NP_liver) - (kcat_plasma*E_plasma*C_NP_plasma)/(Km_plasma+C_NP_plasma);

%%% Equation S2 %%%
dC_NP_tumordt = k_NP_tumor*B_tumor*(C_NP_plasma-C_NP_tumor) - (kcat_tumor*E_tumor*C_NP_tumor)/(Km_tumor+C_NP_tumor);

%%% Equation S3 %%%
dC_NP_liverdt = k_NP_liver*B_liver*(C_NP_plasma-C_NP_liver) - (kcat_liver*E_liver*C_NP_liver)/(Km_liver+C_NP_liver);

%%% Equation S4 %%%
dC_R_tumordt = (kcat_tumor*E_tumor*C_NP_tumor)/(Km_tumor+C_NP_tumor) - k_R_tumor*(C_R_tumor-C_R_plasma);

%%% Equation S5 %%%
dC_R_liverdt = (kcat_liver*E_liver*C_NP_liver)/(Km_liver+C_NP_liver) - k_R_liver*(C_R_liver-C_R_plasma);

%%% Equation S6 %%%
dC_R_plasmadt = k_R_tumor*(C_R_tumor-C_R_plasma) + k_R_liver*(C_R_liver-C_R_plasma) - l_R_kiClear*C_R_plasma + (kcat_plasma*E_plasma*C_NP_plasma)/(Km_plasma+C_NP_plasma);

%%% Equation S7 %%%
dC_R_urinedt = l_R_kiClear*C_R_plasma - a_urine/V_urine*C_R_urine;

dydt = [dC_NP_plasmadt; dC_NP_tumordt; dC_NP_liverdt; dC_R_tumordt; dC_R_liverdt; dC_R_plasmadt; dC_R_urinedt];
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = SEconc(x)                      %Correspondence between macrophage membrane content and sucrose ester content
    y = 21.55/0.9*(1-x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = memratio(x)                    %Computational fit to in vivo experimental results
yb = 0.34841;
A1 = 0.2258;
A2 = 0.2258;
Tau1 = 0.10485;
Tau2 = 0.10485;
y = yb+A1*(1-exp(-x/Tau1))++A2*(1-exp(-x/Tau2));
return 
end