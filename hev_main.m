% load driving cycle ### JN1015 load
load NEDC_input

% create grid
clear grd
grd.Nx{1}    = 1001; 
grd.Xn{1}.hi = 0.9; 
grd.Xn{1}.lo = 0.2;

grd.Nu{1}    = 1001;
grd.Un{1}.hi = 1;
grd.Un{1}.lo = -1;	% Att: Lower bound may vary with engine size.

%engine start-stop function     %increase the computation effort    %even it can be done with conditions defining the engine-on or -off, which was done in my last code
grd.Nu{2} = 2;
grd.Un{2}.hi = 1;
grd.Un{2}.lo = 0;

% set initial state
grd.X0{1} = 0.50;           %intial battery charge

% final state constraints
grd.XN{1}.hi = 0.51;
grd.XN{1}.lo = 0.50;

% define problem
clear prb
prb.W{1} = w_MGB_NEDC; % (661 elements)
prb.W{2} = dw_MGB_NEDC; % (661 elements)
prb.W{3} = T_MGB_NEDC; % (661 elements)

prb.Ts = 1;
prb.N  = 1220*1/prb.Ts + 1;

% set options
options = dpm();
options.MyInf = 1e3;                        %simulated with different Inf value   %result do not change   %Affect if the model is not correct or cost value upbound this Inf value
options.BoundaryMethod = 'Line'; % also possible: 'none' or 'LevelSet';
if strcmp(options.BoundaryMethod,'Line') 
    %these options are only needed if 'Line' is used
    options.Iter = 5;
    options.Tol = 1e-8;
    options.FixedGrid = 0;
end
[res, dyn] = dpm(@hev,[],grd,prb,options);

%% saving the results for NEDC cycle
save NEDC_outputs.mat res dyn

%% Torque split ratio
% Run this section to see the results in simulink

load NEDC_outputs.mat
load NEDC_input

% calculate the Split-ratio
u = res.Tm./(res.Tm + res.Te);
u(isnan(u)) = 0;            % In case of zero torque, 'u' sets to zero   

% save the split ratio
save split_NEDC_3.mat u

% Below given parameters helps to pass the split-ratio in 
% qss_hybrid_electric_vehicle_example.mdl file
u_in = u;
time = 0:1:1220; 
table1 = table(time', u_in');
file_name= 'split_done.xlsx';
writetable(table1, file_name)
save split_done.mat u_in time;  % Split-ratio stored in file with respected time
