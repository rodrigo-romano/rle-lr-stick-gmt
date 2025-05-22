%
% load_rle_sim_stickGMT.m
% Load data required to run simulations based on the telescope stick model.
%

% clear
simulink_fname = 'rle_sim_stickGMT';
fprintf("\n|- Loading End-to-end simulation variables -|\n|-\t  Version: %s \t-|\n",simulink_fname);

%% General settings
%%
% Flag to run all RLE case
run_all_rle = true;
% Flag to save the simulation results
save_res = true;

% % RLE case id
% rle_id = 7;

% Set telescope structural dynamics damping
rayleigh_z = true;
if ~rayleigh_z
    sysDamp = 0.02; 
end

% Set structural dynamics model sampling period
FEM_Ts = 1/1e3;

% Zenith angle string (affects the control model for now)
sZa = "30";

% - - - - - Simulation setting flags (0:disables) - - - - -
clear osim
osim.reduce_model = 0;  % DO NOT ENABLE UNTIL TODO IS DONE![bool] model reduction feature
osim.dc_mm_comp = 0;    % [bool] DC mismatch compensation


%% Load telescope structural dynamics model
%%


% ModelFolder = "20241125_0627_pier_ODC_202411_stickTelescope_largeMass";
ModelFolder = "20250331_0753_pier_ODC_202411_stickTelescope_SIIDOM_HBS_largeMass";
% ModelFolder = "20250331_0825_pier_ODC_202411_stickTelescope_SIIDOM_EDS_largeMass";
% ModelFolder = "20250519_1258_pier_ODC_202411_stickTelescope_IDOM_verification";
FileName = "modal_state_space_model_2ndOrder.mat";
    
if(~exist('inputTable','var') || 1)
    load(fullfile(im.lfFolder,ModelFolder,FileName),'inputs2ModalF','modalDisp2Outputs',...
        'eigenfrequencies','inputTable','outputTable');

    fprintf('Model %s loaded from\n%s\n', FileName, ModelFolder);
    if osim.dc_mm_comp
        % Static solution gain matrix
        staticSolFile = fullfile(im.lfFolder,ModelFolder,"static_reduction_model.mat");
        try
            load(staticSolFile,'gainMatrixMountControlled');
            gainMatrix = gainMatrixMountControlled;
        catch
            load(staticSolFile,'gainMatrix');
            warning("Loading static solution model gainMatrix.");
        end
        fprintf('Static gain matrix loaded from\n%s\n', staticSolFile);
    end
end


%% Pick model IOs according to inputTable and outputTables
%%

% INPUTS
desiredInputLabels = ["Pier_6F"; "OSS00_Ground_Acc"];
isDesired = zeros(size(inputs2ModalF,2),1);
modelMuxDims = zeros(numel(desiredInputLabels),1);

for i1 = 1:numel(desiredInputLabels)
    aux = inputTable{desiredInputLabels{i1},"indices"}{1}(:);
    isDesired(aux) = 1;
    modelMuxDims(i1) = length(aux);
end
indDesInputs = find(isDesired ~= 0);
modelMuxDims(modelMuxDims == 0) = [];

% OUTPUTS
desiredOutputLabels = ["Pier_6D";...
    "Stick_telescope_6D";...
    "OSS00_Ground_6D";...
    "Pier_SI_6D"];

isDesired = zeros(size(modalDisp2Outputs,1),1);
modelDemuxDims = zeros(numel(desiredOutputLabels),1);

for i1 = 1:numel(desiredOutputLabels)
    aux = outputTable{desiredOutputLabels(i1),"indices"}{1}(:);
    isDesired(aux) = 1;
    modelDemuxDims(i1) = length(aux);
end
indDesOutputs = find(isDesired ~= 0);
modelDemuxDims(modelDemuxDims == 0) = [];



%% Structural model discretization
%%

if osim.reduce_model
    % Compute the approximate Hankel singular values
    % <-- TO DO: Always keep the first 3 modes
    [gamma,~] = utils.approxhsv(sqrt(om2), sysDamp, phiB, phiC);
    [~,si] = sort(gamma,'descend');
    th1 = 1e-7;
    gammamax = max(gamma(:,1+3));
    
    nr = length(find((gamma./gammamax) >= th1));
    warning('\n-> The number of modes for TH=%.2g is %d\n',th1,nr);
    mode_ind_vec = si(1:nr);
else
    nr = length(eigenfrequencies(:));
    mode_ind_vec = 1:nr;
end

% om^2 and 2*zeta*om vectors 
om2 = (2*pi*eigenfrequencies(mode_ind_vec)).^2;

if rayleigh_z
    alpha_z = 0.1555;
    beta_z = 6.303e-4;
    sysDamp = (alpha_z/2)./sqrt(om2) + 0.5*(beta_z)*sqrt(om2);
    sysDamp(1:3) = sysDamp(4);
    fprintf('** Structural model damping set according to Rayleigh curve.\n');
    if false
        figure(222);
        semilogx(2*pi./sqrt(om2),sysDamp,'.'); hold on;
        om_ = 2*pi*linspace(1/10,1/0.01,1001);
        ray_zeta = (alpha_z/2)./om_ + (beta_z/2)*om_;
        semilogx(2*pi./om_, ray_zeta,'--'); hold off;
        xlim([2*pi/max(om_),2*pi/min(om_)]);
        grid on;
        xlabel('Period (s)');
        ylabel('Damping (adim)');
    end
else
    sysDamp = sysDamp(1)* ones(size(om2)); %#ok<*UNRCH>
    fprintf('Plant model damping ratio set to:%.2g\n', sysDamp(1));
end

% Perform discretization and provide the 2nd order form DT simulation parameters
twice_zom = 2*sysDamp .* sqrt(om2);
PG = zeros(length(om2),6);
for i = 1:length(om2)
    PhiGamma = expm([0 1 0; -om2(i) -twice_zom(i) 1; 0 0 0]*FEM_Ts);
    PG(i,:) = [PhiGamma(1,1:3), PhiGamma(2,1:3)];
end

phiB = inputs2ModalF(mode_ind_vec,indDesInputs);
phiC = modalDisp2Outputs(indDesOutputs,mode_ind_vec);

s_rate_msg = 'Telescope structural model sampling rate set to %gkHz\n';
if((FEM_Ts ~= 1/8e3) && (FEM_Ts ~= 1e-3)), warning(s_rate_msg,1/FEM_Ts/1e3);
else, fprintf(s_rate_msg,1/FEM_Ts/1e3);
end

%% Open Simulink model
%%
open_system(simulink_fname);


%% Static gain mismatch compensation
%%

StaticModelFolder = fullfile(im.lfFolder,ModelFolder);
staticFileName = fullfile(StaticModelFolder,"static_reduction_model.mat");

struct_dyn_label = "/Telescope model/Structural Dynamics GMT";
memory_label = simulink_fname+struct_dyn_label+"/Psi_ss_memory";
matmult_label = simulink_fname+struct_dyn_label+"/Psi_ss";
zoh_label = simulink_fname+struct_dyn_label+"/Psi_ssZOH";
rt_label = simulink_fname+struct_dyn_label+"/Psi_ssRT";

if osim.dc_mm_comp
    try
        K_ss = phiC(:,4:end)* diag(1./((2*pi*eigenfrequencies(4:end)).^2))* phiB(4:end,:);        
        Psi_ss = gainMatrix(indDesOutputs,indDesInputs) - K_ss;
        %
        num_mnt_ax_in = contains(desiredInputLabels,...
            ["OSS_ElDrive_Torque"; "OSS_AzDrive_Torque"; "OSS_RotDrive_Torque"]);
        v_in = [1; 1; 1; zeros(numel(desiredInputLabels)-3,1)];
        num_mnt_ax_out = contains(desiredOutputLabels,...
            ["OSS_ElEncoder_Angle"; "OSS_AzEncoder_Angle"; "OSS_RotEncoder_Angle"]);
        v_out = [1; 1; 1; zeros(numel(desiredOutputLabels)-3,1)];
        
        if(all(num_mnt_ax_in == v_in) && all(num_mnt_ax_out == v_out))
%             Psi_ss(1:sum(outputTable.size(1:3)),1:sum(inputTable.size(1:3))) = 0;
            Psi_ss(1:sum(outputTable.size(1:3)),:) = 0;
            Psi_ss(:,1:sum(inputTable.size(1:3))) = 0;
        else
            warning("No entries of Psi_ss set to zero!"+...
                "\nCheck for the chosen MNT IOs.\n");
        end
        
        Psi_ssTs = FEM_Ts;        
        set_param(memory_label,'Commented','off');
        set_param(matmult_label,'Commented','off');
        set_param(zoh_label,'Commented','off');
        set_param(rt_label,'Commented','off');
    catch
        warning('Unable to compute static compensation matrix.');
        warning('Disabling static compensation matrix.');
        set_param(memory_label,'Commented','on');
        set_param(matmult_label,'Commented','on');
        set_param(zoh_label,'Commented','on');
        set_param(rt_label,'Commented','on');
        StaticModelFolder = [];
    end
else
    set_param(memory_label,'Commented','on');
    set_param(matmult_label,'Commented','on');
    set_param(zoh_label,'Commented','on');
    set_param(rt_label,'Commented','on');
    StaticModelFolder = [];
end

%% Run all RLE cases
%%

if(run_all_rle)
    rle_Hz = [200, 100, 200, 80, 100, 100, 100];
    
%     rle_path_str = '/Users/rromano/Workspace/grsim/rle_acc/data/RLE0%d_1kHz_dt.mat';
    for i_rle = 1:1%7
        rle_path_str = sprintf(...
            '/Users/rromano/Workspace/grsim/rle_acc/data/RLE_2018_%d_%dHz.mat',...
            i_rle, rle_Hz(i_rle));
        % For upsampled RLE dataset
%         load(sprintf(rle_path_str,i_rle),sprintf('RLE0%d_1kHz_dt',i_rle));
%         eval(sprintf("RLE_dt = RLE0%d_1kHz_dt;",i_rle));
%         sssha_acc_dt.time = 1e-3*(0:length(RLE_dt)-1)';
        % For RLE dataset with original sample
        load(sprintf(rle_path_str,i_rle,rle_Hz(i_rle)),'acc_dt','time');
        eval("RLE_dt = acc_dt;");
        sssha_acc_dt.time = time;
        sssha_acc_dt.signals.values = [RLE_dt(:,1), RLE_dt(:,2), RLE_dt(:,3)];
        sssha_acc_dt.signals.dimensions = 3;
    
        % Run RLE simulation
        fprintf('Running RLE-%02d.\n',i_rle);
        sim(simulink_fname);
        
        if(save_res)
            % Save simulation data
            if(rayleigh_z), fname_str  = "rle0%d_%s_stickGMT_rayZ_FOH_sim_dt";
            else, fname_str  = "rle0%d_%s_stickGMT_sim_dt";
            end
            fname = sprintf(fname_str,i_rle, extractBefore(ModelFolder,14));
            save(fname,'pier6D_dt','ground6D_dt','sssha_acc_dt');
            fprintf("Done! Input date duration is %gs.\nSimulation saved in %s\n",...
                sssha_acc_dt.time(end), fname);
        end
    end
    
end

    

%% Auxiliary functions
%%
%eof

