%
% Post-process pier displacement response to RLE simulation data 
%
% May 2025
%

clearvars

%% Script settings & flags
%%
rle_id = 1;
plot_sssha_acc = true;
print_figs = false;
save_xlsx = false;

useLogLog = false;
% Create function handle based on useLogLog
if useLogLog, compplotFunc = @loglog;
else, compplotFunc = @semilogx;
end

% PSD settings (Welch algorithm)
psd_opts.M = 1024;
psd_opts.nFFT = 2^14;
psd_opts.detrend = 'none';

% Spectral acceleration seetings
srs.zeta = .02; % payload damping ratio
srs.steps = 4; % number of frequency steps per FWHM, 4 yields <3.6% scalloping
srs.fmin = 1; %  minimum frequency [Hz]

% % --- S3 link path
% s3_path = "/home/rromano/mnt";

model_set = ["20250331_0753";...
    "20250331_0825";...
    "20250519_1258"];

mset_info = [...
    "Stick_Tel/HBS config - IDOM SIS Lat stiffness";...
    "Stick_Tel/EDS config - IDOM SIS Lat stiffness";...
    "Stick_Tel/IDOM config - IDOM SIS Lat stiffness"];

mset_info_red = ["StickHBS-IDOM";"StickEDS-IDOM";"StickIDOM"];


sel_model_idx = 1:3;
model_set = model_set(sel_model_idx);
mset_info = mset_info(sel_model_idx);
mset_info_red = mset_info_red(sel_model_idx);

%%
% Pre allocate data structures
pier_acc_dt = cell(numel(model_set),1);
pier_b_Txyz_ps = cell(numel(model_set),1);
pier_sa_dt = cell(numel(model_set),1);
dt_file_label = cell(numel(model_set),1);

for i_ = 1:numel(model_set)
    dtfile_str = ['/Users/rromano/Workspace/gmt-im-studies/',...
        'seismic_pproc/rle_stickGMT/rle%02d_%s_stickGMT_rayZ_sim_dt'];
    dt_file = sprintf(dtfile_str,rle_id,model_set{i_});
    fprintf("Post-processing (%s) data from \n%s\n",mset_info_red(i_),dt_file);
    if(i_ == 1)
        [pier_acc_dt{i_}, time, gnd_acc] = load_pier_acc(dt_file);
        Ts = diff(time(1:2));
        fSRS = get_SRS_fSRS(srs.zeta, srs.steps, srs.fmin, Ts);
    else        
        [pier_acc_dt{i_}, time_] = load_pier_acc(dt_file);
        assert(length(time)==length(time_),...
            'Incompatible data length for model %d',i_);
    end
    label_struct = replace(...
        extractBetween(dt_file,sprintf("rle%02d_",rle_id),"_stickGMT_rayZ"),...
        '_','\_');
    dt_file_label{i_} = label_struct{1};
    
    % Acc PSD
%     pier_b_Txyz_ps{i_} = zeros(psd_opts.nFFT/2+1, 3);
    for iu=1:3
        [pier_b_Txyz_ps{i_}(:,iu),freqP] = utils.pwelch(pier_acc_dt{i_}(:,iu),...
            2*psd_opts.M,[],psd_opts.nFFT,1/Ts,'onesided',psd_opts.detrend);
    end

    % SRS
    pier_sa_dt{i_}(:,1) = fSRS;
    pier_sa_dt{i_}(:,2) = SpectralA04(pier_acc_dt{i_}(:,1), fSRS, Ts, srs.zeta);
    pier_sa_dt{i_}(:,3) = SpectralA04(pier_acc_dt{i_}(:,2), fSRS, Ts, srs.zeta);
    pier_sa_dt{i_}(:,4) = SpectralA04(pier_acc_dt{i_}(:,3), fSRS, Ts, srs.zeta);
end

%% Save data in a spreadsheet
%%

if(save_xlsx)
    for i_ = 1:numel(model_set) %#ok<*UNRCH> 
        header = {sprintf('%s - model_id:%s',mset_info{i_},model_set(i_))};
        dt_fname = sprintf('pierAcc_compDT_RLE0%d_IMsim2025Apr.xlsx',rle_id);
        writecell(header,dt_fname,'Sheet',mset_info_red{i_},'Range','B2');
        fprintf('Saving data from %s\nto spreadsheet %s\n',header{1},dt_fname);

        description = {'Pier acceleration (time-domain data)',[],[],[];...
            'Time (s)','H1 (m/s^2)','H2 (m/s^2)','V (m/s^2)'};
        writecell(description,dt_fname,'Sheet',mset_info_red{i_},'Range','B4');
        writematrix([time, pier_acc_dt{i_}],dt_fname,'Sheet',mset_info_red{i_},'Range','B6');
        %
        description = {sprintf(...
            'Spectral Acceleration (SA) - Payload damping:%g',srs.zeta),[],[],[];...
            'Frequency (Hz)','H1 (m/s^2)','H2 (m/s^2)','V (m/s^2)'};
        writecell(description,dt_fname,'Sheet',mset_info_red{i_},'Range','G4');
        writematrix(pier_sa_dt{i_},dt_fname,'Sheet',mset_info_red{i_},'Range','G6');
        %
        description = {'Pier ACC Power Spectrum Density (PSD)',[],[],[];...
            'Frequency (Hz)','H1 (m^2/s^4/Hz)','H2 (m^2/s^4/Hz)','V (m^2/s^4/Hz)'};
        writecell(description,dt_fname,'Sheet',mset_info_red{i_},'Range','M4');
        writematrix([freqP, pier_b_Txyz_ps{i_}],dt_fname,'Sheet',mset_info_red{i_},'Range','M6');
    end
end

%% IDOM SA data
%%

load(fullfile(...
    '/Users/rromano/Workspace/gmt-im-studies/seismic_pproc/rle_stickGMT',...
    'P102482-LIS-1R1-RLE_Acc.mat'),'pier_rle');
% pier_rle{i_rle}.time
% pier_rle{i_rle}.acc ::: in [g]

idom_sa_data = zeros(length(fSRS), 3);
idom_dT = diff(pier_rle{rle_id}.time(1:2));
for j = 1:3
    idom_sa_data(:,j) = SpectralA04(...
        pier_rle{rle_id}.acc(:,j), fSRS, idom_dT, srs.zeta);
end

%% Fixed damping pier Acc data
%%

if(~exist('pier_acc_cteZ','var') || 0)
    xls_fname = sprintf('pierAcc_compDT_RLE0%d_IMsim2025Apr.xlsx',rle_id);
    sheet = 'StickHBS-IDOM';
    
    xls_cols = 'G:J';
    
    fprintf('Importing RLE%0d pier motion data from\n%s\n', rle_id, xls_fname);
    xls_opts = detectImportOptions(...
        xls_fname,'Sheet',sheet,'Range',xls_cols);
    M = readmatrix(xls_fname,xls_opts);
    dt = M(isfinite(M(:, 1)),:);
    % fprintf('Time series sampled at %gHz.\n', 1/diff(dt(1:2,1)));
    pier_acc_cteZ.fSRS = dt(:,1);
    pier_acc_cteZ.acc = dt(:,2:4);
end



%% Ground ACC
%%
% gndacc_ps = zeros(psd_opts.nFFT/2+1, 3);
% for ik = 1:3
%     [gndacc_ps(:,ik),freqP] = utils.pwelch(gnd_acc(ik,:)'...
%         ,2*psd_opts.M,[],psd_opts.nFFT,1/Ts,'onesided',psd_opts.detrend);
% end
% gnd_sa_data = zeros(length(fSRS), min(size(pier_acc_dt))+1);
% for j = 1:3
%     gnd_sa_data(:,j) = SpectralA04(gnd_acc(j,:)', fSRS, Ts, srs.zeta);
% end


%% Analysis Plots
%%

t_range = [0.005, max(time)];%[];%[4.19,4.208];%[];%
if(isempty(t_range))
    t_idx  = [1, length(time)];
else
    t_idx = [find(time >= t_range(1),1,"first"),find(time < t_range(2),1,"last")];
end

ylabel_str = ["H1","H2","V"];
lw_vec = [1.2, 1.2, 1.4, 1.4, 1, 1];
ltype_vec = ['-', '--', '-', '--', '-', '--'];

% Spectral Acceleration Plot  
fig1 = figure(20000-rle_id);
set(gcf,'position',[423   350   740   400])
for ik = 1:3    
    subplot(3,1,ik)
%     plot(fSRS,gnd_sa_data(:,ik),':','LineWidth',1.5,'Color',[.3 .3 .3]);
%     hold on;
    for i_ = 1:numel(model_set)
        plot(fSRS,pier_sa_dt{i_}(:,ik+1),ltype_vec(sel_model_idx(i_)),...
            'Linewidth',lw_vec(sel_model_idx(i_)));
        hold on;
    end
    plot(fSRS, 9.80665*idom_sa_data(:,ik),'-','LineWidth',1.4);
    plot(pier_acc_cteZ.fSRS, pier_acc_cteZ.acc(:,ik),'--','LineWidth',1.2);
    xlim([0, 60])
    grid on; ylabel(ylabel_str{ik}+" (m/s^2)"); hold off;
end
xlabel("Frequency (Hz)")
subplot(3,1,1);
title(sprintf("RLE%d - SA response (payload damping ratio:%g)",...
    rle_id, srs.zeta));
% legend(["GND Acc"; mset_info_red(:)],'NumColumns',2);
legend([mset_info_red(:); "IDOM 1R1 Acc"; "HBS-IDOM(fixed zeta=0.02)"],'NumColumns',1);


% Acc PSD Plot
fig2 = figure(30000-rle_id);
set(gcf,'position',[423   150   740   400])
ylabel_str = ["H1 PSD","H2 PSD","V PSD"];
for ik = 1:3    
    subplot(3,1,ik)
%     compplotFunc(freqP, gndacc_ps(:,ik),':','LineWidth',1.5,'Color',[.3 .3 .3]);
%     hold on;
    for i_ = 1:numel(model_set)
        compplotFunc(freqP, pier_b_Txyz_ps{i_}(:,ik),ltype_vec(sel_model_idx(i_)),...
            'Linewidth',lw_vec(sel_model_idx(i_)));
        hold on;
    end
    
    xlim([0.5, 50])
    grid on; ylabel(ylabel_str{ik}+" (m^2/s^4/Hz)"); hold off;
end
xlabel("Frequency (Hz)")
subplot(3,1,1);
title(sprintf("RLE%d - Bottom pier Acc PSD", rle_id));
% legend(["GND Acc"; mset_info_red(:)],'Location','northwest','NumColumns',2);
legend(mset_info_red(:),'Location','northwest','NumColumns',1);
legend box off

if(print_figs)
    print(fig1, '-dpng',...
        sprintf('RLE0%d_comparison.pdf',100*srs.zeta,rle_id));
    print(fig2, '-dpng',...
        sprintf('RLE0%d_comparison.pdf',rle_id));
end

return

%% Verification plots
%%
% Pier ACC Plot
fid_offset = 100;%200;
pier_out_dt = pier_acc_dt{1};
y_label1 = 'Pier Acc time response (%s/s^2)';
y_label2 = 'Pier Acc PSD (%s^2/s^4/Hz)';
plot_title_str = "Bottom Pier Node Acceleration";

legend_str = {'X^{\rightarrow}','Y^{\rightarrow}','Z^{\rightarrow}'};
% TRANSLATIONS
fig3 = figure(1000-fid_offset+rle_id);
set(gcf,'position',[123   80   740   400])
subplot(2,1,1)
plot(time(t_idx(1):t_idx(2)),pier_out_dt(t_idx(1):t_idx(2),1:3));
ylabel(sprintf(y_label1,'m'));
xlabel('Time (s)'); grid on; axis tight;
legend(legend_str);
title(plot_title_str)
subplot(2,1,2)
semilogx(freqP, pier_b_Txyz_ps{1});
xlabel('Frequency (Hz)'); ylabel(sprintf(y_label2,'m')); grid on; axis tight;

if(print_figs)
    print(fig3, '-dpdf', '-bestfit',...
        sprintf('PierAcc_RLE0%d_comparison.pdf',rle_id));
end

% ROTATIONS    
%     figure(000)
%     set(gcf,'position',[423   150   740   400])
%     subplot(2,1,1)
%     plot(t(t_idx(1):t_idx(2)),pier_out_dt(:,4:6));
%     ylabel(sprintf(y_label1,'rad'));
%     xlabel('Time (s)'); grid on; axis tight;
%     legend(legend_str);
%     title(plot_title_str)
%     subplot(2,1,2)
%     semilogx(freqP, pier_b_Rxyz_ps);
%     xlabel('Frequency (Hz)'); ylabel(sprintf(y_label2,'rad')); grid on; axis tight;

if(plot_sssha_acc)
    figure(rle_id)
    set(gcf,'position',[123   230   400   400])
    subplot(2,1,1)
    plot(time(t_idx(1):t_idx(2)),gnd_acc(:, t_idx(1):t_idx(2))');
    set(gca,'ColorOrderIndex',1); hold on;
    plot(time([t_idx(1),t_idx(2)]), kron([1 1],mean(gnd_acc,2)), '--');
    ylabel('Accelerations (m/s^2)');
    xlabel('Time (s)');
    legend('H1','H2','V');
    title(sprintf('%s - RLE0%d',dt_file_label{1},rle_id))
    grid on; axis tight; hold off

    subplot(2,1,2)
    semilogx(freqP, gndacc_ps);
    xlabel('Frequency (Hz)');
    ylabel('GND Acc ((m/s^2)^2/Hz)');
    grid on; axis tight;
end

gnd_num_x = cumsum(Ts*cumsum(Ts*gnd_acc(:, t_idx(1):t_idx(2))'));
parquetINFO = parquetinfo(dt_file);
if (any(contains(parquetINFO.VariableNames,"OSS00Ground6D")) && 1)
    sssha_data = parquetread(dt_file,"SampleRate",1e3,...
        "SelectedVariableNames",parquetINFO.VariableNames);
    gnd_D = reshape(cell2mat(sssha_data.OSS00Ground6D),6,[]);
    account4lmk = false;
    if(account4lmk)
        gnd_num_x = gnd_num_x - 8.5294e9/2.1605e12*...
            cumsum(Ts*cumsum(Ts*gnd_D(1:3, t_idx(1):t_idx(2))'));
        warning("Accounting for the large mass spring on the GND position calculation!");
    end
else
    gnd_D = [];
end

% GND motion verification plot
figure(10+rle_id)
set(gcf,'position',[423   250   740   400])
plot(time(t_idx(1):t_idx(2)),gnd_num_x, '-');
hold on;
legend_str = {'\int\int H1 acc','\int\int H2 acc','\int\int V acc'};

if(~isempty(gnd_D))
    set(gca,'ColorOrderIndex',1);
    plot(time(t_idx(1):t_idx(2)), gnd_D(1:3, t_idx(1):t_idx(2))',':','Linewidth',2);
    legend_str = [legend_str(:)',{'GND x motion'},{'GND y motion'},{'GND z motion'}];
    %         plot(t(t_idx(1):t_idx(2)), gnd_D(4:6, t_idx(1):t_idx(2))', '-'); % ZERO
end
if(any(contains(parquetINFO.VariableNames,"Pier6D")) && 1)
    pier_D = reshape(cell2mat(sssha_data.Pier6D),12,[]);
    plot(time(t_idx(1):t_idx(2)), pier_D(1:3, t_idx(1):t_idx(2))', '-.');
    legend_str = [legend_str(:)',{'Pier(B) X^{\rightarrow}'},...
        {'Pier(B) Y^{\rightarrow}'},{'Pier(B) Z^{\rightarrow}'}];
    set(gca,'ColorOrderIndex',4);
    plot(time(t_idx(1):t_idx(2)), pier_D((1:3)+6, t_idx(1):t_idx(2))',...
        '--','Linewidth',1.5);
    legend_str = [legend_str(:)',{'Pier(T) X^{\rightarrow}'},...
        {'Pier(T) Y^{\rightarrow}'},{'Pier(T) Z^{\rightarrow}'}];
end

ylabel('Motion (m)');
xlabel('Time (s)');
legend(legend_str,'Location','Southeast');
title(sprintf('%s - RLE0%d',dt_file_label{1},rle_id))
grid on; axis tight; hold off;


%% Get pier acceleration from displacements in a parquet file 
%%
function [pier_acc_dt,time,gnd_acc] = load_pier_acc(dt_file)
% try
%     parquetINFO = parquetinfo(dt_file);
%     sssha_data = parquetread(dt_file,"SampleRate",1e3,...
%         "SelectedVariableNames",parquetINFO.VariableNames);
%     % "OSSPayloads6D";"OSS00GroundAcc";"OSSHardpointD";"OSSM1Lcl";"MountEncoders"
% catch
%     warning('Unable to run parquetread(). Try Matlab 2022b, or later.');
% end

load(dt_file, 'pier6D_dt');
Ts = 1e-3;

time = (0:(size(pier6D_dt,1)-1))* Ts;
 

pier_D = [zeros(2,12); pier6D_dt];
pier_acc_dt = 1/Ts*diff(1/Ts*diff(pier_D(:,1:3)));

gnd_acc = [];


end

%% Function to get SRS frequency points
%%
function [fSRS] = get_SRS_fSRS(zeta, STEPS, FMIN, Ts)
% Source
% [SRS_STRUCT] = compute_response_spectra (EQ_STRUCT, ZETA)
% B. Smith, 18 Sept 2024
% - zeta: payload damping ratio
% - STEPS: number of frequency steps per FWHM, 4 yields <3.6% scalloping
% - FMIN:  minimum frequency [Hz]

fs = 1/Ts; % sampling frequency
fmax = fs/2;
q = 1/(2*zeta); % resonance Q factor, Q = f/FWHM

nFreqs = round(log(fmax/FMIN)/log(1+1/(q*STEPS)));
fSRS = logspace(log10(FMIN),log10(fmax),nFreqs); % frequency vector [Hz]

end

%% Spectral Acceleration response
%%
function [Sa] = SpectralA04(InputSignal,f,dt,d)
%   SpectralA04(InputSignal,f,dt,d)
%   delivers a vector of the spectral acceleration response of the signals
%   represented in the InputSignal variable.
%   f, corresponds to the frequency vector at which the Spectral [Hz]
%   acceleration is evaluated.
%   dt corresponds to the time step of the input signal [s]
%   d corresponds to the damping ratio at which the spectral response need
%   to be evaluated
%
% AO4 - updated to max abs
    ns=min(size(InputSignal)); % number of signals
    Sa=zeros(length(f),ns); % preallocate result variable
    for i1=1:length(f)

      % continuous transfer function for harmonic oscillator
      filt1=tf([f(i1)*4*pi*d (f(i1)*2*pi)^2],[1 d*4*f(i1)*pi (f(i1)*2*pi)^2]);

      % convert continuous TF to z domain
      % FOH seems to work best per from SpectralA04_check.m
%     filt1d=c2d(filt1,dt,'zoh');
%     filt1d=c2d(filt1,dt,'prewarp',f(i1)*2*pi);
      filt1d=c2d(filt1,dt,'foh');
%     temp=get(filt1d);  %lists structure fieldnames
%      b=filt1d.num{1}; % extract coeffs (Octave compatible)
%      a=filt1d.den{1};

      [b, a] = tfdata(filt1d,'v'); % extract coeffs (Octave/Matlab compatible)

      for i2=1:ns  % filter signals
        Sa(i1,i2)=max(abs(filter(b,a,InputSignal(:,i2))));
      end

    end
end








