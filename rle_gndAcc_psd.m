%
% 
%
% Apr 2025
%

clearvars

%% Script settings & flags
%%
print_figs = false;


useLogLog = false;
% Create function handle based on useLogLog
if useLogLog, compplotFunc = @loglog; %#ok<*UNRCH> 
else, compplotFunc = @semilogx;
end

% PSD settings (Welch algorithm)
psd_opts.M = 1024;
psd_opts.nFFT = 2^14;
psd_opts.detrend = 'none';


%%
% Pre allocate data structures
gnd_acc_dt = cell(7,1);
gndacc_psd = cell(7,1);

for i_rle = 1:7
    dt_file = sprintf('./data/model-20250408_1231-RLE0%d_.parquet',i_rle);
    [gnd_acc_dt{i_rle}, time] = load_gnd_acc(dt_file);

    Ts = diff(time(1:2));
    gndacc_psd{i_rle} = zeros(psd_opts.nFFT/2+1, 3);
    for ik = 1:3
        [gndacc_psd{i_rle}(:,ik),freqP] = utils.pwelch(gnd_acc_dt{i_rle}(ik,:)'...
            ,2*psd_opts.M,[],psd_opts.nFFT,1/Ts,'onesided',psd_opts.detrend);
    end
    
end




% if(print_figs)
%     print(fig1, '-dpdf', '-bestfit',...
%         sprintf('PierAcc%dpctSA_RLE0%d_comparison.pdf',100*srs.zeta,rle_id));
%     print(fig2, '-dpdf', '-bestfit',...
%         sprintf('PierAccPSD_RLE0%d_comparison.pdf',rle_id));
% end



%% GND Acc plot
%%
f_ = figure(3);
% set(f_,'defaultAxesColorOrder',1);
set(gcf,'Position',[421, 267, 750, 600]);%[421, 267, 750, 420]
try axes(ha_gnd_acc(1));
catch, clear ha_gnd_acc; end
if(~exist('ha_gnd_acc','var'))
    [ha_gnd_acc, ~] = utils.tight_subplot(3,1,[.06 .07],[.09 .06],[.07 .02]);
end

plotFontSize = 10;
dir_v = ["X^{\rightarrow}","Y^{\rightarrow}","Z^{\rightarrow}"];

for i_ = 1:3
    axes(ha_gnd_acc(i_)); %#ok<*LAXES>
    for i_rle = 1:7
        semilogx(freqP, gndacc_psd{i_rle}(:,i_));
        if (i_rle == 1), hold on; end
    end
    
    axis tight;
    ylabel(sprintf('Acc {%s} PSD (m^2/s^4/Hz)',dir_v(i_)),'FontSize',plotFontSize);
    if(i_ == 1), title('GND Acc PSD'); end
    if(i_ == 3), xlabel('Frequency (Hz)'); end
    grid on; xlim([0.5,min(90,max(freqP))]);
    hold off;
end


%% Get pier acceleration from displacements in a parquet file 
%%
function [gnd_acc, time] = load_gnd_acc(dt_file)
try
    parquetINFO = parquetinfo(dt_file);
    sssha_data = parquetread(dt_file,"SampleRate",1e3,...
        "SelectedVariableNames",parquetINFO.VariableNames);
    % "OSSPayloads6D";"OSS00GroundAcc";"OSSHardpointD";"OSSM1Lcl";"MountEncoders"

    time = seconds(sssha_data.Time);

    if (any(contains(parquetINFO.VariableNames,"OSS00GroundAcc")) && 1)
        gnd_acc = reshape(cell2mat(sssha_data.OSS00GroundAcc),3,[]);
    else, gnd_acc = [];
    end
catch
    warning('Unable to run parquetread(). Try Matlab 2022b, or later.');
end

end









