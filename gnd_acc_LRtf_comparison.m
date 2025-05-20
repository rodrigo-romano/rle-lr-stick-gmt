% - - - - - -
%
% gnd_acc_LRtf_comparisons
%
% Script to compare transfer functions of different models.
%
% Rodrigo A. Romano
%
% - - - - - -

rad2mas = (180/pi * 3600 * 1000);

%% General analysis settings
%%

% clear
% close all
axisLabels = {'AZ','EL','GIR'};
% Plot font size
plotFontSize = 12;
% Bode and Nichols plot option structures
[hbode, hnichols] = get_fr_plot_opts(plotFontSize);
% Figure number to avoid overlapping
figNumber = 100;
% Flag to include requirements in the plots
plot_reqs = true;
% Flag to save plots as EPS files
save_plots = false;

% Number of frequency response points
Nom = 1000;
% Maximum strucutral dynamics model frequency
max_model_Fr = [];
% Frequency range for analysis
wrange = [0.1,100];      
% Frequency points (rad/s)
omega = logspace(log10(wrange(1)),log10(wrange(2)),Nom)'*2*pi;


%% Telescope structural models
%%
clear model1 model2 model3

% model1 = "20250306_2131_zen_30_M1_202110_FSM_202305_Mount_202305_pier_202411_heavy_SIReq_movingGround_lockedDrives";
model1 = "20250306_2132_zen_30_M1_202110_FSM_202305_Mount_202305_pier_202411_heavy_SIIDOM_movingGround_lockedDrives";
model2 = "20250331_0753_pier_ODC_202411_stickTelescope_SIIDOM_HBS_largeMass";
% model2 = "20230828_2322_zen_30_M1_202110_FSM_202305_Mount_202305_IDOM_concreteAndFoundation_largeMass_softSI";
% model2 = "20230817_1808_zen_30_M1_202110_FSM_202305_Mount_202305_IDOM_concreteAndFoundation_largeMass";
% model3 = "20240408_1535_zen_30_M1_202110_FSM_202305_Mount_202305_IDOM_concreteAndFoundation_finalSI_largeMass";
model3 = "20250331_0825_pier_ODC_202411_stickTelescope_SIIDOM_EDS_largeMass";


modelLabels = {...
    'Locked rotor (SIIDOM)',...
    'Stick HBS-SIIDOM','Stick EDS-SIIDOM'};

% Plot legend labels
FileFolder_char{1} = char(strrep(model1,'_','\_'));
FileFolder_char{2} = char(strrep(model2,'_','\_'));
try
    FileFolder_char{3} = char(strrep(model3,'_','\_'));
catch
    FileFolder_char{3} = char('XXXXYYZZ\_WWWW');
end
legend_txt = {sprintf("%s:%s",modelLabels{1},FileFolder_char{1}(1:14)),...
    sprintf("%s:%s",modelLabels{2},FileFolder_char{2}(1:14)),...
    sprintf("%s:%s",modelLabels{3},FileFolder_char{3}(1:14))};

max_f_Hz = 90;


%% Compute plant frequency responses
%%
if(~exist('G1_jw','var') || true)
    om = cell(3,1);
    % Model 1
    [om{1}, G1_jw] = compute_freq_resp(model1, max_model_Fr, omega);
    % Model 2
    [om{2}, G2_jw] = compute_freq_resp(model2, max_model_Fr, omega);
    % Model 3    
    [om{3}, G3_jw] = compute_freq_resp(model3, max_model_Fr, omega);
%     % Model 4    
%     [om{4}, G4_jw] = compute_freq_resp(model4, max_model_Fr, omega);
    
end
n_u = size(G1_jw.u,1);
n_y = size(G1_jw.y,1);

%% Incorporate mount axis feedback loop effect
%%

% Load ODC controllers
o = {load_ODC_control('30', 'OO', '', ''),...
    load_ODC_control('30', 'OO', '', ''),...
    load_ODC_control('30', 'OO', '', '')};%load_ODC_control('30', '20', '11', '02')};

%%

if(exist('G3_jw','var')), n_models = 3;
% Just for pier TF comparison
else, n_models = 2;
end

H = cell(n_models,1);
Hacc = cell(n_models,1);
S_mnt = H; L_mnt = H; Gmf_mnt = H;
for ii=1:numel(H)
    % Controller+filter frequency response (FR)
    if(0)
        azF17 = notchF(16.9, 3.5, 1.8);
        C_(1,1) = azF17; C_(2,2) = tf(1,1); C_(3,3) = tf(1,1);
        mntCfb_fr = frd(C_*o{ii}.Hp, omega);
        warning('Updated FB controller!') %#ok<*UNRCH>
    else, mntCfb_fr = frd(o{ii}.Hp, omega);
    end
    
    Mratio = eye(3);
    if(strcmp(eval(['FileFolder_char{',num2str(ii),'}(1:14)']),'20230522\_1804') && (ij==3))
        warning('Using *G_adjust* for %s TF',axisLabels{ij});
        Mratio(3,3) = 1.4791;
    end
    mntCff_fr = frd(o{ii}.Hff, omega);
    drive_fr = frd(o{ii}.Hdrive, omega);
    delay_fr = eye* frd(pade(tf(1,1,'IOdelay',4.0e-3), 10), omega);
    
    % Sensitivity function
    L_mnt{ii} = eval(['G',num2str(ii),'_jw((1:3)+(n_y-3),(1:3)+(n_u-3))']) * Mratio *...
        drive_fr*delay_fr*mntCfb_fr;
    S_mnt{ii} = eye(3)/(eye(3)+L_mnt{ii});
    
    if(true)    % Locked-rotor model?
        H{ii} = eval(['G',num2str(ii),'_jw']);
    else
        H{ii} = feedback(eval(['G',num2str(ii),'_jw']),...
            Mratio *drive_fr * delay_fr * mntCfb_fr, (1:3)+(n_u-3), (1:3)+(n_y-3));
    end
    Hacc{ii} = frd(tf([1 0 0],1),omega)* eval(['G',num2str(ii),'_jw']);
end


%% Reorder colormap to match the profile used in the analysis with grounded models
%%
hf_ = figure(99);
map_c = get(gca,'ColorOrder');
map_c = map_c(1:end,:); % Discard 0 first colors.
% set(gca, 'ColorOrder',map_c, 'NextPlot','ReplaceChildren');
close(hf_);

%% Model eigenfrequencies
%%
f_eigenfreq = figure(1);
set(f_eigenfreq,'defaultAxesColorOrder',map_c);
set(gcf,'Position',[421+40*0, 267, 750, 420]);

semilogx(om{1}/2/pi,'+'); hold on; grid on;
semilogx(om{2}/2/pi,'o');
semilogx(om{3}/2/pi,'.');
hold off

max_mode = 200; % Number of modes for the comparison
ylim([0,100]);
xlim([1,max_mode])
title(sprintf('Comparison of the first %d eigenfrequencies',max_mode));
ylabel('Mode # eigenfrequency (Hz)','FontSize',plotFontSize);
xlabel('Model mode #','FontSize',plotFontSize);
legend(legend_txt,'Location','northeast','FontSize',10); legend box off;    
%
m_range = 3:18;
x1 = m_range(1);
y1 = om{2}(m_range(1))/2/pi;
dx = m_range(end)-m_range(1);
dy = abs(om{2}(m_range(end))-om{2}(m_range(1)))/2/pi;

r = rectangle('Position',[x1,y1,dx,dy]);
r.LineStyle = '--';
% create a new pair of axes inside current figure
axes('position',[.18 .42 .33 .48])
box on % put box around new pair of axes
semilogx(m_range, om{1}(m_range)/2/pi,'+'); hold on; grid on
semilogx(m_range, om{2}(m_range)/2/pi,'o');
semilogx(m_range, om{3}(m_range)/2/pi,'.');
xlim([x1, x1+dx])
ylim([y1, y1+dy]);
hold off; 


%%
% f1 = figure(figNumber+1);
% set(f1,'defaultAxesColorOrder',map_c);
% set(gcf,'Position',[421+40*1, 267-60*1, 750, 340]);
% % addpath('/Users/rromano/Workspace/im-tools');
% addpath('/home/rromano/Workspace/gmt-im-tools');
% try axes(ha_3ax(1));
% catch, clear ha_3ax; end
% if(~exist('ha_3ax','var'))
%     [ha_3ax, ~] = utils.tight_subplot(1,3,[.05 .07],[.11 .06],[.07 .02]);
% end
% 
% for ii=1:3
%     axes(ha_3ax(ii)); %#ok<*LAXES> 
%     
%     G1_i = squeeze(G1_jw.ResponseData(ii+(n_y-3),ii+(n_u-3),:));
%     semilogx(omega/2/pi, 20*log10(abs(G1_i)),'-','Linewidth',1.5);
%     hold on;
%     G2_i = squeeze(G2_jw.ResponseData(ii+(n_y-3),ii+(n_u-3),:));
%     if(strcmp(FileFolder_char{2}(1:14),'20230522\_1804') && (ii==3))
%         warning('Using *G_adjust* for %s TF',axisLabels{ii});
%         G2_i = G2_i*1.4791;
%     end
%     semilogx(omega/2/pi, 20*log10(abs(G2_i)),'--','Linewidth',1);
%     if(n_models > 2)
%         G3_i = squeeze(G3_jw.ResponseData(ii+(n_y-3),ii+(n_u-3),:));
%         semilogx(omega/2/pi, 20*log10(abs(G3_i)),'-.','Linewidth',1);
%     end
%     axis tight;
%     
%     if(ii == 2), legend(legend_txt,'Location','Northeast'); legend box off; end 
%     xlim([0.8,max_f_Hz]);
%     title(sprintf('%s Axis',axisLabels{ii}),'FontSize',plotFontSize);
%     ylabel(sprintf('Avg %s ENC (dB)',axisLabels{ii}),'FontSize',plotFontSize);
%     xlabel('Frequency (Hz)','FontSize',plotFontSize);
%     
%     grid on; hold off;
% 
% end


%% GND_ACC -> Top pier node displacements
%%
f2b = figure(figNumber+20);
set(f2b,'defaultAxesColorOrder',map_c);
set(gcf,'Position',[421, 267, 750, 420]);
try axes(ha_gnd_acc_pier6D(1));
catch, clear ha_gnd_acc_pier6D; end
if(~exist('ha_gnd_acc_pier6D','var'))
    [ha_gnd_acc_pier6D, ~] = utils.tight_subplot(2,2,[.06 .07],[.09 .06],[.07 .02]);
end

ind_tpier_Txyz = (7:9);
ind_tpier_Rot = (10:12);
ind_gnd_acc = (1:3);

axes(ha_gnd_acc_pier6D(1)); %#ok<*LAXES>
H1_i = 1e6*squeeze(sqrt(sum(H{1}.ResponseData(ind_tpier_Txyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H1_i),'-','Linewidth',1.5);
hold on;
H2_i = 1e6*squeeze(sqrt(sum(H{2}.ResponseData(ind_tpier_Txyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H2_i),'--','Linewidth',1);
H3_i = 1e6*squeeze(sqrt(sum(H{3}.ResponseData(ind_tpier_Txyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H3_i),'-.','Linewidth',1);
axis tight;
legend(legend_txt,'Location','Southwest'); legend box off;
xlim([0.1,max_f_Hz]);
ylabel('Top pier T_{xyz} RSS (um/m/s^2)','FontSize',plotFontSize);
title('Top pier node translations induced by GND Accelerations');
grid on; hold off;

axes(ha_gnd_acc_pier6D(3));
cH1_i=cumsum(diff(omega).*abs(H1_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH1_i,'-','Linewidth',1.5); hold on;
cH2_i=cumsum(diff(omega).*abs(H2_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH2_i,'--','Linewidth',1);
cH3_i=cumsum(diff(omega).*abs(H3_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH3_i,'-.','Linewidth',1);
xlim([0.1,max_f_Hz]);
ylabel('Cumulative Pier T_{xyz} RSS (um/m/s^2)','FontSize',plotFontSize);
xlabel('Frequency (Hz)','FontSize',plotFontSize);
grid on; hold off;

% Zoom f > XXHz 
idx_min = find(omega/2/pi > 80, 1, 'first');
idx_max = length(omega)-1;
om_range = idx_min:idx_max;
x1 = omega(idx_min)/2/pi;
y1 = min([cH1_i(om_range); cH2_i(om_range)]);
dx = omega(idx_max)/2/pi - omega(idx_min)/2/pi;
dy = abs(max([cH1_i(om_range); cH2_i(om_range)]) - y1);
ax_box_pos_size = [.26 .17 .18 .15];

r = rectangle('Position',[x1-15,0.98*y1,dx+15,8*dy]); r.LineStyle = '--';
% create a new pair of axes inside current figure
hh_ = axes('position',ax_box_pos_size);
box on % put box around new pair of axes
semilogx(hh_, omega(om_range)/2/pi, cH1_i(om_range),'-','Linewidth',1.5);
hold on; grid on
semilogx(hh_, omega(om_range)/2/pi, cH2_i(om_range),'--','Linewidth',1);
semilogx(hh_, omega(om_range)/2/pi, cH3_i(om_range),'-.','Linewidth',1);
xlim([0.99*x1, (x1+dx)])
ylim([0.98*y1, 1.02*(y1+dy)]);
hold off;


% Rot
axes(ha_gnd_acc_pier6D(2)); %#ok<*LAXES>
H1_i = 1e6*squeeze(sqrt(sum(H{1}.ResponseData(ind_tpier_Rot,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H1_i),'-','Linewidth',1.5);
hold on;
H2_i = 1e6*squeeze(sqrt(sum(H{2}.ResponseData(ind_tpier_Rot,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H2_i),'--','Linewidth',1);
H3_i = 1e6*squeeze(sqrt(sum(H{3}.ResponseData(ind_tpier_Rot,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H3_i),'-.','Linewidth',1);
axis tight;
legend(legend_txt,'Location','Southwest'); legend box off;
xlim([0.1,max_f_Hz]);
ylabel('Top pier R_{xyz} RSS (urad/m/s^2)','FontSize',plotFontSize);
title('Top pier node rotations induced by GND Accelerations');
grid on; hold off;
%
axes(ha_gnd_acc_pier6D(4));
cH1_i=cumsum(diff(omega).*abs(H1_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH1_i,'-','Linewidth',1.5); hold on;
cH2_i=cumsum(diff(omega).*abs(H2_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH2_i,'--','Linewidth',1);
cH3_i=cumsum(diff(omega).*abs(H3_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH3_i,'-.','Linewidth',1);
xlim([0.1,max_f_Hz]);
ylabel('Cumulative pier R_{xyz} RSS (urad/m/s^2)','FontSize',plotFontSize);
xlabel('Frequency (Hz)','FontSize',plotFontSize);
grid on; hold off

% Zoom f > XXHz 
idx_min = find(omega/2/pi > 85, 1, 'first');
idx_max = length(omega)-1;
om_range = idx_min:idx_max;
x1 = omega(idx_min)/2/pi;
y1 = min([cH1_i(om_range); cH2_i(om_range)]);
dx = omega(idx_max)/2/pi - omega(idx_min)/2/pi;
dy = abs(max([cH1_i(om_range); cH2_i(om_range)]) - y1);
ax_box_pos_size = [.62 .25 .18 .2];

r = rectangle('Position',[x1-10,0.98*y1,dx,1.2*dy]); r.LineStyle = '--';
% create a new pair of axes inside current figure
hh_ = axes('position',ax_box_pos_size);
box on % put box around new pair of axes
semilogx(hh_, omega(om_range)/2/pi, cH1_i(om_range),'-','Linewidth',1.5);
hold on; grid on
semilogx(hh_, omega(om_range)/2/pi, cH2_i(om_range),'--','Linewidth',1);
semilogx(hh_, omega(om_range)/2/pi, cH3_i(om_range),'-.','Linewidth',1);
xlim([0.99*x1, (x1+dx)])
ylim([0.98*y1, 1.02*(y1+dy)]);
hold off;



%% GND_ACC -> Bottom pier node accelerations
%%
f2c = figure(figNumber+21);
set(f2c,'defaultAxesColorOrder',map_c);
set(gcf,'Position',[421, 267, 750, 420]);
try axes(ha_gnd_acc_pierBacc(1));
catch, clear ha_gnd_acc_pierBacc; end
if(~exist('ha_gnd_acc_pierBacc','var'))
    [ha_gnd_acc_pierBacc, ~] = utils.tight_subplot(2,2,[.06 .07],[.09 .06],[.07 .02]);
end

ind_tpier_Txyz = (1:3);
ind_tpier_Rot = (4:6);
ind_gnd_acc = (1:3);

axes(ha_gnd_acc_pierBacc(1)); %#ok<*LAXES>
H1_i = 1e6*squeeze(sqrt(sum(Hacc{1}.ResponseData(ind_tpier_Txyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H1_i),'-','Linewidth',1.5);
hold on;
H2_i = 1e6*squeeze(sqrt(sum(Hacc{2}.ResponseData(ind_tpier_Txyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H2_i),'--','Linewidth',1);
H3_i = 1e6*squeeze(sqrt(sum(Hacc{3}.ResponseData(ind_tpier_Txyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H3_i),'-.','Linewidth',1);
axis tight;
legend(legend_txt,'Location','Southwest'); legend box off;
xlim([0.1,max_f_Hz]);
ylabel('Pier LinAcc_{xyz} RSS (adim)','FontSize',plotFontSize);
title('Bottom pier node linear Acc induced by GND Acc');
grid on; hold off;

axes(ha_gnd_acc_pierBacc(3));
cH1_i=cumsum(diff(omega).*abs(H1_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH1_i,'-','Linewidth',1.5); hold on;
cH2_i=cumsum(diff(omega).*abs(H2_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH2_i,'--','Linewidth',1);
cH3_i=cumsum(diff(omega).*abs(H3_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH3_i,'-.','Linewidth',1);
xlim([0.1,max_f_Hz]);
ylabel('Cumulative lin Acc RSS (adim)','FontSize',plotFontSize);
xlabel('Frequency (Hz)','FontSize',plotFontSize);
grid on; hold off;

% Zoom f > XXHz 
idx_min = find(omega/2/pi > 75, 1, 'first');
idx_max = length(omega)-1;
om_range = idx_min:idx_max;
x1 = omega(idx_min)/2/pi;
y1 = min([cH1_i(om_range); cH2_i(om_range); cH3_i(om_range)]);
dx = omega(idx_max)/2/pi - omega(idx_min)/2/pi;
dy = abs(max([cH1_i(om_range); cH2_i(om_range); cH3_i(om_range)]) - y1);
ylim([0, 1.02*(y1+dy)]);
ax_box_pos_size = [.13 .19 .18 .2];

r = rectangle('Position',[x1-15,0.98*y1,dx+15,1.02*dy]); r.LineStyle = '--';
% create a new pair of axes inside current figure
hh_ = axes('position',ax_box_pos_size);
box on % put box around new pair of axes
semilogx(hh_, omega(om_range)/2/pi, cH1_i(om_range),'-','Linewidth',1.5);
hold on; grid on
semilogx(hh_, omega(om_range)/2/pi, cH2_i(om_range),'--','Linewidth',1);
semilogx(hh_, omega(om_range)/2/pi, cH3_i(om_range),'-.','Linewidth',1);
xlim([0.99*x1, (x1+dx)])
ylim([0.98*y1, 1.02*(y1+dy)]);
hold off;


% Angular
axes(ha_gnd_acc_pierBacc(2)); %#ok<*LAXES>
H1_i = 1e6*squeeze(sqrt(sum(Hacc{1}.ResponseData(ind_tpier_Rot,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H1_i),'-','Linewidth',1.5);
hold on;
H2_i = 1e6*squeeze(sqrt(sum(Hacc{2}.ResponseData(ind_tpier_Rot,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H2_i),'--','Linewidth',1);
H3_i = 1e6*squeeze(sqrt(sum(Hacc{3}.ResponseData(ind_tpier_Rot,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H3_i),'-.','Linewidth',1);
axis tight;
legend(legend_txt,'Location','Southeast'); legend box off;
xlim([0.1,max_f_Hz]);
ylabel('AngAcc_{xyz} RSS (adim)','FontSize',plotFontSize);
title('Bottom pier node angular Acc induced by GND Acc');
grid on; hold off;
%
axes(ha_gnd_acc_pierBacc(4));
cH1_i=cumsum(diff(omega).*abs(H1_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH1_i,'-','Linewidth',1.5); hold on;
cH2_i=cumsum(diff(omega).*abs(H2_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH2_i,'--','Linewidth',1);
cH3_i=cumsum(diff(omega).*abs(H3_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH3_i,'-.','Linewidth',1);
xlim([0.1,max_f_Hz]);
ylabel('Cumulative AngAcc RSS (adim)','FontSize',plotFontSize);
xlabel('Frequency (Hz)','FontSize',plotFontSize);
grid on; hold off

% Zoom f > XXHz 
idx_min = find(omega/2/pi > 75, 1, 'first');
idx_max = length(omega)-1;
om_range = idx_min:idx_max;
x1 = omega(idx_min)/2/pi;
y1 = min([cH1_i(om_range); cH2_i(om_range); cH3_i(om_range)]);
dx = omega(idx_max)/2/pi - omega(idx_min)/2/pi;
dy = abs(max([cH1_i(om_range); cH2_i(om_range); cH3_i(om_range)]) - y1);
ylim([0, 1.02*(y1+dy)]);
ax_box_pos_size = [.62 .19 .18 .2];

r = rectangle('Position',[x1-10,0.98*y1,dx,1.02*dy]); r.LineStyle = '--';
% create a new pair of axes inside current figure
hh_ = axes('position',ax_box_pos_size);
box on % put box around new pair of axes
semilogx(hh_, omega(om_range)/2/pi, cH1_i(om_range),'-','Linewidth',1.5);
hold on; grid on
semilogx(hh_, omega(om_range)/2/pi, cH2_i(om_range),'--','Linewidth',1);
semilogx(hh_, omega(om_range)/2/pi, cH3_i(om_range),'-.','Linewidth',1);
xlim([0.99*x1, (x1+dx)])
ylim([0.98*y1, 1.02*(y1+dy)]);
hold off;


%% Bottom pier 6F -> Bottom pier node displacements
%%
f2d = figure(figNumber+22);
set(f2d,'defaultAxesColorOrder',map_c);
set(gcf,'Position',[401, 367, 750, 420]);
try axes(ha_pier6D_tfs(1));
catch, clear ha_pier6D_tfs; end
if(~exist('ha_pier6D_tfs','var'))
    [ha_pier6D_tfs, ~] = utils.tight_subplot(2,2,[.06 .07],[.09 .06],[.07 .02]);
end

ind_pierB_Txy = [1,2];
ind_pierB_Tz = 3;
ind_U1 = (1:2);
ind_U2 = (3);

axes(ha_pier6D_tfs(1)); %#ok<*LAXES>
H1_i = 1e6*squeeze(sqrt(sum(H{1}.ResponseData(ind_pierB_Txy,ind_U1,:).^2,[1,2])));
loglog(omega/2/pi, abs(H1_i),'-','Linewidth',1.5);
hold on;
H2_i = 1e6*squeeze(sqrt(sum(H{2}.ResponseData(ind_pierB_Txy,ind_U1,:).^2,[1,2])));
loglog(omega/2/pi, abs(H2_i),'--','Linewidth',1);
H3_i = 1e6*squeeze(sqrt(sum(H{3}.ResponseData(ind_pierB_Txy,ind_U1,:).^2,[1,2])));
loglog(omega/2/pi, abs(H3_i),'-.','Linewidth',1);
axis tight;
legend(legend_txt,'Location','Southwest'); legend box off;
xlim([0.1,max_f_Hz]);
ylabel('Pier T_{xy} RSS (um/N)','FontSize',plotFontSize);
title('Bottom pier Lat translations induced by forces');
grid on; hold off;

axes(ha_pier6D_tfs(3));
cH1_i=cumsum(diff(omega).*abs(H1_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH1_i,'-','Linewidth',1.5); hold on;
cH2_i=cumsum(diff(omega).*abs(H2_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH2_i,'--','Linewidth',1);
cH3_i=cumsum(diff(omega).*abs(H3_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH3_i,'-.','Linewidth',1);
xlim([0.1,max_f_Hz]);
ylabel('Cumulative T_{xy} RSS (um/N)','FontSize',plotFontSize);
xlabel('Frequency (Hz)','FontSize',plotFontSize);
grid on; hold off;

% Zoom f > XXHz 
idx_min = find(omega/2/pi > 80, 1, 'first');
idx_max = length(omega)-1;
om_range = idx_min:idx_max;
x1 = omega(idx_min)/2/pi;
y1 = min([cH1_i(om_range); cH2_i(om_range)]);
dx = omega(idx_max)/2/pi - omega(idx_min)/2/pi;
dy = abs(max([cH1_i(om_range); cH2_i(om_range)]) - y1);
% ax_box_pos_size = [.26 .17 .18 .15];
ax_box_pos_size = [.12 .23 .18 .15];

r = rectangle('Position',[x1-15,0.98*y1,dx+15,1.1*dy]); r.LineStyle = '--';
% create a new pair of axes inside current figure
hh_ = axes('position',ax_box_pos_size);
box on % put box around new pair of axes
semilogx(hh_, omega(om_range)/2/pi, cH1_i(om_range),'-','Linewidth',1.5);
hold on; grid on
semilogx(hh_, omega(om_range)/2/pi, cH2_i(om_range),'--','Linewidth',1);
semilogx(hh_, omega(om_range)/2/pi, cH3_i(om_range),'-.','Linewidth',1);
xlim([0.99*x1, (x1+dx)])
ylim([0.98*y1, 1.02*(y1+dy)]);
hold off;


% V
axes(ha_pier6D_tfs(2)); %#ok<*LAXES>
H1_i = 1e6*squeeze(sqrt(sum(H{1}.ResponseData(ind_pierB_Tz,ind_U2,:).^2,[1,2])));
loglog(omega/2/pi, abs(H1_i),'-','Linewidth',1.5);
hold on;
H2_i = 1e6*squeeze(sqrt(sum(H{2}.ResponseData(ind_pierB_Tz,ind_U2,:).^2,[1,2])));
loglog(omega/2/pi, abs(H2_i),'--','Linewidth',1);
H3_i = 1e6*squeeze(sqrt(sum(H{3}.ResponseData(ind_pierB_Tz,ind_U2,:).^2,[1,2])));
loglog(omega/2/pi, abs(H3_i),'-.','Linewidth',1);
axis tight;
legend(legend_txt,'Location','Southwest'); legend box off;
xlim([0.1,max_f_Hz]);
ylabel('Pier T_{z} RSS (um/N)','FontSize',plotFontSize);
title('Bottom pier Vert translations induced by F_{z}');
grid on; hold off;
%
axes(ha_pier6D_tfs(4));
cH1_i=cumsum(diff(omega).*abs(H1_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH1_i,'-','Linewidth',1.5); hold on;
cH2_i=cumsum(diff(omega).*abs(H2_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH2_i,'--','Linewidth',1);
cH3_i=cumsum(diff(omega).*abs(H3_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH3_i,'-.','Linewidth',1);
xlim([0.1,max_f_Hz]);
ylabel('Cumulative T_{z} (um/N)','FontSize',plotFontSize);
xlabel('Frequency (Hz)','FontSize',plotFontSize);
grid on; hold off

% Zoom f > XXHz 
idx_min = find(omega/2/pi > 80, 1, 'first');
idx_max = length(omega)-1;
om_range = idx_min:idx_max;
x1 = omega(idx_min)/2/pi;
y1 = min([cH1_i(om_range); cH2_i(om_range)]);
dx = omega(idx_max)/2/pi - omega(idx_min)/2/pi;
dy = abs(max([cH1_i(om_range); cH2_i(om_range)]) - y1);
% ax_box_pos_size = [.77 .15 .18 .2];
ax_box_pos_size = [.63 .22 .18 .2];

r = rectangle('Position',[x1-10,0.98*y1,dx,2*dy]); r.LineStyle = '--';
% create a new pair of axes inside current figure
hh_ = axes('position',ax_box_pos_size);
box on % put box around new pair of axes
semilogx(hh_, omega(om_range)/2/pi, cH1_i(om_range),'-','Linewidth',1.5);
hold on; grid on
semilogx(hh_, omega(om_range)/2/pi, cH2_i(om_range),'--','Linewidth',1);
semilogx(hh_, omega(om_range)/2/pi, cH3_i(om_range),'-.','Linewidth',1);
xlim([0.99*x1, (x1+dx)])
ylim([0.98*y1, 1.02*(y1+dy)]);
hold off;



%% Save analysis plots
%%

plot_file_label = strcat('./figures/acc_tf_comp',...
    strrep(FileFolder_char{n_models}(1:14),'\_','_'));

if(save_plots || 0)
%     print(f_eigenfreq,sprintf('%s_eig_freq_comparison', plot_file_label),'-depsc');
%     print(f0,sprintf('%s_Gmimo_fr', plot_file_label),'-depsc');
%     print(f1,sprintf('%s_Gii_fr', plot_file_label),'-depsc');
    print(f2,sprintf('%s_acc2pttWFE_fr', plot_file_label),'-depsc');
    print(f2b,sprintf('%s_acc2TopPier_fr', plot_file_label),'-depsc');    
end

return

%% Other plots
%%

f2b_ = figure(figNumber+21);
set(f2b_,'defaultAxesColorOrder',map_c);
set(gcf,'Position',[421, 267, 750, 420]);
addpath('/Users/rromano/Workspace/im-tools');
try axes(ha_gnd_acc_gir6D(1));
catch, clear ha_gnd_acc_gir6D; end
if(~exist('ha_gnd_acc_gir6D','var'))
    [ha_gnd_acc_gir6D, ~] = utils.tight_subplot(2,2,[.06 .07],[.09 .06],[.07 .02]);
end

ind_gir_Txyz = (1:3)+12;
ind_gir_Rxyz = (4:6)+12;
ind_gnd_acc = (1:3);

axes(ha_gnd_acc_gir6D(1)); %#ok<*LAXES>
H1_i = 1e6*squeeze(sqrt(sum(H{1}.ResponseData(ind_gir_Txyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H1_i),'-','Linewidth',1.5);
hold on;
H2_i = 1e6*squeeze(sqrt(sum(H{2}.ResponseData(ind_gir_Txyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H2_i),'--','Linewidth',1);
H3_i = 1e6*squeeze(sqrt(sum(H{3}.ResponseData(ind_gir_Txyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H3_i),'-.','Linewidth',1);
axis tight;
legend(legend_txt,'Location','Southwest'); legend box off;
xlim([0.5,max_f_Hz]);
ylabel('GIR node T_{xyz} RSS (um/m/s^2)','FontSize',plotFontSize);
title('GIR node translations induced by GND Accelerations');
grid on; hold off;

axes(ha_gnd_acc_gir6D(3));
cH1_i=cumsum(diff(omega).*abs(H1_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH1_i,'-','Linewidth',1.5); hold on;
cH2_i=cumsum(diff(omega).*abs(H2_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH2_i,'--','Linewidth',1);
cH3_i=cumsum(diff(omega).*abs(H3_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH3_i,'-.','Linewidth',1);
xlim([0.5,max_f_Hz]);
ylabel('Cumulative GIR T_{xyz} RSS (um/m/s^2)','FontSize',plotFontSize);
xlabel('Frequency (Hz)','FontSize',plotFontSize);
grid on; hold off;


% R_xyz
axes(ha_gnd_acc_gir6D(2)); %#ok<*LAXES>
H1_i = 1e6*squeeze(sqrt(sum(H{1}.ResponseData(ind_gir_Rxyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H1_i),'-','Linewidth',1.5);
hold on;
H2_i = 1e6*squeeze(sqrt(sum(H{2}.ResponseData(ind_gir_Rxyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H2_i),'--','Linewidth',1);
H3_i = 1e6*squeeze(sqrt(sum(H{3}.ResponseData(ind_gir_Rxyz,ind_gnd_acc,:).^2,[1,2])));
loglog(omega/2/pi, abs(H3_i),'-.','Linewidth',1);

axis tight;
legend(legend_txt,'Location','Southwest'); legend box off;
xlim([0.5,max_f_Hz]);
ylabel('GIR R_{xyz} RSS (urad/m/s^2)','FontSize',plotFontSize);
title('GIR node rotations induced by GND Accelerations');
grid on; hold off;
%
axes(ha_gnd_acc_gir6D(4));
cH1_i=cumsum(diff(omega).*abs(H1_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH1_i,'-','Linewidth',1.5); hold on;
cH2_i=cumsum(diff(omega).*abs(H2_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH2_i,'--','Linewidth',1);
cH3_i=cumsum(diff(omega).*abs(H3_i(1:end-1)));
semilogx(omega(1:end-1)/2/pi, cH3_i,'-.','Linewidth',1);

xlim([0.5,max_f_Hz]);
ylabel('Cumulative GIR R_{xyz} RSS (urad/m/s^2)','FontSize',plotFontSize);
xlabel('Frequency (Hz)','FontSize',plotFontSize);
grid on; hold off


%% AUXILIAR FUNCTIONS
%%
%% Function to create structures with Bode and Nichols plot options
%%
function [hbode, hnichols] = get_fr_plot_opts(font_size)
% - - -
% Bode plot
hbode = bodeoptions;
hbode.TickLabel.FontSize = font_size;
hbode.YLabel.FontSize = font_size;
hbode.XLabel.FontSize = font_size;
hbode.FreqUnits = 'Hz';
hbode.Grid = 'on';
hbode.PhaseMatching = 'on';
hbode.PhaseMatchingValue = -180;
% - - -
% Nichols plot
hnichols = nicholsoptions;
hnichols.TickLabel.FontSize = font_size;
hnichols.YLabel.FontSize = font_size;
hnichols.XLabel.FontSize = font_size;
hnichols.FreqUnits = 'Hz';
hnichols.Grid = 'on';
hnichols.PhaseMatching = 'on';
hnichols.PhaseMatchingValue = -180;
end
 
%% Function to compute the FR used to assess the effect of the pier properties
%%
function [om0, G_jw] = compute_freq_resp(FileFolder, max_model_Fr, omega)

RootFolder = '/home/rromano/Workspace/gmt-data';
FileName = "modal_state_space_model_2ndOrder.mat";
fprintf('Loading model %s\nfrom folder\n%s\n',FileName,FileFolder)
load(fullfile(RootFolder,FileFolder,FileName),...
    'inputs2ModalF','modalDisp2Outputs',...
    'eigenfrequencies','inputTable','outputTable');

include_dcmc = false;
StaticFileName = fullfile(RootFolder,FileFolder,"static_reduction_model_dontUSE.mat");
if(exist(StaticFileName,'file'))
    try
        load(StaticFileName,'gainMatrixMountControlled');
        gainMatrix = gainMatrixMountControlled;
    catch, load(StaticFileName,'gainMatrix');
    end
    fprintf('Static gain matrix loaded from \n%s\n', FileFolder);
    include_dcmc = true;
end

% Handle modal parameters
om0 = eigenfrequencies(:)*2*pi;
fprintf("The model maximum eigenfrequency is %.5fHz\n", om0(end)/2/pi);

if nargout == 1
    G_jw = [];  % Return just the vector of eigenfrequencies
    return
end

% First model mode
first_mode = 4;%1;

% Keep modes just up to frequency <max_model_Fr>
if(isempty(max_model_Fr))
    max_model_Fr = max(eigenfrequencies);
end
[~,i_trim] = min(abs(eigenfrequencies - max_model_Fr));
% 2nd order modal model parameters
phiB = inputs2ModalF(first_mode:i_trim,:);
phiC = modalDisp2Outputs(:,first_mode:i_trim);
% Eigenfrequencies and structural damping ratio vectors
om0 = om0(first_mode:i_trim);
if(first_mode  ~= 1)
    warning('Discarding the first %d modes of the FRs computations.',first_mode);
end
fprintf('\nFirst 8 eigenfrequencies:\n%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\n',...
    eigenfrequencies(1:8));

% Set structural damping
sysDamp = 0.02; %0.005;%
damp = sysDamp*ones(size(om0));
fprintf('\nPlant model damping ratio set to:%.3g\n',sysDamp);

% Pick selected mount subsystem and wrap as a CMM object
% Inputs & output indexes

try
    % Azimuth
    in1 = inputTable{'OSS_AzDrive_Torque','indices'}{1}';
    out1 = outputTable{'OSS_AzEncoder_Angle','indices'}{1}';
    % Elevation
    in2 = inputTable{'OSS_ElDrive_Torque','indices'}{1}';
    out2 = outputTable{'OSS_ElEncoder_Angle','indices'}{1}';
    % GIR
    in3 = inputTable{'OSS_RotDrive_Torque','indices'}{1}';
    out3 = outputTable{'OSS_RotEncoder_Angle','indices'}{1}';

    mnt_phiB = [mean(phiB(:,in1),2), mean(phiB(:,in2),2), mean(phiB(:,in3),2)];
    mnt_phiC = [mean(phiC(out1,:),1); mean(phiC(out2,:),1); mean(phiC(out3,:),1)];
catch
    warning('Mount IOs not available in inputTable/outputTable!')
    mnt_phiB = zeros(size(phiB,1),3);
    mnt_phiC = zeros(3,size(phiC,2));
end

% Ground acceleration inputs
try in4 = inputTable{'OSS00_Ground_Acc','indices'}{1}';
catch
    in4 = 1:3;
    warning('Model %s\ndoes not provide OSS00_Ground_Acc inputs!',FileFolder);
end
% Top of the pier node inputs
try in5 = inputTable{'Pier_6F','indices'}{1}(1:6)'; % (1:6 Bottom)|(7:12 Top node)
catch
    in5 = 1:6;
    warning('Model %s\ndoes not provide Pier_6F inputs!',FileFolder);
end

% Pier 6D
pier_top_out = outputTable{'Pier_6D','indices'}{1}(1:12)';
% OSS GIR 6D
try gir_or_stick_out = outputTable{'OSS_GIR_6d','indices'}{1}(1:6)';
catch
    gir_or_stick_out = outputTable{'Stick_telescope_6D','indices'}{1}(1:6)'; % (1:6 Top)|(7:12 Bottom node)
    warning('Model %s\ndoes not provide GIR_6D outputs!',FileFolder);
end
% Subset of matrices B and C
phiB_ = [phiB(:,in4), phiB(:,in5), mnt_phiB];
n_u = size(phiB_,2);
phiC_ = [phiC(pier_top_out,:); phiC(gir_or_stick_out,:); mnt_phiC];
n_y = size(phiC_,1);

% Compute frequency responses
fprintf('Computing G(jw) with %d modes and w=2*pi*[%.3g ~ %.3g]\n',...
    length(om0), min(omega)/2/pi ,max(omega)/2/pi);

G_jw = frd(...
    freqresp_2ndorder_model(...
    phiB_,phiC_,omega,om0,damp,1:n_u,1:n_y), omega);

if(include_dcmc)
    invOM2 = diag(1./((2*pi*eigenfrequencies(first_mode:i_trim)).^2));
    G_ = phiC_(1:end-3,:) * invOM2 * phiB_(:,1:end-3);
    Ks = gainMatrix([pier_top_out,gir_or_stick_out], [in4,in5]);
    D = [[Ks - G_, zeros(size(Ks,1),3)]; zeros(3, size(Ks,2)+3)];
    dcmc_delay = 1e-3;
    G_jw = G_jw + frd(D*tf(1,1,'IOdelay',dcmc_delay), omega);
end

end

%% Function to load ODC controller object
%%
function o = load_ODC_control(sZa, sVer, sSubVer, sDamping)

% oTest.sZa: elevation zenith angle (ZA) as string e.g. '00','30','60'
oTest.sZa = sZa;
% oTest.sVer: FEM version as string e.g. '19'
oTest.sVer = sVer;
% oTest.sSubVer: FEM subversion as string e.g. '1'
oTest.sSubVer = sSubVer; %'11'; %'2'; %'9';
% oTest.sDamping: now '02' means 2% structural dumping
oTest.sDamping = sDamping;

if(~strcmp('02',sDamping))
    warning('Unable to load controller for 0.02 damping');
    % Open-loop
    o.Hp = tf({0, 0, 0; 0, 0, 0; 0, 0, 0}, 1);
    o.Hff = tf({0, 0, 0; 0, 0, 0; 0, 0, 0}, 1);
    o.Hdrive = tf(1,1);
    return
end

switch sVer
    case '20'
        % Folder with ODC end-to-end functions
        odc_e2efile_folder = fullfile('/Users/rromano/Workspace/mnt-odc',...
            "fdr2023/MatlabFilesE2E_2023-05-10");
        
        odc_base_util_folder = fullfile(odc_e2efile_folder,'base/util');
        odc_base_conf_folder = fullfile(odc_e2efile_folder,'base/conf');
        addpath(odc_base_util_folder,odc_base_conf_folder);
        fprintf('\nIncluding folder\n %s\ninto MatLab path.\n', odc_e2efile_folder);
        fprintf('Getting model parameters ...\n');
      
        % ODC Simulink model used (located in ../base)
        oTest.sRoot = 'root';
        % oTest.sHbsConf: name of HBS configuration: e.g. 'HbTp19'
        oTest.sHbsConf = 'HaTp19'; %'HbTp19'
        % oTest.sViscFrCase: one of 3 cases w.r.t. viscosity: ['ViscFrLow', 'ViscFrMedium', 'ViscFrHigh']  see fun_mueByTempAzEl
        oTest.sViscFrCase = 'ViscFrLow'; %lowest viscosity=>lowest damping
        % oTest.sModelDirIn: directory relative to ../ss_model/ where the state space models are located
        % returns structure [o] with all configuration parameters
        oTest.sModelDirIn = 'v20.11/n100HzR800';
        oTest.bUseReducedModel = true; %false;
        o = fun_confBase(oTest);
        
        % Remove folders from Matlab path
        fprintf('\nRemoving folders\n%s\n%s\nfrom MatLab path.\n',...
            odc_base_util_folder, odc_base_conf_folder);
        rmpath(odc_base_util_folder,odc_base_conf_folder);
        
    case 'OO'
        % Open-loop
        o.Hp = tf({0, 0, 0; 0, 0, 0; 0, 0, 0}, 1);
        o.Hff = tf({0, 0, 0; 0, 0, 0; 0, 0, 0}, 1);
        o.Hdrive = tf(1,1);
        
    case '19'        
        % Folder with ODC End-2-end functions
        odc_base_folder = '/Users/rromano/Workspace/mnt-odc/pdr2021 - E2E_Matlab_Model/base/conf';
        addpath(odc_base_folder);
        fprintf('\nSetting parameters ...\n');
        o = fun_confBase('root', sZa, sVer, sSubVer, sDamping, 0, true);
        % Clear ODC model variable
        rmpath(odc_base_folder);
end

end

%% Function to plot the envelope of the mount rejection transfer function (RTF)
%%
function hlr = plot_mountRTF()
bcolor = [0.3 0.3 0.3];
ax = findobj(gcf,'type','axes');

try xlim_ = get(ax(2),'XLim');
catch, xlim_ = get(ax(end),'XLim'); end

S_env_w = [0.1, 1, 1.8, 30, xlim_(2)]';
S_env_mag = [-60, 0, 7.2, 0, 0]';
hlr = semilogx(S_env_w,S_env_mag,'-.','color',bcolor,'Linewidth',2);

end
 
%% Function to plot closed-loop response envelope
%%
function hlr = plot_cl_OAD_req(CL_Bw0,hifi_w, hifi_mag)

if(nargin < 2), hifi_w = 2*CL_Bw0; hifi_mag = 6; end

bcolor = [0.3 0.3 0.3];
bcolor2 = [0.5 0.5 0.5];
ax = findobj(gcf,'type','axes');

xlim_ = get(ax(end),'XLim');
ylim_ = get(ax(end),'YLim');

semilogx([xlim_(1) CL_Bw0],[-3 -3],'-.','color',bcolor,'Linewidth',2);
hlr = semilogx([CL_Bw0 CL_Bw0],[ylim_(1) -3],'-.','color',bcolor,'Linewidth',2);
semilogx([xlim_(1) CL_Bw0],[6 6],'-.','color',bcolor2,'Linewidth',1.5);
semilogx([hifi_w xlim_(2)],[hifi_mag hifi_mag],'-.','color',bcolor,'Linewidth',2);
semilogx([CL_Bw0 hifi_w],[6 hifi_mag],'-.','color',bcolor,'Linewidth',1.5);

text(0.9*CL_Bw0,-10,sprintf('%.2g Hz\nbandwidth\nlimit',CL_Bw0),...
    'HorizontalAlignment','right',...
    'VerticalAlignment','top',...
    'color',bcolor,'Fontsize',12);
text(18,-4,sprintf('+6dB\nhigh-freq\nbound'),...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom',...
    'color',bcolor,'Fontsize',12);
end

%% Function to plot nichols robustness boundaries
%%
function plot_nichols_ellipse(VM)
bcolor = [0.3 0.3 0.3];%'k';%[0.6350 0.0780 0.1840];

thetaS = linspace(-0.9*pi/2,0,501)';
y = (1/VM)*exp(1i*thetaS);

S = y./(1-y);
ph=angle(S)*(180/pi); ph=ph-(360*ceil(ph/360)); mag=20*log10(abs(S));
line(ph,mag,'color',bcolor,'LineStyle','--','Linewidth',2);
line(ph,-mag,'color',bcolor,'LineStyle','--','Linewidth',2);

% thetaT = linspace(0,0.9*pi/2,501)';
% y = (1/VM)*exp(1i*thetaT);
% T = (1-y)./y;
% ph=angle(T)*(180/pi); ph=ph-(360*ceil(ph/360)); mag=20*log10(abs(T));
% line(ph,mag,'color',bcolor,'LineStyle','--','Linewidth',2);

if 1
    GM = 6; %PM = 30;
    line([-450 -180],GM*[-1 -1],'color',bcolor,'LineStyle','--','Linewidth',2.5)
    line([-180 -180],[GM 20],'color',bcolor,'LineStyle','--','Linewidth',2.5)
end

end

%% Function to compute mode frequency response
%%
function G_jw = freqresp_2ndorder_model(phiB,phiC,om,om0,damp,ind_in,ind_out)
% G_jw = freqresp_2ndorder_model(phiB,phiC,om,om0,damp,ind_in,ind_out)
% * om is vector of frequencies to apply bode at (rad/s)
% * om0 is a vector of natural frequencies (rad/s)
% * damp is a vector of damping ratios (0-1)
% * ind_in and ind_out are input and output indices requested.
% * phiB and phiC are the nonzero partitions of the input and output
% matrices of the state-space model.
%

% Trim input/output matrices as appropriate.
phiB = phiB(:,ind_in);
phiC = phiC(ind_out,:);
%D = D(ind_out,ind_in);

% Define useful dimensions.
[n_m, n_u] = size(phiB);
n_y = size(phiC,1);
n_om = length(om);

s = sqrt(-1)*om;
% Initialize output
G_jw = zeros(n_y, n_u, n_om);

for i =1:n_y %loop through outpus
    for j = 1:n_u % loop through inputs
        for k = 1:n_m %loop through modes
            G_jw(i,j,:) = G_jw(i,j,:)+...
                phiC(i,k)*...
                reshape(1./(s.^2+2*damp(k)*om0(k).*s+om0(k).^2),1,1,n_om)*...
                phiB(k,j);
        end
        %G_jw(i,j,:) = G_jw(i,j,:) + D(i,j);
    end
end

end


%% eof
%%