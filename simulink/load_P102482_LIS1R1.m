%
% This script loads IDOM's 1R1 pier motion results from an Excel
% spreadsheet and saves them into a MatLab data file.
%

matlab_dt_fname = 'P102482-LIS-1R1-RLE_Acc';
if(exist([matlab_dt_fname,'.mat'],'file'))
    warning('Matlab data file %s.mat already exists!\n',matlab_dt_fname);
    return
end

xls_fname = 'P102482-LIS-183-1R1-GMT Pier RLE_OLE_SLE Accelerations to the Telescope.xlsx';
sheet = 'RLE OUTPUTS (Isolated T-H)';

xls_cols = {'A:D'; 'F:I'; 'K:N'; 'P:S'; 'U:X'; 'Z:AC'; 'AE:AH'};

pier_rle = cell(7,1);
for i_rle=1:7
    fprintf('Importing RLE%0d pier motion data. ', i_rle);
    xls_opts = detectImportOptions(...
        xls_fname,'Sheet',sheet,'Range',xls_cols{i_rle});
    M = readmatrix(xls_fname,xls_opts);
    dt = M(isfinite(M(:, 1)),:);
    fprintf('Time series sampled at %gHz.\n', 1/diff(dt(1:2,1)));
    pier_rle{i_rle}.time = dt(:,1);
    pier_rle{i_rle}.acc = (1/9.80665)*dt(:,2:4);
end

if(true)
    save(matlab_dt_fname,'pier_rle','xls_fname','sheet');
end