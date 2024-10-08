function mag_sts2mat(filename)

load("paths.mat", "paths", 'slash');
writepath = paths.mag_ss_mat;

fileID=fopen(filename,'r');
tline = fgetl(fileID);
headersNumber=0;
while ~strcmp(tline(1:4),'  20')
    headersNumber=headersNumber+1;
    tline = fgetl(fileID);
end
%frewind(fileID);   % ðåøèëè ñ Ñåðåêîì, ÷òî íå íóæíî
%tline = fgetl(fileID);   %  ðåøèëè ñ Ñåðåêîì, ÷òî íå íóæíî
fclose(fileID);
% mf_data = dlmread(mf_filename, '', 145);
mf_data = dlmread(filename, '', headersNumber,0);

mf_epoch = mf_data(:, 7) + datenum(['00-Jan-', num2str(mf_data(1, 1)), ' 00:00:00']); %#ok<NASGU>
%     choose_ind = find(mf_epoch>=timefrom & mf_epoch<=timeto);
%     mf_data2 = mf_data(choose_ind, :);
Bx = mf_data(:, 8);
By = mf_data(:, 9);
Bz = mf_data(:, 10);
B = sqrt(Bx.^2 + By.^2 + Bz.^2); %#ok<NASGU>
pos_sc_mso = mf_data(:, [12 13 14]); %#ok<NASGU>
% new_filename = strcat('mf_', datestr(mf_epoch(1),'yyyy-mm-dd'),'.mat');
% save(new_filename,'mf_epoch','Bx','By','Bz','B')

if ~exist(writepath, 'dir')
    mkdir(writepath)
end

save([writepath, slash, filename(end-40:end),'.mat'],'mf_epoch','Bx','By','Bz','B', 'pos_sc_mso')
end