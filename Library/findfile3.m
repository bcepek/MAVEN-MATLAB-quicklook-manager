function [products_H, products_O, products_O2, newDatefrom] = findfile3(date)
%find data where oxygen channels are contaminated by protons

date = datestr(date, 'yyyymmdd');
root_STATIC = '\\193.232.6.100\общие\Владимир\MAVEN_data_in_MatLab_format';
root_MAG = '\\193.232.6.100\общие\Владимир\MAVEN_MAG_data_in_MatLab_format';
% root = '\\193.232.6.100\Data\Maven\Static\d1_32e4d16a8m';

isfound_mag = 0;
isfound_H = 0;
isfound_O = 0;
isfound_O2 = 0;



monthpath_STA = root_STATIC;  %[root '\' cur_year '\' cur_month '\'];
month_filelist_STA_H =  dir([monthpath_STA '\' '*H.mat']);
month_filelist_STA_O =  dir([monthpath_STA '\' '*O.mat']);
month_filelist_STA_O2 =  dir([monthpath_STA '\' '*O2.mat']);

monthpath_MAG = root_MAG;  %[root '\' cur_year '\' cur_month '\'];
month_filelist_MAG =  dir([monthpath_MAG '\' '*sts.mat']);

%{
for i = 1:length(month_filelist_MAG)
    isfound_mag = false;
    if (length(month_filelist_MAG(i).name)<10)
        continue
    end
    if(all(month_filelist_MAG(i).name(22:29) == datefrom)) % for MAG data
        filename = [monthpath_MAG '\' month_filelist_MAG(i).name];
        isfound_mag = true;
        products_mag = load (filename);
        break
    elseif i == length(month_filelist_MAG)
        isfound_mag = false;
        disp([datefrom ' for MAG is not found'])
        products_mag = [];
    end
end
%}
for i = 1:length(month_filelist_STA_H)
    if (length(month_filelist_STA_H(i).name)<10)
        continue
    end
    if(all(month_filelist_STA_H(i).name(26:33) == date)) % for H data
        filename = [monthpath_STA '\' month_filelist_STA_H(i).name];
        isfound_H = true;
        products_H = load (filename);
        break
    elseif i == length(month_filelist_STA_H)
        isfound_H = false;
        disp([date ' for STA is not found'])
        products_H = [];
    end
end
%BEGIN O
for i = 1:length(month_filelist_STA_O)
    if (length(month_filelist_STA_O(i).name)<10)
        continue
    end
    if(all(month_filelist_STA_O(i).name(26:33) == date)) % for O data
        filename = [monthpath_STA '\' month_filelist_STA_O(i).name];
        isfound_O = true;
        products_O = load (filename);
        break
    elseif i == length(month_filelist_STA_O)
        isfound_O = false;
%         disp([datefrom ' for STA O is not found'])
        products_O = [];
    end
end
%END O

%BEGIN O2
for i = 1:length(month_filelist_STA_O2)
    if (length(month_filelist_STA_O2(i).name)<10)
        continue
    end
    if(all(month_filelist_STA_O2(i).name(26:33) == date)) % for O2 data
        filename = [monthpath_STA '\' month_filelist_STA_O2(i).name];
        isfound_O2 = true;
        products_O2 = load (filename);
        break
    elseif i == length(month_filelist_STA_O2)
        isfound_O2 = false;
%         disp([datefrom ' for STA O2 is not found'])
        products_O2 = [];
    end
end


% i = i + 1;
% if i ~= length(month_filelist_STA_O2)
%     newDatefrom = month_filelist_STA_O2(i).name(26:33);
% end