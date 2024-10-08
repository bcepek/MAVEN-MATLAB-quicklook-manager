function [products_H, products_O, products_O2] = find_STA_mat_file_cleaned(date)
% find cleaned data

date = datestr(date, 'yyyymmdd');
load('paths.mat', 'paths', 'slash')
root_STATIC = paths.sta_moments;
%root_STATIC = '/Users/sesh2112/data/Matlab_MAVEN/pre-calculated/STATIC_moments';
%root_STATIC = 'E:\MatLab';



monthpath_STA = root_STATIC;  %[root '\' cur_year '\' cur_month '\'];
month_filelist_STA_H =  dir([monthpath_STA slash '*H.mat']);
month_filelist_STA_O =  dir([monthpath_STA slash '*O.mat']);
month_filelist_STA_O2 =  dir([monthpath_STA slash '*O2.mat']);


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
        filename = [monthpath_STA slash month_filelist_STA_H(i).name];
        products_H = load (filename);
        break
    elseif i == length(month_filelist_STA_H)
%        disp([date ' for STA H is not found'])
        products_H = [];
    end
end
%BEGIN O
for i = 1:length(month_filelist_STA_O)
    if (length(month_filelist_STA_O(i).name)<10)
        continue
    end
    if(all(month_filelist_STA_O(i).name(26:33) == date)) % for O data
        filename = [monthpath_STA slash month_filelist_STA_O(i).name];
        products_O = load (filename);
        break
    elseif i == length(month_filelist_STA_O)
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
        filename = [monthpath_STA slash month_filelist_STA_O2(i).name];
        products_O2 = load (filename);
        break
    elseif i == length(month_filelist_STA_O2)
%         disp([datefrom ' for STA O2 is not found'])
        products_O2 = [];
    end
end

