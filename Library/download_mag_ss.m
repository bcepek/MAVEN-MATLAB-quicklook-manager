function fname = download_mag_ss(time)

%time = datenum('2019-12-14 06:43:00', 'yyyy-mm-dd HH:MM:SS');

mvn_year = num2str(year(time));
mvn_month = num2str(month(time));
if(length(mvn_month)<2)
    mvn_month = ['0' mvn_month];
end
mvn_day = num2str(day(time));
if(length(mvn_day)<2)
    mvn_day = ['0' mvn_day];
end

load("paths.mat", 'paths', 'slash')
download_path = [paths.mag_ss_sts slash mvn_year slash mvn_month slash];
if ~exist(download_path, 'dir')
    mkdir(download_path)
end

url = ['https://lasp.colorado.edu/maven/sdc/public/data/sci/mag/l2/' mvn_year '/' mvn_month '/'];
web_page = webread(url);
pattern = 'mvn_mag_l2_' + digitsPattern(7) + 'ss_' + mvn_year+mvn_month+mvn_day + '_v' + digitsPattern(2) + '_r' + digitsPattern(2) + '.sts';
fname_ind = strfind(web_page, pattern);
fname = web_page(fname_ind(1):fname_ind(1)+40);

websave([download_path, fname], [url, fname]);
fname = [download_path, fname];
end