function filename = download_swe_svyspec(time)
mvn_year = num2str(year(time));
mvn_month = num2str(month(time));
if(length(mvn_month)<2)
    mvn_month = ['0' mvn_month];
end
mvn_day = num2str(day(time));
if(length(mvn_day)<2)
    mvn_day = ['0' mvn_day];
end

load("paths.mat", 'paths', 'slash');
download_path = [paths.swe_svyspec slash mvn_year slash mvn_month slash];
if ~exist(download_path, 'dir')
    mkdir(download_path)
end

url = ['https://lasp.colorado.edu/maven/sdc/public/data/sci/swe/l2/' mvn_year '/' mvn_month '/'];
web_page = webread(url);

pattern = ['mvn_swe_l2_svyspec_' mvn_year mvn_month mvn_day '_v'] + digitsPattern(2) + '_r' + digitsPattern(2) + '.cdf';
fname_ind = strfind(web_page, pattern);
fname = web_page(fname_ind(1):fname_ind(1)+38);
websave([download_path, fname], [url, fname]);

filename = [download_path, fname];
end