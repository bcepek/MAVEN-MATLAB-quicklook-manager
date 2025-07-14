function filename = download_lpw_bursthf(time)

server = "berkeley";
%server = "lasp";
if server == "berkeley"
    server_url = 'https://sprg.ssl.berkeley.edu/data/maven/data/sci/lpw/l2/';
elseif server == "lasp"
    server_url = 'https://lasp.colorado.edu/maven/sdc/public/data/sci/lpw/l2/';
end

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
download_path = [paths.lpw_bursthf slash mvn_year slash mvn_month slash];
if ~exist(download_path, 'dir')
    mkdir(download_path)
end

url = [server_url, mvn_year, '/', mvn_month, '/'];
web_page = webread(url);

pattern = ['mvn_lpw_l2_we12bursthf_' mvn_year mvn_month mvn_day '_v'] + digitsPattern(2) + '_r' + digitsPattern(2) + '.cdf';
fname_ind = strfind(web_page, pattern);
if isempty(fname_ind)
    filename = -1;
    return
end
fname = web_page(fname_ind(end):fname_ind(end)+42);
websave([download_path, fname], [url, fname]);

filename = [download_path, fname];
end