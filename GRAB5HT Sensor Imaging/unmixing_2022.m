function [] = unmixing_2022
%%
% Read raw mixed spectra time series from Oceanview output, and unmix the spectra time series one by one using user defined spectra reference.
% Tzu-Hao Harry Chao 2020/02/14
%%

clc


Data=dir('*.txt'); % extract the names of all txt files
path_data=pwd; % get current folder path

% [dataID,path_data] = uigetfile('*.txt','Select data');
% cd(path_data)

[refID,path_ref] = uigetfile('*.csv','Select reference');
ref=csvread([path_ref refID],1,1);

for ii = 1:length(Data)

    i=0; test = {{'a'}}; % skipping headers
    while isnan(str2double(test{1,1})) == 1
        file = fopen([path_data,'/',Data(ii).name],'r');
%         file = fopen([path_data dataID],'r');
        test = textscan(file, '%s',1,'HeaderLines',i);
        fclose(file);
        i=i+1;
    end
    
    % to fix additional gap line before the wavelength column
    if i == 2
        i = i+1; 
    else
        i = i;
    end


    class=[];
    for j=1:1044
        class=[class '%f '];
    end

    file = fopen([path_data,'/',Data(ii).name],'r');
% file = fopen([path_data dataID],'r');   
    data = textscan(file, ['%s' '%s' '%s' class],'HeaderLines',i); 
    data = cell2mat(data(4:end));  % remove the first 3 columns
    data = data';
    fclose(file);

    coef=zeros(size(ref,2),size(data,2));

    for j=1:size(data,2)
        coef(:,j)=max(0,lsqnonneg(ref(140:500,:), data(140:500,j)));
        clc
        [num2str(j/size(data,2)*100) '%']
    end

    figure
    for j=1:size(coef,1)
        subplot(size(coef,1),1,j)
        plot(0.1:0.1:size(data,2)/10,coef(j,:)/mean(coef(j,:),2)*100-100)
%title('GCaMP time course','FontWeight','bold','FontSize',12)
        xlabel('Time (s)','FontWeight','bold','FontSize',12)
        ylabel('dF/F (%)','FontWeight','bold','FontSize',12)
    end

    save([Data(ii).name(1:length(Data(ii).name)-4) '.mat'],'coef')
end
