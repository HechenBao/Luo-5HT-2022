function [] = unmixing
clc

%[dataID,path_data] = uigetfile('*.txt','Select data');
[refID,path_ref] = uigetfile('*.csv','Select reference');
ref=csvread([path_ref refID],1,1);

Data=dir('*.txt'); % extract the names of all txt files
path_data=pwd; % get current folder path

for i = 1:length(Data)
file = fopen([path_data,'/',Data(i).name],'r');
data = textscan(file, '%s','HeaderLines',15); %orignal 13 with left, or 15 with right desktop
data=reshape([0;0;data{1,1}],1046,(length(data{1,1})+2)/1046);
data=data(3:end,2:end);fclose(file);

coef=zeros(size(ref,2),size(data,2));

    for j=1:size(data,2)
        coef(:,j)=max(0,lsqnonneg(ref(140:500,:),str2double(data(140:500,j))));
        clc
        [num2str(j/size(data,2)*100) '%']
    end

save([Data(i).name(1:length(Data(i).name)-4) '.mat'],'coef')

data=[]; % clear up memory
coef=[]; % clear up memory

end



%save([dataID(1:length(dataID)-4) '.mat'],'coef')