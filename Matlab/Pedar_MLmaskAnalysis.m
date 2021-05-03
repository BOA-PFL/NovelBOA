clear
clc

addpath('C:\Users\kate.harrison\Documents\GitHub\NovelBOA\NovelBOA\Matlab')
input_dir = 'C:\Users\kate.harrison\Dropbox (Boa)\AgilityPerformance\BOA_Basketball_Mar2021\PressureData';

cd(input_dir)
files = dir('*.mva');
dataList = {files.name};
[f,~] = listdlg('PromptString','Select data files','SelectionMode','multiple','ListString',dataList);

NumbTrials = length(f);

Subject = cell(2, 1);
Shoe = cell(2, 1);
Movement = cell(2, 1);

sdHeel = zeros(2,1);
meanHeel = zeros(2,1);
sdLatFF = zeros(2,1);
meanLatFF = zeros(2,1);
sdMedFF = zeros(2,1);
meanMedFF = zeros(2,1);
sdLatMF = zeros(2,1);
meanLatMF = zeros(2,1);
sdMedMF = zeros(2,1);
meanMedMF = zeros(2,1);
sdToes = zeros(2,1);
meanToes = zeros(2,1);

maxHeel = zeros(2,1);
maxLatFF = zeros(2,1);
maxMedFF = zeros(2,1);
maxLatMF = zeros(2,1);
maxMedMF = zeros(2,1);
maxToes = zeros(2,1);

MAXmaxHeel = zeros(2,1);
MAXmaxLatFF = zeros(2,1);
MAXmaxMedFF = zeros(2,1);
MAXmaxLatMF = zeros(2,1);
MAXmaxMedMF = zeros(2,1);
MAXmaxToes = zeros(2,1);

contact_time = zeros(2,1);
r = 0;

for i = 1:NumbTrials
    
    FileName = dataList(f(i));
    FileLoc = char(strcat(input_dir,'\', FileName));
    data = dlmread(FileLoc, '\t', 17 , 0);
    
    names = split(FileName, ["_", " ", "."]);
    sub = names{1};
    shoe = names{2};
    move = names{3};
        
    time = data(:,1);
    Heel_Force = data(:,2);
    Heel_MaxP = data(:,3);
    Heel_MeanP = data(:,4);
    Heel_Pct = data(:,5);
    MedMF_Force = data(:,6);
    MedMF_MaxP = data(:,7);
    MedMF_MeanP = data(:,8);
    MedMF_Pct = data(:,9);
    LatMF_Force = data(:,10);
    LatMF_MaxP = data(:,11);
    LatMF_MeanP = data(:,12);
    LatMF_Pct = data(:,13);
    MedFF_Force = data(:,14);
    MedFF_MaxP = data(:,15);
    MedFF_MeanP = data(:,16);
    MedFF_Pct = data(:,17);
    LatFF_Force = data(:,18);
    LatFF_MaxP = data(:,19);
    LatFF_MeanP = data(:,20);
    LatFF_Pct = data(:,21);
    Toes_Force = data(:,22);
    Toes_MaxP = data(:,23);
    Toes_MeanP = data(:,24);
    Toes_Pct = data(:,25);
    
    Force_Tot = Heel_Force + MedMF_Force + LatMF_Force + MedFF_Force + LatFF_Force + Toes_Force; % sum force under each mask to find total force
    
    hold off
    plot (Force_Tot)
    title(FileName)
    
    [start, fThresh] = ginput(1); %select start of usable data, at y-value of desired minimum force threshold
    start = round(start);
    
    [finish, ~] = ginput(1); %select end of useable data (y-value doesn't matter)
    finish = round(finish);
    
    Force_Tot(Force_Tot<fThresh) = 0; % zero all force values below threshold
    Force_Tot(1:start) = 0; % elminate data before start of usable data
    Force_Tot(finish:end) = 0; % eliminate data after end of usable data
    
    R_steps = zeros(1,1);
    counter = 1;

    for j = 1:length(Force_Tot)-1
        if Force_Tot(j) == 0 && Force_Tot(j+1) > 0
           R_steps(counter) = j;
           counter = counter + 1;
        end
    end

    R_ends = zeros(1,1);
    counter = 1;

    for j = 1:length(Force_Tot)-1
        if Force_Tot(j) > 0 && Force_Tot(j+1) == 0 
           R_ends(counter) = j;
           counter = counter + 1;
        end
    end

    if R_ends(1)<R_steps(1)
        R_ends = R_ends(2:end);
    end

    if length(R_steps) > length(R_ends)
            R_steps = R_steps(1:length(R_ends)+1);
    elseif length(R_steps) < length(R_ends)
            R_ends = R_ends(1:length(R_steps));
    end
        
    steps = length(R_ends);
    
    for j = 1:steps
        Subject{r+j} = sub;
        Shoe{r+j} = shoe;
        Movement{r+j} = move;
        
        contact_time(r+j)= R_ends(j)-R_steps(j);
        
        meanHeel(r+j) = mean((Heel_MeanP(R_steps(j):R_ends(j))));
        sdHeel(r+j) = std((Heel_MeanP(R_steps(j):R_ends(j))));
        maxHeel(r+j) = max((Heel_MeanP(R_steps(j):R_ends(j))));
        MAXmaxHeel(r+j) = max((Heel_MaxP(R_steps(j):R_ends(j))));
        
        meanMedMF(r+j) = mean((MedMF_MeanP(R_steps(j):R_ends(j))));
        sdMedMF(r+j) = std((MedMF_MeanP(R_steps(j):R_ends(j))));
        maxMedMF(r+j) = max((MedMF_MeanP(R_steps(j):R_ends(j))));
        MAXmaxMedMF(r+j) = max((MedMF_MaxP(R_steps(j):R_ends(j))));
        
        meanLatMF(r+j) = mean((LatMF_MeanP(R_steps(j):R_ends(j))));
        sdLatMF(r+j) = std((LatMF_MeanP(R_steps(j):R_ends(j))));
        maxLatMF(r+j) = max((LatMF_MeanP(R_steps(j):R_ends(j))));
        MAXmaxLatMF(r+j) = max((LatMF_MaxP(R_steps(j):R_ends(j))));
        
        meanMedFF(r+j) = mean((MedFF_MeanP(R_steps(j):R_ends(j))));
        sdMedFF(r+j) = std((MedFF_MeanP(R_steps(j):R_ends(j))));
        maxMedFF(r+j) = max((MedFF_MeanP(R_steps(j):R_ends(j))));
        MAXmaxMedFF(r+j) = max((MedFF_MaxP(R_steps(j):R_ends(j))));
        
        meanLatFF(r+j) = mean((LatFF_MeanP(R_steps(j):R_ends(j))));
        sdLatFF(r+j) = std((LatFF_MeanP(R_steps(j):R_ends(j))));
        maxLatFF(r+j) = max((LatFF_MeanP(R_steps(j):R_ends(j))));
        MAXmaxLatFF(r+j) = max((LatFF_MaxP(R_steps(j):R_ends(j))));
        
        meanToes(r+j) = mean((Toes_MeanP(R_steps(j):R_ends(j))));
        sdToes(r+j) = std((Toes_MeanP(R_steps(j):R_ends(j))));
        maxToes(r+j) = max((Toes_MeanP(R_steps(j):R_ends(j))));
        MAXmaxToes(r+j) = max((Toes_MaxP(R_steps(j):R_ends(j))));
        
    end 
    
    r=length(contact_time);
    
end 
   

Titles = {'Subject', 'Shoe', 'Movement', 'ContactTime', 'meanHeel', 'sdHeel', 'maxHeel', 'MAXmaxHeel', 'meanMedMF', 'sdMedMF', 'maxMedMF', 'MAXmaxMedMF','meanLatMF', 'sdLatMF', 'maxLatMF', 'MAXmaxLatMF','meanMedFF', 'sdMedFF', 'maxMedFF', 'MAXmaxMedFF', 'meanLatFF', 'sdLatFF', 'maxLatFF', 'MAXmaxLatFF', 'meanToes', 'sdToes', 'maxToes', 'MAXmaxToes'};
data = [contact_time, meanHeel, sdHeel, maxHeel, MAXmaxHeel, meanMedMF, sdMedMF, maxMedMF, MAXmaxMedMF, meanLatMF, sdLatMF, maxLatMF, MAXmaxLatMF, meanMedFF, sdMedFF, maxMedFF, MAXmaxMedFF, meanLatFF, sdLatFF, maxLatFF, MAXmaxLatFF, meanToes, sdToes, maxToes, MAXmaxToes];
data = num2cell(data);
data = horzcat(Subject, Shoe, Movement, data);
data = vertcat(Titles, data);

writecell(data, 'CompiledPressureData.csv')
    
disp ('done!')

   
    
    
    