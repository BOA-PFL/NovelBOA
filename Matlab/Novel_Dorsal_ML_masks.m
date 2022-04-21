clear
clc

addpath('C:\Users\kate.harrison\Documents\GitHub\NovelBOA\NovelBOA\Matlab')
input_dir = 'C:\Users\kate.harrison\Dropbox (Boa)\EndurancePerformance\TNF_Scrambler_Apr_21\Novel_Data';

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
maxHeel = zeros(2,1);

sdLatFF = zeros(2,1);
meanLatFF = zeros(2,1);
maxLatFF = zeros(2,1);

sdMedFF = zeros(2,1);
meanMedFF = zeros(2,1);
maxMedFF = zeros(2,1);

sdToes = zeros(2,1);
meanToes = zeros(2,1);
maxToes = zeros(2,1);

sdDorsMedFF = zeros(2,1);
meanDorsMedFF = zeros(2,1);
maxDorsMedFF = zeros(2,1);

sdDorsLatFF = zeros(2,1);
meanDorsLatFF = zeros(2,1);
maxDorsLatFF = zeros(2,1);

sdDorsMedMF = zeros(2,1);
meanDorsMedMF = zeros(2,1);
maxDorsMedMF = zeros(2,1);

sdDorsLatMF = zeros(2,1);
meanDorsLatMF = zeros(2,1);
maxDorsLatMF = zeros(2,1);

contact_time = zeros(2,1);
r = 0;

for i = 1:NumbTrials
    
    FileName = dataList(f(i));
    FileLoc = char(strcat(input_dir,'\', FileName));
    data = dlmread(FileLoc, '\t', 18 , 0);
    
    names = split(FileName, ["_", " ", "."]);
    sub = names{1};
    shoe = names{2};
    move = names{3};
        
    time = data(:,1);
    
    Heel_Force = data(:,2);
%     Heel_Force(Heel_Force == 0) = NaN;
%     Heel_Force = interp1gap(Heel_Force, 10, 'spline','interpval', 0, 'extrap', 0);
    
    Heel_MeanP = data(:,4);
%     Heel_MeanP(Heel_MeanP == 0) = NaN;
%     Heel_MeanP = interp1gap(Heel_MeanP, 10, 'spline','interpval', 0, 'extrap', 0);
    
    LatFF_Force = data(:,6);
%     LatFF_Force(LatFF_Force == 0) = NaN;
%     LatFF_Force = interp1gap(LatFF_Force, 10, 'spline','interpval', 0, 'extrap', 0);
   
    LatFF_MeanP = data(:,8);
%     LatFF_MeanP(LatFF_MeanP == 0) = NaN;
%     LatFF_MeanP = interp1gap(LatFF_MeanP, 10, 'spline','interpval', 0, 'extrap', 0);
        
    MedFF_Force = data(:,10);
%     MedFF_Force(MedFF_Force == 0) = NaN;
%     MedFF_Force = interp1gap(MedFF_Force, 10, 'spline','interpval', 0, 'extrap', 0);
    
    MedFF_MeanP = data(:,12);
%     MedFF_MeanP(MedFF_MeanP == 0) = NaN;
%     MedFF_MeanP = interp1gap(MedFF_MeanP, 10, 'spline','interpval', 0, 'extrap', 0);
    
    Toes_Force = data(:,14);
%     Toes_Force(Toes_Force == 0) = NaN;
%     Toes_Force = interp1gap(Toes_Force, 10, 'spline','interpval', 0, 'extrap', 0);
    
    Toes_MeanP = data(:,16);
%     Toes_MeanP(Toes_MeanP == 0) = NaN;
%     Toes_MeanP = interp1gap(Toes_MeanP, 10, 'spline','interpval', 0, 'extrap', 0);
    
    DorsMedFF_Force = data(:,18);
%     DorsMedFF_Force(DorsMedFF_Force == 0) = NaN;
%     DorsMedFF_Force = interp1gap(DorsMedFF_Force, 10, 'spline','interpval', 0, 'extrap', 0);
    
    DorsMedFF_MeanP = data(:,20);
%     DorsMedFF_MeanP(DorsMedFF_MeanP == 0) = NaN;
%     DorsMedFF_MeanP = interp1gap(DorsMedFF_MeanP, 10, 'spline','interpval', 0, 'extrap', 0);
    
    DorsLatFF_Force = data(:,22);
%     DorsLatFF_Force(DorsLatFF_Force == 0) = NaN;
%     DorsLatFF_Force = interp1gap(DorsLatFF_Force, 10, 'spline','interpval', 0, 'extrap', 0);
    
    DorsLatFF_MeanP = data(:,24);
%     DorsLatFF_MeanP(DorsLatFF_MeanP == 0) = NaN;
%     DorsLatFF_MeanP = interp1gap(DorsLatFF_MeanP, 10, 'spline','interpval', 0, 'extrap', 0);
    
    DorsMedMF_Force = data(:,26);
%     DorsMedMF_Force(DorsMedMF_Force == 0) = NaN;
%     DorsMedMF_Force = interp1gap(DorsMedMF_Force, 10, 'spline','interpval', 0, 'extrap', 0);
    
    DorsMedMF_MeanP = data(:,28);
%     DorsMedMF_MeanP(DorsMedMF_MeanP == 0) = NaN;
%     DorsMedMF_MeanP = interp1gap(DorsMedMF_MeanP, 10, 'spline','interpval', 0, 'extrap', 0);
    
    DorsLatMF_Force = data(:,30);
%     DorsLatMF_Force(DorsLatMF_Force == 0) = NaN;
%     DorsLatMF_Force = interp1gap(DorsLatMF_Force, 10, 'spline','interpval', 0, 'extrap', 0);
    
    DorsLatMF_MeanP = data(:,32);
%     DorsLatMF_MeanP(DorsLatMF_MeanP == 0) = NaN;
%     DorsLatMF_MeanP = interp1gap(DorsLatMF_MeanP, 10, 'spline','interpval', 0, 'extrap', 0);
    
    Force_Tot = Heel_Force + MedFF_Force + LatFF_Force + Toes_Force; % sum force under each mask to find total force
    
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
                             
        meanMedFF(r+j) = mean((MedFF_MeanP(R_steps(j):R_ends(j))));
        sdMedFF(r+j) = std((MedFF_MeanP(R_steps(j):R_ends(j))));
        maxMedFF(r+j) = max((MedFF_MeanP(R_steps(j):R_ends(j))));
                
        meanLatFF(r+j) = mean((LatFF_MeanP(R_steps(j):R_ends(j))));
        sdLatFF(r+j) = std((LatFF_MeanP(R_steps(j):R_ends(j))));
        maxLatFF(r+j) = max((LatFF_MeanP(R_steps(j):R_ends(j))));
                
        meanToes(r+j) = mean((Toes_MeanP(R_steps(j):R_ends(j))));
        sdToes(r+j) = std((Toes_MeanP(R_steps(j):R_ends(j))));
        maxToes(r+j) = max((Toes_MeanP(R_steps(j):R_ends(j))));
                
        meanDorsMedFF(r+j) = mean((DorsMedFF_MeanP(R_steps(j):R_ends(j))));
        sdDorsMedFF(r+j) = std((DorsMedFF_MeanP(R_steps(j):R_ends(j))));
        maxDorsMedFF(r+j) = max((DorsMedFF_MeanP(R_steps(j):R_ends(j))));
                
        meanDorsLatFF(r+j) = mean((DorsLatFF_MeanP(R_steps(j):R_ends(j))));
        sdDorsLatFF(r+j) = std((DorsLatFF_MeanP(R_steps(j):R_ends(j))));
        maxDorsLatFF(r+j) = max((DorsLatFF_MeanP(R_steps(j):R_ends(j))));
        
        meanDorsMedMF(r+j) = mean((DorsMedMF_MeanP(R_steps(j):R_ends(j))));
        sdDorsMedMF(r+j) = std((DorsMedMF_MeanP(R_steps(j):R_ends(j))));
        maxDorsMedMF(r+j) = max((DorsMedMF_MeanP(R_steps(j):R_ends(j))));
                
        meanDorsLatMF(r+j) = mean((DorsLatMF_MeanP(R_steps(j):R_ends(j))));
        sdDorsLatMF(r+j) = std((DorsLatMF_MeanP(R_steps(j):R_ends(j))));
        maxDorsLatMF(r+j) = max((DorsLatMF_MeanP(R_steps(j):R_ends(j))));
        
    end 
    
    r=length(contact_time);
    
end 
   

Titles = {'Subject', 'Shoe', 'Movement', 'ContactTime',...
    'meanHeel', 'sdHeel', 'maxHeel',...
    'meanMedFF', 'sdMedFF', 'maxMedFF',...
    'meanLatFF', 'sdLatFF', 'maxLatFF',...
    'meanToes', 'sdToes', 'maxToes',...
    'meanDorsMedFF', 'sdDorsMedFF', 'maxDorsMedFF',...
    'meanDorsLatFF', 'sdDorsLatFF', 'maxDorsLatFF',...
    'meanDorsMedMF', 'sdDorsMedMF', 'maxDorsMedMF',...
    'meanDorsLatMF', 'sdDorsLatMF', 'maxDorsLatMF'};
data = [contact_time,...
    meanHeel, sdHeel, maxHeel,...
    meanMedFF, sdMedFF, maxMedFF, ...
    meanLatFF, sdLatFF, maxLatFF, ...
    meanToes, sdToes, maxToes,...
    meanDorsMedFF, sdDorsMedFF, maxDorsMedFF,...
    meanDorsLatFF, sdDorsLatFF, maxDorsLatFF, ...
    meanDorsMedMF, sdDorsMedMF, maxDorsMedMF,...
    meanDorsLatMF, sdDorsLatMF, maxDorsLatMF
    ];
data = num2cell(data);
data = horzcat(Subject, Shoe, Movement, data);
data = vertcat(Titles, data);

writecell(data, 'CompiledPressureData.csv')
    
disp ('done!')