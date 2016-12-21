%==========================================================================
%
%   WRITTEN BY:  
%   Mohamed Amgad TagEl-Din
%   OIST Graduate University, Japan
%   Faculty of Medicine, Cairo University, Egypt   
%   
%   LAST EDITED: 
%   Nov 17th 2015
%
%   PURPOSE:
%   This script calculates K~ values along with the critical quantiles for
%   grayscale images in a pre-specified input directory and extension. 
%   The output is an individual K~ plot for each image and a combined excel
%   file containing the K profile at the full range of radii for each
%   image.
% 
%   SAMPLE INPUT:
%   A GUI is used that should makes it easier to input pre-specified
%   parameters.
%   
%   SAMPLE OUTPUT:
%   See Output directory after running a test run of the script.
%
%==========================================================================

clear all ; close all ; clc ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   GUI SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prompt = {'Specify the input directory  :  ' , ... 
          'Specify the output directory  :  ' , ...
          'Specify image extension  : ' ...
          'Start radius (pixels)  : ' ...
          'Radius increment (Pixels)  : ' ...
          'End radius (pixels)  : ' ...
          'Ignore Zero pixels? (1=Yes,0=No)'};
          
dlg_title = 'Calculating K~ for grayscale images (zero pixels included in field of view)';
num_lines = 1;

def = { 'ProjectInput\' , ...
        'ProjectOutput\' , ...
        '.tiff' , ...
        '3' , ...
        '1' , ...
        '9' ,...
        '0'};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='none';
    
answer = inputdlg(prompt,dlg_title,num_lines,def,options);

% for i = 2:length(answer(:,1))
%     answer{i,1} = str2num(answer{i,1});
% end
    
%%

InputDir = answer{1,1};
OutputDir = answer{2,1};
imExten = answer{3,1};
StartRadius = answer{4,1};
StepRadius = answer{5,1};
EndRadius = answer{6,1};
IgnoreZero = answer{7,1};

StartRadius = str2num(StartRadius);
StepRadius = str2num(StepRadius);
EndRadius = str2num(EndRadius);
IgnoreZero = str2num(IgnoreZero);

KRadius = StartRadius:StepRadius:EndRadius;

%%

AllIms = [InputDir,'*',imExten];
d1 = dir(fullfile(AllIms));

for i=1:length(d1)
    fpath1{i,1}=InputDir; fpath1{i,2}=d1(i).name;
end
 
fpathfull1=[strcat(fpath1(:,1),fpath1(:,2))]; 

%%

for imNo = 1:length(fpathfull1)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculating K~ and Critical Quantiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    im = double(imread(fpathfull1{imNo,1}));
    %im = dlmread(fpathfull1{imNo,1}); %If image was a floating point .txt file
    imName = fpath1{imNo,2};
    allK = zeros(length(KRadius),4);
    
    for Krad = KRadius
            
        clc
        fprintf(1,'Current K radius =  %d \n', Krad)
        fprintf(1,'            Out of  %d \n', EndRadius)
        fprintf(1,' %d \n', [])
        fprintf(1,['Current image : ',imName,' %d \n'], [])
        fprintf(1,['                (',num2str(imNo),'   Out of  ', num2str(length(fpathfull1(:,1))),') %d \n'], [])
        fprintf(1,' %d \n', [])
        
        if IgnoreZero == 0
        [K,KNC,Q01,Q99] = RipleysK1 (im,Krad);
        else
        [K,KNC,Q01,Q99] = RipleysK2 (im,Krad);
        end
        
        allK(Krad-StartRadius+1,1) = Krad;
        allK(Krad-StartRadius+1,2) = KNC;
        allK(Krad-StartRadius+1,3) = Q01;
        allK(Krad-StartRadius+1,4) = Q99; 
        
    end
    
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Saving K profile numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[k_m,k_n] = size(allK);
KCell = zeros(k_m+1,k_n);
KCell(2:k_m+1,1:k_n)=allK;
KCell=num2cell(KCell);
KCell{1,1}='Radius';KCell{1,2}='K~';KCell{1,3}='Q01';KCell{1,4}='Q99';

xlswrite([OutputDir,'KProfiles','.xlsx'],KCell,imName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Saving K profile figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%h=figure('visible','on');
h=figure('visible','off');

plot (allK(:,2), 'b - .')
grid off
title (['K~ Profile for ',imName])
xlabel('Radius')
labels = allK(:,1)';
set(gca, 'XTick', 1:length(labels)); % Change x-axis ticks
set(gca, 'XTickLabel', labels); % Change x-axis ticks labels
ylabel('K~')

hold on
plot (allK(:,3:4), 'r --')
hold off

saveas(h,[OutputDir,'KProfile_',imName,'.tiff'],'tiff');

end