%==========================================================================
%
%   WRITTEN BY:  
%   Mohamed Amgad TagEl-Din
%   OIST Graduate University, Japan
%   Faculty of Medicine, Cairo University, Egypt   
%   
%   LAST EDITED: 
%   Aug 30th 2015
%
%   PURPOSE:
%   This script creates randomly-located, variable-sized disk-shaped 
%   clusters (without overlap).
%   This algorithm modifies the number of aggregates based on the
%   input parameters, including the aggregate-to-diffuse ratio (ADR) and 
%   the signal-to-background ratio (SBR). Note that the signal is added to
%   the background, so at SBR = 1.5, the pixels belonging to aggregates 
%   will become background pixels plus 1.5X background pixels (i.e. 
%   overall, 2.5X background pixels).
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

prompt = {'Specify the output directory  :  ' , ... 
          'Number of images per set  : ' , ...
          'Minimum ADR  : ' ...
          'ADR increment  : ' ...
          'Maximum ADR : ' ... 
          'Image length (pixels) : ' ...
          'Event density (intensity units per pixel) : ' ...
          'SBR : ' ...
          'Minimum cluster offset : ' ...
          'Minimum cluster radius : ' ...
          'Cluster radius offset : ' ...
          'Maximum cluster radius : ' ...
          'PSF - diameter (diameter of Gaussian blur) : ' ...
          'PSF - diameter (sigma of Gaussian blur) : ' ...
          'SNR : ' };
          
dlg_title = 'Simulating Clustering (SBR constant)';
num_lines = 1;

def = { 'ProjectOutput\' , ...
        '3' , ...
        '0.02' , ...
        '0.02' , ...
        '0.16' , ...
        '256' , ...
        '10' , ...
        '0.4' , ...
        '1' , ...
        '2' , ...
        '2' , ...
        '4' , ...
        '3' , ...
        '1' , ...
        '10' };

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='none';
    
answer = inputdlg(prompt,dlg_title,num_lines,def,options);

for i = 2:length(answer(:,1))
    answer{i,1} = str2num(answer{i,1});
end
    
%%

OutputDir = answer{1,1};

NoOfImages = answer{2,1}; %No of images per set

StartADR = answer{3,1}; %Aggregate-to-diffuse ratio
StepADR = answer{4,1}; 
EndADR = answer{5,1};

imLength = answer{6,1};
ImIntens = answer{7,1}; %Average intensity per pixel
 
SBR = answer{8,1}; %Signal-to-background ratio
 
AggOffset = answer{9,1}; %minimum distance between edge of one aggregate and another

StartRadius = answer{10,1}; %Aggregate radius
StepRadius = answer{11,1};
EndRadius = answer{12,1};
 
PSF_d = answer{13,1}; %Point spread function - diameter of gaussian blur filter
PSF_s = answer{14,1}; %Point spread function - sigma of gaussian blur filter
SNR = answer{15,1}; %Signal-to-noise ratio


%%
TotInt = ImIntens * (imLength^2);

OutputDir = {OutputDir};
mkdir(OutputDir{1,1},'Images');
mkdir(OutputDir{1,1},'GroundTruth(Combined)');
mkdir(OutputDir{1,1},'GroundTruth(PerRadius)');
 
%%
 
Parameters {1,1} = 'Image No';

Parameters {1,2} = 'Aggregate-to-diffuse ratio (ADR)';

Parameters {1,3} = 'Image length';
Parameters {1,4} = 'Average no of particles per pixel';

Parameters {1,5} = 'Aggregation offset (min distance between aggregates)';
Parameters {1,6} = 'Radius of aggregates';
Parameters {1,7} = 'Frequency of aggregates';

Parameters {1,8} = 'Signal-to-background ratio (SBR)';
Parameters {1,9} = 'Post-Hoc (Empirical) SBR';

Parameters {1,10} = 'PSF - diameter of gaussian blur';
Parameters {1,11} = 'PSF - sigma of gaussian blur';
Parameters {1,12} = 'Signal-to-Noise ratio (SNR)'; 


 
%%
 
ADR_All = StartADR:StepADR:EndADR;
 
 
for SetNo = 1:length(ADR_All)
 
ADR = ADR_All(1,SetNo);    
 
TotAggInt = (TotInt * ADR) / (1+ADR) ; %Total aggregate intensity
NpAgg = floor( (ADR*(imLength^2)) / (SBR+ADR) ); %Number of pixels assigned to aggregates
    
radii = StartRadius:StepRadius:EndRadius;
 
%%
%
%   Creating structuring elements for aggregates at different sizes 
%
 
StrelSums = zeros(length(radii(:)),1);
 
for s = 1:length(radii)

strel_radius = radii(1,s);
    
%
%   Aggregate "template"
%

% IMPORTANT NOTE:
%
% This section was modified from 
% Peter Kovesi's function: CIRCULARSTRUCT
% Function to construct a circular structuring element
% for morphological operations.
%
% Peter Kovesi   March 2000

dia = ceil(2*strel_radius);  % Diameter of structuring element
 
if mod(dia,2) == 0     % If diameter is a odd value
 dia = dia + 1;        % add 1 to generate a `centre pixel'
end
 
r = fix(dia/2);
[x,y] = meshgrid(-r:r);
rad = sqrt(x.^2 + y.^2);  
strel = rad <= r;
strel = double(strel);

Strels{s,1} = strel;
StrelSums(s,1) = sum(strel(:));
 
end
 
%%
%
%   Assigning each aggregate radius a frequency
%
 
AggRemaining = NpAgg;
StrelFreq = zeros(size(StrelSums));
 
for t = 0:length(StrelSums(:,1))-1
 
if sum(StrelSums(1:length(StrelSums(:,1))-t,1)) > AggRemaining
continue
end
 
a = floor( AggRemaining / sum(StrelSums(1:length(StrelSums(:,1))-t,1)) );
AggRemaining = AggRemaining - a* sum(StrelSums(1:length(StrelSums(:,1))-t,1));
 
StrelFreq(1:length(StrelSums(:,1))-t,1) = StrelFreq(1:length(StrelSums(:,1))-t,1) + a;
 
end
 
%% 
for ImNo = 1:NoOfImages
    
 clc
fprintf(1,'Current image no. %d \n', ImNo)
fprintf(1,'           Out of %d \n', NoOfImages )
fprintf(1,' %d \n', [] )
fprintf(1,'Current set no. %d \n', SetNo)
fprintf(1,'         Out of %d \n', length(ADR_All) )
fprintf(1,' %d \n', [] )    
 
%%
 
%
%   Initializing
%
 
%Main field
N = zeros(imLength,imLength);
[m,n] = size(N);
[X,Y] = meshgrid(1:imLength); %X- and Y- Coordinates
 
%Initializing the Inhibition map at for aggregates at different radii
 
for i = 1:length(radii)
    
    MaxAggR = max(radii(:));
    InhD = ceil(2*(MaxAggR+radii(1,i))) + 2*AggOffset;
    
    if mod(InhD,2) == 0     % If diameter is a odd value
    InhD = InhD + 1;        % add 1 to generate a `centre pixel'
    end
 
    InhR = fix(InhD/2);
    
    %Adding an artificial margin of AggR around the inhibition maps (beside the AggR zone towards the inside)
    %to allow for marginal aggregates (NOTE: this essentially means you'll have to
    %adjust the coordinates of the aggregates)
    InhibitionMap{i,1} = ones((m+2*InhR),(n+2*InhR));
    
    InhibitionMap{i,1}(2*InhR+1:m,2*InhR+1:n) = 0;
    
    InhibitionMap{i,2} = InhR;
   
end
 
%%
 
for ii = 1:length(radii)
 
radius = radii(1,ii);
 
AggNo = StrelFreq(ii,1);
 
%Radius-specific ground truth image
RGTruth = zeros(m,n);
 
%
% Aggregation "template"
%
strel = Strels{ii,1};
 
for k = 1:length(radii)
 
%
%   Aggregate Inhibition zones
%

Inhdia = ceil(2*(radius+radii(1,k))) + 2*AggOffset;  % Diameter of inhibition structuring element
 
 
if mod(Inhdia,2) == 0     % If diameter is a odd value
 Inhdia = Inhdia + 1;        % add 1 to generate a `centre pixel'
end
 
InhRad = fix(Inhdia/2);
[InhX,InhY] = meshgrid(-InhRad:InhRad);
InhRad2 = sqrt(InhX.^2 + InhY.^2);  
InhStrel = InhRad2 <= (InhRad);
InhStrel = double(InhStrel);
 
%%
 
InhibitionMap{k,3} = InhStrel; %Inhibition template
InhibitionMap{k,4} = InhRad; %Inhibition radius
 
end
 
 
%%
for j = 1:AggNo
 
%Listing the set of potential locations for center of aggregate
clear PotentialLocations
PotentialLocations = 1:length(N(:));
PotentialLocations = reshape(PotentialLocations,m,n);
 
%CHECK THIS!!!!!
InhMargin = InhibitionMap{ii,2};
InhTemp = InhibitionMap{ii,1}(InhMargin+1:m+InhMargin,InhMargin+1:n+InhMargin);
 
PotentialLocations(InhTemp==1) = 0;
 
PotentialLocations = PotentialLocations(:);
PotentialLocations(PotentialLocations==0)=[];
 
%Randomly choosing one of those potential locations to place our aggregate
clear Loc
Loc = randi(length(PotentialLocations));
Loc = PotentialLocations(Loc);
%Getting the X- Y- coordinate of corresponding location
Loc_X = X(Loc);
Loc_Y = Y(Loc);
 
%
%   Placing the aggregate at that location
%
 
 
Window = N((Loc_Y-radius:Loc_Y+radius),(Loc_X-radius:Loc_X+radius));
Window_temp = Window;
        
        
%Only considering elements within SE
Window = ones(size(Window));
Window(strel==0)=0;
        
Window = Window+Window_temp;                
%Replacing the values of M with the new values
N((Loc_Y-radius:Loc_Y+radius),(Loc_X-radius:Loc_X+radius)) = Window;
 
%Adding the aggregate at the corresponding location in the radius-specific
%ground truth image
RGTruth((Loc_Y-radius:Loc_Y+radius),(Loc_X-radius:Loc_X+radius)) = Window;
 
%%
%
%   Updating the inhibition maps
%
 
for i = 1:length(radii)
 
strel1 = InhibitionMap{i,3};
InhR1 = InhibitionMap{i,4};
 
%Correcting the X,Y Coordinates of shuffled locations so that they relate
%to their new positions in the inhibition map (which has an extra margin added)
Loc_X1 = Loc_X + InhibitionMap{i,2};
Loc_Y1 = Loc_Y + InhibitionMap{i,2};
 
Window1 = InhibitionMap{i,1}((Loc_Y1-InhR1:Loc_Y1+InhR1),(Loc_X1-InhR1:Loc_X1+InhR1));
Window_temp1 = Window1;
 
%Only considering elements within SE
Window1 = ones(size(Window1));
Window1(strel1==0)=0;
        
Window1 = Window1+Window_temp1;                
%Replacing the values of Inhibition map with the new values
InhibitionMap{i,1}((Loc_Y1-InhR1:Loc_Y1+InhR1),(Loc_X1-InhR1:Loc_X1+InhR1)) = Window1;
 
InhibitionMap{i,1}(InhibitionMap{i,1}>0) = 1;
 
end
 
end

%
%   Writing per-radius ground truth
%
imwrite(RGTruth,[OutputDir{1,1},'GroundTruth(PerRadius)/Set',num2str(SetNo),'Im',num2str(ImNo),'GTruthR',num2str(radius),'.tiff'])

end
 
 
%%
% 
% Randomly allocating particles to either 
% aggregates or background at pre-defined 
% densities
%
 
%
%   Creating a background to superimpose aggregates on
%
 
Positions = randi([1 imLength^2],1,floor(TotInt-TotAggInt));                               
Positions(1,1)=1;
Positions(1,2)=imLength^2;
 
% Getting the position "count" by getting a frequency histogram
M = tabulate(Positions(:));
M=M(:,2);
M=reshape(M,m,n);
 
 
%
%   Assigning particles assigned to aggregates to their respective
%   locations randomly
%
 
% Coordinates of all points
AllCoordinates = 1:(imLength^2);
AllCoordinates = reshape(AllCoordinates,m,n);
 
% Coordinates of aggregates
AggCoordinates = AllCoordinates;
AggCoordinates(N==0) = 0;
AggCoordinates(AggCoordinates==0) = [];
 
AggPositions = randi([1 length(AggCoordinates)],1,floor(TotAggInt));                   
AggPositions(1,1)=1;                                                                                    
AggPositions(1,2)=length(AggCoordinates);                                                     
 
% Getting the position "count" by getting a frequency histogram
AllAggs = tabulate(AggPositions(:));
AllAggs(:,1) = AggCoordinates';
 
% Replacing pixels at aggregate locations with their assigned particle
% count values
MAgg = changem(AllCoordinates,AllAggs(:,2),AllAggs(:,1));
MAgg(N==0)=0;
 
%
% Calculating the POST-HOC SBR after the "quantization" made earlier to
% "squeeze" aggregate intensity within aggregated of pre-specified radii
% 

Agg_PostHoc = sum(MAgg(:))./sum(N(:));                                     
Bck_PostHoc = sum(M(:))./((m*n)-sum(N(:)));                                
 
SBR_PostHoc = Agg_PostHoc ./ Bck_PostHoc;                                  
 
% Adding aggregates to background of "diffuse" protein
M = MAgg+M;
 
%%
%
%   Multiplying the "base" ("ideal") image by a point spread function (PSF)
%   as happens in actual confocal microscopy
%
 
%  
H = fspecial('gaussian',PSF_d,PSF_s);
M = imfilter(M,H,'symmetric');
% %Rescaling to same total intensity
M = M ./ sum(M(:));
M = M .* TotInt;
 
%%
%
%   Adding poisson noise
%
 
%Adding poisson noise according to prespecified SNR
[M,noise,SNR_measured]= genSignalForSNR(M,SNR);
 
%%
 
if max(M(:)) > 255
    error('Some pixels have passed the saturation limit, please do one or more of the following: a- Reduce ImIntens (preferable); b- Reduce ADR; c-Reduce SBR.');
end
 
%Rescaling so that the total intensity in the normalized image (where
%intensity ranges from 0 to 1) remains the same so that different images
%could be compared as long as ImIntens remains the same. 
MMat = M ./ 255;
 
%%
 
Parameters {ImNo+1,1} = ImNo;

Parameters {ImNo+1,2} = ADR;

Parameters {ImNo+1,3} = imLength;
Parameters {ImNo+1,4} = ImIntens;
 
Parameters {ImNo+1,5} = AggOffset;
Parameters {ImNo+1,6} = num2str(radii);
Parameters {ImNo+1,7} = num2str(StrelFreq');
 

Parameters {ImNo+1,8} = SBR;
Parameters {ImNo+1,9} = SBR_PostHoc;
 
Parameters {ImNo+1,10} = PSF_d;
Parameters {ImNo+1,11} = PSF_s;
Parameters {ImNo+1,12} = SNR;
 
%%
%
%   Writing output and image-specific parameters and ground truth
%

imwrite(MMat,[OutputDir{1,1},'Images/Set',num2str(SetNo),'Im',num2str(ImNo),'.tiff'])
imwrite(N,[OutputDir{1,1},'GroundTruth(Combined)/Set',num2str(SetNo),'Im',num2str(ImNo),'GTruth.tiff'])
 
xlswrite([OutputDir{1,1},'Parameters.xlsx'],Parameters,['SetNo = ',num2str(SetNo)])
 
end
end