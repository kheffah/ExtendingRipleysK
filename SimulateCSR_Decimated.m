function CSR_Decimated = SimulateCSR_Decimated (imLength,imIntens,NAgg,AggOffset,StartRadius,StepRadius,EndRadius)

%==========================================================================
%
%   WRITTEN BY:  
%   Mohamed Amgad TagEl-Din
%   OIST Graduate University, Japan
%   Faculty of Medicine, Cairo University, Egypt   
%   
%   LAST EDITED: 
%   Aug 25th 2015
%
%   PURPOSE:
%   This function generates a homogenous poisson distribution of
%   events in a field of view (Complete Spatial Randomness), allowing 
%   multiple particles to co-localize at the same pixel. 
%   This modified function decimates the field of view to create an 
%   irregular field to test the shape dependency of the K function.
% 
%   SAMPLE INPUT:
%
%   imLength = 2^8 ;  %(length of the square field of view)
%   imIntens = 10 ;   %(intensity of particles)
%   NAgg = 3 ;        %(number of decimations of each radius)
%   AggOffset = 1;    %(minimum offset between decimations)
%   StartRadius = 10; %(minimum radius of decimations)
%   StepRadius = 4;   %(increment by which decimation radius increases)
%   EndRadius = 18;   %(maximum radius of decimations)
%   
%   CSR_Decimated = SimulateCSR_Decimated
%   (imLength,imIntens,NAgg,AggOffset,StartRadius,StepRadius,EndRadius);
%
%   SAMPLE OUTPUT:
%
%   whos CSR_Decimated
%   Name               Size        Bytes     Class     Attributes
%   CSR_Decimated      256x256     524288    double     
%
%==========================================================================




%%
TotParticles = ceil((imLength^2)*imIntens) ; 
radii = StartRadius:StepRadius:EndRadius;

%%
%
%   Creating structuring elements for decimations at different sizes 
%
 
StrelSums = zeros(length(radii(:)),1);
 
for s = 1:length(radii)
 
strel_radius = radii(1,s);
    
%%
%
% Decimation "template"
%

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
 
%%
 
Strels{s,1} = strel;
StrelSums(s,1) = sum(strel(:));
 
end
 
%%
%
%   Assigning each decimation radius a frequency
%

NpAgg = NAgg .* sum(StrelSums(:)); %Number of pixels assigned to decimations   
 
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
 
%
%   Initializing
%
 
%Main field
N = zeros(imLength,imLength);
[m,n] = size(N);
[X,Y] = meshgrid(1:imLength); %X- and Y- Coordinates
 
%Initializing the Inhibition map at for decimations at different radii
 
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
% Decimation "template"
%
strel = Strels{ii,1};
 
for k = 1:length(radii)
 
%
%   Decimation Inhibition zones
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
 
%====================================
 
InhibitionMap{k,3} = InhStrel; %Inhibition template
InhibitionMap{k,4} = InhRad; %Inhibition radius
 
end
 
 
%%
for j = 1:AggNo
 
%Listing the set of potential locations for center of decimation
clear PotentialLocations
PotentialLocations = 1:length(N(:));
PotentialLocations = reshape(PotentialLocations,m,n);
 

InhMargin = InhibitionMap{ii,2};
InhTemp = InhibitionMap{ii,1}(InhMargin+1:m+InhMargin,InhMargin+1:n+InhMargin);
 
PotentialLocations(InhTemp==1) = 0;
 
PotentialLocations = PotentialLocations(:);
PotentialLocations(PotentialLocations==0)=[];
 
%Randomly choosing one of those potential locations to place our decimation
clear Loc
Loc = randi(length(PotentialLocations));
Loc = PotentialLocations(Loc);
%Getting the X- Y- coordinate of corresponding location
Loc_X = X(Loc);
Loc_Y = Y(Loc);
 
%
%   Placing the decimation at that location
%
 
Window = N((Loc_Y-radius:Loc_Y+radius),(Loc_X-radius:Loc_X+radius));
Window_temp = Window;
                
%Only considering elements within SE
Window = ones(size(Window));
Window(strel==0)=0;
        
Window = Window+Window_temp;                
%Replacing the values of M with the new values
N((Loc_Y-radius:Loc_Y+radius),(Loc_X-radius:Loc_X+radius)) = Window;
 
%Adding the decimation at the corresponding location in the radius-specific
%ground truth image
RGTruth((Loc_Y-radius:Loc_Y+radius),(Loc_X-radius:Loc_X+radius)) = Window;
 
 
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
  
end

%%

%
%   Assigning particles to their respective
%   locations randomly
%

% Coordinates of all points
AllCoordinates = 1:(imLength^2);
AllCoordinates = reshape(AllCoordinates,m,n);
 
% Coordinates of everything EXCEPT decimations
AggCoordinates = AllCoordinates;
AggCoordinates(N==1) = 0;
AggCoordinates(AggCoordinates==0) = [];
 
AggPositions = randi([1 length(AggCoordinates)],1,floor(TotParticles));                   
AggPositions(1,1)=1;                                                                                     
AggPositions(1,2)=length(AggCoordinates);                                                    
 
% Getting the position "count" by getting a frequency histogram
AllAggs = tabulate(AggPositions(:));
AllAggs(:,1) = AggCoordinates';
 
% Replacing pixels at aggregate locations with their assigned particle
% count values
MAgg = changem(AllCoordinates,AllAggs(:,2),AllAggs(:,1));
MAgg(N==1)=0;


CSR_Decimated = MAgg;

end