function [K,K_NC,Q01,Q99] = RipleysK1 (M,radius)

%==========================================================================
%
%   WRITTEN BY:  
%   Mohamed Amgad TagEl-Din
%   OIST Graduate University, Japan
%   Faculty of Medicine, Cairo University, Egypt   
%   
%   LAST EDITED: 
%   Aug 29th 2015
%
%   PURPOSE:
%   This function calculates Ripley's K-function for a rectangular
%   grayscale field of view. It does NOT ignore zero pixels. That is, it
%   considers zero pixels to be regions devoid of events rather than
%   excluding them from the analysis. It also calculated the critical
%   quantiles using the Cornishe-Fisher expansion, as proposed and
%   validated by Lagache et al [1].
%
%   1. Analysis of the spatial organization of molecules with robust 
%      statistics.
%      Lagache T, Lang G, Sauvonnet N, Olivo-Marin JC.
%      PLoS One. 2013 Dec 4;8(12):e80914. doi: 10.1371/
%      journal.pone.0080914. eCollection 2013.
%
%   SAMPLE INPUT:
%   
%   im = double(imread('TestImage.tiff'));
%   radius = 3; %(radius at which K-NC is calculated)
%
%   [K,K_NC,Q01,Q99] = RipleysK1 (im,radius);
%
%   SAMPLE OUTPUT:
% 
% K =               %(Raw K-function, not normalized or centralized)
% 
%    28.5510
% 
% 
% K_NC =            %(Normalized and centralized K-function)
% 
%   102.3934
% 
% 
% Q01 =             %(1st quantile)
% 
%    -2.2983
% 
% 
% Q99 =             %(99th quantile)
% 
%     2.3588    
%
%==========================================================================


%%
%
% IMPORTANT NOTE:
%
% This section was modified from 
% Peter Kovesi's function: CIRCULARSTRUCT
% Function to construct a circular structuring element
% for morphological operations.
%
% Peter Kovesi   March 2000

if radius < 1
  error('radius must be >= 1');
end

dia = ceil(2*radius);  % Diameter of structuring element

if mod(dia,2) == 0     % If diameter is a odd value
 dia = dia + 1;        % add 1 to generate a `centre pixel'
end

r = fix(dia/2);
[x,y] = meshgrid(-r:r);
rad = sqrt(x.^2 + y.^2);  
strel = rad <= radius;


%%
%
%   Adding an artificial margin around the image (later, edge correction):
%

[m, n] = size(M);

N = nan((m+2*radius),(n+2*radius));
N ((radius+1:m+radius),(radius+1:n+radius)) = M;

%%
%
%   Getting an "aggregation map":
%

strel = double(strel);
strel(strel==0) = nan;

AggMap = zeros(size(N));


for i = radius+1:m+radius
    for j = radius+1:n+radius
        
        Window = N((i-radius):(i+radius),(j-radius):(j+radius));
        Window(isnan(strel)==1) = 0; %only considering elements within CIRCULAR SE
 
        %
        % Specifying center and surround
        %
         
        [mS, nS] = size(strel);
        equator_X = (mS/2)+0.5;
        equator_Y = (nS/2)+0.5;
        
        Center = Window(equator_X,equator_Y);
        
        Surround = Window;
        Surround(equator_X,equator_Y)= 0; %removing central pixel
        
        Surround(isnan(Surround)==1) = 0;
        
        
        %
        % Edge correction: Besag's correction, where we divide by 
        % proportion of area included within ROI.
        %
        
        Window_Binary = Window;
        Window_Binary = Window_Binary + 1 ; 
        Window_Binary(isnan(Window_Binary)==1) = 0;
        Window_Binary(Window_Binary>0) = 1;
        
        b = Window_Binary+strel;
        b(b<2) = 0;
        b(b>1) = 1;
        
        b(isnan(b)==1)=0;
        strel1 = strel;
        strel1(isnan(strel1)==1) = 0;
        
        EdgeCorrection = sum(b(:)) / (pi*(radius^2));
        
        %
        % Building an aggregation map
        %
        
        if Center >0
      
        AggMap(i,j) = (Center*(Center - 1) + Center*(sum(Surround(:))))...
           ./ EdgeCorrection;
        
        end
        
    end
end

AggMap = AggMap((radius+1:m+radius),(radius+1:n+radius)); %trimming off xs

%%

%
%   Calculating K statistics:
%

M_Sum = sum(M(:));
M_Area = length(M(:));

K = sum(AggMap(:)) * (M_Area / (M_Sum*(M_Sum-1)));

%%

%
% Calculating the perimeter of the ROI
% 
 
[mM,nM] = size(M);

M_Perimeter = 2*(mM+nM);

% 
% 
% Calculating Variance of K
% 
% 
BetaR = (pi*(radius^2)) / M_Area;
GammaR = (M_Perimeter * radius) / M_Area;

Var_K = ( (2*(M_Area^2) * BetaR)/(M_Sum^2) ) * ...
            ( 1 + 0.305*GammaR + BetaR*( -1+0.0132*M_Sum*GammaR ) );

% 
% 
% Calculating Normalized and Centered K
% 

K_NC = (K - (pi*(radius^2))) / (sqrt(Var_K));

% 
% 
% Estimating the critical quantiles
% 
% 
Skew_K = ( (4*(M_Area^3)*BetaR) / ((M_Sum^4)*(Var_K^(3/2))) ) * ...
            (1 + 0.76*GammaR + ...
            M_Sum*BetaR*(1.173+0.414*GammaR)+ ...
            M_Sum*(BetaR^2)*(-2+0.012*M_Sum*GammaR));
                
Kurt_K = ( ((M_Area^4)*BetaR) / ((M_Sum^6)*(Var_K^2)) ) * ...
            (8 + 11.52*GammaR + ... 
            M_Sum*BetaR*((104.3+12*M_Sum)+ ...
            (78.7+7.32*M_Sum)*GammaR+1.116*M_Sum*(GammaR^2))+ ...
            M_Sum*(BetaR^2)* ((-304.3-1.92*M_Sum)+ ...
            (-97.9+2.69*M_Sum+0.317*(M_Sum^2))*GammaR+ ...
            0.0966*(M_Sum^2)*(GammaR^2))+ ...
            (M_Sum^2)*(BetaR^3)*(-36+0.0021*(M_Sum^2)*(GammaR^2)));

Z01 = norminv(0.01,0,1);
Z99 = norminv(0.99,0,1);

Q01 = Z01 + (1/6)*(Z01^2-1)*Skew_K + ...
    (1/24)*((Z01^3)-3*Z01)*(Kurt_K-3) - ...
    (1/36)*(2*(Z01^3)-5*Z01)*(Skew_K^2);

Q99 = Z99 + (1/6)*(Z99^2-1)*Skew_K + ...
    (1/24)*((Z99^3)-3*Z99)*(Kurt_K-3) - ...
    (1/36)*(2*(Z99^3)-5*Z99)*(Skew_K^2);

end