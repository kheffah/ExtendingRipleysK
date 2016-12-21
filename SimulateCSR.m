function CSR = SimulateCSR (imLength,imIntens)

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
% 
%   SAMPLE INPUT:
%
%   imLength = 2^8 ; %(length of the square field of view)
%   imIntens = 10 ;  %(intensity of particles)
%   imCSR = SimulateCSR (imLength,imIntens) ;
%
%   SAMPLE OUTPUT:
%
%   whos imCSR
%   Name       Size        Bytes     Class     Attributes
%   imCSR      256x256     524288    double     
%
%==========================================================================



%%
TotParticles = ceil((imLength^2)*imIntens) ; 

%%


Positions = randi([1 imLength^2],floor(sqrt(TotParticles)));
Positions(1,1)=1;
Positions(1,2)=imLength^2;

M = tabulate(Positions(:));
M=M(:,2);
M=reshape(M,imLength,imLength);

CSR = M;

end