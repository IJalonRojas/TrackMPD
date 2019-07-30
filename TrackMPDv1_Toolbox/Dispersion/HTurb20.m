function [varargout]=HTurb20(n,dt,varargin)

% ============================================================================
% $RCSfile$
% $Source$
% $Revision$
% $Date$
% $Author$
% $Name$
%
% USAGE: [TurbHx,TurbHy]=HTurb20(1000,2)
% DESCRIPTION: 
% The Random Displacement Model (RDM) zero order assumes that the
% turbulence at each point is isotropic in the horizontal, then turbulence
% is characterized by the horizontal diffusivity Kh=K11=K22, and the 
% vertical diffusivity K33 (Rodean, 1996). 
%
% Reference:
% Rodean, H. C. 1996. Stochastic Lagrangian Models of Turbulent Diffusion. American Meteorological
% Society, Boston, MA, USA. pp 84.
%
% ============================================================================
%
% Author: Erick Fredj
%
% Copyright @ 2018
% ============================================================================


% [TurbHx,TurbHy]=HTurb0(1000,2,1)
% Apply Random Walk Model
% default options
if ~exist('n','var')
    n = 0;
end

if ~exist('dt','var')
    dt = 0;
end
Kh = 1.0;  % horizontal diffusivity m^2/s

% look at any remaining options
k=1;
while k<length(varargin)
    optn=lower(varargin{k});
    switch optn
        case 'kh'
            Kh = varargin{k+1};
    end
    k=k+2;
end

devX=2*rand(n,1)-1;     % the random deviate in the X direction
devY=2*rand(n,1)-1;     % the random deviate in the Y direction

% Apply random walk model to calculate horizontal turbulent particle
% displacement
TurbHx= devX*sqrt(2*Kh*dt);
TurbHy= devY*sqrt(2*Kh*dt);

varargout{1} = TurbHx;
varargout{2} = TurbHy;

end