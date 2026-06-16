function [TurbHx,TurbHy,TurbV] = HTurb_drw3D(dt,Kh_shiftx,Kh_shifty,Kv_shiftz,dKh_dX, dKh_dY, dKv_dZ)
% ============================================================================
% $RCSfile$
% $Source$
% $Revision$
% $Date$
% $Author$
% $Name$
%
% USAGE: [TurbHx,TurbHy,TurbV]=HTurb_drw3D(dt,Kh_shiftx,Kh_shifty,Kv_shiftz,dKh_dX,dKh_dY,dKv_dZ)
% DESCRIPTION: 
%
% References:Visser, Andy. (1997). Using random walk models to simulate the vertical distribution of particles in a turbulent water column.
% Marine Ecology-progress Series - MAR ECOL-PROGR SER. 158. 275-281. 10.3354/meps158275. 
% 
% Pilechi, A., Mohammadian, A., and Murphy, E.: A numerical framework for modeling fate and transport of microplastics in inland and
% coastal waters, Marine P ollution Bulletin, 184, 114 119, 2022.
% 
% ============================================================================

% Apply Random Walk Model
% default options

devX=2*rand(1,1)-1;     % the random deviate in the X direction
devY=2*rand(1,1)-1;     % the random deviate in the Y direction
devZ=2*rand(1,1)-1;     % the random deviate in the Z direction

% Apply the diffusive random walk model to calculate horizontal and vertical turbulent 
% particle  displacement


% Compute particle displacement (extended diffusion equation)
TurbHx = (dKh_dX)*dt + ... % dKp/dX term in the original equation
           sqrt(3)*devX .* ... % Random walk component
           (sqrt(2 * Kh_shiftx * dt));

TurbHy= (dKh_dY)*dt + ... % dKp/dY term in the original equation
           sqrt(3)*devY .* ... % Random walk component
           (sqrt(2 * Kh_shifty * dt));

TurbV= (dKv_dZ)*dt + ... % dKp/dZ term in the original equation
           sqrt(3)*devZ .* ... % Random walk component
           (sqrt(2 * Kv_shiftz * dt));

end