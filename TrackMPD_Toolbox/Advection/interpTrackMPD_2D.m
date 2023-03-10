function K=interpTrackMPD_2D(X,Y,TT,k,posPart,time)

ki = interp3( X, Y, TT, squeeze(k(:,:,1,:)), posPart(1,1), posPart(1,2),time, 'linear' );

% Removing Nan values before interpolation
% Nan values can be at bottom in shallow waters with hybrid layers ?
ki(isnan(ki))=0;
%Ind = find(~isnan(ki));
%ki = ki(Ind);
%hi = hi(Ind);


% Comment this bloc for testing 2D %IJR 2D
% % Settling => posPart(1,3)<hi(end)
% % Out of domain => isnan(sum(hi))
% % Inland => sum(hi)>0 (Z(Land) set to 1 in transform*inputs.m) 
% if isnan(sum(hi)) || posPart(1,3)<hi(end) || sum(hi)>=0
%   K=0;
% else
%   % if out of water, value at the surface
%   if posPart(1,3)>=0   
%     K=ki(1);
%   else
%     K=interp1(hi,ki,posPart(1,3));
%   end
% end

K=ki; %IJR 2D

end