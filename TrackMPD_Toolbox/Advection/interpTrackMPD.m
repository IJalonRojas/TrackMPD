function K=interpTrackMPD(X,Y,ZZ,TT,k,posPart,time)

for i=1:size(k,3)
  ki(i) = interp3( X, Y, TT, squeeze(k(:,:,i,:)), posPart(1,1), posPart(1,2),time, 'linear' );
  hi(i) = interp3( X, Y, TT, squeeze(ZZ(:,:,i,:)), posPart(1,1), posPart(1,2),time, 'linear' );
end

% Removing Nan values before interpolation
% Nan values can be at bottom in shallow waters with hybrid layers ?
ki(isnan(ki))=0;
%Ind = find(~isnan(ki));
%ki = ki(Ind);
%hi = hi(Ind);

% Settling => posPart(1,3)<hi(end)
% Out of domain => isnan(sum(hi))
% Inland => sum(hi)>0 (Z(Land) set to 1 in transform*inputs.m) 
if isnan(sum(hi)) || posPart(1,3)<hi(end) || sum(hi)>0
  K=0;
else
  % if out of water, value at the surface
  if posPart(1,3)>=0   
    K=ki(1);
  else
    K=interp1(hi,ki,posPart(1,3));
  end
end

end
