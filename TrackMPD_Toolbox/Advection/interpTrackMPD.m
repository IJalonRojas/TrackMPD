function K=interpTrackMPD(X,Y,ZZ,TT,k,posPart,time,mask_water)

for i=1:size(k,3)
  ki(i) = interp3( X, Y, TT, squeeze(k(:,:,i,:)), posPart(1,1), posPart(1,2),time, 'linear' );
  hi(i) = interp3( X, Y, TT, squeeze(ZZ(:,:,i,:)), posPart(1,1), posPart(1,2),time, 'linear' );
end

% Mask Water interpolated at particle position
MaskWaterInterp = interp2( X, Y, mask_water, posPart(1,1), posPart(1,2), 'linear' );

% Removing Nan values before interpolation
% Nan values can be at bottom in shallow waters with hybrid layers ?
ki(isnan(ki))=0;
%Ind = find(~isnan(ki));
%ki = ki(Ind);
%hi = hi(Ind);

% Settling => posPart(1,3)<hi(end)
% Out of domain => isnan(sum(hi))
% Inland => MaskWateri<0.2 MODIFIED BY MARIEU 2025/01 to take into account the new level reference system 
if isnan(sum(hi)) || posPart(1,3)<hi(end) || MaskWaterInterp<0.2
  K=0;
else
  % if out of water, value at the surface
  if posPart(1,3)>=hi(1)   
    K=ki(1);
  else
    K=interp1(hi,ki,posPart(1,3));
  end
end

end
