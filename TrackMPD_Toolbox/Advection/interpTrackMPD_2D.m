function K=interpTrackMPD_2D(X,Y,TT,k,posPart,time,mask_water)

K = interp3( X, Y, TT, squeeze(k(:,:,1,:)), posPart(1,1), posPart(1,2),time, 'linear' );

% Mask Water interpolated at particle position
MaskWaterInterp = interp2( X, Y, mask_water, posPart(1,1), posPart(1,2), 'linear' );

% Removing Nan values before interpolation
K(isnan(K))=0;

% Values Inland
if MaskWaterInterp<0.2
  K=0;
end

end
