function K=interpTrackMPD_2D(X,Y,TT,k,posPart,time,mask_water)

K = interp3( X, Y, TT, squeeze(k(:,:,1,:)), posPart(1,1), posPart(1,2),time, 'linear' );

% Mask Water interpolated at particle position. Could be NaN if out of domain
try 
  MaskWaterInterp = interp2( X, Y, mask_water, posPart(1,1), posPart(1,2), 'linear' );
catch
  MaskWaterInterp = 0;
end  

% Removing Nan values before interpolation
K(isnan(K))=0;

% Values Inland
if MaskWaterInterp<0.2
  K=0;
end

end
