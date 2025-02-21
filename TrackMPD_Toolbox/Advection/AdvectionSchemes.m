function [ln_ADV,lt_ADV,h_ADV] = AdvectionSchemes(Scheme,longitude,latitude,Depths,TTs,Us,Vs,Ws,LL1,LL2,LL3,Time,TimeStepCalc,Direction,mask_water)

% Advection Euler 1st order
if strcmpi(Scheme,'Euler')
  UU=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[LL1 LL2 LL3],Time,mask_water);
  VV=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[LL1 LL2 LL3],Time,mask_water);
  WW=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[LL1 LL2 LL3],Time,mask_water);
  ln_ADV=LL1+Direction*UU*TimeStepCalc;
  lt_ADV=LL2+Direction*VV*TimeStepCalc;
  h_ADV=LL3+Direction*WW*TimeStepCalc;
  
  % Advection Runge-Kutta 2nd order
elseif strcmpi(Scheme,'RK2')
  UU1=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[LL1 LL2 LL3],Time,mask_water);
  VV1=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[LL1 LL2 LL3],Time,mask_water);
  WW1=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[LL1 LL2 LL3],Time,mask_water);
  
  ln_ADV=LL1+Direction*UU1*TimeStepCalc;
  lt_ADV=LL2+Direction*VV1*TimeStepCalc;
  h_ADV=LL3+Direction*WW1*TimeStepCalc;
  UU2=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[ln_ADV,lt_ADV,h_ADV],Time+Direction*TimeStepCalc,mask_water);
  VV2=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[ln_ADV,lt_ADV,h_ADV],Time+Direction*TimeStepCalc,mask_water);
  WW2=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[ln_ADV,lt_ADV,h_ADV],Time+Direction*TimeStepCalc,mask_water);
  
  ln_ADV=LL1+Direction*(UU1+UU2)*0.5*TimeStepCalc;
  lt_ADV=LL2+Direction*(VV1+VV2)*0.5*TimeStepCalc;
  h_ADV=LL3+Direction*(WW1+WW2)*0.5*TimeStepCalc;
  
  % Advection Runge-Kutta 4th order
elseif strcmpi(Scheme,'RK4')
  UU1=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[LL1 LL2 LL3],Time,mask_water);
  VV1=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[LL1 LL2 LL3],Time,mask_water);
  WW1=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[LL1 LL2 LL3],Time,mask_water);
  
  ln_ADV=LL1+Direction*UU1*0.5*TimeStepCalc;
  lt_ADV=LL2+Direction*VV1*0.5*TimeStepCalc;
  h_ADV=LL3+Direction*WW1*0.5*TimeStepCalc;
  UU2=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[ln_ADV,lt_ADV,h_ADV],Time+Direction*0.5*TimeStepCalc,mask_water);
  VV2=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[ln_ADV,lt_ADV,h_ADV],Time+Direction*0.5*TimeStepCalc,mask_water);
  WW2=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[ln_ADV,lt_ADV,h_ADV],Time+Direction*0.5*TimeStepCalc,mask_water);
  
  ln_ADV=LL1+Direction*UU2*0.5*TimeStepCalc;
  lt_ADV=LL2+Direction*VV2*0.5*TimeStepCalc;
  h_ADV=LL3+Direction*WW2*0.5*TimeStepCalc;
  UU3=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[ln_ADV,lt_ADV,h_ADV],Time+Direction*0.5*TimeStepCalc,mask_water);
  VV3=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[ln_ADV,lt_ADV,h_ADV],Time+Direction*0.5*TimeStepCalc,mask_water);
  WW3=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[ln_ADV,lt_ADV,h_ADV],Time+Direction*0.5*TimeStepCalc,mask_water);
  
  ln_ADV=LL1+Direction*UU3*TimeStepCalc;
  lt_ADV=LL2+Direction*VV3*TimeStepCalc;
  h_ADV=LL3+Direction*WW3*TimeStepCalc;
  UU4=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[ln_ADV,lt_ADV,h_ADV],Time+Direction*TimeStepCalc,mask_water);
  VV4=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[ln_ADV,lt_ADV,h_ADV],Time+Direction*TimeStepCalc,mask_water);
  WW4=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[ln_ADV,lt_ADV,h_ADV],Time+Direction*TimeStepCalc,mask_water);
  
  ln_ADV=LL1+Direction*(UU1+2*UU2+2*UU3+UU4)*TimeStepCalc/6;
  lt_ADV=LL2+Direction*(VV1+2*VV2+2*VV3+VV4)*TimeStepCalc/6;
  h_ADV=LL3+Direction*(WW1+2*WW2+2*WW3+WW4)*TimeStepCalc/6;

end
