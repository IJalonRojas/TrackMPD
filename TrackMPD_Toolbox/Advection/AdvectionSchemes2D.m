function [ln_ADV,lt_ADV] = AdvectionSchemes2D(Scheme,longitude,latitude,TTs,Us,Vs,LL1,LL2,Time,TimeStepCalc,Direction,mask_water)


% Advection Euler 1st order
if strcmpi(Scheme,'Euler')
  UU=interpTrackMPD_2D(longitude,latitude,TTs,Us,[LL1 LL2],Time,mask_water);
  VV=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[LL1 LL2],Time,mask_water);
  ln_ADV=LL1+Direction*UU*TimeStepCalc;
  lt_ADV=LL2+Direction*VV*TimeStepCalc;
  
  % Advection Runge-Kutta 2nd order
elseif strcmpi(Scheme,'RK2')
  UU1=interpTrackMPD_2D(longitude,latitude,TTs,Us,[LL1 LL2],Time,mask_water);
  VV1=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[LL1 LL2],Time,mask_water);
  
  ln_ADV=LL1+Direction*UU1*TimeStepCalc;
  lt_ADV=LL2+Direction*VV1*TimeStepCalc;
  UU2=interpTrackMPD_2D(longitude,latitude,TTs,Us,[ln_ADV,lt_ADV],Time+Direction*TimeStepCalc,mask_water);
  VV2=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[ln_ADV,lt_ADV],Time+Direction*TimeStepCalc,mask_water);
  
  ln_ADV=LL1+Direction*(UU1+UU2)*0.5*TimeStepCalc;
  lt_ADV=LL2+Direction*(VV1+VV2)*0.5*TimeStepCalc;
  
  % Advection Runge-Kutta 4th order
elseif strcmpi(Scheme,'RK4')
  UU1=interpTrackMPD_2D(longitude,latitude,TTs,Us,[LL1 LL2],Time,mask_water);
  VV1=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[LL1 LL2],Time,mask_water);
  
  ln_ADV=LL1+Direction*UU1*0.5*TimeStepCalc;
  lt_ADV=LL2+Direction*VV1*0.5*TimeStepCalc;
  UU2=interpTrackMPD_2D(longitude,latitude,TTs,Us,[ln_ADV,lt_ADV],Time+Direction*0.5*TimeStepCalc,mask_water);
  VV2=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[ln_ADV,lt_ADV],Time+Direction*0.5*TimeStepCalc,mask_water);
  
  ln_ADV=LL1+Direction*UU2*0.5*TimeStepCalc;
  lt_ADV=LL2+Direction*VV2*0.5*TimeStepCalc;
  UU3=interpTrackMPD_2D(longitude,latitude,TTs,Us,[ln_ADV,lt_ADV],Time+0.5*Direction*TimeStepCalc,mask_water);
  VV3=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[ln_ADV,lt_ADV],Time+0.5*Direction*TimeStepCalc,mask_water);
  
  ln_ADV=LL1+Direction*UU3*TimeStepCalc;
  lt_ADV=LL2+Direction*VV3*TimeStepCalc;
  UU4=interpTrackMPD_2D(longitude,latitude,TTs,Us,[ln_ADV,lt_ADV],Time+Direction*TimeStepCalc,mask_water);
  VV4=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[ln_ADV,lt_ADV],Time+Direction*TimeStepCalc,mask_water);
  
  ln_ADV=LL1+Direction*(UU1+2*UU2+2*UU3+UU4)*TimeStepCalc/6;
  lt_ADV=LL2+Direction*(VV1+2*VV2+2*VV3+VV4)*TimeStepCalc/6;


end
