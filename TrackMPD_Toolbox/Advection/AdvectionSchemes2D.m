function [ln_ADV,lt_ADV] = AdvectionSchemes2D(Scheme,longitude,latitude,TTs,Us,Vs,LL1,LL2,Time,TimeStepCalc)

% Advection Euler 1st order
if strcmpi(Scheme,'Euler')
  UU=interpTrackMPD_2D(longitude,latitude,TTs,Us,[LL1 LL2],Time);
  VV=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[LL1 LL2],Time);
  ln_ADV=LL1+UU*TimeStepCalc;
  lt_ADV=LL2+VV*TimeStepCalc;
  
  % Advection Runge-Kutta 2nd order
elseif strcmpi(Scheme,'RK2')
  UU1=interpTrackMPD_2D(longitude,latitude,TTs,Us,[LL1 LL2],Time);
  VV1=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[LL1 LL2],Time);
  
  ln_ADV=LL1+UU1*TimeStepCalc;
  lt_ADV=LL2+VV1*TimeStepCalc;
  UU2=interpTrackMPD_2D(longitude,latitude,TTs,Us,[ln_ADV,lt_ADV],Time+TimeStepCalc);
  VV2=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[ln_ADV,lt_ADV],Time+TimeStepCalc);
  
  ln_ADV=LL1+(UU1+UU2)*0.5*TimeStepCalc;
  lt_ADV=LL2+(VV1+VV2)*0.5*TimeStepCalc;
  
  % Advection Runge-Kutta 4th order
elseif strcmpi(Scheme,'RK4')
  UU1=interpTrackMPD_2D(longitude,latitude,TTs,Us,[LL1 LL2],Time);
  VV1=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[LL1 LL2],Time);
  
  ln_ADV=LL1+UU1*0.5*TimeStepCalc;
  lt_ADV=LL2+VV1*0.5*TimeStepCalc;
  UU2=interpTrackMPD_2D(longitude,latitude,TTs,Us,[ln_ADV,lt_ADV],Time+0.5*TimeStepCalc);
  VV2=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[ln_ADV,lt_ADV],Time+0.5*TimeStepCalc);
  
  ln_ADV=LL1+UU2*0.5*TimeStepCalc;
  lt_ADV=LL2+VV2*0.5*TimeStepCalc;
  UU3=interpTrackMPD_2D(longitude,latitude,TTs,Us,[ln_ADV,lt_ADV],Time+0.5*TimeStepCalc);
  VV3=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[ln_ADV,lt_ADV],Time+0.5*TimeStepCalc);
  
  ln_ADV=LL1+UU3*TimeStepCalc;
  lt_ADV=LL2+VV3*TimeStepCalc;
  UU4=interpTrackMPD_2D(longitude,latitude,TTs,Us,[ln_ADV,lt_ADV],Time+TimeStepCalc);
  VV4=interpTrackMPD_2D(longitude,latitude,TTs,Vs,[ln_ADV,lt_ADV],Time+TimeStepCalc);
  
  ln_ADV=LL1+(UU1+2*UU2+2*UU3+UU4)*TimeStepCalc/6;
  lt_ADV=LL2+(VV1+2*VV2+2*VV3+VV4)*TimeStepCalc/6;
end