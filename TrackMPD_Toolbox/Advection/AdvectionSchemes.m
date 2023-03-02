function [ln_ADV,lt_ADV,h_ADV] = AdvectionSchemes(Scheme,longitude,latitude,Depths,TTs,Us,Vs,Ws,LL1,LL2,LL3,Time,TimeStepCalc)

% Advection Euler 1st order
if strcmpi(Scheme,'Euler')
  UU=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[LL1 LL2 LL3],Time);
  VV=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[LL1 LL2 LL3],Time);
  WW=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[LL1 LL2 LL3],Time);
  ln_ADV=LL1+UU*TimeStepCalc;
  lt_ADV=LL2+VV*TimeStepCalc;
  h_ADV=LL3+WW*TimeStepCalc;
  
  % Advection Runge-Kutta 2nd order
elseif strcmpi(Scheme,'RK2')
  UU1=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[LL1 LL2 LL3],Time);
  VV1=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[LL1 LL2 LL3],Time);
  WW1=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[LL1 LL2 LL3],Time);
  
  ln_ADV=LL1+UU1*TimeStepCalc;
  lt_ADV=LL2+VV1*TimeStepCalc;
  h_ADV=LL3+WW1*TimeStepCalc;
  UU2=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[ln_ADV,lt_ADV,h_ADV],Time+TimeStepCalc);
  VV2=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[ln_ADV,lt_ADV,h_ADV],Time+TimeStepCalc);
  WW2=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[ln_ADV,lt_ADV,h_ADV],Time+TimeStepCalc);
  
  ln_ADV=LL1+(UU1+UU2)*0.5*TimeStepCalc;
  lt_ADV=LL2+(VV1+VV2)*0.5*TimeStepCalc;
  h_ADV=LL3+(WW1+WW2)*0.5*TimeStepCalc;
  
  % Advection Runge-Kutta 4th order
elseif strcmpi(Scheme,'RK4')
  UU1=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[LL1 LL2 LL3],Time);
  VV1=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[LL1 LL2 LL3],Time);
  WW1=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[LL1 LL2 LL3],Time);
  
  ln_ADV=LL1+UU1*0.5*TimeStepCalc;
  lt_ADV=LL2+VV1*0.5*TimeStepCalc;
  h_ADV=LL3+WW1*0.5*TimeStepCalc;
  UU2=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[ln_ADV,lt_ADV,h_ADV],Time+0.5*TimeStepCalc);
  VV2=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[ln_ADV,lt_ADV,h_ADV],Time+0.5*TimeStepCalc);
  WW2=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[ln_ADV,lt_ADV,h_ADV],Time+0.5*TimeStepCalc);
  
  ln_ADV=LL1+UU2*0.5*TimeStepCalc;
  lt_ADV=LL2+VV2*0.5*TimeStepCalc;
  h_ADV=LL3+WW2*0.5*TimeStepCalc;
  UU3=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[ln_ADV,lt_ADV,h_ADV],Time+0.5*TimeStepCalc);
  VV3=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[ln_ADV,lt_ADV,h_ADV],Time+0.5*TimeStepCalc);
  WW3=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[ln_ADV,lt_ADV,h_ADV],Time+0.5*TimeStepCalc);
  
  ln_ADV=LL1+UU3*TimeStepCalc;
  lt_ADV=LL2+VV3*TimeStepCalc;
  h_ADV=LL3+WW3*TimeStepCalc;
  UU4=interpTrackMPD(longitude,latitude,Depths,TTs,Us,[ln_ADV,lt_ADV,h_ADV],Time+TimeStepCalc);
  VV4=interpTrackMPD(longitude,latitude,Depths,TTs,Vs,[ln_ADV,lt_ADV,h_ADV],Time+TimeStepCalc);
  WW4=interpTrackMPD(longitude,latitude,Depths,TTs,Ws,[ln_ADV,lt_ADV,h_ADV],Time+TimeStepCalc);
  
  ln_ADV=LL1+(UU1+2*UU2+2*UU3+UU4)*TimeStepCalc/6;
  lt_ADV=LL2+(VV1+2*VV2+2*VV3+VV4)*TimeStepCalc/6;
  h_ADV=LL3+(WW1+2*WW2+2*WW3+WW4)*TimeStepCalc/6;
end