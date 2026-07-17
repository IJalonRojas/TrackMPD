function conf = read_input_TELEMAC(InputFile)
fn = fopen(InputFile);

Ligne = fgetl(fn);
while ischar(Ligne)

  % Removing lignes that does not have '='
  Eq =  strfind(Ligne,'=');
  if ~isempty(Eq)

    % Taking everything before first comment
    Com = strfind(Ligne,'%');
    if ~isempty(Com)
      Ligne = Ligne(1:max(1,Com(1)-1));
    end

    % Taking everything before first ;
    Com = strfind(Ligne,';');
    if ~isempty(Com)
      Ligne = Ligne(1:max(1,Com(1)-1));
    end

    % Removing lignes that does not have '=' before first comment
    Eq =  strfind(Ligne,'=');
    if ~isempty(Eq)


      % Separation between key and Value
      Sep =  strfind(Ligne,'=');
      Key = Ligne(1:Sep(1)-1);
      Value = Ligne(Sep(1)+1:end);

      % Spaces and tabs removal
      while strcmp(Key(1),' ')==1 || strcmp(Key(1),'	')==1
        Key = Key(2:end);
      end
      while strcmp(Key(end),' ')==1 || strcmp(Key(end),'	')==1
        Key = Key(1:end-1);
      end
      while strcmp(Value(1),' ')==1 || strcmp(Value(1),'	')==1
        Value = Value(2:end);
      end
      while strcmp(Value(end),' ')==1 || strcmp(Value(end),'	')==1
        Value = Value(1:end-1);
      end

      % Key and value recording
      switch Key


        % GENERAL PARAMETERS

        case {'OGCM.TELEMACFile','OGCM.TelemacFile','OGCM.telemacfile','ogcm.telemacfile'}
          conf.OGCM.TELEMACFile = eval(Value);
        case {'OGCM.T0','OGCM.t0','ogcm.t0'}
          conf.OGCM.t0 = eval(Value);

        % 1 GRID CONFIGURATION CLASSIC

        case {'OGCM.NumLonGrid','OGCM.numlongrid','ogcm.numlongrid'}
          conf.OGCM.NumLonGrid = eval(Value);
        case {'OGCM.NumLatGrid','OGCM.numlatgrid','ogcm.numlatgrid'}
          conf.OGCM.NumLatGrid = eval(Value);

        % 2 GRIDS CONFIGURATION VARIABLE RESOLUTION

        case {'OGCM.NumLonGrid1','OGCM.numlongrid1','ogcm.numlongrid1'}
          conf.OGCM.NumLonGrid1 = eval(Value);
        case {'OGCM.NumLonGrid2','OGCM.numlongrid2','ogcm.numlongrid2'}
          conf.OGCM.NumLonGrid2 = eval(Value);
        case {'OGCM.NumLatGrid1','OGCM.numlatgrid1','ogcm.numlatgrid1'}
          conf.OGCM.NumLatGrid1 = eval(Value);
        case {'OGCM.NumLatGrid2','OGCM.numlatgrid2','ogcm.numlatgrid2'}
          conf.OGCM.NumLatGrid2 = eval(Value);
          
        % INTERFACE COORDINATES BETWEEN THE TWO GRIDS
        
        case {'OGCM.intLat','OGCM.intlat','ogcm.intlat'}
          conf.OGCM.intLat = eval(Value);
        case {'OGCM.intLon','OGCM.intlon','ogcm.intlon'}
          conf.OGCM.intLon = eval(Value);

        otherwise
          error(['This key is not available in parameter file: ' Key]);
      end

      % Display
      disp([Key ': ' Value])
    end
  end

  Ligne = fgetl(fn);
end
fclose(fn);



%%%
