function conf = read_input_FVCOM(InputFile)
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

        % OGCM
        case {'OGCM.FVCOMFile','OGCM.FvcomFile','OGCM.fvcomfile','ogcm.fvcomfile'}
          conf.OGCM.FVCOMFile = eval(Value);
        case {'OGCM.FVCOMGrid','OGCM.FvcomGrid','OGCM.fvcomgrid','ogcm.fvcomgrid'}
          conf.OGCM.FVCOMGrid = eval(Value);
        case {'OGCM.NumLonGrid','OGCM.numlongrid','ogcm.numlongrid'}
          conf.OGCM.NumLonGrid = eval(Value);
        case {'OGCM.NumLatGrid','OGCM.numlatgrid','ogcm.numlatgrid'}
          conf.OGCM.NumLatGrid = eval(Value);
        case {'OGCM.NameTime','OGCM.nametime','ogcm.nametime'}
          conf.OGCM.NameTime = eval(Value);
        case {'OGCM.NameW','OGCM.namew','ogcm.namew'}
          conf.OGCM.NameW = eval(Value);

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
