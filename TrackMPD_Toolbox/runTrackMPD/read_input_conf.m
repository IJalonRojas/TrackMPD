function conf=read_input_conf(InputFile)

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

        % Data
        case {'Data.BaseDir','data.BaseDir','Data.basedir','data.dasedir'}
          conf.Data.BaseDir = eval(Value) ;
        case {'Data.Domain','Data.domain','data.Domain','data.domain'}
          conf.Data.Domain = eval(Value);
        case {'Data.ParticlesFile','Data.particlesfile','data.ParticlesFile','data.particlesfile'}
          conf.Data.ParticlesFile = eval(Value);

          % OGCM
        case {'OGCM.Model_Name','OGCM.Model_name','OGCM.model_name','OGCM.ModelName','ogcm.model_name'}
          conf.OGCM.Model_name = eval(Value);
        case {'OGCM.DomainName','OGCM.domainname','OGCM.Domain_Name','OGCM.domain_name','ogcm.domain_name'}
          conf.OGCM.DomainName = eval(Value);
        case {'OGCM.BaseDir','OGCM.basedir'}
          conf.OGCM.BaseDir = eval(Value);
        case {'OGCM.TimeStep','OGCM.timestep','ogcm.timestep'}
          conf.OGCM.TimeStep = eval(Value);
        case {'OGCM.Coordinates','OGCM.coordinates','ogcm.coordinates'}
          conf.OGCM.Coordinates = eval(Value);
        case {'OGCM.VerticalLayer','OGCM.verticalayer','ogcm.verticallayer'}
          conf.OGCM.VerticalLayer = eval(Value);
        case {'OGCM.BottomType','OGCM.bottomtype','ogcm.bottomtype'}
          conf.OGCM.BottomType = eval(Value);
        case {'OGCM.Cut','OGCM.cut','ogcm.cut'}
          conf.OGCM.cut = eval(Value);
        case {'OGCM.minLon','OCGM.minlon','ogcm.minlon'}
          conf.OGCM.minLon = eval(Value);
        case {'OGCM.maxLon','OCGM.maxlon','ogcm.maxlon'}
          conf.OGCM.maxLon = eval(Value);
        case {'OGCM.minLat','OCGM.minlat','ogcm.minlat'}
          conf.OGCM.minLat = eval(Value);
        case {'OGCM.maxLat','OCGM.maxlat','ogcm.maxlat'}
          conf.OGCM.maxLat = eval(Value);
	  

          % Trajectory
        case {'Traj.Verbose','Traj.verbose','traj.verbose'}
          conf.Traj.Verbose = eval(Value);
        case {'Traj.Mode','Traj.Mode','traj.Mode'}
          conf.Traj.Mode = eval(Value);
        case {'Traj.NbCores','Traj.nbcores','traj.cbcores'}
          conf.Traj.NbCores = eval(Value);
        case {'Traj.NbNodes','Traj.nbnodes','traj.nbnodes'}
          conf.Traj.NbNodes = eval(Value);
        case {'Traj.NbCoresPerNode','Traj.nbcorespernode','traj.nbcorespernode'}
          conf.Traj.NbCoresPerNode = eval(Value);
        case {'Traj.ReleaseTime','Traj.releasetime','Traj.releasetime'}
          conf.Traj.ReleaseTime = eval(Value);
        case {'Traj.TrajectoryDuration','Traj.trajectoryduration','traj.trajectoryduration'}
          conf.Traj.TrajectoryDuration = eval(Value) ;
        case {'Traj.TimeStepCalc','Traj.timestepcalc','traj.timestepcalc'}
          conf.Traj.TimeStepCalc = eval(Value) ;
        case {'Traj.Direction','Traj.direction','traj.direction'}
          conf.Traj.Direction = eval(Value);
        case {'Traj.Scheme','Traj.scheme','traj.scheme'}
          conf.Traj.Scheme = eval(Value);
        case {'Traj.TimeStepOut','Traj.timestepout','traj.timestepout'}
          conf.Traj.TimeStepOut = eval(Value);
        case {'Traj.ScenarioName','Traj.scenarioname','traj.scenarioname'}
          conf.Traj.ScenarioName = eval(Value);
        case {'Traj.BaseDir','Traj.basedir','traj.basedir'}
          conf.Traj.BaseDir = eval(Value);
        case {'Traj.Chunklen','Traj.chunklen','traj.chunklen'}
          conf.Traj.chunklen = eval(Value);
        case  {'Traj.KvOption','Traj.kvoption','traj.kvoption'}
          conf.Traj.KvOption = eval(Value);
        case  {'Traj.KhOption','Traj.khoption','traj.khoption'}
          conf.Traj.KhOption = eval(Value);
        case  {'Traj.Kh','Traj.kh','traj.kh'}
          conf.Traj.Kh = eval(Value);
        case  {'Traj.Kv','Traj.kv','traj.kv'}
          conf.Traj.Kv = eval(Value);
        case  {'Traj.Tau0Formula','Traj.tau0formula','traj.tau0formula'}
          conf.Traj.Tau0Formula = eval(Value);
        case  {'Traj.z0','Traj.Z0','traj.z0'}
          conf.Traj.z0 = eval(Value);
        case  {'Traj.alpha','Traj.Alpha','traj.alpha'}
          conf.Traj.alpha = eval(Value);
        case  {'Traj.Beaching','Traj.beaching','traj.beaching'}
          conf.Traj.Beaching = eval(Value);
        case  {'Traj.Refloating','Traj.refloating','traj.refloating'}
          conf.Traj.Refloating = eval(Value);
        case  {'Traj.Tw','Traj.tw','traj.tw'}
          conf.Traj.Tw = eval(Value);
        case  {'Traj.RefloatingAtHighTide','Traj.refloatingathightide','traj.refloatingathightide'}
          conf.Traj.RefloatingAtHighTide = eval(Value);
        case  {'Traj.Xcoord_Elev','Traj.xcoord_elev','traj.xcoord_elev'}
          conf.Traj.xcoord_elev = eval(Value);
        case  {'Traj.Ycoord_Elev','Traj.ycoord_elev','traj.ycoord_elev'}
          conf.Traj.ycoord_elev = eval(Value);
        case  {'Traj.Deposition','Traj.deposition','traj.deposition'}
          conf.Traj.Deposition = eval(Value);
        case  {'Traj.Resuspension','Traj.resuspension','traj.resuspension'}
          conf.Traj.Resuspension = eval(Value);
        case  {'Traj.Sliding','Traj.sliding','traj.sliding'}
          conf.Traj.Sliding = eval(Value);
        case  {'Traj.HresMax','Traj.hresmax','traj.hresmax'}
          conf.Traj.HresMax = eval(Value);
        case  {'Traj.ResOption','Traj.resoption','traj.resoption'}
          conf.Traj.ResOption = eval(Value);        
        case  {'Traj.Tauc1','Traj.tauc1','traj.tauc1','Traj.Taucr1','Traj.taucr1','traj.taucr1','Traj.Tau1','Traj.tau1','traj.tau1'}
          conf.Traj.taucr1 = eval(Value);
        case  {'Traj.Tauc2','Traj.tauc2','traj.tauc2','Traj.Taucr2','Traj.taucr2','traj.taucr2','Traj.Tau2','Traj.tau2','traj.tau2'}
          conf.Traj.taucr2 = eval(Value);
          
          % Behavior
        case  {'Beh.WaterDensity','Beh.waterdensity','Beh.waterdensity'}
          conf.Beh.WaterDensity = eval(Value);
        case  {'Beh.WaterViscosity','Beh.waterviscosity','Beh.waterviscosity'}
          conf.Beh.WaterViscosity = eval(Value);
        case  {'Beh.SediD50','Beh.sedid50','Beh.sedid50'}
          conf.Beh.SediD50 = eval(Value);
        case  {'Beh.Cm','Beh.cm','Beh.cm'}
          conf.Beh.Cm = eval(Value);
        case  {'Beh.Cd','Beh.cd','Beh.cd'}
          conf.Beh.Cd = eval(Value);
        case  {'Beh.PolymerAcronym','Beh.polymeracronym','Beh.polymeracronym'}
          conf.Beh.PolymerAcronym = eval(Value);
        case  {'Beh.ParticleDensity','Beh.particledensity','Beh.particledensity'}
          conf.Beh.ParticleDensity = eval(Value);
        case  {'Beh.ParticleSize','Beh.particlesize','Beh.particlesize'}
          conf.Beh.ParticleSize = eval(Value);
        case  {'Beh.ParticleDequi','Beh.particledequi','Beh.particledequi'}
          conf.Beh.ParticleDequi = eval(Value);
        case  {'Beh.ParticleShape','Beh.particleshape','Beh.particleshape'}
          conf.Beh.ParticleShape = eval(Value);
        case  {'Beh.InitialWsOption','Beh.Initialwsoption','Beh.Initialwsoption'}
          conf.Beh.InitialWsOption = eval(Value);
        case  {'Beh.Ws0','Beh.ws0','Beh.ws0'}
          conf.Beh.Ws0 = eval(Value);
        case  {'Beh.WsForm','Beh.wsform','Beh.wsform'}
          conf.Beh.WsForm = eval(Value);
        case  {'Beh.CSF','Beh.csf','Beh.csf'}
          conf.Beh.CSF = eval(Value);
        case  {'Beh.P','Beh.p','Beh.p'}
          conf.Beh.P = eval(Value);
        case  {'Beh.Diameter','Beh.diameter','Beh.diameter'}
          conf.Beh.Diameter = eval(Value);
        case  {'Beh.Phi','Beh.phi','Beh.phi'}
          conf.Beh.Phi = eval(Value);
        case  {'Beh.WsRate','Beh.wsrate','Beh.wsrate'}
          conf.Beh.WsRate = eval(Value);
        case  {'Beh.Biofouling','Beh.biofouling','Beh.biofouling'}
          conf.Beh.Biofouling = eval(Value);
        case  {'Beh.BiofoulingThickness','Beh.biofoulingthickness','Beh.biofoulingthickness'}
          conf.Beh.BiofoulingThickness = eval(Value);
        case  {'Beh.FilmDensity','Beh.filmdensity','Beh.filmdensity'}
          conf.Beh.FilmDensity = eval(Value);
        case  {'Beh.BiofoulingRate','Beh.biofoulingrate','Beh.biofoulingrate'}
          conf.Beh.BiofoulingRate = eval(Value);
        case  {'Beh.Degradation','Beh.degradation','Beh.degradation'}
          conf.Beh.Degradation = eval(Value);
        case  {'Beh.DegradationRate','Beh.degradationrate','Beh.degradationrate'}
          conf.Beh.DegradationRate = eval(Value);

        otherwise
          error(['This key is not available in parameter file: ' Key]);
      end

      % Display
      %disp([Key ': ' Value])
    end
  end

  Ligne = fgetl(fn);
end
fclose(fn);

% Full files (Directory + File)
conf.Data.Domain = [conf.Data.BaseDir '/' conf.Data.Domain];
conf.Data.ParticlesFile = [conf.Data.BaseDir '/' conf.Data.ParticlesFile];