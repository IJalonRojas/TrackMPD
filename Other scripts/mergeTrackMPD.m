function mergeTrackMPD(conf_file)
% MERGETRACKMPD  Merge the per-node output files of a parallel run into one file.
%
%   After a multi-node run launched with runTrackMPDParallel, each node writes
%   <ScenarioName>_Node<k>.mat in the output directory conf.Traj.BaseDir. This
%   function concatenates them into a single <ScenarioName>_Merged.mat, ready
%   for mapping and analysis.
%
%   Usage:
%     mergeTrackMPD(conf_file)   % conf_file is a .par file or a config .m
%
%   The number of nodes, the output directory and the scenario name are read
%   from the same configuration file used for the run, so the merge is always
%   consistent with what runTrackMPDParallel produced. Node files that are
%   missing are skipped with a warning.

    % --- Read the configuration file (same logic as runTrackMPDParallel) ---
    [~, name_only, ext] = fileparts(conf_file);
    if strcmpi(ext, '.par')
        conf = read_input_conf(conf_file);
    else
        conf = feval(name_only);
    end

    if ~isfield(conf.Traj, 'NbNodes')
        error('Traj.NbNodes is not defined in the configuration file.');
    end

    NbNodes      = conf.Traj.NbNodes;
    outputDir    = conf.Traj.BaseDir;
    scenarioName = conf.Traj.ScenarioName;
    mergedFile   = fullfile(outputDir, [scenarioName '_Merged.mat']);

    % Per-particle fields, concatenated along dimension 1 (particles)
    particleFields = {'Lon', 'Lat', 'Depth', 'FateType', ...
                      'FinalLonLatDepth', 'TrajectoryDuration'};
    % Global fate counters, summed across nodes
    counterFields  = {'WaterParticles', 'BeachedParticles', ...
                      'BottomParticles', 'OutDomainParticles'};

    fprintf('Starting merge of %d node files from %s\n', NbNodes, outputDir);

    merged_TRAJ = struct();
    first   = true;   % becomes false once the first existing file is loaded
    nMerged = 0;

    for n = 1:NbNodes
        filename = fullfile(outputDir, sprintf('%s_Node%d.mat', scenarioName, n));

        if ~isfile(filename)
            warning('File %s not found. Skipping node %d.', filename, n);
            continue;
        end

        % Load only the TRAJ structure to save memory
        data = load(filename, 'TRAJ');

        if first
            % Initialise on the first file actually found (not necessarily node 1).
            % Shared metadata (TimeStamp, conf, Grid, Domain) is inherited here.
            merged_TRAJ = data.TRAJ;
            first = false;
        else
            % Concatenate particle-specific arrays
            for f = 1:numel(particleFields)
                fn = particleFields{f};
                if isfield(data.TRAJ, fn) && isfield(merged_TRAJ, fn)
                    merged_TRAJ.(fn) = cat(1, merged_TRAJ.(fn), data.TRAJ.(fn));
                end
            end
            % Sum global counters
            for f = 1:numel(counterFields)
                fn = counterFields{f};
                if isfield(data.TRAJ, fn) && isfield(merged_TRAJ, fn)
                    merged_TRAJ.(fn) = merged_TRAJ.(fn) + data.TRAJ.(fn);
                end
            end
        end

        nMerged = nMerged + 1;
    end

    if first
        error(['No node file found for scenario "%s" in %s. ' ...
               'Nothing to merge.'], scenarioName, outputDir);
    end

    % --- Save the unified dataset ---
    % TRAJ is the variable name expected by the mapping/analysis tools.
    TRAJ = merged_TRAJ;
    fprintf('Merged %d of %d node files ; %d particles total.\n', ...
            nMerged, NbNodes, size(TRAJ.Lon, 1));
    fprintf('Saving merged dataset to %s\n', mergedFile);

    % '-v7.3' is required to save files larger than 2 GB
    save(mergedFile, 'TRAJ', '-v7.3');

    disp('Merge completed successfully. Data is ready for mapping.');
end