function runTrackMPDParallel(conf_file, telemac_file, NodeID)
    % RUNTRACKMPDPARALLEL Distributes computations across multiple nodes.
    % NbNodes is read from the configuration file (Traj.NbNodes).
    % If SLURM_ARRAY_TASK_MAX is set, it is compared to Traj.NbNodes
    % and a WARNING is raised if they differ.
    %
    % Usage:
    %   runTrackMPDParallel(conf_file, telemac_file, NodeID)  

    % Conversion for compiled code compatibility (string to number)
    if nargin >= 3
        if ischar(NodeID) || isstring(NodeID)
            NodeID = str2double(NodeID);
        end
    end

    % Read the main configuration file
    [~, name_only, ext] = fileparts(conf_file);
    if strcmpi(ext, '.par')
        conf = read_input_conf(conf_file);
    else
        conf = feval(name_only);
    end

    % Read NbNodes from configuration file
    if isfield(conf.Traj, 'NbNodes')
        NbNodes = conf.Traj.NbNodes;
    else
        error('Traj.NbNodes is not defined in the configuration file.');
    end

    % SLURM consistency check 
    slurm_nb_nodes_str = getenv('SLURM_ARRAY_TASK_MAX');
    if ~isempty(slurm_nb_nodes_str)
        NbNodes_SLURM = str2double(slurm_nb_nodes_str);
        if NbNodes_SLURM ~= NbNodes
            fprintf('\n');
            fprintf('WARNING - MISMATCH BETWEEN SLURM AND PARAMETER FILE \n');
            fprintf('-> Traj.NbNodes in .par file  : %d\n', NbNodes);
            fprintf('-> >>> Using SLURM value (%d nodes) <<<\n', NbNodes_SLURM);
            fprintf('-> Please synchronize Traj.NbNodes in your .par file.\n');
            fprintf('\n');
            NbNodes = NbNodes_SLURM;  % SLURM is authoritative
        else
            fprintf('[OK] Traj.NbNodes (%d) matches SLURM --array (%d).\n', NbNodes, NbNodes_SLURM);
        end
    end


    % Count total number of particles
    try
        M_total = dlmread(conf.Data.ParticlesFile, ',', 0, 0);
        NbTotalParticles = size(M_total, 1);
    catch
        error('Unable to read the particles file.');
    end

    % Workload distribution (Load Balancing)
    BaseSize = floor(NbTotalParticles / NbNodes);
    Remainder = mod(NbTotalParticles, NbNodes);

    fprintf('Running Node %d out of %d\n', NodeID, NbNodes);
    node_list = NodeID;

    % Main execution loop
    for n = node_list
        % Calculate exact start index for the current node
        if n <= Remainder
            NodeSize = BaseSize + 1;
            idx_start = (n - 1) * NodeSize + 1;
        else
            NodeSize = BaseSize;
            idx_start = Remainder * (BaseSize + 1) + (n - 1 - Remainder) * BaseSize + 1;
        end

        if NodeSize == 0
            fprintf('Node %d has no particles to process. Skipping.\n', n);
            continue;
        end

        idx_end = idx_start + NodeSize - 1;
        indices_to_compute = idx_start:idx_end;

        fprintf('\n--- Starting Node %d (%d particles, indices %d to %d) ---\n', ...
                n, length(indices_to_compute), idx_start, idx_end);

        % Actual launch of the trajectory computation
        runTrackMPD(conf_file, telemac_file, indices_to_compute, n);
    end

    disp('++++++ Parallel distribution completed successfully. ++++++');
end