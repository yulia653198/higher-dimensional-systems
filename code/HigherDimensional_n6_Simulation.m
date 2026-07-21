function HigherDimensional_n6_Simulation
%% HIGHERDIMENSIONAL_N6_SIMULATION
% Additional sixth-order validation for the residual-gated,
% event-triggered, adaptive synchronization protocol.
%
% State of agent i:
%   x_i = [p_x, p_y, p_z, v_x, v_y, v_z]^T in R^6.



clc;
close all;

cfg = make_configuration();

if exist('care','file') ~= 2
    error('The function care is required (Control System Toolbox).');
end

if rank(ctrb(cfg.A,cfg.B)) < cfg.n
    error('The selected sixth-order pair (A,B) is not controllable.');
end

[cfg.P,~,~] = care(cfg.A,cfg.B,cfg.Q,cfg.R);
cfg.K = cfg.R\(cfg.B'*cfg.P);
[cfg.Ad,cfg.Bd] = zoh_discretization(cfg.A,cfg.B,cfg.dt);

fprintf('Sixth-order system: n=%d, m=%d, rank(ctrb(A,B))=%d\n', ...
    cfg.n,cfg.m,rank(ctrb(cfg.A,cfg.B)));
fprintf('Maximum real part of eig(A): %.4f\n',max(real(eig(cfg.A))));

%% Residual-scale calibration for the sixth-order model
fprintf('[1/5] Estimating attack-free residual scales for n=6...\n');
[sigmaFF,sigmaLF,calibration] = estimate_residual_scales(cfg);

fprintf('  sigmaFF=%.6e, sigmaLF=%.6e\n',sigmaFF,sigmaLF);

%% Dimension-specific candidate evaluation using Algorithm 1
fprintf('[2/5] Selecting the n=6 threshold pair by closed-loop cost...\n');
[cfg,thresholdTable] = select_threshold_pair(cfg,sigmaFF,sigmaLF,calibration);

fprintf('  selected gammaFF=%.6e (gammaFFn=%.2f)\n', ...
    cfg.gammaFF,cfg.gammaFFn);
fprintf('  selected gammaLF=%.6e (gammaLFn=%.2f)\n', ...
    cfg.gammaLF,cfg.gammaLFn);

%% Monte Carlo evaluation
fprintf('[3/5] Running %d Monte Carlo simulations...\n',cfg.numMonteCarlo);

metrics = repmat(empty_metrics(),cfg.numMonteCarlo,1);

for run = 1:cfg.numMonteCarlo
    seed = cfg.baseSeed + run;

    tic;
    [~,runMetrics] = simulate_once(cfg,seed,true,false,false);
    runMetrics.cpuTime = toc;
    metrics(run) = runMetrics;

    if mod(run,10)==0 || run==cfg.numMonteCarlo
        fprintf('  completed %d/%d runs\n',run,cfg.numMonteCarlo);
    end
end

% Select the run closest to the median steady-state MSE. This avoids
% presenting either a particularly favorable or unfavorable realization.
allSteadyMSE = [metrics.totalMSE];
medianSteadyMSE = median(allSteadyMSE);
[~,representativeRun] = min(abs(allSteadyMSE-medianSteadyMSE));
cfg.representativeRun = representativeRun;
representativeSeed = cfg.baseSeed+representativeRun;
[representative,~] = simulate_once(cfg,representativeSeed,true,true,false);
fprintf('  representative run: %d (closest to median totalMSE)\n', ...
    representativeRun);

%% Export the numerical results
fprintf('[4/5] Exporting summary tables...\n');

if ~exist(cfg.outputDirectory,'dir')
    mkdir(cfg.outputDirectory);
end

runTable = metrics_to_table(metrics);
writetable(runTable,fullfile(cfg.outputDirectory,'n6_metrics_all_runs.csv'));
writetable(thresholdTable,fullfile(cfg.outputDirectory, ...
    'n6_threshold_candidates.csv'));

summaryTable = summarize_metrics(runTable);
writetable(summaryTable,fullfile(cfg.outputDirectory,'n6_metrics_summary.csv'));

systemTable = table(cfg.n,cfg.m,cfg.N,cfg.dt,cfg.T,cfg.gammaFFn, ...
    cfg.gammaLFn,cfg.gammaFF,cfg.gammaLF,cfg.deltaFF,cfg.deltaLF, ...
    cfg.varrhoFF,cfg.varrhoLF, ...
    'VariableNames',{'n','m','N','dt','T','gammaFFn','gammaLFn', ...
    'gammaFF','gammaLF','deltaFF','deltaLF','varrhoFF','varrhoLF'});
writetable(systemTable,fullfile(cfg.outputDirectory,'n6_system_settings.csv'));

save(fullfile(cfg.outputDirectory,'n6_complete_results.mat'), ...
    'cfg','runTable','summaryTable','representative');

disp('--- Sixth-order Monte Carlo summary ---');
disp(summaryTable);

%% Generate figures
fprintf('[5/5] Generating figures...\n');
make_adaptive_gain_figure(cfg,representative);
make_tracking_error_figure(cfg,representative);
make_control_input_figure(cfg,representative);
make_trigger_error_figure(cfg,representative);
make_trajectory_figure(cfg,representative);
make_mse_figure(cfg,representative);
make_residual_gate_figure(cfg,representative);

fprintf('\nResults saved in: %s\n',cfg.outputDirectory);
fprintf('Recommended files for the supplementary material:\n');
fprintf('  n6_adaptive_gains.pdf\n');
fprintf('  n6_tracking_errors.pdf\n');
fprintf('  n6_control_inputs.pdf\n');
fprintf('  n6_triggering_errors.pdf\n');
fprintf('  n6_3d_trajectories.pdf\n');
fprintf('  n6_mean_square_errors.pdf\n');
fprintf('  n6_residual_gates.pdf\n');
fprintf('  n6_threshold_candidates.csv\n');
fprintf('  n6_metrics_summary.csv\n');

end

%% =====================================================================
function cfg = make_configuration()

cfg.n = 6;
cfg.m = 3;
cfg.N = 6;

% A genuinely sixth-order translational model. Omega introduces weak
% cross-axis velocity coupling, so this is not three independently executed
% copies of an n=2 simulation.
damping = 0.25;
Omega = [ 0     0.08 -0.04;
         -0.08  0     0.06;
          0.04 -0.06  0   ];

cfg.A = [zeros(3), eye(3);
         zeros(3), Omega-damping*eye(3)];
cfg.B = [zeros(3);
         eye(3)];

cfg.Q = eye(6);
cfg.R = eye(3);

cfg.T  = 12;
cfg.dt = 0.01;
cfg.time = 0:cfg.dt:cfg.T;

% The adjacency convention is a_ij>0 for information j -> i.
cfg.a = [0 1 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 1 0 0 0;
         0 0 0 1 0 1;
         0 0 0 0 1 0];
cfg.g = ones(cfg.N,1);

% Event-triggering rule in the manuscript:
%   t_{k+1}^i = inf{t>t_k^i: ||e_i(t)|| >= c1*exp(-alpha*t)}.
cfg.c1Event    = 0.50;
cfg.alphaEvent = 0.30;

% Leader packets arrive intermittently and are processed in batches.
cfg.leaderPacketPeriod = 0.05;
cfg.leaderPacketStride = max(1,round(cfg.leaderPacketPeriod/cfg.dt));

% Bounded scalar adaptive coupling gain used in the revised manuscript:
%   dot{c_hat}_i = c_hat_i^{-2}*bar_xi_i'*bar_xi_i
% in the interior of the prescribed interval, with zero outward
% derivative at the upper projection boundary.
cfg.chatInitial = ones(cfg.N,1);
cfg.chatLower   = 0.50*ones(cfg.N,1);
cfg.chatUpper   = 12.0*ones(cfg.N,1);

% Bernoulli MITM model. At an attacked packet, the additive tampering
% vector satisfies ||epsilon|| <= varrho*||x_transmitted||.
cfg.deltaFF  = 0.10;
cfg.deltaLF  = 0.08;
cfg.varrhoFF = 0.80;
cfg.varrhoLF = 0.60;

% The original normalized pair is retained as one candidate. Additional
% n=6 candidates are generated from attack-free residual quantiles and are
% evaluated using the same closed-loop cost as Algorithm 1.
cfg.gammaFFn = 0.80;
cfg.gammaLFn = 0.50;
cfg.gammaFF  = Inf;
cfg.gammaLF  = Inf;
cfg.sigmaFloor = 1e-6;

cfg.numCalibrationRuns = 15;
cfg.numThresholdRuns = 10;
cfg.numMonteCarlo = 50;
cfg.baseSeed = 246810;
cfg.steadyStateFraction = 0.25;

cfg.weightError = 0.50;
cfg.weightTrigger = 0.20;
cfg.weightLeakage = 0.30;

cfg.outputDirectory = fullfile(pwd,'n6_results');
end

%% =====================================================================
function [Ad,Bd] = zoh_discretization(A,B,dt)
n = size(A,1);
m = size(B,2);
Md = expm([A B; zeros(m,n+m)]*dt);
Ad = Md(1:n,1:n);
Bd = Md(1:n,n+1:n+m);
end

%% =====================================================================
function [sigmaFF,sigmaLF,calibration] = estimate_residual_scales(cfg)
rFF = [];
rLF = [];

calCfg = cfg;
calCfg.gammaFF = Inf;
calCfg.gammaLF = Inf;

for run = 1:cfg.numCalibrationRuns
    seed = cfg.baseSeed - 1000 + run;
    [out,~] = simulate_once(calCfg,seed,false,false,true);
    rFF = [rFF; out.normalResidualFF(:)]; %#ok<AGROW>
    rLF = [rLF; out.normalResidualLF(:)]; %#ok<AGROW>
end

if isempty(rFF)
    error('No FF residual was collected during calibration.');
end
if isempty(rLF)
    error('No LF residual was collected during calibration.');
end

sigmaFF = std(rFF,0);
sigmaLF = std(rLF,0);

% With an exact leader model, attack-free LF residuals can be zero up to
% numerical precision. The floor only prevents a zero numerical threshold.
sigmaFF = max(sigmaFF,cfg.sigmaFloor);
sigmaLF = max(sigmaLF,cfg.sigmaFloor);

calibration.rFF = rFF;
calibration.rLF = rLF;
calibration.qFF = empirical_quantiles(rFF,[0.80 0.90 0.95 0.99]);
calibration.qLF = empirical_quantiles(rLF,[0.80 0.90 0.95 0.99]);

fprintf('  FF residual quantiles [80%% 90%% 95%% 99%%]:\n');
disp(calibration.qFF);
fprintf('  LF residual quantiles [80%% 90%% 95%% 99%%]:\n');
disp(calibration.qLF);
end

%% =====================================================================
function [cfg,selectionTable] = select_threshold_pair(cfg,sigmaFF,sigmaLF,calibration)
% Algorithm-1 evaluation for the sixth-order residual distribution. The
% normalized grid contains the original candidate 0.8 and dimension-aware
% candidates obtained by mapping attack-free residual quantiles to gamma_n.

ffRatios = calibration.qFF/max(sigmaFF,cfg.sigmaFloor);
q99Ratio = ffRatios(end);
ffList = [0.80,ffRatios,1.5*q99Ratio,2.0*q99Ratio, ...
    3.0*q99Ratio,4.0*q99Ratio];
ffList = ffList(isfinite(ffList) & ffList>0);
ffList = unique(round(100*ffList)/100,'stable');

% The attack-free LF residual is zero up to numerical precision because the
% leader and its predictor use the same known A. Hence, changing gamma_Ln
% only rescales a numerical floor; retain the manuscript value 0.5.
lfList = 0.50;

numberCandidates = numel(ffList)*numel(lfList);
rows = zeros(numberCandidates,10);
row = 0;

for i = 1:numel(ffList)
    for j = 1:numel(lfList)
        testCfg = cfg;
        testCfg.gammaFFn = ffList(i);
        testCfg.gammaLFn = lfList(j);
        testCfg.gammaFF = testCfg.gammaFFn*max(sigmaFF,cfg.sigmaFloor);
        testCfg.gammaLF = testCfg.gammaLFn*max(sigmaLF,cfg.sigmaFloor);

        candidateMetrics = repmat(empty_metrics(),cfg.numThresholdRuns,1);
        for run = 1:cfg.numThresholdRuns
            seed = cfg.baseSeed+5000+run;
            [~,candidateMetrics(run)] = simulate_once( ...
                testCfg,seed,true,false,false);
        end

        Ebar = mean([candidateMetrics.horizonMSE]);
        taubar = mean([candidateMetrics.triggerRate]);
        ellFF = mean_omitnan([candidateMetrics.ffLeakageRate]);
        ellLF = mean_omitnan([candidateMetrics.lfLeakageRate]);
        ellbar = mean_omitnan([ellFF ellLF]);
        falseRejectionFF = mean_omitnan( ...
            [candidateMetrics.ffFalseRejectionRate]);
        acceptanceFF = mean_omitnan([candidateMetrics.ffAcceptanceRate]);

        J = cfg.weightError*Ebar+cfg.weightTrigger*taubar+ ...
            cfg.weightLeakage*ellbar;

        row = row+1;
        rows(row,:) = [testCfg.gammaFFn,testCfg.gammaLFn, ...
            testCfg.gammaFF,testCfg.gammaLF,Ebar,taubar,ellFF,ellLF, ...
            falseRejectionFF,J];

        fprintf(['  gammaFFn=%5.2f, gammaLFn=%4.2f: ', ...
            'Ebar=%.4g, tau=%.4g, ellFF=%.4g, FR_FF=%.4g, J=%.4g\n'], ...
            testCfg.gammaFFn,testCfg.gammaLFn,Ebar,taubar,ellFF, ...
            falseRejectionFF,J);

        if acceptanceFF==0
            fprintf('    warning: this candidate rejects all FF packets.\n');
        end
    end
end

selectionTable = array2table(rows,'VariableNames', ...
    {'gammaFFn','gammaLFn','gammaFF','gammaLF','E_bar','tau_bar', ...
     'ell_FF','ell_LF','FR_FF','J'});

% Exclude the all-rejection deadlock whenever at least one feasible
% candidate has an FF false-rejection rate below one.
feasible = selectionTable.FR_FF < 1-1e-12;
if any(feasible)
    feasibleRows = find(feasible);
    [~,localIndex] = min(selectionTable.J(feasible));
    best = feasibleRows(localIndex);
else
    [~,best] = min(selectionTable.J);
end

cfg.gammaFFn = selectionTable.gammaFFn(best);
cfg.gammaLFn = selectionTable.gammaLFn(best);
cfg.gammaFF = selectionTable.gammaFF(best);
cfg.gammaLF = selectionTable.gammaLF(best);

selectionTable.Selected = false(height(selectionTable),1);
selectionTable.Selected(best) = true;
end

%% =====================================================================
function [out,metric] = simulate_once(cfg,seed,attacksEnabled,recordHistory,collectCalibration)
rng(seed,'twister');

n = cfg.n;
m = cfg.m;
N = cfg.N;
time = cfg.time;
numSteps = numel(time);

% Leader and follower initial conditions.
x0 = [2.0; -1.0; 1.2; 0.70; 0.35; -0.25];
x = repmat(x0,1,N) + [1.50*randn(3,N); 0.45*randn(3,N)];

chat = min(max(cfg.chatInitial,cfg.chatLower),cfg.chatUpper);
net = initialize_network(x,x0,n,N);
zeta = compute_disagreement(x,net,cfg.a,cfg.g);
net.zetaHold = zeta;

stats = initialize_statistics();

if recordHistory
    out.time = time;
    out.x0 = zeros(n,numSteps);
    out.x = zeros(n,N,numSteps);
    out.positionMSE = zeros(1,numSteps);
    out.velocityMSE = zeros(1,numSteps);
    out.totalMSE = zeros(1,numSteps);
    out.componentMaxError = zeros(n,numSteps);
    out.eventError = zeros(N,numSteps);
    out.eventThreshold = zeros(N,numSteps);
    out.couplingGain = zeros(N,numSteps);
    out.control = zeros(m,N,numSteps);
    out.controlNorm = zeros(N,numSteps);
    out.zetaNorm = zeros(N,numSteps);
else
    out = struct();
end

for k = 1:numSteps
    t = time(k);

    % Predictor flow from the previous sampling instant to the current one.
    if k > 1
        net = propagate_predictors(net,cfg.Ad);
    end

    % Intermittent LF packet batch: residual is evaluated before reset.
    if mod(k-1,cfg.leaderPacketStride)==0
        [net,stats] = process_leader_batch(net,x0,cfg,t,stats, ...
            attacksEnabled,collectCalibration);
    end

    % Pre-transmission disagreement and event-triggering error.
    zeta = compute_disagreement(x,net,cfg.a,cfg.g);
    eventError = column_norms(net.zetaHold-zeta);
    eventThreshold = cfg.c1Event*exp(-cfg.alphaEvent*t)*ones(N,1);

    [net,stats] = process_follower_transmissions(net,x,zeta,cfg,t,stats, ...
        eventError,eventThreshold,attacksEnabled,collectCalibration);

    % Recompute the disagreement after all accepted predictor resets.
    zetaAfterPackets = compute_disagreement(x,net,cfg.a,cfg.g);

    % Bounded scalar adaptive law in the revised manuscript.
    for i = 1:N
        zi = net.zetaHold(:,i);
        rawDerivative = (zi'*zi)/(chat(i)^2+1e-12);

        if chat(i) >= cfg.chatUpper(i) && rawDerivative > 0
            chatDot = 0;
        elseif chat(i) <= cfg.chatLower(i) && rawDerivative < 0
            chatDot = 0;
        else
            chatDot = rawDerivative;
        end

        chat(i) = min(cfg.chatUpper(i), ...
            max(cfg.chatLower(i),chat(i)+cfg.dt*chatDot));
    end

    % Sampled-data control input.
    u = zeros(m,N);
    for i = 1:N
        u(:,i) = -chat(i)*cfg.K*net.zetaHold(:,i);
    end

    % Record trajectories and performance before the next flow step.
    d = x-repmat(x0,1,N);
    positionMSE = mean(sum(d(1:3,:).^2,1));
    velocityMSE = mean(sum(d(4:6,:).^2,1));
    totalMSE = mean(sum(d.^2,1));

    stats.positionMSE(k) = positionMSE;
    stats.velocityMSE(k) = velocityMSE;
    stats.totalMSE(k) = totalMSE;

    if recordHistory
        out.x0(:,k) = x0;
        out.x(:,:,k) = x;
        out.positionMSE(k) = positionMSE;
        out.velocityMSE(k) = velocityMSE;
        out.totalMSE(k) = totalMSE;
        out.componentMaxError(:,k) = max(abs(d),[],2);
        out.eventError(:,k) = eventError;
        out.eventThreshold(:,k) = eventThreshold;
        out.couplingGain(:,k) = chat;
        out.control(:,:,k) = u;
        out.controlNorm(:,k) = column_norms(u);
        out.zetaNorm(:,k) = column_norms(zetaAfterPackets);
    end

    % Exact zero-order-hold discretization of follower dynamics.
    if k < numSteps
        x = cfg.Ad*x+cfg.Bd*u;
        x0 = cfg.Ad*x0;
    end
end

if collectCalibration
    out.normalResidualFF = stats.normalResidualFF;
    out.normalResidualLF = stats.normalResidualLF;
end

if recordHistory
    out.ffLog = stats.ffLog;
    out.lfLog = stats.lfLog;
end

steadyStart = max(1,floor((1-cfg.steadyStateFraction)*numSteps));
idx = steadyStart:numSteps;

metric = empty_metrics();
metric.horizonMSE = mean(stats.totalMSE);
metric.positionMSE = mean(stats.positionMSE(idx));
metric.velocityMSE = mean(stats.velocityMSE(idx));
metric.totalMSE = mean(stats.totalMSE(idx));
finalError = abs(x-repmat(x0,1,N));
metric.finalMaxError = max(finalError(:));
metric.triggerRate = stats.triggerCount/(N*numSteps);
metric.ffAcceptanceRate = safe_ratio(stats.ffAccepted,stats.ffTotal);
metric.lfAcceptanceRate = safe_ratio(stats.lfAccepted,stats.lfTotal);
metric.ffFalseRejectionRate = safe_ratio(stats.ffNormalRejected,stats.ffNormalTotal);
metric.lfFalseRejectionRate = safe_ratio(stats.lfNormalRejected,stats.lfNormalTotal);
metric.ffLeakageRate = safe_ratio(stats.ffAttackAccepted,stats.ffAttackTotal);
metric.lfLeakageRate = safe_ratio(stats.lfAttackAccepted,stats.lfAttackTotal);
metric.maxCouplingGain = max(chat);
metric.cpuTime = NaN;
end

%% =====================================================================
function net = initialize_network(x,x0,n,N)
net.xhat = zeros(n,N,N);
for receiver = 1:N
    for sender = 1:N
        net.xhat(:,receiver,sender) = x(:,sender);
    end
end
net.xhat0 = repmat(x0,1,N);
net.zetaHold = zeros(n,N);
net.lastTransmission = -inf(N,1);
end

function net = propagate_predictors(net,Ad)
[~,N,~] = size(net.xhat);
for receiver = 1:N
    for sender = 1:N
        net.xhat(:,receiver,sender) = Ad*net.xhat(:,receiver,sender);
    end
    net.xhat0(:,receiver) = Ad*net.xhat0(:,receiver);
end
end

function zeta = compute_disagreement(x,net,a,g)
[n,N] = size(x);
zeta = zeros(n,N);
for i = 1:N
    for j = 1:N
        if a(i,j) ~= 0
            zeta(:,i) = zeta(:,i)+a(i,j)*(x(:,i)-net.xhat(:,i,j));
        end
    end
    zeta(:,i) = zeta(:,i)+g(i)*(x(:,i)-net.xhat0(:,i));
end
end

%% =====================================================================
function [net,stats] = process_leader_batch(net,x0,cfg,t,stats,attacksEnabled,collectCalibration)
for i = 1:cfg.N
    isAttack = attacksEnabled && (rand < cfg.deltaLF);
    received = x0;

    if isAttack
        epsilonLF = tampering_component(x0,cfg.varrhoLF);
        received = x0+epsilonLF;
    end

    residual = norm(received-net.xhat0(:,i));
    accepted = residual <= cfg.gammaLF;

    stats.lfTotal = stats.lfTotal+1;
    stats.lfAccepted = stats.lfAccepted+double(accepted);

    if isAttack
        stats.lfAttackTotal = stats.lfAttackTotal+1;
        stats.lfAttackAccepted = stats.lfAttackAccepted+double(accepted);
    else
        stats.lfNormalTotal = stats.lfNormalTotal+1;
        stats.lfNormalRejected = stats.lfNormalRejected+double(~accepted);
        if collectCalibration
            stats.normalResidualLF(end+1,1) = residual; 
        end
    end

    if i == 1
        stats.lfLog.time(end+1,1) = t;
        stats.lfLog.residual(end+1,1) = residual;
        stats.lfLog.isAttack(end+1,1) = isAttack;
        stats.lfLog.accepted(end+1,1) = accepted;
    end

    if accepted
        net.xhat0(:,i) = received;
    end
end
end

function [net,stats] = process_follower_transmissions(net,x,zeta,cfg,t,stats, ...
    eventError,eventThreshold,attacksEnabled,collectCalibration)

for sender = 1:cfg.N
    triggered = eventError(sender) >= eventThreshold(sender);

    if ~triggered
        continue;
    end

    stats.triggerCount = stats.triggerCount+1;
    net.zetaHold(:,sender) = zeta(:,sender);
    net.lastTransmission(sender) = t;

    transmitted = x(:,sender);

    for receiver = 1:cfg.N
        if cfg.a(receiver,sender)==0
            continue;
        end

        isAttack = attacksEnabled && (rand < cfg.deltaFF);
        received = transmitted;

        if isAttack
            epsilonFF = tampering_component(transmitted,cfg.varrhoFF);
            received = transmitted+epsilonFF;
        end

        residual = norm(received-net.xhat(:,receiver,sender));
        accepted = residual <= cfg.gammaFF;

        stats.ffTotal = stats.ffTotal+1;
        stats.ffAccepted = stats.ffAccepted+double(accepted);

        if isAttack
            stats.ffAttackTotal = stats.ffAttackTotal+1;
            stats.ffAttackAccepted = stats.ffAttackAccepted+double(accepted);
        else
            stats.ffNormalTotal = stats.ffNormalTotal+1;
            stats.ffNormalRejected = stats.ffNormalRejected+double(~accepted);
            if collectCalibration
                stats.normalResidualFF(end+1,1) = residual; %#ok<AGROW>
            end
        end

        % Representative FF link 2 -> 1 for residual visualization.
        if receiver==1 && sender==2
            stats.ffLog.time(end+1,1) = t;
            stats.ffLog.residual(end+1,1) = residual;
            stats.ffLog.isAttack(end+1,1) = isAttack;
            stats.ffLog.accepted(end+1,1) = accepted;
        end

        if accepted
            net.xhat(:,receiver,sender) = received;
        end
    end
end
end

function epsilon = tampering_component(reference,varrho)
direction = randn(size(reference));
direction = direction/(norm(direction)+1e-12);
% A declared attack event produces a non-negligible tampering component,
% while the sector bound ||epsilon|| <= varrho*||x|| remains satisfied.
relativeMagnitude = varrho*(0.5+0.5*rand);
epsilon = relativeMagnitude*norm(reference)*direction;
end

%% =====================================================================
function stats = initialize_statistics()
stats.triggerCount = 0;

stats.ffTotal = 0;
stats.ffAccepted = 0;
stats.ffNormalTotal = 0;
stats.ffNormalRejected = 0;
stats.ffAttackTotal = 0;
stats.ffAttackAccepted = 0;

stats.lfTotal = 0;
stats.lfAccepted = 0;
stats.lfNormalTotal = 0;
stats.lfNormalRejected = 0;
stats.lfAttackTotal = 0;
stats.lfAttackAccepted = 0;

stats.normalResidualFF = [];
stats.normalResidualLF = [];

stats.ffLog = empty_packet_log();
stats.lfLog = empty_packet_log();

stats.positionMSE = [];
stats.velocityMSE = [];
stats.totalMSE = [];
end

function logData = empty_packet_log()
logData.time = [];
logData.residual = [];
logData.isAttack = [];
logData.accepted = [];
end

function metric = empty_metrics()
metric.horizonMSE = NaN;
metric.positionMSE = NaN;
metric.velocityMSE = NaN;
metric.totalMSE = NaN;
metric.finalMaxError = NaN;
metric.triggerRate = NaN;
metric.ffAcceptanceRate = NaN;
metric.lfAcceptanceRate = NaN;
metric.ffFalseRejectionRate = NaN;
metric.lfFalseRejectionRate = NaN;
metric.ffLeakageRate = NaN;
metric.lfLeakageRate = NaN;
metric.maxCouplingGain = NaN;
metric.cpuTime = NaN;
end

function values = column_norms(X)
values = sqrt(sum(X.^2,1))';
end

function value = safe_ratio(numerator,denominator)
if denominator==0
    value = NaN;
else
    value = numerator/denominator;
end
end

function value = mean_omitnan(values)
values = values(~isnan(values));
if isempty(values)
    value = 0;
else
    value = mean(values);
end
end

function q = empirical_quantiles(samples,probabilities)
samples = sort(samples(:));
numberSamples = numel(samples);
q = zeros(size(probabilities));

for k = 1:numel(probabilities)
    probability = min(1,max(0,probabilities(k)));
    index = 1+(numberSamples-1)*probability;
    lowerIndex = floor(index);
    upperIndex = ceil(index);

    if lowerIndex==upperIndex
        q(k) = samples(lowerIndex);
    else
        fraction = index-lowerIndex;
        q(k) = (1-fraction)*samples(lowerIndex)+ ...
            fraction*samples(upperIndex);
    end
end
end

%% =====================================================================
function runTable = metrics_to_table(metrics)
runTable = struct2table(metrics);
runTable.Run = (1:height(runTable))';
runTable = movevars(runTable,'Run','Before',1);
end

function summaryTable = summarize_metrics(runTable)
names = runTable.Properties.VariableNames;
names(strcmp(names,'Run')) = [];

metricName = strings(numel(names),1);
meanValue = zeros(numel(names),1);
stdValue = zeros(numel(names),1);
medianValue = zeros(numel(names),1);
percentile95 = zeros(numel(names),1);
maximumValue = zeros(numel(names),1);

for k = 1:numel(names)
    metricName(k) = string(names{k});
    values = runTable.(names{k});
    meanValue(k) = mean(values,'omitnan');
    stdValue(k) = std(values,'omitnan');
    validValues = values(~isnan(values));
    medianValue(k) = median(validValues);
    percentile95(k) = empirical_quantiles(validValues,0.95);
    maximumValue(k) = max(validValues);
end

summaryTable = table(metricName,meanValue,stdValue,medianValue, ...
    percentile95,maximumValue, ...
    'VariableNames',{'Metric','Mean','StandardDeviation','Median', ...
    'Percentile95','Maximum'});
end

%% =====================================================================
function make_adaptive_gain_figure(cfg,out)
fig = figure('Color','w','Position',[120 120 1000 520]);
ax = axes(fig);
hold(ax,'on');

colors = [0.0000 0.4470 0.7410;
          0.8500 0.3250 0.0980;
          0.9290 0.6940 0.1250;
          0.4940 0.1840 0.5560;
          0.4660 0.6740 0.1880;
          0.3010 0.7450 0.9330];
markers = {'o','>','s','d','^','v'};
markerStep = max(1,round(numel(out.time)/32));

for i = 1:cfg.N
    plot(ax,out.time,out.couplingGain(i,:), ...
        'Color',colors(i,:),'LineWidth',1.5,'Marker',markers{i}, ...
        'MarkerIndices',1:markerStep:numel(out.time),'MarkerSize',4, ...
        'DisplayName',sprintf('Agent %d',i));
end

xlabel(ax,'Time (s)','Interpreter','latex');
ylabel(ax,'$\hat c_i(t)$','Interpreter','latex');
xlim(ax,[0 cfg.T]);
legend(ax,'Location','north','Orientation','horizontal', ...
    'NumColumns',cfg.N,'Interpreter','latex');
apply_axes_style(ax);
save_publication_figure(fig,fullfile(cfg.outputDirectory,'n6_adaptive_gains'));
end

function make_tracking_error_figure(cfg,out)
fig = figure('Color','w','Position',[40 40 1450 940]);
tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
colors = lines(cfg.n);

for i = 1:cfg.N
    ax = nexttile;
    Xi = squeeze(out.x(:,i,:));
    trackingError = Xi-out.x0;
    hold(ax,'on');
    for stateIndex = 1:cfg.n
        plot(ax,out.time,trackingError(stateIndex,:), ...
            'Color',colors(stateIndex,:),'LineWidth',1.15, ...
            'DisplayName',sprintf('$e_{%d%d}$',i,stateIndex));
    end
    xlabel(ax,'Time (s)','Interpreter','latex');
    ylabel(ax,'Tracking error','Interpreter','latex');
    xlim(ax,[0 cfg.T]);
    legend(ax,'Location','north','Orientation','horizontal', ...
        'NumColumns',3,'Interpreter','latex','FontSize',8);
    apply_axes_style(ax);
end

save_publication_figure(fig,fullfile(cfg.outputDirectory,'n6_tracking_errors'));
end

function make_control_input_figure(cfg,out)
fig = figure('Color','w','Position',[45 45 1450 940]);
tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
controlColors = [0 0.20 0.90; 0.90 0.05 0.05; 0.10 0.60 0.20];

for i = 1:cfg.N
    ax = nexttile;
    Ui = squeeze(out.control(:,i,:));
    hold(ax,'on');
    for inputIndex = 1:cfg.m
        plot(ax,out.time,Ui(inputIndex,:), ...
            'Color',controlColors(inputIndex,:),'LineWidth',1.15, ...
            'DisplayName',sprintf('$\\nu_{%d%d}$',i,inputIndex));
    end
    xlabel(ax,'Time (s)','Interpreter','latex');
    ylabel(ax,'Control input','Interpreter','latex');
    xlim(ax,[0 cfg.T]);
    legend(ax,'Location','north','Orientation','horizontal', ...
        'NumColumns',cfg.m,'Interpreter','latex','FontSize',9);
    apply_axes_style(ax);
end

save_publication_figure(fig,fullfile(cfg.outputDirectory,'n6_control_inputs'));
end

function make_trigger_error_figure(cfg,out)
fig = figure('Color','w','Position',[50 50 1400 930]);
tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');

for i = 1:cfg.N
    ax = nexttile;
    plot(ax,out.time,out.eventError(i,:),'r-','LineWidth',1.05, ...
        'DisplayName',sprintf('$\\|e_%d(t)\\|$',i));
    hold(ax,'on');
    plot(ax,out.time,out.eventThreshold(i,:),'k--','LineWidth',1.25, ...
        'DisplayName','Threshold');
    xlabel(ax,'Time (s)','Interpreter','latex');
    ylabel(ax,sprintf('$\\|e_%d(t)\\|$',i),'Interpreter','latex');
    xlim(ax,[0 cfg.T]);
    ylim(ax,[0,1.08*max([out.eventError(i,:),out.eventThreshold(i,:)])]);
    legend(ax,'Location','northeast','Interpreter','latex','FontSize',9);
    apply_axes_style(ax);
end

save_publication_figure(fig,fullfile(cfg.outputDirectory,'n6_triggering_errors'));
end

function make_trajectory_figure(cfg,out)
fig = figure('Color','w','Position',[120 80 900 720]);
ax = axes(fig);
hold(ax,'on');
colors = lines(cfg.N);

for i = 1:cfg.N
    Xi = squeeze(out.x(:,i,:));
    plot3(ax,Xi(1,:),Xi(2,:),Xi(3,:),'Color',colors(i,:), ...
        'LineWidth',1.25,'DisplayName',sprintf('Agent %d',i));
    plot3(ax,Xi(1,1),Xi(2,1),Xi(3,1),'o','Color',colors(i,:), ...
        'MarkerFaceColor',colors(i,:),'MarkerSize',5, ...
        'HandleVisibility','off');
end
plot3(ax,out.x0(1,:),out.x0(2,:),out.x0(3,:),'k--', ...
    'LineWidth',2.0,'DisplayName','Leader');

xlabel(ax,'$p_x$','Interpreter','latex');
ylabel(ax,'$p_y$','Interpreter','latex');
zlabel(ax,'$p_z$','Interpreter','latex');
legend(ax,'Location','best','Interpreter','latex');
axis(ax,'equal');
view(ax,38,24);
apply_axes_style(ax);
save_publication_figure(fig,fullfile(cfg.outputDirectory,'n6_3d_trajectories'));
end

function make_mse_figure(cfg,out)
fig = figure('Color','w','Position',[130 130 1000 530]);
ax = axes(fig);
semilogy(ax,out.time,max(out.positionMSE,1e-14),'b-', ...
    'LineWidth',1.5,'DisplayName','$E_p(t)$');
hold(ax,'on');
semilogy(ax,out.time,max(out.velocityMSE,1e-14),'r-', ...
    'LineWidth',1.5,'DisplayName','$E_v(t)$');
semilogy(ax,out.time,max(out.totalMSE,1e-14),'k--', ...
    'LineWidth',1.5,'DisplayName','$E_x(t)$');
xlabel(ax,'Time (s)','Interpreter','latex');
ylabel(ax,'Mean-square synchronization error','Interpreter','latex');
xlim(ax,[0 cfg.T]);
legend(ax,'Location','best','Interpreter','latex');
apply_axes_style(ax);
save_publication_figure(fig,fullfile(cfg.outputDirectory,'n6_mean_square_errors'));
end

function make_residual_gate_figure(cfg,out)
fig = figure('Color','w','Position',[80 180 1400 500]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
plot_packet_residual(out.ffLog,cfg.gammaFF,'FF link $2\to1$');
apply_axes_style(ax1);

ax2 = nexttile;
plot_packet_residual(out.lfLog,cfg.gammaLF,'LF link $0\to1$');
apply_axes_style(ax2);

save_publication_figure(fig,fullfile(cfg.outputDirectory,'n6_residual_gates'));
end

function plot_packet_residual(logData,gammaValue,titleText)
if isempty(logData.time)
    axis off;
    text(0.5,0.5,'No packet was recorded','HorizontalAlignment','center');
    return;
end

normalIndex = ~logical(logData.isAttack);
attackRejected = logical(logData.isAttack) & ~logical(logData.accepted);
attackAccepted = logical(logData.isAttack) & logical(logData.accepted);

hold on;
scatter(logData.time(normalIndex),logData.residual(normalIndex),14, ...
    [0.2 0.45 0.8],'filled','DisplayName','Normal packet');
scatter(logData.time(attackRejected),logData.residual(attackRejected),24, ...
    [0.85 0.2 0.2],'x','LineWidth',1.2,'DisplayName','Rejected attack');
scatter(logData.time(attackAccepted),logData.residual(attackAccepted),28, ...
    [0.55 0.1 0.65],'d','filled','DisplayName','Admitted attack');
plot([logData.time(1),logData.time(end)],[gammaValue,gammaValue],'k--', ...
    'LineWidth',1.2,'DisplayName','Gate threshold');
grid on;
xlabel('Time (s)','Interpreter','latex');
ylabel('Residual','Interpreter','latex');
title(titleText,'Interpreter','latex');
legend('Location','best');
end

function apply_axes_style(ax)
set(ax,'FontName','Times New Roman','FontSize',11, ...
    'LineWidth',0.8,'Box','on','TickDir','in', ...
    'XMinorTick','on','YMinorTick','on');
grid(ax,'on');
ax.GridColor = [0.84 0.84 0.84];
ax.GridAlpha = 0.75;
ax.MinorGridAlpha = 0.25;
end

function save_publication_figure(fig,baseName)
savefig(fig,[baseName '.fig']);

if exist('exportgraphics','file')==2
    exportgraphics(fig,[baseName '.pdf'],'ContentType','vector');
    exportgraphics(fig,[baseName '.png'],'Resolution',300);
else
    print(fig,[baseName '.pdf'],'-dpdf','-painters');
    print(fig,[baseName '.png'],'-dpng','-r300');
end
end
