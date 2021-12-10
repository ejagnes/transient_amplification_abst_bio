function [output] = initialise_rate_dynamics_hetero_n( W, gradient, r0, rmax, ...
    fr_flag, pl_noise, tau, tfinal, icnoise, initial_cond, noise )

%% Grab initial conditions for dynamics
if icnoise ~= 0
%     InitialCond = eigenvectors(:, 1) + 0.01*norm(eigenvectors(:,1))*randn(size(eigenvectors(:,1)));  %Grab initial condition 0.008
    InitialCond = awgn(initial_cond(:,1),icnoise, 'measured');
    %InitialCond=awgn(zeros(200,1),icnoise); %spontaneous
else
    InitialCond = initial_cond(:, 1);
end

%% Initialise all parameters
param.over_tau_d = 1/tau;
param.const = gradient;        %Const for the linear gain function
param.r0 = r0;
param.rmax = rmax;
param.a_g = rmax-r0;
param.b_g = gradient/(rmax-r0);
%options = odeset('RelTol', 1e-6);
param.tfinal = tfinal;

% Firing rate function to use
if strcmp(fr_flag, 'L')
    param.g = 'g_linear';
    param.gf = 'g_linear';
elseif strcmp(fr_flag, 'NL')
    param.g = 'g_non_linear';
    param.gf = 'g_final_nl';
else
    warning('Incorrect firing rate function flag given, using linear firing rate');
    param.g = 'g_linear';
end

%% Solve the ODEs

[output.t, output.X] = ode45(@rate_dynamics_inh_ode, ...
    linspace(0, tfinal, tfinal), InitialCond, [], W, param, noise);

param.const = repmat((param.const)', tfinal,1); %HERE replace tfinal with 400 for linspace of 400
output.X = feval(param.gf, output.X, param);
output.param = param;

end

function k = rate_dynamics_inh_ode(t, X, W, param, noise)

noise_current = interp1(linspace(0,param.tfinal,param.tfinal),noise,t);

%Network dynamics
k = -X + feval(param.g, W*X, param) + noise_current'; 
k = param.over_tau_d*k;

end

function out = g_non_linear(X, param)

out = zeros(size(X));
I = X < 0;
out(I) = param.r0*tanh(param.const(I).*X(I)/param.r0);
I2 = logical(1 - I);
out(I2) = ...
    (param.a_g)*tanh(param.b_g(I2).*X(I2));

end

function out = g_final_nl(X, param)

out = zeros(size(X));
I = X < 0;
out(I) = param.r0*tanh(param.const(I).*X(I)/param.r0);
I2 = logical(1 - I);
out(I2) = ...
    (param.rmax - param.r0)*tanh(param.const(I2).*X(I2)/(param.rmax-param.r0));

end

function out = g_linear(X, param)

out = param.const.*X;

end