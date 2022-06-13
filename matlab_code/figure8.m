gamma_inh=40;
rate=5;
threshold=0.5;

%Initial Wsoc and dynamics:


    
W=initialnet_gamma(200, 0.1, gamma_inh); %initial random matrix
[Wsoc, e] = create_inh_soc_gamma(W, rate, threshold,gamma_inh);  %SOC network
Wsoc=100/norm(Wsoc, 'fro') *Wsoc; %normalise norm
    
%[eigenvectors, trial(times).abs_eigvals, ~] = MaximiseIC(Wsoc, 1);
[abs_eigvals, percentage_amplified]=max_norm_analytical(Wsoc);

energy=sort(abs_eigvals, 'descend');
    
    
   
 NN = length(Wsoc);
[eigenvectors, ~] = MaximiseIC(Wsoc, 1);
initial_cond = 20*eigenvectors;



%% initialise parameters
n_exc = NN/2; r0 = 20; rmax = 100; tau = 200;
tfinal = 1000; gain_fct_flag = 'NL'; plasticity_noise = 0; %noise=Noisefn(NN,tfinal); 
noise=zeros(NN,tfinal)';
gains0 = ones(NN,1); %Initial gains;


%% Run dynamics
dyncurrent = initialise_rate_dynamics_hetero_n( Wsoc, gains0, r0, rmax, ...
    gain_fct_flag, plasticity_noise, tau, tfinal, 0, initial_cond(:,1), noise);

 
 




    %PLOTS

figure;
subplot(2,2,1);
plot(dyncurrent.X)
ylabel('r (Hz)')
xlabel('time (ms)')
set(gca, 'TickDir', 'out');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
ylim([-20, 40])
box off



subplot(2,2,3);
plot(energy)
ylabel('max norm')
xlabel('condition (%)')
set(gca, 'TickDir', 'out');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
box off
ylim([0, 5])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Now scale inhibitory weights down by 40:%%%%%


Wsoc_new=Wsoc; %SOC network
Wsoc_new(:,101:end)= 1/40 * Wsoc(:, 101:end); %scale inhibitory weights

    

[abs_eigvals, percentage_amplified]=max_norm_analytical(Wsoc_new);

energy_new=sort(abs_eigvals, 'descend');
    
    
   
 NN = length(Wsoc);
[eigenvectors, ~] = MaximiseIC(Wsoc_new, 1);
initial_cond = 20*eigenvectors;



%% initialise parameters
n_exc = NN/2; r0 = 20; rmax = 100; tau = 200;
tfinal = 1000; gain_fct_flag = 'NL'; plasticity_noise = 0; %noise=Noisefn(NN,tfinal); 
noise=zeros(NN,tfinal)';
gains0 = ones(NN,1); %Initial gains;


%% Run dynamics
dyncurrent_new = initialise_rate_dynamics_hetero_n( Wsoc_new, gains0, r0, rmax, ...
    gain_fct_flag, plasticity_noise, tau, tfinal, 0, initial_cond(:,1), noise);

 
 




    %PLOTS

subplot(2,2,2);
plot(dyncurrent_new.X)
ylabel('r (Hz)')
xlabel('time (ms)')
set(gca, 'TickDir', 'out');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
ylim([-20, 40])
box off


subplot(2,2,4);
plot(energy_new)
ylabel('max norm')
xlabel('condition (%)')
set(gca, 'TickDir', 'out');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
box off
ylim([0, 5])





    