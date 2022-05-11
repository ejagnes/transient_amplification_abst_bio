% real distribution is random uniform. 
str=75;
%possible_r=[0,0.01,0.1,1,10, 100, 1000, 10000];
   % pos = randi(length(possible_r));
    %r = possible_r(pos);
    r=100;
    
    
    M=create_matrix(r);
    [M,T1,T2]=ffd_random_and_parts(M,str);
    
    NN = length(M);
    
    
    %case 1, weak:
    
    x=700; %spectral norm
    y=100;  %ffd norm 
    
    
M_weak= x/norm(T1,'fro') * T1 + y/norm(T2, 'fro')*T2;
    
    [~, vals]=eig(M_weak);
    vals=diag(vals);
    vals=real(vals);
    abscissa=max(vals);
    
  if abscissa>0.5
      abscissa=abscissa-0.5;
      M_weak=M_weak-abscissa*eye(NN);
  end
  
  cond=1;
 
 [eigenvectors, abs_eigvals] = MaximiseIC(M_weak, 1);
 %initial_conditions = 20*eigenvectors;
  initial_conditions = eigenvectors;



NN=200;
%% initialise parameters
n_exc = NN/2; r0 = 1; rmax = 5; tau = 200;
tfinal = 10000; gain_fct_flag = 'NL'; plasticity_noise = 0; %noise=Noisefn(NN,tfinal); 
noise=zeros(NN,tfinal)';

gains0 = ones(NN,1); %Initial gains;

dyncurrent_weak = initialise_rate_dynamics_hetero_n(M_weak, gains0, r0, rmax, ...
    gain_fct_flag, plasticity_noise, tau, tfinal, 0, initial_conditions(:,cond), noise);


response_norm_weak=zeros(1,tfinal);  
for t=1:tfinal
    response_norm_weak(t)=norm(dyncurrent_weak.X(t,:));
end

normalised_norm_weak=response_norm_weak/sqrt(NN);


 %case 2,  short transient:
 
  x=500; %spectral norm
    y=500;  %ffd norm 
    
    
M_short= x/norm(T1,'fro') * T1 + y/norm(T2, 'fro')*T2;
    
    [~, vals]=eig(M_short);
    vals=diag(vals);
    vals=real(vals);
    abscissa=max(vals);
    
  if abscissa>0.5
      abscissa=abscissa-0.5;
      M_short=M_short-abscissa*eye(NN);
  end
  
  cond=1;
 
 [eigenvectors, abs_eigvals] = MaximiseIC(M_short, 1);
 %initial_conditions = 20*eigenvectors;
  initial_conditions = eigenvectors;



NN=200;
%% initialise parameters
n_exc = NN/2; r0 = 1; rmax = 5; tau = 200;
tfinal = 10000; gain_fct_flag = 'NL'; plasticity_noise = 0; %noise=Noisefn(NN,tfinal); 
noise=zeros(NN,tfinal)';

gains0 = ones(NN,1); %Initial gains;

dyncurrent_short = initialise_rate_dynamics_hetero_n(M_short, gains0, r0, rmax, ...
    gain_fct_flag, plasticity_noise, tau, tfinal, 0, initial_conditions(:,cond), noise);


response_norm_short=zeros(1,tfinal);  
for t=1:tfinal
    response_norm_short(t)=norm(dyncurrent_short.X(t,:));
end

normalised_norm_short=response_norm_short/sqrt(NN);


 %case 2,  long transient:
 
  x=100; %spectral norm
    y=700;  %ffd norm 
    
    
M_long= x/norm(T1,'fro') * T1 + y/norm(T2, 'fro')*T2;
    
    [~, vals]=eig(M_long);
    vals=diag(vals);
    vals=real(vals);
    abscissa=max(vals);
    
  if abscissa>0.5
      abscissa=abscissa-0.5;
      M_long=M_long-abscissa*eye(NN);
  end
  
  cond=1;
 
 [eigenvectors, abs_eigvals] = MaximiseIC(M_long, 1);
 %initial_conditions = 20*eigenvectors;
  initial_conditions = eigenvectors;



NN=200;
%% initialise parameters
n_exc = NN/2; r0 = 1; rmax = 5; tau = 200;
tfinal = 10000; gain_fct_flag = 'NL'; plasticity_noise = 0; %noise=Noisefn(NN,tfinal); 
noise=zeros(NN,tfinal)';

gains0 = ones(NN,1); %Initial gains;

dyncurrent_long = initialise_rate_dynamics_hetero_n(M_long, gains0, r0, rmax, ...
    gain_fct_flag, plasticity_noise, tau, tfinal, 0, initial_conditions(:,cond), noise);


response_norm_long=zeros(1,tfinal);  
for t=1:tfinal
    response_norm_long(t)=norm(dyncurrent_long.X(t,:));
end

normalised_norm_long=response_norm_long/sqrt(NN);



    
   



 
 
 
 
 
 
figure;
a1=plot(normalised_norm_weak); M1='weak';
hold on;
a2=plot(normalised_norm_short); M2='short transient';
hold on;
a3=plot(normalised_norm_long); M3='long transient';
hold on;
yline(1/sqrt(200))
legend([a1,a2,a3],M1,M2,M3)



    
    