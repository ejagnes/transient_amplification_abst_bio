function [linear_percentage, nonlinear_percentage, eff_rank_linear, eff_rank_nonlinear]=dim_analysis(M)
%function computing the effective rank for the dynamical matrix used in the
%paper as well as the percentage of amplified conditions
%the function returns results for the linear and the nonlinear case


NN = length(M);
[eigenvectors, abs_eigvals] = MaximiseIC(M, 1);
initial_conditions = eigenvectors;
  
pc_matrix_linear=[];
pc_matrix_nonlinear=[];
linear_amplified=0;
nonlinear_amplified=0;


  
for cond=1:NN
    
    
    linear_logical=0;
    nonlinear_logical=0;


%% initialise parameters
n_exc = NN/2; r0 = 1; rmax = 5; tau = 200;
tfinal = 10000; gain_fct_flag = 'L'; plasticity_noise = 0; %noise=Noisefn(NN,tfinal); 
noise=zeros(NN,tfinal)';

gains0 = ones(NN,1); %Initial gains;

dyncurrent = initialise_rate_dynamics_hetero_n(M, gains0, r0, rmax, ...
    gain_fct_flag, plasticity_noise, tau, tfinal, 0, initial_conditions(:,cond), noise);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for time=1:tfinal
     norm_vector(time)=norm(dyncurrent.X(time,:));
 end
%    
 m(cond)=max(norm_vector);
%    
      if m(cond) > 2
         linear_amplified=linear_amplified+1;
         linear_logical=1;
      end




if linear_logical==1

[coeff,score,latent,tsquared,explained,mu] = pca(dyncurrent.X);

pc_matrix_linear=[pc_matrix_linear, coeff(:,1)];


end


        
%now same analysis for nonlinear case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dyncurrent1 = initialise_rate_dynamics_hetero_n(M, gains0, r0, rmax, ...
    'NL', plasticity_noise, tau, tfinal, 0, initial_conditions(:,cond), noise);

for time=1:tfinal
     norm_vector(time)=norm(dyncurrent1.X(time,:));
 end
%    
 m1(cond)=max(norm_vector);
%    
      if m1(cond) > 2
         nonlinear_amplified=nonlinear_amplified+1;
         nonlinear_logical=1;
      end



if nonlinear_logical==1

[coeff1,score,latent,tsquared,explained1,mu] = pca(dyncurrent1.X);
pc_matrix_nonlinear=[pc_matrix_nonlinear, coeff1(:,1)];

end



end

%now compute effective rank of the matrix constructed from the Principal
%Components of the dynamics of the amplified conditions

s_linear = (svd(pc_matrix_linear))';
if length(s_linear)>0
eff_rank_linear=effective_rank(s_linear);
else
    eff_rank_linear=0;
end


s_nonlinear=(svd(pc_matrix_nonlinear))';

if length(s_nonlinear)>0
eff_rank_nonlinear=effective_rank(s_nonlinear);
else
    eff_rank_nonlinear=0;
end


%now find them as percentages
linear_percentage=linear_amplified*100/NN;
nonlinear_percentage=nonlinear_amplified*100/NN;
eff_rank_linear=eff_rank_linear*100/NN;
eff_rank_nonlinear=eff_rank_nonlinear*100/NN;





