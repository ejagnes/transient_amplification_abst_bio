%linear, analytical solution to max firing rate

function [m, percentage_amplified]=max_norm_analytical(W)

NN = length(W);
[eigenvectors, ~] = MaximiseIC(W, 1);
initial_cond=eigenvectors;


%% initialise parameters
 tau = 200;
tfinal = 1000;



for t=1:tfinal
    
 P_t=expm(t*(W-eye(NN))/tau);   
 
 for cond=1:NN
rate_norm(t, cond)=norm(P_t*initial_cond(:,cond));
 end
end

for cond=1:NN
m(cond)=max(rate_norm(:, cond));
end

amplified=0;

for k=1:NN
    if m(k)>1.5
        amplified=amplified+1;
    end
end

percentage_amplified=amplified*100/length(W);


