inhibitory_strengths=[3,5,10,20,30,40];
rate=5;
threshold=0.5;
TIMES=1;

index=0;

for gamma_inh=inhibitory_strengths
    index=index+1;
    
    for times=1:TIMES
    
    W=initialnet_gamma(200, 0.1, gamma_inh);
    [Wsoc, e] = create_inh_soc_gamma(W, rate, threshold,gamma_inh);
    Wsoc=100/norm(Wsoc, 'fro') *Wsoc;
   
    
    
[linear_percentage(times), nonlinear_percentage(times), eff_rank_linear(times), eff_rank_nonlinear(times)]=dim_analysis(Wsoc);
   
    
end
    
perc_ampl_linear(index)=mean(linear_percentage);
perc_ampl_nonlinear(index)=mean(nonlinear_percentage);

effect_rank_linear(index)=mean(eff_rank_linear);
effect_rank_nonlinear(index)=mean(eff_rank_nonlinear);
    
    
    
    
    
end



%PLOT

    figure;
    a1=plot(perc_ampl_linear); M1='amplified linear';
    hold on;
    a2=plot(perc_ampl_nonlinear); M2='amplified nonlinear';
    hold on;
    a3=plot(effect_rank_linear); M3='effective directions, linear';
    hold on;
    a4=plot(effect_rank_nonlinear); M4='effective directions, nonlinear';
    legend([a1,a2,a3,a4],M1,M2,M3,M4)






