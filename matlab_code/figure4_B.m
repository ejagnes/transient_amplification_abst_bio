% script for figure 4.B You can use the same script for both panels by
% commenting/uncommenting the relevant part. You can either run those as a
% function of the imaginary vector or as a function of the str_vector.

index=0;
TIMES=1;

%str_vector=[10:50:510];
imaginary_vector=[0,10,20,30,40,50,60,70,80,90,100];


str=75;
for r=imaginary_vector
%r=20;
%for str=str_vector
  index=index+1;
    
    



    
  for times=1:TIMES
        
        

     M=create_matrix(r);    %create spectrum

     M=ffd_random(M, str);  %create random feedforward structure

    [linear_percentage(times), nonlinear_percentage(times), eff_rank_linear(times), eff_rank_nonlinear(times)]=dim_analysis(M);

  end



perc_ampl_linear(index)=mean(linear_percentage);
perc_ampl_nonlinear(index)=mean(nonlinear_percentage);

effect_rank_linear(index)=mean(eff_rank_linear);
effect_rank_nonlinear(index)=mean(eff_rank_nonlinear);


end



    figure;
    a1=plot(perc_ampl_linear); M1='% amplified linear';
    hold on;
    a2=plot(perc_ampl_nonlinear); M2='% amplified nonlinear';
    hold on;
    a3=plot(effect_rank_linear); M3='effective rank, linear';
    hold on;
    a4=plot(effect_rank_nonlinear); M4='effective rank, nonlinear';
    legend([a1,a2,a3,a4],M1,M2,M3,M4)
    ylabel('% in neuornal space')
    %xlabel('ffd norm')
    xlabel('imaginary diameter')
    set(gca, 'TickDir', 'out')
    set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
    set(gca,'LineWidth',3);
    set(gca,'fontsize', 20);
    box off;
    
    
   
    
    
   