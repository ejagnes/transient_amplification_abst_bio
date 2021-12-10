% script for figure 3.A and potentially 3.B if you run for str% 
index=0;
TIMES=1; %number of repetitions
%for str=[10:10:10000]  %either compute as function of ffd strength
    %r=20;
    
for r=[0,0.01,0.1,1,10, 100, 1000, 10000]  %is as function of imaginary diameter
str=75;


    
   
        
    index=index+1;
    for times=1:TIMES
        M=create_matrix(r);
        M=ffd_random(M,str);
        
        
        
        %compute percentage of eigenvector overlaps
        [hist_overlap, ~]=eigvecs_overlaps(M);

        D_overlap=0;
        for k=1:length(hist_overlap)
           if hist_overlap(k)<=45
              D_overlap=D_overlap+1;
           end
        end

        D_overlap=D_overlap*100/length(hist_overlap); %percentage of overlaps less than 45 degrees
        over(times)=D_overlap;
        
        %compute maximum norm and percentage of amplified conditions
        
        [m_all, percentage_ampl(times)]=max_norm_analytical(M);

        maximum_norm(times)=m_all(1);
        
        

 end
        
    
    overlaps(index)=mean(over);
    maxim_norm(index)=mean(maximum_norm);
    percent_ampl(index)=mean(percentage_ampl);

    
end





figure; 
plot(overlaps(1,:))
ylabel('% of overlap angles < 45')
xlabel('imaginary radius')
set(gca, 'TickDir', 'out')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
box off;
xlim([1 8]);
xticks([1  3  5  7 ])
xticklabels({'0','0.1','10', '1000'})


figure; 
plot(percent_ampl(1,:))
ylabel('% amplified directions from norm')
xlabel('imaginary radius')
set(gca, 'TickDir', 'out')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
box off;
xlim([1 8]);
xticks([1  3  5  7 ])
xticklabels({'0','0.1','10', '1000'})


figure; 
plot(log10(maxim_norm(1,:)))
ylabel('log of maximum firing norm')
xlabel('imaginary radius')
set(gca, 'TickDir', 'out')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
box off;
xlim([1 8]);
xticks([1  3  5  7 ])
xticklabels({'0','0.1','10', '1000'})







