%computing effective rank of the eigenvector matrix not of the principal
%components either as a function of the imaginary diameter or as a function
%of the ffd strength (comment/uncomment accordingly)

index=0;
TIMES=20;

%str_vector=[0:10:500];
imaginary_vector=[0:10:500];

str=75;
for r=imaginary_vector
%r=20;
%for str=str_vector

    index=index+1;
    
    
    
for times=1:TIMES
    
    M=create_matrix(r);
    M=ffd_random(M,str);
    
   
    [eigvecs,~]=eig(M);

   S=svd(eigvecs);
   S_for_M=svd(M);

   ef_rank_reps(times)=effective_rank(S);
   ef_rank_M_reps(times)=effective_rank(S_for_M);

end
    
   
    ef_rank(index)=mean(ef_rank_reps);
    ef_rank_M(index)=mean(ef_rank_M_reps);
    
end

NN=length(M);
ef_rank=ef_rank*100/NN;
ef_rank_M=ef_rank*100/NN;

    
    figure;
    plot(ef_rank)
    xlabel('imaginary diameter')
    %xlabel('ffd norm')
    ylabel('effective rank of V (%)')
    set(gca, 'TickDir', 'out');
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    set(gca,'LineWidth',3);
    set(gca,'fontsize', 20);
    box off
    ylim([0 100])
    xlim([0 50])
   
    
    figure;
    plot(ef_rank_M)
    %xlabel('imaginary diameter')
    xlabel('ffd norm')
    ylabel('effective rank of M (%)')
    set(gca, 'TickDir', 'out');
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    set(gca,'LineWidth',3);
    set(gca,'fontsize', 20);
    box off
    ylim([0 100])
    xlim([0 50])
   
    
    
    
   