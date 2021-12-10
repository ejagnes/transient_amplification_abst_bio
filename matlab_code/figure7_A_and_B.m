inhibitory_strengths=[3,5,10,20,30,40];
rate=5;
threshold=0.5;
TIMES=10; %tdifferent simulations to then average over

norm_spectrum=zeros(1,length(inhibitory_strengths)); 
norm_ffd=zeros(1,length(inhibitory_strengths));
energies=zeros(200,length(inhibitory_strengths));



column_index=0;

for gamma_inh=inhibitory_strengths
    column_index=column_index+1;
    
    for times=1:TIMES
    
    W=initialnet_gamma(200, 0.1, gamma_inh); %initial random matrix
    [Wsoc, e] = create_inh_soc_gamma(W, rate, threshold,gamma_inh);  %SOC network
    Wsoc=100/norm(Wsoc, 'fro') *Wsoc; %normalise norm
    
    %[eigenvectors, trial(times).abs_eigvals, ~] = MaximiseIC(Wsoc, 1);
    [trial(times).abs_eigvals, percentage_amplified]=max_norm_analytical(Wsoc);
    
    trial(times).abs_eigvals=trial(times).abs_eigvals';
   
    size_matrix=length(Wsoc);
 
 %finding the ffd and diagonal components of Schur decomposition
 
[U,T]=schur(Wsoc);
T2=T;  

[vecs,vals]=eig(Wsoc);

[~,D2] = cdf2rdf(vecs,vals);

[rId, cId] = find(D2);

indices=[rId,cId];

for i=1:length(indices)
    T2(indices(i,1), indices(i,2))=0;
end

T1=D2;

%computing the norm of the two components
trial(times).normSpectrum=norm(T1, 'fro');
trial(times).normffd=norm(T2, 'fro');

    
 end
    
    
  %now compute average over simulations for every inhibitory dominance value  
    
    matrix_energy=zeros(200,1);
    normSpectrum=0;
    normffd=0;
    for k=1:TIMES
        matrix_energy=matrix_energy+trial(k).abs_eigvals;
        normSpectrum=normSpectrum+trial(k).normSpectrum;
        normffd=normffd+trial(k).normffd;
     end
    
    matrix_energy=matrix_energy/TIMES;
    average_norm_spectrum=normSpectrum/TIMES;
    average_norm_ffd=normffd/TIMES;
    
    
    
    energies(:, column_index)=matrix_energy;
    norm_spectrum(column_index)=average_norm_spectrum;
    norm_ffd(column_index)=average_norm_ffd;
   
    
end


%PLOTS

figure;
a1=plot(energies(:,1)); M1='ratio=3';
hold on;
a2=plot(energies(:,2)); M2='ratio=5';
hold on;
a3=plot(energies(:,3)); M3='ratio=10';
hold on;
a4=plot(energies(:,4)); M4='ratio=20';
hold on;
a5=plot(energies(:,5)); M5='ratio=30';
hold on;
a6=plot(energies(:,6)); M6='ratio=40';
legend([a1,a2,a3,a4,a5,a6],M1,M2,M3,M4,M5,M6)
xlabel('condition')
ylabel('max norm')
set(gca, 'TickDir', 'out');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
box off


figure;
a1=plot(norm_spectrum); M1='spectrum norm';
hold on;
a2=plot(norm_ffd); M2='ffd norm';
legend([a1,a2],M1,M2)
xlabel('inh strength')
ylabel('Frobenius norm')
set(gca, 'TickDir', 'out');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
box off


    