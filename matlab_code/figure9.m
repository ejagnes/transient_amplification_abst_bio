
load('standard_20_stash.mat')

scales=[0.4,0.6,0.8,1,1.2,1.4,1.6,1.8];
TIMES=20;
norm_spectrum=zeros(1,length(scales)); 
norm_ffd=zeros(1,length(scales));
energies=zeros(200,length(scales));
imaginary_diameter=zeros(1,length(scales));



column_index=0;
for kappa=scales
    column_index=column_index+1;
    
    for times=1:TIMES
        Wsoc=stash(times).Wsoc;
        
        new_Wsoc=Wsoc;
        new_Wsoc(101:end,1:100)=kappa*new_Wsoc(101:end,1:100);
        new_Wsoc=create_inh_soc_gamma(new_Wsoc, 5, 0.5,3);
new_Wsoc=100/norm(new_Wsoc, 'fro') * new_Wsoc; %normalise norm

        
        
        



[trial(times).max_norm, percentage_amplified]=max_norm_analytical(new_Wsoc);
%[~, trial(i).abs_eigvals, ~] = MaximiseIC(new_EI(i).Wsoc, 1);
    
trial(times).max_norm=trial(times).max_norm';

   [~,vals]=eig(new_Wsoc);
    vals=diag(vals);
    vals_imag=imag(vals);
    B=sort(vals_imag);
    im_diameter(times)=B(end)-B(1); 
    
    
    
     size_matrix=length(Wsoc);
 
 %finding the ffd and diagonal components of Schur decomposition
 
[U,T]=schur(new_Wsoc);
T2=T;  

[vecs,vals]=eig(new_Wsoc);

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
    imaginary=0;
    for k=1:TIMES
        matrix_energy=matrix_energy+trial(k).max_norm;
        normSpectrum=normSpectrum+trial(k).normSpectrum;
        normffd=normffd+trial(k).normffd;
        imaginary=imaginary+im_diameter(k);
     end
    
    matrix_energy=matrix_energy/TIMES;
    average_norm_spectrum=normSpectrum/TIMES;
    average_norm_ffd=normffd/TIMES;
    average_diameter=imaginary/TIMES;
    
    
    
    energies(:, column_index)=matrix_energy;
    norm_spectrum(column_index)=average_norm_spectrum;
    norm_ffd(column_index)=average_norm_ffd;
    imaginary_diameter(column_index)=average_diameter;
    
   
    
end








%% 


figure;
set(gcf, 'Position',  [100, 100, 500, 1200])
subplot(3,1,3);
a1=plot(energies(:,1)); M1='0.4';
hold on;
a2=plot(energies(:,2)); M2='0.6';
hold on;
a3=plot(energies(:,3)); M3='0.8';
hold on;
a4=plot(energies(:,4)); M4='1';
hold on;
a5=plot(energies(:,5)); M5='1.2';
hold on;
a6=plot(energies(:,6)); M6='1.4';
hold on;
a7=plot(energies(:,7)); M7='1.6';
hold on;
a8=plot(energies(:,8)); M8='1.8';
legend([a1,a2,a3,a4,a5,a6,a7,a8],M1,M2,M3,M4,M5,M6,M7,M8)
xlabel('condition (%)')
ylabel('max norm')
set(gca, 'TickDir', 'out');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
box off
xlim([0,200])
xticks([0 100 200])
xticklabels({'0','50','100'})

subplot(3,1,2);
a1=plot(norm_spectrum, '-o'); M1='spectrum';
hold on;
a2=plot(norm_ffd, '-o'); M2='feedforward';
legend([a1,a2],M1,M2)
xlabel('IE/EE')
ylabel('Frobenius norm')
set(gca, 'TickDir', 'out');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
box off
xlim([0,8])
xticks([0 2 4 6 8])
xticklabels({'0.2','0.6','1','1.4','1.8'})

subplot(3,1,1);
plot(imaginary_diameter, '-o')
xlabel('IE/EE')
ylabel('imaginary diameter')
set(gca, 'TickDir', 'out');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3);
set(gca,'fontsize', 20);
box off
xlim([0,8])
xticks([0 2 4 6 8])
xticklabels({'0.2','0.6','1','1.4','1.8'})















