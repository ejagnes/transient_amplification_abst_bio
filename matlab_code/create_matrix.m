%create matrix given imaginary variance r
%the real distribution can be either of the three from figure 2.A
%comment/uncomment according to which real distribution you would like to
%run

function M=create_matrix(r)

load('Wsoc.mat')
[U,S]=schur(Wsoc);
S_upper=upper_triang_Schur(Wsoc);
[vecs, vals]=eig(S);
vals=diag(vals);
re_vals=real(vals);
im_vals=imag(vals);

new_imaginary=(nnz(im_vals))/2;


vec=-r/2+r*rand(new_imaginary,1);

J=im_vals~=0;
I=im_vals==0;

im_vals_nonz=im_vals(J);
real_vals_nonz=re_vals(J);

im_vals_zero=im_vals(I);
re_vals_zero=re_vals(I);


k=1;
for i=1:2:length(im_vals_nonz)
   
        im_vals_nonz(i)=vec(k);
        im_vals_nonz(i+1)=-vec(k);
        k=k+1;
        
end




%% REAL PART OF SPECTRUM% %

%Choose only one of the three options below, comment the others

% uniform between -0.5 and 0.5 %
vec=-0.5+rand(new_imaginary,1); %for uniform between -0.5 and 0.5

k=1;
for i=1:2:length(real_vals_nonz)
   
        real_vals_nonz(i)=vec(k);
        real_vals_nonz(i+1)=vec(k);
        k=k+1;
        
end

new_vals_nonz=real_vals_nonz+1j*im_vals_nonz;
new_vals_rest=re_vals_zero;


new_vals_rest(:)=-0.5+rand(6,1);

new_vals=[new_vals_nonz; new_vals_rest];








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixed real at 0 %


% fixed_real=0;
% vec=fixed_real * ones(new_imaginary,1);
% 
% 
% 
% k=1;
% for i=1:2:length(real_vals_nonz)
%    
%         real_vals_nonz(i)=vec(k);
%         real_vals_nonz(i+1)=vec(k);
%         k=k+1;
%         
% end
% 
% new_vals_nonz=real_vals_nonz+1j*im_vals_nonz;
% new_vals_rest=re_vals_zero;
% 
% 
% new_vals_rest(:)=fixed_real;
% 
% new_vals=[new_vals_nonz; new_vals_rest];
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% OR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % a negative outlier at -20 %
% 
% new_vals_nonz=zeros(length(im_vals_nonz),1);
% dom=-20;
% new_vals_nonz(1:end-40)=fixed_real+1j*im_vals_nonz(1:end-40);
% new_vals_nonz(end-40+1:end)=0.5+1j*im_vals_nonz(end-40+1:end);
% new_vals_rest=re_vals_zero;
% new_vals_rest(:)=fixed_real;
% new_vals_rest(end)=-20; %global inhibitory balance
% new_vals=[new_vals_nonz; new_vals_rest];
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D=diag(new_vals);
[~,D2] = cdf2rdf(vecs,D);
M=D2+S_upper;


