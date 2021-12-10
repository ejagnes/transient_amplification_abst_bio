% script taking a matrix M and putting random connections to the upper triangular
% part, with norm specified by str

function M_new=ffd_random(M, str)


N=length(M); %fetch size of matrix


A=-1+2*rand(N); %create a random matrix to use for the new upper triangular part
[U,T]=schur(M);  %now do the Schur transform of M

[vecs, vals]=eig(M);  %find the eigenvalues 
vals=diag(vals);
re_vals=real(vals);
im_vals=imag(vals);
I1=find(im_vals)';  %where the imaginary eigenvalues are
I2=find(im_vals==0)'; %where the real eigenvalues are

new_imaginary=(nnz(im_vals))/2;  %find how many eigenvalues are complex, and hence need
                                 %to occupy a 2x2 block on the diagonal.

%We need to separate the diagonal from the feedforward part, in order to assign the
%strength only at the feedforward part of the matrix and not at the spectrum. 
%We need to be careful because this is the real Schur transform and 
%we have 2 x 2 blocks along the diagonal when an eigenvalue is complex so
%we need to make sure this is associated to the spectrum, rather than the
%ffd part.

T2=triu(T,2);  %initial separation of the stictly upper triangular

T1=T-T2;  %T1 is the part of the upper triangular part that might belong to the
          %spectrum 

T2=triu(A,2); %now we assign to T2 the random values from A.

for row=I1
    nonzero_index=0;
    keep_col=[];
    for col=1:N
        if T1(row,col)~=0
            nonzero_index=nonzero_index+1;
            keep_col=[keep_col,col];
        end
    
    
    if nonzero_index >2
        
        T2(row, keep_col(3))=T1(row, keep_col(3));
        T1(row, keep_col(3))=0;
        break;
    end
    end
end


for row=I2
    nonzero_index=0;
    keep_col=[];
    for col=1:N
        if T1(row,col)~=0
            nonzero_index=nonzero_index+1;
            keep_col=[keep_col,col];
        end
    
    
       if nonzero_index ==2
        T2(row, keep_col(2))=T1(row, keep_col(2));
        T1(row, keep_col(2))=0;
        break;
       end
    end
end

    
ffd_str=norm(T2, 'fro'); %current norm of the ffd part
new_strength=str;  %the norm we want to assign


S_aux=(new_strength/ffd_str)*T2;

D=diag(vals);
[~,D2] = cdf2rdf(vecs,D);
M_new=D2+S_aux;


%now attach the new strictly upper triangular part:



%M_new=U*M_new*U';  %return back to the initial basis
