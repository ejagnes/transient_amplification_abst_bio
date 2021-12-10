%script computing the distribution of the inner products of eigenvectors 
%also if the angle phi is more thn 90, takes the angle 180-phi

function [D,arg_diff]=eigvecs_overlaps(M)

[vecs, vals]=eig(M);

vals=diag(vals);
argument=imag(vals);

N=length(M);

%normalise eienvectors
for i=1:N
    vecs(:,i)=vecs(:,i)/norm(vecs(:,i));
end

%now compute overlaps:
k=0;
for i=1:N
    for j=i+1:N
        k=k+1;
        A(k)=real(dot(vecs(:,i), vecs(:,j)));
        D(k)=real(rad2deg(acos(A(k))));
        arg_diff(k)=abs((argument(i))-(argument(j)));
    end
end

for index=1:length(D)
    if D(index)>90
        D(index)=180-D(index);
    end
end


        