%separate eigenvalues from ffd in real Schur and then write spectrum with
%real and imaginary parts on the blocks
function T2=upper_triang_Schur(M)



[U,T]=schur(M);
T2=T;  

[vecs,vals]=eig(M);

[~,D2] = cdf2rdf(vecs,vals);

[rId, cId] = find(D2);

indices=[rId,cId];

for i=1:length(indices)
    T2(indices(i,1), indices(i,2))=0;
end


M_new=D2+T2;



