%function takes as input singular values and returns effective rank
%according to the definition given in this paper: THE EFFECTIVE RANK:
%A MEASURE OF EFFECTIVE DIMENSIONALITY , Olivier Roy and Martin Vetterli?

function ef_rank=effective_rank(Sing_vals)

N=length(Sing_vals);


one_norm=sum(abs(Sing_vals));

for k=1:N
    p(k)=Sing_vals(k)/one_norm;
    
end

log_matrix=zeros(1, length(p));

for k=1:N
    if p(k)==0
        log_matrix(k)=0;
    else
        log_matrix(k)=log(p(k));
    end
end



entropy=(-1)* sum(p.*log_matrix);

ef_rank=exp(entropy);

