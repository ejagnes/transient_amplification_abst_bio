function [eigenvectors, abs_eigvals, Q] = MaximiseIC(W, flag)



if flag
    
    Q = lyap( (W-eye(size(W)))', 2*eye(size(W)));
    
else
    
    [ev, D] = eig(W);
    
    [~,I] = sort(abs(imag(diag((D)))),'ascend');
    
    C = [(real(ev(:,I(end))))'; (imag(ev(:,I(end))))'];
    
    Q = lyap((W-eye(size(W)))', 2*(C')*C);
    
end

% Eigenvalue decomposition of Q gives initial conditions that maximise
% energy
[eigenvectors, eigenvalues] = eig(Q);

abs_eigvals = abs(diag(eigenvalues));

% Sort eigenvalues by magnitude
[~,index] = sort(abs_eigvals, 'descend');

eigenvectors = eigenvectors(:,index);

abs_eigvals = abs_eigvals(index);

%figure;
%plot([1:200], abs_eigvals) 