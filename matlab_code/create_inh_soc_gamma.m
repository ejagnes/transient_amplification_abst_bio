%make soc but keep overall strength of inhibition gamma times stronger than
%excitation

function [Wsoc, e] = create_inh_soc_gamma(W, rate, threshold,gamma)

%Inputs: W (initial matrix), rate (rate of descent - see Guillaume's
%paper), threshold (point at which code stops - see below).

%Outputs: Wsoc (SOC matrix after optimisation), e (largest real eigenvalue
%at each iteration), Hist (history of the matrix over optimisation)

Wsoc = W;
% figure
% subplot(2,1,1)
% plot(eig(Wsoc),'*');
% xlim([-5 5]);
% ylim([-5 5]);
% xlabel('real part');
% ylabel('imaginary part');

i = 2;
e(1) = 1000;
e(2) = 100;
% Hist = zeros(length(W)^2,1000);
% Hist(:,1) = W(:);
while e(i) > threshold
    
    i = i+1 %Output counter
    
    [Wsoc] = ssaCode(1.5,Wsoc,rate, gamma);
%     Hist(:,i-1) = Wsoc(:);
    e(i) = max(real(eig(Wsoc)));
    
    e(i) %Output max eigenvalue
    
%     if mod(i,20) == 0 %Plot every 20 iterations
%         
%         pause(0.01)
%         
%         subplot(2,1,1)
%         plot(eig(Wsoc),'*');
%         xlim([-5 5]);
%         ylim([-5 5]);
%         xlabel('real part');
%         ylabel('imaginary part');
%         
%         subplot(2,1,2)
%         plot(e(3:end))
%         
%     end
    
    
    
end


end

function [Wo, Emax] = ssaCode(C,Wi,rate, gamma)

Wr = length(Wi);

end_exc = Wr/2;
start_inh = Wr/2 + 1;

Emax = max(real(eig(Wi)));

if Emax < 0 %Stop if max eigenvalue is below 0
    
    Wo = Wi;
    
else
    
    s = max(Emax*C, Emax + 0.2);
    
    A = Wi - s*eye(size(Wi));
    
    X = 2*eye(size(Wi));
    
    Q = lyap(A',X); %Solve the lyapunov equations
    
    P = lyap(A,X);
    
    temp = Q*P;
    grad = temp/(trace(temp)); %calculate gradient
    
    Wo = Wi;
    Wo(:,start_inh:end) = Wi(:,start_inh:end) - (rate * grad(:,start_inh:end));
    
    %Now perform all constraints
    I = Wo > 0; %Clip positive I connections
    I(:,1:end_exc) = 0;
    Wo(I) = 0;
    
    %Make all inhibitory weights on average 3 times stronger than the
    %excitatory ones.
    
    exc_to_exc = Wo(1:end_exc,1:end_exc);
    meanEE = mean(exc_to_exc(:));
    
    exc_to_inh = Wo(start_inh:end,1:end_exc);
    meanEI = mean(exc_to_inh(:));
    
    inh_to_exc = Wo(1:end_exc,start_inh:end);
    meanIE = mean(inh_to_exc(:));
    
    inh_to_inh = Wo(start_inh:end,start_inh:end);
    meanII = mean(inh_to_inh(:));
    
    Wo(1:end_exc,start_inh:end) = -gamma*(meanEE/meanIE) * Wo(1:end_exc,start_inh:end);
    Wo(start_inh:end,start_inh:end) = -gamma*(meanEI/meanII) * Wo(start_inh:end,start_inh:end);
    
    %Make only 40% of the inh weights non-zero by clipping the weakest to
    %zero
    inh_list = Wo(:,start_inh:end);
    inh_list = inh_list(:);
    [~,I] = sort(inh_list,'ascend');
    thres = round(0.4*length(inh_list));
    inh_list(I(thres:end)) = 0;
    Wo(:,start_inh:end) = reshape(inh_list,end_exc*2,end_exc);
    
    %d=diag(Wo);
    %d(1:150)=zeros(150,1);
    %Wo=Wo-diag(d);
    
    Wo = Wo - diag(diag(Wo));
    
end

end