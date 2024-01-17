function [s,nc]=IAAFT(t,num_surr,maxiter)
%IAAFT generates # of surrogates (num_surr) from template time series (t)
%credit to algorithm by victor.venema@uni-bonn.de available at website:
%https://uni-bonn.viven.org//themes/surrogates/iaaft/iaaft_algorithm.html
%INPUTS
%   t: template data (vector)
%   num_surr: # of surrogates to be generated from template (value)
%   maxiter: stop after this many iterations whether or not converged (value)
%OUTPUTS
%   s: surrogate timeserie(s) generated (vector or matrix)
%   nc: number of surrogates where maxiter is reached (value)

%shuffle template to generate initial surrogate data
s = zeros(length(t), num_surr);
for col=1:num_surr
    s(:, col) = t(randperm(length(t)));
end
mag=abs(fft(t)); [t,i]=sort(t); iter = zeros(1, num_surr);
for n=1:num_surr
    c = 1;
    old = i;
    converge = 0;
    while c<=maxiter && converge == 0 
        % Phase-randomize surrogates
        phase=angle(fft(s(:,n))); s(:,n)=mag.*exp(phase.*1i);    
        s(:,n)=real(ifft(s(:,n)));      
        % Match distrubition of surrogates to template
        [~,d]=sort(s(:,n)); [~,new]=sort(d); s(:,n)=t(new);                                        
        % Check if converged
        if new==old
            converge=1;
        else
            old=new;
            c=c+1;
        end
    end
    %record iteration count
    iter(n)=c;
end
nc = sum(sum(iter>=maxiter));
end