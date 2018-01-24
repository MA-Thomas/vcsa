function [covfunc,hyp2] = set_Kernel(cov_iter,Dim)
   
        R = Dim;
        
        if cov_iter == 1
        covfunc = {@covMaternard,3};hyp2.cov = rand(Dim+1,1);hyp2.lik = log(0.1);
        elseif cov_iter == 2
        covfunc = {@covMaternard,5};hyp2.cov = rand(Dim+1,1);hyp2.lik = log(0.1);
        elseif cov_iter == 3
        covfunc = @covRQard; hyp2.cov = rand(Dim+2,1); hyp2.lik = log(0.1);
        elseif cov_iter == 4
        covfunc = @covRQiso; hyp2.cov = [0;0;0]; hyp2.lik = log(0.1);
        elseif cov_iter == 5
        covfunc = @covGaborard; hyp2.cov = rand(Dim+Dim,1); hyp2.lik = log(0.1);
        elseif cov_iter == 6
        covfunc = @covNNone; hyp2.cov = [0; 0]; hyp2.lik = log(0.1);
        elseif cov_iter == 7
        covfunc = @covSEard; hyp2.cov = rand(Dim+1,1); hyp2.lik = log(0.1);
        
        
        elseif cov_iter == 8
        covfunc = {'covADD',{1:R,'covNNone'} }; 
        hyp2.cov = [ log(ones(1,2*Dim)), log(ones(1,R))]; hyp2.lik = log(0.1);        
        elseif cov_iter == 9
        covfunc = { 'covADD',{1:R,'covSEiso'} };
        hyp2.cov = [ log(ones(1,2*Dim)), log(ones(1,R))]; hyp2.lik = log(0.1);
        end
        
%         D = Dim;
%         R = Dim;
%         if cov_iter == 1
% 
%         elseif cov_iter == 5
%         covfunc = { 'covADD',{1:R,'covGaboriso'} }; hyp2.cov = [ log(ones(1,2*D)), log(ones(1,R))]; hyp2.lik = log(0.1);
%         elseif cov_iter == 6
%         covfunc = { 'covADD',{1:R,'covNNone'} }; hyp2.cov = [ log(ones(1,2*D)), log(ones(1,R))]; hyp2.lik = log(0.1);
%         elseif cov_iter == 7
%         covfunc = { 'covADD',{1:R,'covSEiso'} }; hyp2.cov = [ log(ones(1,2*D)), log(ones(1,R))]; hyp2.lik = log(0.1);
%         end       
end
