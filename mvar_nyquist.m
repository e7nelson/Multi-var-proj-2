%
% multivariable nyquist plot
%

function lambdaG=mvar_nyquist(G,w)

for i=1:length(w)
    lambdaG(:,i)=eig(evalfr(G,j*w(i)));
end

end