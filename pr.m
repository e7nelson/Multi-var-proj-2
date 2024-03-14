function y=pr(sys,w)
  
  y=zeros(size(w));
  for i=1:length(w)
    F=evalfr(sys,j*w(i));
    y(i)=.5*min(real(eig(F+F')));
  end
  
