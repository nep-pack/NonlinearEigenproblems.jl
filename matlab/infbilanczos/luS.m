function y = luS(x,L,U,p,q2,R)

y = U\(L\(R(:,p)\x)); 
y = y(q2); 