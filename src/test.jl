#x=[-2 1 0.1 2]
#x=[.1 .2 -.1 .5]
x=0.1
m=4
yy=cos.((0:m)'.*acos(x));

y=function TT(x,m)
    y=zeros(1,m)
    II=(0:m)'
    if abs(x)<1
        y=cos.(II.*acos(x));
    elseif x>1
        y=cosh.(II.*acosh(x));
    else
        y=cosh.((II.^(-1)).*acosh(x));
    end
end


y=TT(x,m)
