function [I V]=whichvalue(x,v)

if(length(v)==1)
    [V I]=min(abs(x-v));
else
    V=zeros(1,length(v));
    I=zeros(1,length(v));
    for i=1:length(v)
        [V(i) I(i)]=min(abs(x-v(i)));
    end
end

