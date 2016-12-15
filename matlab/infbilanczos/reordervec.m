function [Itest,lvv]=reordervec(ref,  lv)
%  reorder lv such according to v0 and ref
ii=1:length(lv);
I=[];
vs=lv;
for i=1:length(ref)
    if (~isnan(ref(i)))
        [Y,I0]=min(abs(ref(i)-vs));
%         disp([Y,I0])
%         keyboard
        if isnan(Y)
            I(end+1)=NaN;
            lvv(i)=NaN;
        else
            ii(I0)=NaN;
            vs(I0)=NaN;
            I(end+1)=I0;
            lvv(i)=lv(I0);
        end
    else
        I(end+1)=NaN;
        lvv(i)=NaN;
    end
end
Itest = I;
% Itest = I(find(isnan(ii))); 
% Itest = ii;
% Itest = Itest(find(~isnan(Itest))); 
% I=[I,ii(~isnan(find(ii)))];
