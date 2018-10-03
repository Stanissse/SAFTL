function [V_out,BG,V] = EverBG(V)

[x,y,sz]=size(V);

Vf=mean(mean(V));
Vmin=min(Vf);
VvF=Vf/Vmin;
Vr=mean(Vf/Vmin);

for l=1:x 
   for m=1:y
   Vtempmin(l,m)=min(V(l,m,:));
   end
end

A=0.9433*exp(0.03336*Vr)-186100*exp(-15.89*Vr);
B=-600.6*exp(-4.131*Vr)-22.9*exp(-0.2063*Vr);

BG=Vf/Vmin.*(Vtempmin-B)/A;


for h=1:sz
V_out(:,:,h) = V(:,:,h)-BG(:,:,h);
end

end