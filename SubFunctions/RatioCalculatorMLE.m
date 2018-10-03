function [X, Y, Z, Int_, BG,ZEst,Start]=RatioCalculatorMLE(imagedata,rx ,ry)
global pxSize DAT I

%%

N_Images = max(DAT(:,1));

Start=[];
Z=[];
X=[];
Y=[];
Int_=[];
iii=[];
BG=[];
ZEst=[];

for j = 1:25%N_Images-1
    j
    %% Bild aufteilen
    data = imagedata(:,:,j);
    [dx dy dz] =size(data);

    
    [Dat_TOT,Dat_UAF]=DataHandlerNew(j);

%     for k=1:length(Dat_UAF)
    if Dat_TOT(j)~=0
        %%
        clear Verh_P

        l=Dat_TOT(j);
        if l ==0
            Verh_P = []; 
        else
            for i=1:l
                
                if j==1
                    posX = ceil(DAT(i,2)/pxSize+.5);
                    posY = ceil(DAT(i,3)/pxSize+.5);
                    Int = DAT(i,5);             
                else
                    posX = ceil(DAT(Dat_UAF(j-1)+i,2)/pxSize+.5);
                    posY = ceil(DAT(Dat_UAF(j-1)+i,3)/pxSize+.5);
                    Int = DAT(Dat_UAF(j-1)+i,5);
                end

                % condition if fluoropühore lies to close to the edges
                if(posY > 5 & posX > 5 & posY < dx-5 & posX < dy-5 & Int > 50000 & Int <100000)

                    % estimate z
                    R=sqrt((posX*115-rx*115).^2+(posY*115-ry*115).^2);
                    ZestAnal=-sqrt(-(R)^2+(1e6)^2)+1e6;

                    I=data(posY-5:posY+4,posX-5:posX+4);

                    [x,y,z,ii,bg_,ZEst_] = MLE_Estimator_TOT(double(I),50 ,Int);%++ZestAnal/5
                    
                    pause(.01)
                    figure(1)
                    imagesc(I)
                    title(ZEst_)
                    
                    Verh_P(i)=z;
                    ZEst_
                    Int_=[Int_; ii];
                    ZEst = [ZEst; ZEst_];
                    Start =[Start; ZestAnal/5];
                    X = [X; DAT(Dat_UAF(j)+i-1,2)];
                    Y = [Y; DAT(Dat_UAF(j)+i-1,3)];
                else
                    Int=[Int; 0];

                    if i==1
                        Verh_P = [];  
                    end
                    Verh_P = [Verh_P ; []]; 

                end

            end
        end

        Z = [Z; Verh_P'];
%     end
    end
%     I = [I; iii'];
%     BG = [BG; bg];
end

 
% % % %%
% % % figure(3)
% % % X_=Dat_TOT(:,2);
% % % Y_=Dat_TOT(:,3);
% % % scatter(X_,Y_,'o');
% % % hold on
% % % X_=Dat_UAF(:,2);
% % % Y_=Dat_UAF(:,3);
% % % scatter(X_,Y_,'x')
% % % axis equal;
% % % hold on

end