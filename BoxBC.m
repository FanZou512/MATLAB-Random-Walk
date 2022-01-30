a=2.5;
x=(-a/2:0.01:a/2)';y=[];
sig=0.4:0.2:1.6;
my=zeros(1,length(sig));
for p=1:length(sig)
%sig=1.6;
sig2=sig(p)^2;
yt=zeros(size(x));
for k=-5:5
    yt=yt+exp(-(x-k*a).^2/sig2);
end
yt=yt/sum(yt);
y=[y,yt];
my(p)=sum(yt.*x.^2);
end
figure(1),hold on,plot([0,sig.^2],[0,my],'c','linewidth',2)
figure,plot(x,y,'linewidth',2)
axis([-1,1,0,max(y(:))])
legend('\sigma=0.4','\sigma=0.6','\sigma=0.8','\sigma=1','\sigma=1.2','\sigma=1.4','\sigma=1.6')
