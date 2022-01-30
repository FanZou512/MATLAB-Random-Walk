function [ varargout ] = RandWalkSim( varargin )
% This function simulate the 3D random walk of a given dot
% Format of use: [ ] = RandWalkSim( N, timed, times, x0 )
% Discription for variables:
%   N:      number of traces
%   timed:	time duration, unit s
%   times:  time resolution, unit s
%   x0:     initial position
%   bc:     boundary condition, bc.type code 1/2/3 are square/sphere/mesh
%           bc.xbc diameter/3D box size

%% set input values and the defaults
switch class(varargin{1})
    case {'char','path'}
        load(varargin{1})
    case {'double','int','long'}
        threshold=0.2;shifts=0;diff_const=0.0012;
        calOP=[0,0,1,1];    % calculation options: [2D MSD, 3D MSD, FPT, cumulative FPT]
        plotOP=[0,0,1,1,1];	% plot options: [trace, 2D MSD, 3D MSD, FPT, cumulative FPT]
        optOP=[0,0,1,1,1];  % output options: [traces, 2D MSD, 3D MSD, FPT, cumulative FPT]
        FPTtype=1;
    switch nargin
    case 0
        %N=10;timed=900;times=10;x0=[];
        %bc.type=2;bc.xbc=[0.5,2,2];
    case 1
        N=varargin{1};timed=900;times=1;x0=[];
        bc.type=2;bc.xbc=[0.5,2,2];
    case 2
        N=varargin{1};timed=varargin{2};times=1;x0=[];
        bc.type=1;bc.xbc=[2,2,2];
    case 3
        N=varargin{1};timed=varargin{2};times=varargin{3};x0=[];
        bc.type=1;bc.xbc=[2,2,2];
    case 4
        N=varargin{1};timed=varargin{2};times=varargin{3};x0=varargin{4};
        bc.type=1;bc.xbc=[2,2,2];
    case 5
        N=varargin{1};timed=varargin{2};times=varargin{3};x0=varargin{4};bc=varargin{5};
    end
end

t=0:times:timed;
tl=length(t)-1;
RW.time=t(2:end);

%% generate initial position
n=size(x0,2);
if n==0
    x0=randn(3,N)*min(bc.xbc)-min(bc.xbc/2);
    for k=1:N
        x0r=sum(x0(:,k).^2);
        while x0r>min(bc.xbc/2)^2
            x0(:,k)=randn(3,1)*min(bc.xbc)-min(bc.xbc/2);
            x0r=sum(x0(:,k).^2);
        end
    end
elseif n<N
    x0=x0(:,1)*ones(1,N);
end
%plot3(x0(1,:),x0(2,:),x0(3,:),'*');% check for the initial positions

%% generate random walk steps
mu=0;sig=sqrt(diff_const*times);	% mu and sigma of the random walk steps
% Diffusion constant from Sussan is about 1.2e-3 um^2/s
x1= mu+sig*randn(tl,N);
x2= mu+sig*randn(tl,N);
x3= mu+sig*randn(tl,N);

% generate final random walk coordinates
x1=[x0(1,:);x1];x2=[x0(2,:);x2];x3=[x0(3,:);x3];
x1=cumsum(x1);x2=cumsum(x2);x3=cumsum(x3);
tl=tl+1;            % length of data including the initial point

% setup boundary conditions
switch bc.type
    case 1          % box condition
        [x1,x2,x3]=BCbox(bc.xbc,x1,x2,x3);
    case 2          % sphere condition
        [x1,x2,x3]=BCsphere(bc.xbc(1)/2,x1,x2,x3);
    case 3          % irregular condition defined by mesh girds
end

%% calculate MSD
if calOP(1)
    disp('Start working on MSD - 2D projection.')
    RW.MSD2D=MnSqDis(x1,x2);	% projection to 2D plane
end
if calOP(2)
    disp('Start working on MSD - 3D.')
    RW.MSD3D=MnSqDis(x1,x2,x3);	% projection to 2D plane
end

%% calculate FPT use pairwise traces
if calOP(3)
    disp('Start working on FPT.')
    [FPT,nFPT]=FPTModified(x1,x2,x3,threshold,FPTtype,shifts);
    FPT=FPT*times;FPT=FPT';
end
if calOP(4)
    RW.FPThist=hist(FPT(:),t(2:end));
    FPTcum=cumsum(RW.FPThist);
end
%% plot results
if plotOP(1)    % plot trace
    n=min(N,5);
    color=rand(n,3)*0.8+0.1;
    figure
    for k=1:n
        plot(x1(:,k),x2(:,k),'color',color(k,:));
        %plot3(x1(:,k),x2(:,k),x3(:,k),'color',color(k,:));
        hold on
    end
    title(['duration ',num2str(timed),'s, time interval ',num2str(times),'s'],'FontSize',16)
    hold off
end
if plotOP(2)    % plot MSD-2D
    figure,plot(t,RW.MSD2D,'k')
    xlabel('time/s'),ylabel('2D MSD um^2/s')
end
if plotOP(3)    % plot MSD-3D
    figure,plot(t,RW.MSD3D,'b')
    xlabel('time/s'),ylabel('3D MSD um^2/s')
end
if plotOP(4)    % plot FPT histogram
    FPTx=FPT(:);FPTx(FPTx<0)=[];
    figure,hist(FPTx(:),30)
end
if plotOP(5)    % plot FPT cumulative curve
    figure,subplot(1,2,1),plot(t(2:end),FPTcum)
    xlabel('time/s'),ylabel('cumulative number of encounters')
    subplot(1,2,2),plot(t(2:end),FPTcum./nFPT)
    xlabel('time/s'),ylabel('percentage of cumulative encounters')
end

%% output traces
RW.FPT=FPT;
RW.FPTcum=FPTcum;
RW.FPTcp=FPTcum./nFPT;
varargout{1}=RW;%varargout{2}=MSD;varargout{3}=x1;varargout{4}=x2;varargout{5}=x3;
end

function [FPT, nFPT]=FPTModified(x1t,x2t,x3t,threshold,option,vargargin)
% This function calculate the FPT based on the coordinates of 3D trajectories
% Notice that this function use all the traces in a pairwise way
% If two traces do not meet, output NaN. Invalid pairs output -1.
% *********** choices ***********
% option 1: regular FPT as default with 3D
% option 2: regular FPT but projection at 2D
% option 3: shifting target trace, here we shift along x-axis, but it can extend
[~,N]=size(x1t);

switch option
    case 2 % FPT of 2D projection
        FPT1=FPT2D(x1t,x2t,threshold);
        FPT2=FPT2D(x1t,x3t,threshold);
        FPT3=FPT2D(x2t,x3t,threshold);
        FPT=[FPT1;FPT2;FPT3];
        nFPT=3*N*(N-1)/2;
    case 3 % FPT with shifts, FPT is an asymetric matrix
        shifts=vargargin;	% shift target trace along x-axis
        angle=randi(90,[1,2]);angle=angle*pi/180; % theta,phi both [0,pi/2] to ensure move in positive x-y-z direction
        shifts=shifts*[cos(angle(1))*sin(angle(2)),sin(angle(1))*sin(angle(2)),cos(angle(2))];
        FPT=NaN(N,N);
        for ko=1:N
            for kt=1:N
                if ko==kt
                    %FPT(ko,kt)=-1;
                    continue;
                end
                DisDot=sqrt((x1t(:,ko)-x1t(:,kt)-shifts(1)).^2+(x2t(:,ko)-x2t(:,kt)-shifts(2)).^2+(x3t(:,ko)-x3t(:,kt)-shifts(3)).^2);
                a=find(DisDot<threshold,1,'first');
                if ~isempty(a)
                    FPT(ko,kt)=a(1);
                end
            end
        end
        nFPT=N*(N-1);
    otherwise % regular FPT in 3D
        FPT=NaN(1,(N-1)*N/2);
        k=0;
        for ko=1:N-1        % original trace
            for kt=ko+1:N   % tartget trace
                k=k+1;
                DisDot=sqrt((x1t(:,ko)-x1t(:,kt)).^2+(x2t(:,ko)-x2t(:,kt)).^2+(x3t(:,ko)-x3t(:,kt)).^2);
                a=find(DisDot<threshold,1,'first');
                if ~isempty(a)
                    FPT(k)=a(1);
                end
            end
        end
        nFPT=N*(N-1)/2;
end
end

function FPT=FPT2D(x1t,x2t,threshold)
% This function calculate the FPT based on the coordinates of 2D trajectories
% Notice that this function use all the traces in a pairwise way
% If two traces do not meet, output NaN. Invalid pairs output -1.
[~,N]=size(x1t);

FPT=NaN(1,(N-1)*N/2);
k=0;
for ko=1:N-1        % original trace
    for kt=ko+1:N   % tartget trace
        k=k+1;
        DisDot=sqrt((x1t(:,ko)-x1t(:,kt)).^2+(x2t(:,ko)-x2t(:,kt)).^2);
        a=find(DisDot<threshold,1,'first');
        if ~isempty(a)
            FPT(k)=a(1);
        end
    end
end
end

function [ MSD ]=MnSqDisModify(varargin)
% This function calculate the MSD of given data
% input should be x,y,z, ... steps (differences) in matrix form, with same size
% each input is one dimension
% ------------ How to use ------------
% [MSD] = MnSqDisModified( xt,yt,... );
[m,n]=size(varargin{1});
MSD=zeros(m,1);

for k=1:m
    xk=zeros(m-k+1,n);
    for p=1:nargin
        xk=xk+(conv2(varargin{p},ones(k,1),'valid')).^2;
    end
    MSD(k)=nanmean(xk(:));
    if mod(sqrt(k),1)==0
        disp([num2str(k),' intervals calculated!'])
    end
end
MSD=[0;MSD];
end

function [varargout]=BCsphere(R,varargin)
% This function check if trace hits boundary, if does then reflect trace
% R is the radius of the sphere
% varargin contains coordinates of every dimension separately, [x, y, ...]
disp('Start working on BCs.')
varargout=varargin;
Rsq=R^2;
[m,n]=size(varargin{1});
r=zeros(m,n);
dim=nargin-1;
for k=1:dim
    r=r+varargout{k}.^2;
end

for kl=1:n
    a=1;
    while 1     % repeat reflection
        ax=find(r(a:m,kl)>Rsq, 1, 'first');
        if isempty(ax)
            break
        else
            a=ax(1)+a-1;
            % find the intersection
            p1=zeros(1,dim);p2=p1;
            for k=1:dim
                p1(k)=varargout{k}(a-1,kl);
                p2(k)=varargout{k}(a,kl);
            end
            p1p2=p2-p1;
            au=nansum(p1p2.^2);
            bu=nansum(p1.*p1p2);    % b=2*op1.*p1p2 in doc, but 2 cancelled
            cu=Rsq-nansum(p1.^2);   % c=op1^2-R^2, here negative cancelled
            u=(sqrt(bu.^2+au*cu)-nansum(bu))/au;
            %if u<0
                %disp('u<0')
            %end
            ps=p1+u.*p1p2;          % intersection point
            shifts=-2*nansum((p2-ps).*ps).*ps/nansum(ps.^2);
            rt=zeros(m-a+1,1);
            for k=1:dim
                varargout{k}(a:m,kl)=varargout{k}(a:m,kl)+shifts(k);
                rt=rt+varargout{k}(a:m,kl).^2;
            end
            r(a:m,kl)=rt;
        end
    end
    if mod(kl,50)==0
        disp([num2str(kl),' traces checked!'])
    end
end
end

function [varargout]=BCbox(xbc,varargin)
% This function check if trace hits boundary, if does then reflect trace
% xbc is the max/min x, i.e. boundary, format [x_length, y_length,...]
% varargin contains coordinates of every dimension separately, [x, y, ...]
[m,n]=size(varargin{1});
varargout=varargin;
for k=1:nargin-1
    x=varargin{k};	% coordinates of the dot
    xbck=xbc(k)/2;  % boundary of this dimension
    for p=1:n
        while 1     % repeat reflection
            a=find(x(:,p)>xbck | x(:,p)<-xbck, 1, 'first');
            if isempty(a)
                break
            else
                a=a(1);
                shifts=sign(x(a,p))*2*xbck-2*x(a,p);
                x(a:m,p)=x(a:m,p)+shifts;
            end
        end
    end
    varargout{k}=x;
end
end