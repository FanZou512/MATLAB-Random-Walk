function varargout = RWsimu(varargin)
% RWSIMU MATLAB code for RWsimu.fig
%      RWSIMU, by itself, creates a new RWSIMU or raises the existing
%      singleton*.
%
%      H = RWSIMU returns the handle to a new RWSIMU or the handle to
%      the existing singleton*.
%
%      RWSIMU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RWSIMU.M with the given input arguments.
%
%      RWSIMU('Property','Value',...) creates a new RWSIMU or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RWsimu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RWsimu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RWsimu

% Last Modified by GUIDE v2.5 09-Oct-2019 13:32:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RWsimu_OpeningFcn, ...
                   'gui_OutputFcn',  @RWsimu_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before RWsimu is made visible.
function RWsimu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RWsimu (see VARARGIN)

% Choose default command line output for RWsimu
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

WaitBar(handles,0,'Not started yet ...')
axes(handles.axes3)
plot(0,0,'or'),hold on, axis off
axes(handles.axes1),axis off
global BCListContent
global BCs
global path
path=['C:\Users\fuz20\Documents\backup\Fan\matlab\Random walk\data',date,'.mat'];
BCListContent={};BCs={};
% UIWAIT makes RWsimu wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = RWsimu_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function SaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global path
[x1,x2] = uiputfile('*.mat','Save your data as ...');
path=[x2,x1];

% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RW
WaitBar(handles,0,'Opening project ...')
[x1,x2] = uigetfile('*.mat','Save your 4-D as ...');
path=[x2,x1];
load(path,'RW')
WaitBar(handles,0,'Project is open!')

% ------------------------refresh-------------------------------------
function WaitBar(handles,x,comments)
% This function updates the waitbar based on x [0,1]
%               updates the string based on comments
axes(handles.WaitBarFig)
barh(x)
axis([0,1,0.75,1.25])
axis off
set(handles.WaitBarTxt,'String',comments);

function listbox1_update(hObject, eventdata, handles)
global BCListContent
set(handles.BClist,'string',BCListContent);

% --- Executes on selection change in BCtype.
function BCtype_Callback(hObject, eventdata, handles)
% hObject    handle to BCtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Value')
    case 1 % sphere
        set(handles.text19,'String','BC radius (um):');
    case 2 % ellipsoid
        set(handles.text19,'String','BC axis radius (um):');
    case 3 % box
        set(handles.text19,'String','BC dimensions (um):');
end

% --- Executes on button press in AddBC.
function AddBC_Callback(hObject, eventdata, handles)
% hObject    handle to AddBC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BCs
global BCListContent
global BCline
nBCl=numel(BCline);
newBC.type= get(handles.BCtype,'Value');
axes(handles.axes3)
switch newBC.type
    case 1 % sphere
        x1=str2num(get(handles.BCsize,'String'));
        newBC.a=x1(1);
        x2=str2num(get(handles.BCcenter,'String'));
        newBC.c=x2(1:3);
        BCline(nBCl+1)=rectangle('Position',[x2(1:2)-newBC.a,2*newBC.a,2*newBC.a],'Curvature',[1,1]);
        info=['Sphere center at ',get(handles.BCcenter,'String'),' um, radius [',get(handles.BCsize,'String'),'] um'];
    case 2 % elliploid
        x1=str2num(get(handles.BCsize,'String'));
        newBC.a=x1(1:3);
        x2=str2num(get(handles.BCcenter,'String'));
        newBC.c=x2(1:3);
        % Add the ellipsoid projection on 2D
        t=-0.1:0.05:2*pi;
        x=newBC.a(1)*cos(t)+newBC.c(1);
        y=newBC.a(2)*sin(t)+newBC.c(2);
        BCline(nBCl+1)=plot(x,y,'-k');
        info=['Elliploid center at ',get(handles.BCcenter,'String'),' um, a/b/c radius [',get(handles.BCsize,'String'),'] um'];
    case 3 % box
        x1=str2num(get(handles.BCsize,'String'));
        if length(x1)<3 && ~isempty(x1)
            newBC.a=x1(1)*ones(1,3);
        else
            newBC.a=x1(1:3);
        end
        x2=str2num(get(handles.BCcenter,'String'));
        newBC.c=x2(1:3);
        BCline(nBCl+1)=rectangle('Position',[x2(1:2)-newBC.a(1:2)/2,newBC.a(1:2)]);
        info=['Box center at ',get(handles.BCcenter,'String'),' um, a*b*c dimensions [',get(handles.BCsize,'String'),'] um'];
end
axis equal
listbox1_update(hObject, eventdata, handles)
l=numel(BCs);
BCs{l+1}=newBC;
BCListContent{l+1}=info;
listbox1_update(hObject, eventdata, handles)

% --- Executes on button press in RemoveBC.
function RemoveBC_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveBC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BCs
global BCListContent
global BCline
x=get(handles.BClist,'value');
axes(handles.axes3)
delete(BCline(x))
BCline(x)=[];
BCs(x)=[];
BCListContent(x)=[];
listbox1_update(hObject, eventdata, handles)

% --- Executes on button press in Generate.
function Generate_Callback(hObject, eventdata, handles)
% hObject    handle to Generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RW
n0=numel(RW);
setup=GetSetup(handles);
n=numel(setup);
RW{n0+n}=[];
for k=1:n
    RW{n0+k}.setup=setup{k};
    [x,y,z,t]=GeneData(setup{k},handles,k);
    RW{n0+k}.x=x(2:end,:);RW{n0+k}.y=y(2:end,:);RW{n0+k}.z=z(2:end,:);RW{n0+k}.t=t;
end
%axes(handles.axes3),hold off

function [ setup ] = GetSetup(handles)
s.N=str2double(get(handles.N_cells,'String'));
s.Td=str2double(get(handles.TDuration,'String'));
s.Ts=str2double(get(handles.TStep,'String'));
s.D=str2num(get(handles.DiffConst,'String'));
s.thr=str2num(get(handles.ThrDis,'String'));
s.sig=str2num(get(handles.InitialPos,'String'));
lD=length(s.D);lthr=length(s.thr);ls=length(s.sig);
if ls>1
    for k=1:ls
        setup{k}=s;
        setup{k}.sig=s.sig(k);
    end
elseif lD>1
    for k=1:lD
        setup{k}=s;
        setup{k}.D=s.D(k);
    end
elseif lthr>1
    for k=1:lthr
        setup{k}=s;
        setup{k}.thr=s.thr(k);
    end
else % customized
    radius=0.35:0.075:1.2;
    for k=1:length(radius)
        setup{k}=s;
        setup{k}.r=radius(k);
    end
    %setup{1}=s;
end

function [x1,x2,x3,t]=GeneData(setup,handles,n_set)
% load data setup
N=setup.N;
timed=setup.Td;
times=setup.Ts;
% start here
t=0:times:timed;
tl=length(t)-1;
% initial pos
WaitBar(handles,0,'Generating initial position ...')
x0=normrnd(0,setup.sig,[3,N]);
WaitBar(handles,1,'Initial position done!')
% generate prelimitary random walk steps
mu=0;sig=sqrt(setup.D*times);	% mu and sigma of the random walk steps
WaitBar(handles,0,'Generating random walk steps ...')
x1= normrnd(mu,sig,[tl,N]);
WaitBar(handles,0.17,'Generating random walk steps ...')
x2= normrnd(mu,sig,[tl,N]);
WaitBar(handles,0.34,'Generating random walk steps ...')
x3= normrnd(mu,sig,[tl,N]);
WaitBar(handles,0.5,'Generating random walk steps ...')
% generate final random walk coordinates
x1=[zeros(1,N);x0(1,:);x1];x1=cumsum(x1);
WaitBar(handles,0.67,'Generating final random walk positions ...')
x2=[zeros(1,N);x0(2,:);x2];x2=cumsum(x2);
WaitBar(handles,0.84,'Generating final random walk positions ...')
x3=[zeros(1,N);x0(3,:);x3];x3=cumsum(x3);
WaitBar(handles,1,'Final random walk positions done!')
% check boundary conditions
WaitBar(handles,0,'Boundary conditions 0% checked!')
global BCs
for k=1:N
    BCs{1}.a=setup.r;
    x=CheckBCs([x1(:,k),x2(:,k),x3(:,k)]);
    x1(:,k)=x(:,1);x2(:,k)=x(:,2);x3(:,k)=x(:,3);
    WaitBar(handles,k/N,['Set No.',num2str(n_set),', boundary conditions ',num2str(k*100/N,3),'% checked!'])
end
WaitBar(handles,1,'Boundary conditions checking is done!')
axes(handles.axes3),hold on,plot(x1(2:end,min(N,5)),x2(2:end,min(N,5)))

function [ x ] =CheckBCs(x)
global BCs
nbc=numel(BCs);
st=1;
[m,~]=size(x);
while 1
    n=[];xs=[]; % where and how much to shift
    for k=1:nbc
        switch BCs{k}.type
            case 1
                [nk,xsk]=CheckSphere(BCs{k},x,st);
            case 2
                [nk,xsk]=CheckEllipsoid(BCs{k},x,st);
            case 3
                [nk,xsk]=CheckBox(BCs{k},x,st);
        end
        n=[n;nk];xs=[xs;xsk];
    end
    if isempty(n)
        break;
    else % find first point to reflect
        st=min(n);
        n_index=find(n==st);
        if length(n_index)>1 % if several BC applies, choice the largest shift
            [~,xs_index]=max(sum((xs(n_index,:)).^2,2));
            n_index=n_index(xs_index);
        end
        xs=xs(n_index,:);
        x(st:m,:)=x(st:m,:)+ones(m-st+1,1)*xs; % add shift, refresh trace
    end
end

function [n,xs]=CheckSphere(BC,x,st)
% find the first outlier after st
% BC - is a single sphere boundary
% x  - is the trace
% st - is a index where start to search for outlier
Rsq=BC.a^2;
[m,~]=size(x);
x=x-ones(m,1)*BC.c;	% make circle center the origin
r=sum((x).^2,2);
n=find(r(st:m)>Rsq, 1, 'first');
if isempty(n)
    xs=[];
else
    n=n(1)+st-1;
    % find the intersection
    p1=x(n-1,:);
    p2=x(n,:);
    p1p2=p2-p1;
    au=sum(p1p2.^2);
    bu=sum(p1.*p1p2);    % b=2*op1.*p1p2, here bu=b/2
    cu=Rsq-sum(p1.^2);   % c=op1^2-R^2, here cu=-c
    u=(sqrt(bu.^2+au*cu)-bu)/au; % (-b+sqrt(b^2-4*a*c))/2a
    ps=p1+u.*p1p2;       % intersection point
    xs=-2*nansum((p2-ps).*ps).*ps/nansum(ps.^2); % shift
end

function [n,xs]=CheckEllipsoid(BC,x,st)
% desscription same as sphere BC
% transform into unit circle center at origin
[m,~]=size(x);
x=(x-ones(m,1)*BC.c)./(ones(m,1)*BC.a);
r=sum((x).^2,2);
% find the reflection as circle
n=find(r(st:m)>1, 1, 'first');
if isempty(n)
    xs=[];
else
    n=n(1)+st-1;
    % find the intersection
    p1=x(n-1,:);
    p2=x(n,:);
    p1p2=p2-p1;
    au=sum(p1p2.^2);
    bu=sum(p1.*p1p2);   % b=2*op1.*p1p2, here bu=b/2
    cu=1-sum(p1.^2);  % c=op1^2-R^2, here cu=-c
    u=(sqrt(bu.^2+au*cu)-bu)/au; % (-b+sqrt(b^2-4*a*c))/2a
    ps=p1+u.*p1p2;      % intersection point
    xs=-2*nansum((p2-ps).*ps).*ps/nansum(ps.^2); % shifts in transform space
    xs=xs.*BC.a;        % reverse transform
end

function [n,xs]=CheckBox(BC,x,st)
% BC - is a single box boundary, others are the same as sphere BC
[m,~]=size(x);
x=x-ones(m,1)*BC.c;	% make circle center the origin
BC.a=BC.a/2;
sigx=abs(x(st:m,:))>(ones(m-st+1,1)*BC.a);
n=find(sum(sigx,2), 1, 'first');
if isempty(n)
    xs=[];
else
    sigx=sigx(n(1),:);
    n=n(1)+st-1;
    xs=2.*sigx.*(sign(x(n,:)).*BC.a-x(n,:));
end

% --- Executes on button press in MSD_2D.
function MSD_2D_Callback(hObject, eventdata, handles)
% hObject    handle to MSD_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RW
n=numel(RW);
for k=1:n
    RW{k}.MSD_2D=MSD_Modify(RW{k}.x,RW{k}.y);
end

% --- Executes on button press in MSD_3D.
function MSD_3D_Callback(hObject, eventdata, handles)
% hObject    handle to MSD_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RW
n=numel(RW);
for k=1:n
    RW{k}.MSD_3D=MSD_Modify(RW{k}.x,RW{k}.y,RW{k}.z);
end

function [ MSD ]=MSD_Modify(varargin)
% This function calculate the MSD of given data
% input should be x,y,z, ... steps (differences) in matrix form, with same size
[m,n]=size(varargin{1});
MSD=zeros(m,1);
for k=1:m
    xk=zeros(m-k+1,n);
    for p=1:nargin
        xk=xk+(conv2(varargin{p},ones(k,1),'valid')).^2;
    end
    MSD(k)=nanmean(xk(:));
    if mod(k,10)==0
        WaitBar(handles,k/m,['MSD ',num2str(k*100/m,3),'% steps done!'])
    end
end
MSD=[0;MSD];

% --- Executes on button press in FPT.
function FPT_Callback(hObject, eventdata, handles)
% hObject    handle to FPT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RW
global path
nRW=numel(RW);
option=get(handles.CalOption,'Value');
xc=str2num(get(handles.CenterDist,'String'));
switch option
    case 1
        for k=1:nRW
            N=floor(RW{k}.setup.N/2);
            xa=RW{k}.x(:,1:N);xb=RW{k}.x(:,N+1:2*N);
            ya=RW{k}.y(:,1:N);yb=RW{k}.y(:,N+1:2*N);
            za=RW{k}.z(:,1:N);zb=RW{k}.z(:,N+1:2*N);
            thr=RW{k}.setup.thr;
            t=RW{k}.t;
            [RW{k}.FPT_num,RW{k}.FPTcump,RW{k}.tau,RW{k}.FPTt]=CalFPT(xa,ya,za,xb,yb,zb,t,xc,thr);
        end
    case 2
        xa=RW{nRW-1}.x;xb=RW{nRW}.x;
        ya=RW{nRW-1}.y;yb=RW{nRW}.y;
        za=RW{nRW-1}.z;zb=RW{nRW}.z;
        thr=min(RW{nRW-1}.setup.thr,RW{nRW}.setup.thr);
        t=RW{nRW-1}.t;
        [RW{nRW}.FPT_num,RW{nRW}.FPTcump,RW{nRW}.tau,RW{nRW}.FPTt]=CalFPT(xa,ya,za,xb,yb,zb,t,xc,thr);
    case 3 % pair
        nRW=nRW/2;
        FPTcump=[];FPT=[];tau=zeros(nRW,1);FPT_num=tau;
        for k=1:nRW
            xa=RW{k}.x;xb=RW{nRW+k}.x;
            ya=RW{k}.y;yb=RW{nRW+k}.y;
            za=RW{k}.z;zb=RW{nRW+k}.z;
            thr=min(RW{k}.setup.thr,RW{nRW+k}.setup.thr);
            t=RW{k}.t;
            [FPT_num(k),FPTcumpt,tau(k),FPTt]=CalFPT(xa,ya,za,xb,yb,zb,t,xc,thr);
            FPTcump=[FPTcump,FPTcumpt];
            FPT=[FPT,FPTt];
            WaitBar(handles,k/nRW,['At ',num2str(k),'-th pair of FPT. ',num2str(k*100/nRW,3),'% done.'])
        end
    case 4 % pairwise
        FPTcump=[];FPT=[];tau=zeros(nRW^2-nRW,1);FPT_num=tau;count=1;
        for k=1:nRW
            for p=k+1:nRW
                xa=RW{k}.x;xb=RW{p}.x;
                ya=RW{k}.y;yb=RW{p}.y;
                za=RW{k}.z;zb=RW{p}.z;
                thr=min(RW{k}.setup.thr,RW{p}.setup.thr);
                t=RW{k}.t;
                [FPT_num(count),FPTcumpt,tau(count),FPTt]=CalFPT(xa,ya,za,xb,yb,zb,t,xc,thr);
                FPTcump=[FPTcump,FPTcumpt];
                FPT=[FPT,FPTt];
                WaitBar(handles,k/nRW,['At ',num2str(k),'-th pair of FPT. ',...
                    num2str(count*200/(nRW^2-nRW),3),'% done.'])
                count=count+1;
            end
        end
    case 5 % split 2 parts, pairwise
        nRW=nRW/2;
        t=RW{nRW-1}.t;FPTcump=[];FPT=[];tau=zeros(nRW,nRW);FPT_num=tau;
        for k=1:nRW
            for p=nRW+1:2*nRW
                xa=RW{k}.x;xb=RW{p}.x;
                ya=RW{k}.y;yb=RW{p}.y;
                za=RW{k}.z;zb=RW{p}.z;
                thr=min(RW{k}.setup.thr,RW{p}.setup.thr);
                [FPT_num(k,p-nRW),FPTcumpt,tau(k,p-nRW),FPTt]=CalFPT(xa,ya,za,xb,yb,zb,t,xc,thr);
                FPTcump=[FPTcump,FPTcumpt];
                FPT=[FPT,FPTt];
                WaitBar(handles,k/nRW,['At ',num2str((k-1)*nRW+p),'-th pair of FPT. ',...
                    num2str(count*200/(nRW^2),3),'% done.'])
            end
        end
    case 6 % customize
        d=0:0.05:0.9;nd=length(d);
        xa=RW{1}.x;xb=RW{2}.x;
        ya=RW{1}.y;yb=RW{2}.y;
        za=RW{1}.z;zb=RW{2}.z;
        thr=min(RW{1}.setup.thr,RW{2}.setup.thr);
        t=RW{nRW-1}.t;
        FPTcump=[];tau=zeros(nd,1);FPT_num=tau;
        for k=1:nd
            xc=[d(k),0,0];
            [FPT_num(k),FPTcumpt,tau(k),FPTt]=CalFPT(xa,ya,za,xb,yb,zb,t,xc,thr);
            FPTcump=[FPTcump(:),FPTcumpt];
            WaitBar(handles,k/nd,['At ',num2str(k),'-th pair of FPT. ',num2str(k*100/nd,3),'% done.'])
        end
end
save(path,'tau','FPTcump','FPT')

function [n,FPTcump,tau,FPT]=CalFPT(xa,ya,za,xb,yb,zb,t,xc,thr)
xb=xb+xc(1);
yb=yb+xc(2);
zb=zb+xc(3);
[~,N1]=size(xa);
[~,N2]=size(ya);
FPT=NaN(N1,N2);
for ko=1:N1      % original trace
    for kt=1:N2	% tartget trace
        DisDot=sqrt(sum(([xa(:,ko),ya(:,ko),za(:,ko)]-[xb(:,kt),yb(:,kt),zb(:,kt)]).^2,2));
        a=find(DisDot<thr,1,'first');
        if ~isempty(a)
            FPT(ko,kt)=a(1);
        end
    end
end
n=N1*N2;
FPT=FPT(:);
FPThist=hist(FPT(:),t);FPThist=FPThist(:);
FPTcum=cumsum(FPThist);FPTcum=FPTcum(:);
FPTcump=FPTcum/n;
[fitresult,~]=DecayFit(t, FPTcump);
tau=fitresult.b;

% --- Executes on button press in FPT_2D.
function [FPT,n]=FPT_2D_Callback(hObject, eventdata, handles)
% hObject    handle to FPT_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RW
nRW=numel(RW);
if nRW==1
    N=floor(RW{1}.setup.N/2);
    xa=RW{1}.x(:,1:N);xb=RW{1}.x(:,N+1:2*N);
    ya=RW{1}.y(:,1:N);yb=RW{1}.y(:,N+1:2*N);
    thr=RW{1}.setup.thr;
elseif numel(RW)>1
    N=RW{nRW-1}.setup.N;
    xa=RW{nRW-1}.x;xb=RW{nRW}.x;
    ya=RW{nRW-1}.y;yb=RW{nRW}.y;
    thr=min(RW{nRW-1}.setup.thr,RW{nRW}.setup.thr);
end
thr=sqrt((thr^2)*2/3);
xc=str2num(get(handles.CenterDist,'String'));
xb=xb+xc(1);
yb=yb+xc(2);
FPT=NaN(N,N);
for ko=1:N      % original trace
    for kt=1:N	% tartget trace
        DisDot=sqrt(sum(([xa(:,ko),ya(:,ko)]-[xb(:,kt),yb(:,kt)]).^2,2));
        a=find(DisDot<thr,1,'first');
        if ~isempty(a)
            FPT(ko,kt)=a(1);
        end
    end
end
n=N^2;

% --- Executes on button press in MakePlot.
function MakePlot_Callback(hObject, eventdata, handles)
% hObject    handle to MakePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RW
global FPT
nRW=numel(RW);
option=get(handles.PlotOption,'Value');
axes(handles.axes1)
switch option
    case 1 % example traces
        for k=1:nRW
            N=RW{k}.setup.N;color=rand(1,3);
            plot(RW{k}.x1(2:end,min(N,5)),RW{k}.x2(2:end,min(N,5)),'color',color),hold on
        end
    case 2 % distribution of dot at t=0
        x=[];y=[];z=[];
        for k=1:nRW
            x=[x,(RW{k}.x(1,:))'];
            y=[y,(RW{k}.y(1,:))'];
            z=[z,(RW{k}.z(1,:))'];
        end
        hist(x,30)
    case 3 % FPT hist
        x=FPT;
        if isempty(FPT)
            for k=1:nRW
                x=[x,RW{k}.FPT(:)];
            end
        end
        hist(x,30)
end
hold off

% --- Executes during object creation, after setting all properties.
function DiffConst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DiffConst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DiffConst_Callback(hObject, eventdata, handles)
% hObject    handle to DiffConst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new DiffConst value
handles.metricdata.density = density;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function N_cells_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function N_cells_Callback(hObject, eventdata, handles)
% hObject    handle to N_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
volume = str2double(get(hObject, 'String'));
if isnan(volume)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new N_cells value
handles.metricdata.volume = volume;
guidata(hObject,handles)

function TDuration_Callback(hObject, eventdata, handles)
% hObject    handle to TDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TDuration as text
%        str2double(get(hObject,'String')) returns contents of TDuration as a double


% --- Executes during object creation, after setting all properties.
function TDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TStep_Callback(hObject, eventdata, handles)
% hObject    handle to TStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function TStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function BCtype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BCtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PlotOption.
function PlotOption_Callback(hObject, eventdata, handles)
% hObject    handle to PlotOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function PlotOption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function BCsize_Callback(hObject, eventdata, handles)
% hObject    handle to BCsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BCsize as text
%        str2double(get(hObject,'String')) returns contents of BCsize as a double


% --- Executes during object creation, after setting all properties.
function BCsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BCsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in BClist.
function BClist_Callback(hObject, eventdata, handles)
% hObject    handle to BClist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BClist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BClist


% --- Executes during object creation, after setting all properties.
function BClist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BClist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function BCcenter_Callback(hObject, eventdata, handles)
% hObject    handle to BCcenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BCcenter as text
%        str2double(get(hObject,'String')) returns contents of BCcenter as a double


% --- Executes during object creation, after setting all properties.
function BCcenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BCcenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ThrDis_Callback(hObject, eventdata, handles)
% hObject    handle to ThrDis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThrDis as text
%        str2double(get(hObject,'String')) returns contents of ThrDis as a double


% --- Executes during object creation, after setting all properties.
function ThrDis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThrDis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InitialPos_Callback(hObject, eventdata, handles)
% hObject    handle to InitialPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InitialPos as text
%        str2double(get(hObject,'String')) returns contents of InitialPos as a double


% --- Executes during object creation, after setting all properties.
function InitialPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InitialPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CenterDist_Callback(hObject, eventdata, handles)
% hObject    handle to CenterDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CenterDist as text
%        str2double(get(hObject,'String')) returns contents of CenterDist as a double


% --- Executes during object creation, after setting all properties.
function CenterDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CenterDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CalOption.
function CalOption_Callback(hObject, eventdata, handles)
% hObject    handle to CalOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CalOption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CalOption


% --- Executes during object creation, after setting all properties.
function CalOption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CalOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
