para.stime = 60;
para.tspan=linspace(0,para.stime,100);

para.dim=3;
para.n=10; % define the population of agents
para.rd = 5;
para.sp = 1/(para.stime*1.2)*2*pi;
para.th1=1;
para.th2=2;
para.r=1.2;
para.la=1;

para.k1=10;
para.k2=1;
para.k3=1;
para.k4=1;
% para.m0=2.5;
para.m0=1;
para.m=1;
para.g=9.8;
para.Gamma=1;

para.ka=3;

n0 = ones(para.dim,1);
[~,th] = f(n0,n0,n0,n0,n0,n0,para);
[Y1,th1] = f1(n0,n0,para);
[Y2,th2] = f2(n0,n0,n0,n0,n0,para);
para.th = [th1;th2];

state = [];
state.x = unifrnd(0,5,[para.dim,para.n]);
state.x(end,:)=ones(size(state.x(end,:)))*(-5);
state.v = unifrnd(-5,5,[para.dim,para.n]);
state.hat_th = zeros([length(th),para.n]);

state.cal_X = zeros(numel(Y1),para.n);
state.cal_Y = zeros(numel(Y2),para.n);
state.cal_T = zeros(para.dim,para.n);

for i=1:para.n
    [Y1,th1] = f1(state.x(:,i),state.v(:,i),para);
    var = (Y1*th1);
    state.cal_X(:,i) = var(:);
end


state.xd=[5;0;0];
state.vd=[0;0.5;0];

data = [];
data.chi_X = '';
data.Y1 = '';
data.chi_Y = '';
data.chi_T = '';

% mem = Memory(data);
% 
% sim = vectorField(state,@(t,s)model(t,s,mem,para));
% sim.setFunctionAfterEachTspan(@(t,state)recording(t,state,mem,para));
% sim.solve('ode45',para.tspan,state);
% [t,result] = sim.result();

figure(1);
clf;
result.q = result.x;
plotFormation_t(t,result,'proposed',[1,para.stime*90/240,para.stime*160/240,para.stime],para)
view([45,25])
title('Flocking trajectories')

axis equal

figure(2)
plotPara(t,result,'a',para)
xlabel('t/s')
ylabel('value')
title('Estimation of model parameters')

figure(3)
plotVel(t,result,'a',para)
xlabel('t/s')
ylabel('ms^{-1}')
title('Velocities of all agents')


function dot_state = model(t,state,mem,para)

x = state.x;
v = state.v;
hat_th = state.hat_th;
xd = state.xd;
vd = state.vd;
cal_X = state.cal_X;
cal_Y = state.cal_Y;
cal_T = state.cal_T;

ad = -0.008*xd;

ddot_q = zeros(size(x));
dot_hat_th = zeros(size(hat_th));

phi = zeros(para.dim,para.n);
% dPhi = zeros(para.dim,para.n);

u_cal_X = state.cal_X.*0;
u_cal_Y = state.cal_Y.*0;
u_cal_T = state.cal_T.*0;

% xd = [   para.rd * cos(t*para.sp);            para.rd * sin(t*para.sp)];
% vd = [  -para.rd * sin(t*para.sp)*para.sp;    para.rd * cos(t*para.sp)*para.sp];
% ad = [  -para.rd * cos(t*para.sp)*para.sp^2; -para.rd * sin(t*para.sp)*para.sp^2];
for i=1:para.n
    xi=x(:,i);
    vi=v(:,i);
    hat_thi=hat_th(:,i);
    for j=1:para.n
        if j==i
            continue
        else
            xj=x(:,j);
            vj=v(:,j);
            dis = norm(xi-xj);
            xij = xi - xj;
            vij = vi - vj;
            phi(:,i)=phi(:,i)+ para.k1*rho(dis^2,para.r^2,0.8)*sig(dis^2,para.la^2,0.2)*xij...
                             + para.k2*rho(dis^2,para.r^2,0.8)*vij;
        end
    end
    phi(:,i)=phi(:,i)+para.k3*(xi-xd);
    
    ri = -vd;
    dot_ri = -ad;
    ei =  ri + vi;
    
    [Y,~]=f(xi,vi,-dot_ri,-ri,xd,vd,para);
    [Y1,th1]=f1(xi,vi,para);
    [Y2,th2]=f2(xi,vi,vi,xd,vd,para);
    ui=-10*ei - 10*phi(:,i) + Y*hat_thi;
    ddot_q(:,i) = 1/para.m*(ui-Y2*th2);
    dot_hat_th(:,i) = -para.Gamma.*Y'*ei;
    
    if mem.n>1
        for k=1:mem.n
            chi_X_k = mem.data.chi_X{k};
            chi_Y_k = mem.data.chi_Y{k};
            chi_T_k = mem.data.chi_T{k};
            Y1_k = mem.data.Y1{k};
            W_k = [reshape(Y1_k(:,i),size(Y1))-reshape(chi_X_k(:,i),size(Y1)),reshape(chi_Y_k(:,i),size(Y2))];
            T_k = chi_T_k(:,i);
            if k==1
                P = W_k'*W_k;
                Q = W_k'*(T_k-W_k*hat_thi);
            else
                P=P+W_k'*W_k;
                Q = Q + W_k'*(T_k-W_k*hat_thi);
            end
            % debug
            if mem.n==length(para.tspan)-2
                th = [th1;th2];
                p1 = W_k'*W_k*(th-hat_thi);
                p2 = W_k'*(T_k-W_k*hat_thi);
                p3 = T_k;
                p4 = W_k*th;
                p5 = P * (th - hat_thi);
                p6 = Q;
                t;
            end
        end
        dot_hat_th(:,i) = dot_hat_th(:,i) + Q;
    end
    
    
    % debug
    
            
    
    % filtered regression
    [Y1,th1] = f1(xi,vi,para);
    [Y2,th2] = f2(xi,vi,vi,xd,vd,para);
    u_cal_X(:,i) = -para.ka * cal_X(:,i) + para.ka * Y1(:);
    u_cal_Y(:,i) = -para.ka * cal_Y(:,i) + Y2(:);
    u_cal_T(:,i) = -para.ka * cal_T(:,i) + ui;
    
    % debug
    if t>3
    p1 = Y1;
    p2 = reshape(cal_X(:,i),size(Y1));
    p3 = reshape(cal_Y(:,i),size(Y2));
    p4 = cal_T(:,i);
    ans1=[p1-p2,p3]*[th1;th2];
    ans2=p4;
    t;
    end
%     assert(norm(p1-p2-p3-p4)<10^-4,'filtered regression fails')
    
    
end
dot_state.x = v;
dot_state.v = ddot_q;
dot_state.hat_th = dot_hat_th;
dot_state.xd = vd;
dot_state.vd = ad;

dot_state.cal_X = u_cal_X;
dot_state.cal_Y = u_cal_Y;
dot_state.cal_T = u_cal_T;
end

function new_state = recording(t,state,mem,para)
new_state = state;
index = mem.n+1;
mem.data.chi_X{index}=state.cal_X;
mem.data.chi_Y{index}=state.cal_Y;

cal_Y1=[];
for i=1:para.n
    [Y1,~]=f1(state.x(:,i),state.v(:,i),para);
    cal_Y1=[cal_Y1,Y1(:)];
end

mem.data.Y1{index}=cal_Y1;
mem.data.chi_T{index}=state.cal_T;
mem.n=mem.n+1;
end

function y=chief(xd,vd,para)
g=para.g;
m0=para.m0;
% y=-1/m0*g/norm(xd)^3*xd;
y=-0.01*xd;
end

function [Y,th] =f(q,dot_q,alpha,beta,xd,vd,para)
[Y1,th1]=f1(q,alpha,para);
[Y2,th2]=f2(q,dot_q,beta,xd,vd,para);
Y = [Y1,Y2];
th= [th1;th2];
end

function [y,a]=f1(x,alpha,para)
y=alpha;
a=[para.m];
end

function [y,a]=f2(x,v,beta,xd,vd,para)
r0=norm(xd);
ri=norm([r0;0;0]+x);
y=zeros(para.dim,para.dim-1);
y(1,1)=-2*beta(2)/r0^(1.5)+3*x(2)/2/r0^(2.5)*xd'*vd/r0;
y(2,1)= 2*beta(1)/r0^(1.5)-3*x(1)/2/r0^(2.5)*xd'*vd/r0;
y(1,2)=-1/r0^3*x(1)+(r0+x(1))/ri^3-1/r0^2;
y(2,2)=-1/r0^3*x(2)+x(2)/ri^3;
y(3,2)=x(3)/ri^3;
a=para.m.*[para.g^0.5;para.g];
end

function [y,a]=f3(x,v,alpha,para)
y=alpha.*0;
a=para.m;
end


function y=rho(x,r,h)
y=smoothLink(x,r*h,1,r,0);
end
function y=dRho(x,r,h)
y=dSmoothLink(x,r*h,1,r,0);
end

function y=sig(x,la,h)
y=smoothLink(x,la-h*la,-1,la+h*la,1);
end
function y=dSig(x,la,h)
y=dSmoothLink(x,la-h*la,-1,la+h*la,1);
end
function y=smoothLink(x,x1,y1,x2,y2)
x_min = min(x1,x2);
if x_min==x1
    y_min = y1;
    y_max = y2; 
else
    y_min = y2;
    y_max = y1; 
end
x_max = max(x1,x2);
if x<x_min
    y=y_min;
elseif x<x_max
    y=((x-x_min)/(x_max-x_min)-sin(2*pi*(x-x_min)/(x_max-x_min))/(2*pi))*(y_max-y_min)+y_min;
else
    y=y_max;
end
end

function y=dSmoothLink(x,x1,y1,x2,y2)
x_min = min(x1,x2);
if x_min==x1
    y_min = y1;
    y_max = y2; 
else
    y_min = y2;
    y_max = y1; 
end
x_max = max(x1,x2);
if x<x_min
    y=0;
elseif x<x_max
    y=(1/(x_max-x_min)-cos(2*pi*(x-x_min)/(x_max-x_min))/(x_max-x_min))*(y_max-y_min);
else
    y=0;
end
end

function plotFormation_t(t,data,methodName,ts,para)

dim=para.dim;



% xd = zeros(length(t),para.dim);
xd = data.xd;
if para.dim==2
plot(data.q(:,1:dim:end),data.q(:,2:dim:end),'b-');hold on
q = data.q;
plot(xd(:,1:dim:end),xd(:,2:dim:end),'r--','linewidth',2);

plot(q(1,1:dim:end),q(1,2:dim:end),'kx');

for k=1:length(ts)
    ts_i = ts(k);
    t_delta = t - ts_i;
    t_delta = abs(t_delta);
    [~,index]=min(t_delta);
    t_i = index;
    plotFormationAtTime(data.q(t_i,:),xd(t_i,:),para)
end

plotFormationAtTime(data.q(end,:),xd(end,:),para)
else
plot3(data.q(:,1:dim:end),data.q(:,2:dim:end),data.q(:,3:dim:end),'b-');hold on
q = data.q;
plot3(xd(:,1:dim:end),xd(:,2:dim:end),xd(:,3:dim:end),'r--','linewidth',2);

plot3(q(1,1:dim:end),q(1,2:dim:end),q(1,3:dim:end),'kx');

for k=1:length(ts)
    ts_i = ts(k);
    t_delta = t - ts_i;
    t_delta = abs(t_delta);
    [~,index]=min(t_delta);
    t_i = index;
    plotFormationAtTime(data.q(t_i,:),xd(t_i,:),para)
end

plotFormationAtTime(data.q(end,:),xd(end,:),para)  
    
end

hold off

grid on
axis equal
end

function plotFormationAtTime(q,xd,para)
dim=para.dim;
n = para.n;

q_mat=reshape(q,[dim,n]);
if dim==2
    
    for i=1:n
        qi=q_mat(:,i);
        for j=1:n
            qj=q_mat(:,j);
            if i==j
                continue
            else
                dis=norm(qi-qj);
                if dis<para.la*0.95
                    line([qi(1);qj(1)],[qi(2);qj(2)],'color',[0.6,0.6,0.6])
                elseif dis<para.r
                    line([qi(1);qj(1)],[qi(2);qj(2)],'color',[0.6,0.6,0.6])
                end
            end
        end
    end
    plot(q(1:dim:end),q(2:dim:end),'o','Color',[1,1,1],'MarkerFaceColor',[0.2,0.6,1]);
    plot(xd(1:dim:end),xd(2:dim:end),'p','Color','r','MarkerFaceColor','r');
else
    for i=1:n
        qi=q_mat(:,i);
        for j=1:n
            qj=q_mat(:,j);
            if i==j
                continue
            else
                dis=norm(qi-qj);
                if dis<para.la*0.95
                    line([qi(1);qj(1)],[qi(2);qj(2)],[qi(3);qj(3)],'color',[0.6,0.6,0.6])
                elseif dis<para.r
                    line([qi(1);qj(1)],[qi(2);qj(2)],[qi(3);qj(3)],'color',[0.6,0.6,0.6])
                end
            end
        end
    end
    plot3(q(1:dim:end),q(2:dim:end),q(3:dim:end),'o','Color',[1,1,1],'MarkerFaceColor',[0.2,0.6,1]);
    plot3(xd(1:dim:end),xd(2:dim:end),xd(3:dim:end),'p','Color','r','MarkerFaceColor','r');
end
end

function plotVel(t,data,methodName,para)
plot(t,data.v-repmat(data.vd,[1,para.n]))

% dim=length(data.xd(end,:));
% plot3(data.q(:,1:dim:end),data.q(:,2:dim:end),data.q(:,3:dim:end));hold on
% plot3(data.xd(:,1:dim:end),data.xd(:,2:dim:end),data.xd(:,3:dim:end),'r--');
% plot3(data.q(end,1:dim:end),data.q(end,2:dim:end),data.q(end,3:dim:end),'o');
% plot3(data.xd(end,1:dim:end),data.xd(end,2:dim:end),data.xd(end,3:dim:end),'rp');
% title(['velocity mismatch (',methodName,')'])
grid on
% axis equal
end

function plotPara(t,data,methodName,para)
plot(t,data.hat_th)
hold on
th_line=[ones(length(t),1).*para.th(1),ones(length(t),1).*para.th(2),ones(length(t),1).*para.th(3)];
plot(t,th_line,'--')
hold off
% title(['parameters estimation (',methodName,')'])
end