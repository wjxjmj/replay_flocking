save_type='.eps';

load('jingyan.mat')
t1=t;
data1=result;
figure(1);
clf;
result.q = result.x;
plotFormation_t(t,result,'proposed',[1,para.stime*90/240,para.stime*160/240,para.stime],para)
view([45,25])

xlabel('x axis')
ylabel('y axis')
set(gcf, 'Color', 'w');
export_fig(['adaptive_formation',save_type])


figure(2)
plotDistanceError(t1,data1,para)
xlabel('t/s')
ylabel('distance')
grid on
set(gcf, 'Color', 'w');
export_fig(['adaptive_distances',save_type])
figure(3)
plotPositionError(t1,data1,para)
xlabel('t/s')
ylabel('position')
grid on
set(gcf, 'Color', 'w');
export_fig(['adaptive_position_error',save_type])
figure(4)
plotSim(t1,data1,'AC',para)
xlabel('t/s')
ylabel('velocity')
grid on
set(gcf, 'Color', 'w');
export_fig(['adaptive_velocity_error',save_type])

figure(5)
plotPara(t,result,'a',para)
xlabel('t/s')
ylabel('estimation')
grid on
set(gcf, 'Color', 'w');
export_fig(['hat_theta',save_type])


% % comparative plot
load('jingyan_noLearning.mat')
t1=t;
data1=result;
figure(6);
clf;
result.q = result.x;
plotFormation_t(t,result,'proposed',[1,para.stime*90/240,para.stime*160/240,para.stime],para)
view([45,25])

xlabel('x axis')
ylabel('y axis')
set(gcf, 'Color', 'w');
export_fig(['adaptive_formation_cmp',save_type])


figure(7)
plotDistanceError(t1,data1,para)
xlabel('t/s')
ylabel('distance')
grid on
set(gcf, 'Color', 'w');
export_fig(['adaptive_distances_cmp',save_type])
figure(8)
plotPositionError(t1,data1,para)
xlabel('t/s')
ylabel('position')
grid on
set(gcf, 'Color', 'w');
export_fig(['adaptive_position_error_cmp',save_type])
figure(9)
plotSim(t1,data1,'AC',para)
xlabel('t/s')
ylabel('velocity')
grid on
set(gcf, 'Color', 'w');
export_fig(['adaptive_velocity_error_cmp',save_type])

figure(10)
plotPara(t,result,'a',para)
xlabel('t/s')
ylabel('estimation')
grid on
set(gcf, 'Color', 'w');
export_fig(['hat_theta_cmp',save_type])


function plotDistanceError(t,data,para)
% xd = [cos(t*1.8*pi/para.stime)*para.rr,...
%       sin(t*1.8*pi/para.stime)*para.rr];
eq = zeros([length(t),para.n]);
for ti=1:length(t)
    q_t = reshape(data.q(ti,:),[para.dim,para.n]);
    for i=1:para.n
        qi = q_t(:,i);
        diss = [];
        for j=1:para.n
            if i==j
                continue
            else
                qj = q_t(:,j);
                diss = [diss,norm(qi-qj)];
            end
        end
        min_dis = min(diss);
        eq(ti,i) = min_dis;
    end
end
plot(t,eq,'-');
end

function plotPositionError(t,data,para)
xd = data.xd;
eq = zeros(size(xd));
for i=1:length(t)
    q_i = reshape(data.q(i,:),[para.dim,para.n]);
    q_e = mean(q_i,2);
    eq(i,:) = (q_e' - xd(i,:));
end
plot(t,eq,'-');
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

function plotSim(t,data,methodName,para)

vd=zeros(length(t),para.dim);
for i=1:length(t)
    vd_t = data.vd(i,:)';
    vd(i,:)=vd_t';
end
data.vd=vd;

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

function [xd,vd,ad] = virtualLeader(t,para)
if para.dim==3
d_angle = 4*pi/para.stime;
angle = t * d_angle;
rs = para.rr;
rh = para.rh;
d_height = para.rh/para.stime;
height = t * d_height;
xd=[cos(angle)*rs;sin(angle)*rs;height];
vd=[-sin(angle)*d_angle*rs;cos(angle)*d_angle*rs;d_height];
ad=[-cos(angle)*(d_angle)^2*rs;-sin(angle)*(d_angle)^2*rs;0];
elseif para.dim==2
rr = para.rr;
xd = [cos(t*1.8*pi/para.stime)*rr;
      sin(t*1.8*pi/para.stime)*rr];
vd = [-sin(t*1.8*pi/para.stime)*rr*1.8*pi/para.stime;
      cos(t*1.8*pi/para.stime)*rr*1.8*pi/para.stime];
ad = [-cos(t*1.8*pi/para.stime)*rr*(1.8*pi/para.stime)^2;
      -sin(t*1.8*pi/para.stime)*rr*(1.8*pi/para.stime)^2];
end
end
function plotPara(t,data,methodName,para)
plot(t,data.hat_th)
hold on
th_line=[ones(length(t),1).*para.th(1),ones(length(t),1).*para.th(2),ones(length(t),1).*para.th(3)];
plot(t,th_line,'--')
hold off
% title(['parameters estimation (',methodName,')'])
end
 