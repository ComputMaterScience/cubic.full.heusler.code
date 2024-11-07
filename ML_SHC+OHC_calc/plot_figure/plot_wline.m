function plot_wline(x,y,w,span,factor,color)
w = smooth(w,span)*factor;
hold on;
r = 0.85:0.01:1;
num_a = linspace(0.3,1,size(r,2));
cnt = 0;
for i = r
    cnt = cnt + 1;
    %h = fill([x;flipud(x)],[(y+(1.01-i)*w);flipud((y-(1.01-i)*w))],color,'EdgeColor','none');
    h = fill([x;flipud(x)],[(y+(1.01-i)*w);flipud((y-(1.01-i)*w))],color,'EdgeColor','none');
    set(h,'FaceAlpha',num_a(cnt));
end
%plot(x,y,'k','LineWidth',3*factor);
end