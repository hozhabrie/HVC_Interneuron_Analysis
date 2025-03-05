function mapEdges()
col=[1,1,1]*.5;
rs=[1.5:2:13.5];
for r=1:length(rs)
   line([0,9],[1,1]*rs(r),'color',col,'linewidth',1)
end
cs=[2.5,4.5,6.5];
for c=1:3
   line([1,1]*cs(c),[0,16],'color',col,'linewidth',1)
end