function mySurfStatColormap( map );
 
%Colormap function for SurfStatView.
%
% Usage: SurfStatColormap( map );
%
% Same as for matlab's colormap function - see help colormap.
if strcmp(map, 'my_rwb')
   load my_rwb
   colormap(my_rwb);
elseif strcmp(map, 'my_rwb1')
   load my_rwb1
   colormap(my_rwb1);
elseif strcmp(map, 'my_rwb2')
   load my_rwb2
   colormap(my_rwb2);
elseif strcmp(map, 'sym_rwb')
   load sym_rwb
   colormap(sym_rwb);
elseif strcmp(map, 'tmap')
   load fourvaluescolor
   colormap(fourvaluescolor);
elseif strcmp(map, 'hot1')
   load hot1
   colormap(hot1);
else
       colormap(map);
end
a=get(gcf,'Children');
k=0;
for i=1:length(a)
    tag=get(a(i),'Tag');
    if strcmp(tag,'Colorbar')
        cb=a(i);
    end
    if length(tag)>12 & strcmp(tag(1:12),'SurfStatView')
        k=k+1;
    end
end
 
% if k>1
    set(cb,'location','South');
    set(cb,'Position',[0.1 0.04 0.7 0.03]);
    set(cb,'XAxisLocation','bottom');
% end
 
return
end